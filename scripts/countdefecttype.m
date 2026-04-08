%% countdefecttype.m
% Classify each detected defect as MF / HDF / Unclassified based on local
% density maps, and export per-frame + total counts and a summary plot.
%
% Expected input directory structure (relative to repo root):
%   Data/<run_name>/
%       defectData.mat
%       segmentation/density_maps_frame*.mat
%       defectimages/director_data_*.mat
%
% Outputs:
%   Data/<run_name>/defect_type_counts.csv
%   Results/defecttype.tif
%
% Copyright (c) 2025 Zhaofei Zheng
% Released under the MIT License. See LICENSE in the repository root.

close all; clear;

% -------------------- Repo-safe paths (do not use pwd) --------------------
repoRoot = fileparts(fileparts(mfilename("fullpath"))); % scripts/ -> repo root

% Input dataset folder (relative to repo root)
folderName = fullfile('Data/sample', 'live_cell_imaging_every30min_10x_4tiles_750_50_03');
baseDir    = fullfile(repoRoot, folderName);

% Output folder for figures (relative to repo root)
resultPath = fullfile(repoRoot, 'Results');
if ~exist(resultPath, 'dir')
    mkdir(resultPath);
end

% --- Load defect information (all frames) ---
defectMatPath = fullfile(baseDir, 'defectData.mat');
D = load(defectMatPath);

% Parameter setting
edgeSize = 0;           % pixels trimmed from each border
classRadius = 50;       % neighborhood radius (px) for classification
classMinDiff = 8e-4;    % confidence threshold |G-R| below this => Unknown
classMethod  = 'mean';

defectX_all      = D.defectX(:);
defectY_all      = D.defectY(:);
defectCharge_all = D.defectCharge(:);
defectFrame_all  = D.defectFrames(:);

% Get total number of frames
vals = unique(defectFrame_all);
numFrame_total = numel(vals);

% Pre-allocate counters
countGreenPos = zeros(1, numFrame_total);
countGreenNeg = zeros(1, numFrame_total);
countPinkPos  = zeros(1, numFrame_total);
countPinkNeg  = zeros(1, numFrame_total);
countUnknown  = zeros(1, numFrame_total);

for frameNum = 80:100
    % --- Load density maps & directors for this frame ---
    densityMatPath = fullfile(baseDir, 'segmentation', sprintf('density_maps_frame%d.mat', frameNum));
    directorsPath  = fullfile(baseDir, 'defectimages', sprintf('director_data_%d.mat', frameNum));

    % Directors (not used for classification here, but kept to match your flow)
    DD = load(directorsPath);
    if ~isfield(DD, 'directors')
        error('director_data_%d.mat does not contain variable "directors".', frameNum);
    end
    directors = DD.directors; 

    % Density maps
    S = load(densityMatPath);   % expects G_cells_per_um2, R_cells_per_um2
    if ~isfield(S, 'G_cells_per_um2') || ~isfield(S, 'R_cells_per_um2')
        error('density_maps_frame%d.mat missing required variables.', frameNum);
    end
    densityMapG = S.G_cells_per_um2;
    densityMapR = S.R_cells_per_um2;

    densityMapG = densityMapG * mean2(densityMapR) / mean2(densityMapG);

    % --- Select defects in this frame ---
    idxFrame     = (defectFrame_all == frameNum);
    defectX      = defectX_all(idxFrame);
    defectY      = defectY_all(idxFrame);
    defectCharge = defectCharge_all(idxFrame);

    % Remove defects near edges
    imgSize = size(densityMapG);
    idxCenter_X = defectX < (imgSize(1) - edgeSize) & defectX > (edgeSize + 1);
    idxCenter_Y = defectY < (imgSize(2) - edgeSize) & defectY > (edgeSize + 1);
    idxCenter = idxCenter_X & idxCenter_Y;

    defectX      = floor(defectX(idxCenter));   % NOTE: using your original row/col convention
    defectY      = floor(defectY(idxCenter));
    defectCharge = defectCharge(idxCenter);

    % --- Classify each defect location by local (G vs R) dominance ---
    defectNum = numel(defectX);
    for t = 1:defectNum
        coord = [defectX(t), defectY(t)];  % [row, col]; if your data are (x,y), swap to [y,x]

        colorClass = classifyPointColor( ...
            densityMapG, densityMapR, coord, ...
            'Radius',   classRadius, ...
            'MinDiff',  classMinDiff, ...
            'Method',   classMethod);

        % Charge class
        if defectCharge(t) == 0.5
            chargeClass = "Positive";
        else
            chargeClass = "Negative";
        end

        % Increment counters
        if colorClass == "Green" && chargeClass == "Positive"
            countGreenPos(frameNum) = countGreenPos(frameNum) + 1;
        elseif colorClass == "Green" && chargeClass == "Negative"
            countGreenNeg(frameNum) = countGreenNeg(frameNum) + 1;
        elseif colorClass == "Pink"  && chargeClass == "Positive"
            countPinkPos(frameNum)  = countPinkPos(frameNum)  + 1;
        elseif colorClass == "Pink"  && chargeClass == "Negative"
            countPinkNeg(frameNum)  = countPinkNeg(frameNum)  + 1;
        else
            countUnknown(frameNum)  = countUnknown(frameNum)  + 1;
        end
    end
end

% --- Per-frame table ---
T = table((1:numel(countGreenPos))', countGreenPos', countGreenNeg', ...
          countPinkPos', countPinkNeg', countUnknown', ...
    'VariableNames', {'Frame','GreenPos','GreenNeg','PinkPos','PinkNeg','Unclassified'});

% --- Totals row ---
totalGreenPos = sum(countGreenPos);
totalGreenNeg = sum(countGreenNeg);
totalPinkPos  = sum(countPinkPos);
totalPinkNeg  = sum(countPinkNeg);
totalUnclassified = sum(countUnknown);

Ttotals = table("Total", totalGreenPos, totalGreenNeg, totalPinkPos, totalPinkNeg, totalUnclassified, ...
    'VariableNames', {'Frame','GreenPos','GreenNeg','PinkPos','PinkNeg','Unclassified'});

% --- Combine & Save ---
Tall   = [T; Ttotals];
outFile = fullfile(baseDir, 'defect_type_counts.csv');
writetable(Tall, outFile);

Tsummary = table(totalGreenPos, totalGreenNeg, totalPinkPos, totalPinkNeg, totalUnclassified);
disp(Tsummary);

%% --- Draw histogram plot ---
vals = table2array(Tsummary);
barLabels = {'MF +1/2', 'MF -1/2', 'HDF +1/2', 'HDF -1/2', 'Unclassified'};

% Custom colors (RGB 0–1)
lightGreen   = [0.6 1.0 0.6];
green        = [0.0 0.6 0.0];
magenta      = [1.0 0.0 1.0];
lightMagenta = [1.0 0.6 1.0];
grayColor    = [0.5 0.5 0.5];

colors = [
    lightGreen;
    green;
    magenta;
    lightMagenta;
    grayColor
];

% -----------------------------
% Plot histogram (bar plot)
% -----------------------------
figure('Position',[100 100 1200 800], 'Color', 'w', 'Visible','off');

b = bar(vals, 'FaceColor', 'flat', 'LineWidth', 2);
pbaspect([1 1 1]);

for i = 1:numel(vals)
    b.CData(i,:) = colors(i,:);
end

set(gca, 'FontSize', 24, 'LineWidth', 2, ...
         'FontName', 'Times New Roman');

set(gca, 'XTick', 1:5, 'XTickLabel', barLabels, ...
         'FontName', 'Times New Roman', 'FontSize', 50, ...
         'LineWidth', 2);

ylabel('Counts', 'FontName', 'Times New Roman', 'FontSize', 50);
box on;

% resultFilename = fullfile(resultPath, 'defecttype.tif');
% exportgraphics(gcf, resultFilename, 'Resolution', 300);

%% Parameter sampling by changing the threshold
threshold = [1e-4, 2e-4, 4e-4, 6e-4, 8e-4, 1e-3];
type1 = [42, 39, 38, 37, 36, 35];
type2 = [68, 64, 61, 57, 57, 51];
type3 = [57, 55, 52, 50, 47, 43];
type4 = [40, 38, 33, 30, 27, 27];
unclassified = [4, 15, 27, 37, 44, 55];

figure('Color', 'w'); 
hold on;

plot(threshold, type1, '-o', ...
    'LineWidth', 2, 'MarkerSize', 8, ...
    'Color', lightGreen, ...
    'MarkerFaceColor', lightGreen);

plot(threshold, type2, '-s', ...
    'LineWidth', 2, 'MarkerSize', 8, ...
    'Color', green, ...
    'MarkerFaceColor', green);

plot(threshold, type3, '-^', ...
    'LineWidth', 2, 'MarkerSize', 8, ...
    'Color', magenta, ...
    'MarkerFaceColor', magenta);

plot(threshold, type4, '-d', ...
    'LineWidth', 2, 'MarkerSize', 8, ...
    'Color', lightMagenta, ...
    'MarkerFaceColor', lightMagenta);

plot(threshold, unclassified, '-d', ...
    'LineWidth', 2, 'MarkerSize', 8, ...
    'Color', grayColor, ...
    'MarkerFaceColor', grayColor);

set(gca, 'XScale', 'log');
set(gca, 'XDir', 'reverse');

xlabel('Threshold', 'FontName', 'Times New Roman', 'FontSize', 28);
ylabel('Count',     'FontName', 'Times New Roman', 'FontSize', 28);

set(gca, 'FontName', 'Times New Roman', ...
         'FontSize', 24, ...
         'LineWidth', 2, ...
         'TickDir', 'out');

legend({'MF +1/2','MF -1/2','HDF +1/2','HDF -1/2','Unclassfied'}, ...
    'Location', 'eastoutside', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 20, ...
    'Box', 'off');

box on;

%% Plot the histogram with error bars
vals = [mean(type1), mean(type2), mean(type3), mean(type4)];
errs = [std(type1), std(type2), std(type3), std(type4)];

barLabels4 = barLabels(1:4);
colors4 = colors(1:4,:);

figure('Position',[100 100 1800 1000], 'Color', 'w');

x = 1:4;

b = bar(x, vals, 0.55, 'FaceColor', 'flat', 'LineWidth', 2);
hold on;

for i = 1:numel(vals)
    b.CData(i,:) = colors4(i,:);
end

errorbar(x, vals, errs, ...
    'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 18);

% -----------------------------
% Set axes size manually
% -----------------------------
ax = gca;
left   = 0.12;
bottom = 0.22;   % give more room for x labels
width  = 0.6;
height = 0.50;   % visual 2:1 style without distortion
ax.Position = [left bottom width height];

% -----------------------------
% Formatting
% -----------------------------
set(gca, 'XTick', x, ...
         'XTickLabel', barLabels4, ...
         'FontName', 'Times New Roman', ...
         'FontSize', 40, ...
         'LineWidth', 4, ...
         'TickDir', 'out', ...
         'Box', 'on' );

ylabel('Counts', 'FontName', 'Times New Roman', 'FontSize', 45);
xlim([0.4 4.6]);
ylim([0 70]);

box on;

resultFilename = fullfile(resultPath, 'defecttype_4bars_with_error.tif');
exportgraphics(gcf, resultFilename, 'Resolution', 300);

% ===================== Helper =====================
function colorClass = classifyPointColor(densityG, densityR, coord, varargin)
% Returns "Green", "Pink", or "Unknown" based on local G vs R around coord.
% Supports distance-weighted aggregation.
%
% Kernel options for weighted methods ("wmean","wtrim","wmedian"):
%   "gaussian" (default), "invdist", "linear", "tricubic"
%   Tunables: 'Sigma' (for gaussian), 'Power' (for invdist), 'TrimQuantile' (for wtrim)

p = inputParser;
p.addParameter('Radius',    5,        @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MinDiff',   5e-4,     @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Method',    'wmean',  @(s)ischar(s) || isstring(s));
p.addParameter('MinPixels', [],       @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% Weighted options
p.addParameter('Kernel',    'gaussian', @(s)ischar(s)||isstring(s));
p.addParameter('Sigma',     [],        @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Power',     2,         @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('TrimQuantile', 0.1,    @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<0.5);
p.parse(varargin{:});

R         = p.Results.Radius;
minDiff   = p.Results.MinDiff;
method    = lower(string(p.Results.Method));
kernel    = lower(string(p.Results.Kernel));
sigma     = p.Results.Sigma;
power     = p.Results.Power;
qtrim     = p.Results.TrimQuantile;

[nRows, nCols] = size(densityG);
if ~isequal(size(densityG), size(densityR))
    error('densityG and densityR size mismatch');
end

row = round(coord(1)); 
col = round(coord(2));
if row < 1 || row > nRows || col < 1 || col > nCols
    colorClass = "Unknown"; 
    return;
end

if isempty(p.Results.MinPixels)
    MinPixels = max(10, ceil(pi*R^2/4));
else
    MinPixels = p.Results.MinPixels;
end
if isempty(sigma)
    sigma = max(R/2, 1);
end

% Window + distances
rmin = max(1, row - R); 
rmax = min(nRows, row + R);
cmin = max(1, col - R); 
cmax = min(nCols, col + R);

[cc, rr] = meshgrid(cmin:cmax, rmin:rmax);
dr = rr - row; 
dc = cc - col;
d2 = dr.^2 + dc.^2;
mask = d2 <= R^2;

Gpatch = densityG(rmin:rmax, cmin:cmax);
Rpatch = densityR(rmin:rmax, cmin:cmax);
valid = mask & ~isnan(Gpatch) & ~isnan(Rpatch);

if nnz(valid) < MinPixels
    colorClass = "Unknown"; 
    return;
end

gVals = Gpatch(valid);
rVals = Rpatch(valid);
dist = sqrt(d2(valid));

% Build distance weights
w = ones(size(dist));
switch kernel
    case "gaussian"
        w = exp(-(dist.^2) / (2*sigma^2));
    case "invdist"
        eps0 = 1e-6;
        w = 1 ./ max(dist, eps0).^power;
    case "linear"
        w = max(0, 1 - dist / max(R, 1));
    case "tricubic"
        rrel = min(dist / max(R,1), 1);
        w = (1 - rrel.^3).^3;
    otherwise
        error('Unknown Kernel: %s', kernel);
end

wsum = sum(w);
if wsum <= 0
    colorClass = "Unknown"; 
    return;
end
w = w / wsum;

% Aggregators
switch method
    case "mean"
        Gagg = mean(gVals); 
        Ragg = mean(rVals);
    case "median"
        Gagg = median(gVals); 
        Ragg = median(rVals);
    case "wmean"
        Gagg = sum(w .* gVals); 
        Ragg = sum(w .* rVals);
    case "wtrim"
        [Gagg, Ragg] = weightedTrimmedMean(gVals, rVals, w, qtrim);
    case "wmedian"
        Gagg = weightedQuantile(gVals, w, 0.5);
        Ragg = weightedQuantile(rVals, w, 0.5);
    otherwise
        error('Method must be "mean", "median", "wmean", "wtrim", or "wmedian".');
end

diffGR = Gagg - Ragg;
if abs(diffGR) < minDiff
    colorClass = "Unknown";
elseif diffGR > 0
    colorClass = "Green";
else
    colorClass = "Pink";
end
end

% --------- helpers ----------
function [mG, mR] = weightedTrimmedMean(g, r, w, q)
loG = weightedQuantile(g, w, q);   
hiG = weightedQuantile(g, w, 1-q);
loR = weightedQuantile(r, w, q);   
hiR = weightedQuantile(r, w, 1-q);

maskG = (g >= loG) & (g <= hiG);
maskR = (r >= loR) & (r <= hiR);

wG = w(maskG); gT = g(maskG); 
wR = w(maskR); rT = r(maskR);

wG = wG / sum(wG); 
wR = wR / sum(wR);

mG = sum(wG .* gT);
mR = sum(wR .* rT);
end

function qv = weightedQuantile(x, w, q)
tmp = sortrows([x(:), w(:)], 1);
x_sorted = tmp(:,1); 
w_sorted = tmp(:,2);

cw = cumsum(w_sorted);
cw = cw / cw(end);

qv = x_sorted(find(cw >= q, 1, 'first'));
end