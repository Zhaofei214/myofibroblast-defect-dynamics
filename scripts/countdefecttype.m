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

clear
clc

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
edgeSize = 0;         % pixels trimmed from each border
classRadius = 50;       % neighborhood radius (px) for classification
classMinDiff = 5e-4;   % confidence threshold |G-R| below this => Unknown
classMethod  = 'mean';

defectX_all      = D.defectX(:);
defectY_all      = D.defectY(:);
defectCharge_all = D.defectCharge(:);
defectFrame_all  = D.defectFrames(:);

% Get total number of frames
vals = unique(defectFrame_all);
numFrame_total = numel(vals);

% Pre-allocate counters
countGreenPos = zeros(1,numFrame_total);
countGreenNeg = zeros(1,numFrame_total);
countPinkPos  = zeros(1,numFrame_total);
countPinkNeg  = zeros(1,numFrame_total);
countUnknown  = zeros(1,numFrame_total);

for frameNum = 80:100
    % --- Load density maps & directors for this frame ---
    densityMatPath = fullfile(baseDir, 'segmentation', sprintf('density_maps_frame%d.mat', frameNum));
    directorsPath  = fullfile(baseDir, 'defectimages', sprintf('director_data_%d.mat', frameNum));
    
    % Directors (not used for classification here, but kept to match your flow)
    DD = load(directorsPath);
    if ~isfield(DD,'directors')
        error('director_data_%d.mat does not contain variable "directors".', frameNum);
    end
    directors = DD.directors; 
    
    % Density maps
    S = load(densityMatPath);   % expects G_cells_per_um2, R_cells_per_um2
    if ~isfield(S,'G_cells_per_um2') || ~isfield(S,'R_cells_per_um2')
        error('density_maps_frame%d.mat missing required variables.', frameNum);
    end
    densityMapG = S.G_cells_per_um2;
    densityMapR = S.R_cells_per_um2;

    densityMapG =  densityMapG*mean2(densityMapR)/mean2(densityMapG);

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

        % Use circular neighborhood classification
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
    'VariableNames', {'Frame','GreenPos','GreenNeg','PinkPos','PinkNeg','Unknown'});

% --- Totals row ---
totalGreenPos = sum(countGreenPos);
totalGreenNeg = sum(countGreenNeg);
totalPinkPos  = sum(countPinkPos);
totalPinkNeg  = sum(countPinkNeg);
totalUnknown  = sum(countUnknown);

Ttotals = table("Total", totalGreenPos, totalGreenNeg, totalPinkPos, totalPinkNeg, totalUnclassified, ...
    'VariableNames', {'Frame','GreenPos','GreenNeg','PinkPos','PinkNeg','Unclassified'});

% --- Combine & Save ---
Tall   = [T; Ttotals];
outFile = fullfile(baseDir, 'defect_type_counts.csv');
writetable(Tall, outFile);
T = table(totalGreenPos, totalGreenNeg, totalPinkPos, totalPinkNeg, totalUnclassified);
disp(T);

%% --- Draw histogram plot ---
% Extract numeric values as an array
vals = table2array(T);
% Labels for each bar
barLabels = {'MF -1/2', 'MF -1/2', 'HDF +1/2', 'HDF +1/2', 'Unclassified'};

% Custom colors (RGB 0–1)
lightGreen   = [0.6 1.0 0.6];
green        = [0.0 0.6 0.0];
magenta      = [1.0 0.0 1.0];
lightMagenta = [1.0 0.6 1.0];
grayColor    = [0.5 0.5 0.5];   % for "unknown" if you want

% Combine into color map
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
% Make a large figure window (width x height in pixels)
figure('Position',[100 100 1200 900]);   % <-- Increase size here

b = bar(vals,'FaceColor','flat','LineWidth',2);
pbaspect([1 1 1]);

% Apply colors to each bar
for i = 1:numel(vals)
    b.CData(i,:) = colors(i,:);
end

% Optional: nicer axes formatting
set(gca, 'FontSize', 24, 'LineWidth', 2, ...
         'FontName', 'Times New Roman');

% -----------------------------
% Axes formatting
% -----------------------------
set(gca, 'XTick', 1:5, 'XTickLabel', barLabels, ...
         'FontName','Times New Roman', 'FontSize',50, ...
         'LineWidth',2);

ylabel('Counts','FontName','Times New Roman','FontSize',50);
%ylim([0 70]);

% Save the figure (PNG, 300 dpi)
resultFilename = fullfile(resultPath, 'defecttype.tif');
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
p.addParameter('Sigma',     [],        @(x)isnumeric(x)&&isscalar(x)&&x>0);  % default set later
p.addParameter('Power',     2,         @(x)isnumeric(x)&&isscalar(x)&&x>0);  % for invdist
p.addParameter('TrimQuantile', 0.1,    @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<0.5);
p.parse(varargin{:});

R         = p.Results.Radius;
minDiff   = p.Results.MinDiff;
method    = lower(string(p.Results.Method));
kernel    = lower(string(p.Results.Kernel));
sigma     = p.Results.Sigma;   % if empty, set to R/2 below
power     = p.Results.Power;
qtrim     = p.Results.TrimQuantile;

[nRows, nCols] = size(densityG);
if ~isequal(size(densityG), size(densityR)), error('densityG and densityR size mismatch'); end

row = round(coord(1)); col = round(coord(2));
if row < 1 || row > nRows || col < 1 || col > nCols
    colorClass = "Unknown"; return;
end

% Default MinPixels ~ quarter of the circle area, but at least 10
if isempty(p.Results.MinPixels)
    MinPixels = max(10, ceil(pi*R^2/4));
else
    MinPixels = p.Results.MinPixels;
end
if isempty(sigma), sigma = max(R/2, 1); end  % reasonable default width

% Window + distances
rmin = max(1, row - R); rmax = min(nRows, row + R);
cmin = max(1, col - R); cmax = min(nCols, col + R);
[cc, rr] = meshgrid(cmin:cmax, rmin:rmax);
dr = rr - row; dc = cc - col;
d2 = dr.^2 + dc.^2;
mask = d2 <= R^2;

Gpatch = densityG(rmin:rmax, cmin:cmax);
Rpatch = densityR(rmin:rmax, cmin:cmax);
valid = mask & ~isnan(Gpatch) & ~isnan(Rpatch);

if nnz(valid) < MinPixels
    colorClass = "Unknown"; return;
end

gVals = Gpatch(valid);
rVals = Rpatch(valid);
dist = sqrt(d2(valid));

% Build distance weights
w = ones(size(dist));
switch kernel
    case "gaussian"
        % center gets largest weight; decays with distance
        w = exp(-(dist.^2) / (2*sigma^2));
    case "invdist"
        % inverse-distance^power; regularize center to avoid Inf
        eps0 = 1e-6;
        w = 1 ./ max(dist, eps0).^power;
    case "linear"
        % linearly decreasing to 0 at radius R
        w = max(0, 1 - dist / max(R, 1));
    case "tricubic"
        % smooth compact support: (1 - (r/R)^3)^3  for r<R
        rrel = min(dist / max(R,1), 1);
        w = (1 - rrel.^3).^3;
    otherwise
        error('Unknown Kernel: %s', kernel);
end

% Normalize weights
wsum = sum(w);
if wsum <= 0, colorClass = "Unknown"; return; end
w = w / wsum;

% Aggregators
switch method
    case "mean"
        Gagg = mean(gVals); Ragg = mean(rVals);
    case "median"
        Gagg = median(gVals); Ragg = median(rVals);
    case "wmean"
        Gagg = sum(w .* gVals); Ragg = sum(w .* rVals);
    case "wtrim"
        % Weighted trimmed mean: drop tails by weighted quantiles
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
% Drop lower q and upper q by weighted quantiles, then re-average
loG = weightedQuantile(g, w, q);   hiG = weightedQuantile(g, w, 1-q);
loR = weightedQuantile(r, w, q);   hiR = weightedQuantile(r, w, 1-q);
maskG = (g >= loG) & (g <= hiG);
maskR = (r >= loR) & (r <= hiR);
wG = w(maskG); gT = g(maskG); 
wR = w(maskR); rT = r(maskR);
wG = wG / sum(wG); wR = wR / sum(wR);
mG = sum(wG .* gT);
mR = sum(wR .* rT);
end

function qv = weightedQuantile(x, w, q)
% x values, weights w (nonnegative), target quantile q in [0,1]
[w, idx] = sortrows([x(:), w(:)], 1); % sort by x
x_sorted = w(:,1); w_sorted = w(:,2);
cw = cumsum(w_sorted);
cw = cw / cw(end);
qv = x_sorted(find(cw >= q, 1, 'first'));
end