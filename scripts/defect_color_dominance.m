%% defect_color_dominance.m
% This script classifies defects by local G/R dominance (from binarized density maps)
% and defect charge (+1/2 or -1/2), then counts each class per frame.
%
% Expected repo structure (relative to repo root):
%   data/sample/defect_color_dominance/
%           defectData.mat
%           segmentation/density_maps_frame*.mat
%           defectimages/director_data_*.mat
%
% Outputs:
%   Results/defect_type_counts_edge.csv
%
% Notes:
% - Repo-safe paths (no pwd / no absolute paths).
% - The output is written to Results/ (not inside data/).
%
% Copyright (c) 2026 Zhaofei Zheng
% Released under the MIT License. See LICENSE in the repository root.

clear; close all

%% -------------------- Repo-safe paths --------------------
repoRoot   = fileparts(fileparts(mfilename("fullpath"))); % scripts/ -> repo root
resultPath = fullfile(repoRoot, 'Results');
if ~exist(resultPath, 'dir'); mkdir(resultPath); end

%% -------------------- Input folder  --------------------
runName    = 'defect_color_dominance';
dataDir    = fullfile(repoRoot, 'data', 'sample', runName);

folderName = 'live_cell_imaging_every30min_10x_4tiles_750_50_03';
baseDir    = fullfile(dataDir, folderName);

%% -------------------- Load defect information --------------------
defectMatPath = fullfile(baseDir, 'defectData.mat');
assert(exist(defectMatPath,'file')==2, 'Missing file: %s', defectMatPath);
D = load(defectMatPath);

% Parameter setting
edgeSize     = 0;        % pixels trimmed from each border
classRadius  = 10;       % neighborhood radius (px) for classification
classMinDiff = 0.0001;   % confidence threshold |(G-R)/(G+R)| below this => Unknown
classMethod  = 'mean';

defectX_all      = D.defectX(:);
defectY_all      = D.defectY(:);
defectCharge_all = D.defectCharge(:);
defectFrame_all  = D.defectFrames(:);

% Frames present in defect table
vals = unique(defectFrame_all);
numFrame_total = numel(vals);

% Pre-allocate counters (index by actual frame number, so use max frame)
maxFrame = max(vals);
countGreenPos = zeros(1, maxFrame);
countGreenNeg = zeros(1, maxFrame);
countPinkPos  = zeros(1, maxFrame);
countPinkNeg  = zeros(1, maxFrame);
countUnknown  = zeros(1, maxFrame);

%% -------------------- Choose frames to analyze --------------------
% EDIT HERE if you want a different frame range:
framesToRun = 60:10:90;

for frameNum = framesToRun

    % --- Load density maps & directors for this frame ---
    densityMatPath = fullfile(baseDir, 'segmentation', sprintf('density_maps_frame%d.mat', frameNum));
    directorsPath  = fullfile(baseDir, 'defectimages', sprintf('director_data_%d.mat', frameNum));

    assert(exist(densityMatPath,'file')==2, 'Missing file: %s', densityMatPath);
    assert(exist(directorsPath,'file')==2,  'Missing file: %s', directorsPath);

    % Directors (not used for classification here, but loaded to match pipeline)
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

    % ---- binarized maps (z-scored then threshold at 0) ----
    densityMapG = (S.G_cells_per_um2 - mean(S.G_cells_per_um2(:))) / std(S.G_cells_per_um2(:)) > 0;
    densityMapR = (S.R_cells_per_um2 - mean(S.R_cells_per_um2(:))) / std(S.R_cells_per_um2(:)) > 0;

    % --- Select defects in this frame ---
    idxFrame     = (defectFrame_all == frameNum);
    defectX      = defectX_all(idxFrame);
    defectY      = defectY_all(idxFrame);
    defectCharge = defectCharge_all(idxFrame);

    % Remove defects near edges
    imgSize = size(densityMapG);
    idxCenter_X = defectX < (imgSize(1) - edgeSize) & defectX > (edgeSize + 1);
    idxCenter_Y = defectY < (imgSize(2) - edgeSize) & defectY > (edgeSize + 1);
    idxCenter   = idxCenter_X & idxCenter_Y;

    defectX      = floor(defectX(idxCenter));  % keep your original row/col convention
    defectY      = floor(defectY(idxCenter));
    defectCharge = defectCharge(idxCenter);

    % Optional visualization (comment out for batch runs)
    figure('Color','w');
    xlim([0 800]);
    ylim([0 800]);
    set(gca, 'YDir', 'reverse');
    pbaspect([1 1 1]);

    % ---- Axis style: 4 sides, no axis numbers ----
    set(gca, ...
        'Box', 'on', ...              % show all 4 sides
        'LineWidth', 2, ...
        'TickDir', 'out', ...
        'XTickLabel', [], ...         % remove x-axis numbers
        'YTickLabel', [], ...         % remove y-axis numbers
        'FontName', 'Times New Roman');
        % --- Classify each defect location by local (G vs R) dominance ---
        defectNum = numel(defectX);
        for t = 1:defectNum
            coord = [defectX(t), defectY(t)];  % [col, row] in this script's convention
    
            colorClass = classifyPointColor( ...
                densityMapG, densityMapR, coord, ...
                'Radius',   classRadius, ...
                'MinDiff',  classMinDiff, ...
                'Method',   classMethod);
    
            text(defectX(t)+2, defectY(t), colorClass, ...
                'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold', ...
                'HorizontalAlignment','left', 'VerticalAlignment','middle');
    
            % Charge class
            if defectCharge(t) == 0.5
                chargeClass = "Positive";
            else
                chargeClass = "Negative";
            end
    
            % Increment counters (indexed by frameNum)
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

%% -------------------- Per-frame table + totals --------------------
frameList = (1:maxFrame).';
T = table(frameList, countGreenPos(:), countGreenNeg(:), countPinkPos(:), countPinkNeg(:), countUnknown(:), ...
    'VariableNames', {'Frame','GreenPos','GreenNeg','PinkPos','PinkNeg','Unknown'});

totalGreenPos = sum(countGreenPos);
totalGreenNeg = sum(countGreenNeg);
totalPinkPos  = sum(countPinkPos);
totalPinkNeg  = sum(countPinkNeg);
totalUnknown  = sum(countUnknown);

fprintf('Green(+):%d  Green(-):%d  Pink(+):%d  Pink(-):%d  Unknown:%d\n', ...
    totalGreenPos,totalGreenNeg,totalPinkPos,totalPinkNeg,totalUnknown);

Ttotals = table("Total", totalGreenPos, totalGreenNeg, totalPinkPos, totalPinkNeg, totalUnknown, ...
    'VariableNames', {'Frame','GreenPos','GreenNeg','PinkPos','PinkNeg','Unknown'});

Tall = [T; Ttotals];

%% -------------------- Save to Results/ --------------------
outFile = fullfile(resultPath, 'defect_type_counts.csv');
writetable(Tall, outFile);
fprintf('Saved: %s\n', outFile);
disp('Counting is finished');

%% ===================== Helper =====================
function colorClass = classifyPointColor(densityG, densityR, coord, varargin)
% Returns "Green", "Pink", or "Unknown" based on local G vs R around coord.
% Uses circular neighborhood with simple aggregation.

p = inputParser;
p.addParameter('Radius',    5,       @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MinDiff',   5e-4,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Method',    'mean',  @(s)ischar(s) || isstring(s));
p.addParameter('MinPixels', [],      @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.parse(varargin{:});

R       = p.Results.Radius;
minDiff = p.Results.MinDiff;

if ~isequal(size(densityG), size(densityR))
    error('densityG and densityR size mismatch');
end

[nRows, nCols] = size(densityG);

% IMPORTANT: in your original code you used row=coord(2), col=coord(1)
row = coord(2);
col = coord(1);

% Default MinPixels ~ quarter of the circle area, but at least 10
if isempty(p.Results.MinPixels)
    MinPixels = max(10, ceil(pi*R^2/4));
else
    MinPixels = p.Results.MinPixels;
end

% Window + distances
rmin = max(1, row - R); rmax = min(nRows, row + R);
cmin = max(1, col - R); cmax = min(nCols, col + R);
[cc, rr] = meshgrid(cmin:cmax, rmin:rmax);
dr = rr - row; dc = cc - col;
d2 = dr.^2 + dc.^2;
mask = d2 <= R^2;

Gpatch = densityG(rmin:rmax, cmin:cmax);
Rpatch = densityR(rmin:rmax, cmin:cmax);
valid  = mask & ~isnan(Gpatch) & ~isnan(Rpatch);

if nnz(valid) < MinPixels
    colorClass = "Unknown"; return;
end

gVals = Gpatch(valid);
rVals = Rpatch(valid);

% Aggregation (binarized images => mean ~ local fraction of 1s)
Gagg = mean(gVals, 'omitnan');
Ragg = mean(rVals, 'omitnan');

% Robust normalized difference
den = (Gagg + Ragg);
if den == 0
    diffGR = 0;
else
    diffGR = (Gagg - Ragg) / den;
end

if abs(diffGR) < minDiff
    colorClass = "Unknown";
elseif diffGR > 0
    colorClass = "Green";
else
    colorClass = "Pink";
end
end