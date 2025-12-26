%% defectvelocity.m
% Compute and visualize defect velocities from time-resolved defect
% trajectories obtained from the director-field analysis pipeline.
%
% Expected input directory structure (relative to repo root):
%   Data/<run_name>/
%       defectData.mat
%
% Outputs:
%   Results/defect_velocity_<run_name>.pdf (or .tif)
%
% This script reproduces defect-velocity statistics shown in the
% main text and Supporting Information.
%
% Copyright (c) 2025 Zhaofei Zheng
% Released under the MIT License. See LICENSE in the repository root.

clear
clc

% -------------------- Repo-safe paths (do not use pwd) --------------------
repoRoot = fileparts(fileparts(mfilename("fullpath"))); % scripts/ → repo root

% --------- User-adjustable dataset name (relative to Data/) --------------
runName = 'defectvelocity';

% Input directory
dataDir = fullfile(repoRoot, 'data/sample', runName);

% Output directory for final figures
resultDir = fullfile(repoRoot, 'Results');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

% Extract filename with patterns
varName    = 'velocity';      % variable name inside each MAT file
pattern    = '*_*.mat';       % file name pattern


% X-axis numerical values for each of the 12 files

xLabels = [8, 9.1, 13.1, 26.3, 19.8, 24, 31.7, 34.7, 43.7, 61.2, 48.2, 55];


files = dir(fullfile(dataDir, pattern));
nFiles = numel(files);

if nFiles == 0
    error('No MAT files found in %s with pattern %s', dataDir, pattern);
end

fileNames = {files.name};

% Preallocate
fileMeans = zeros(nFiles, 1);
fileStds  = zeros(nFiles, 1);
groupVals = zeros(nFiles, 1);   % numeric: 10, 30, 50, 70, ...

for i = 1:nFiles
    fname = fileNames{i};
    fpath = fullfile(dataDir, fname);

    % Extract group from file name, e.g. something_10_01.mat -> 10
    tok = regexp(fname, '_(\d+)_\d+\.mat', 'tokens', 'once');
    if isempty(tok)
        error('Cannot parse group from file name: %s', fname);
    end
    groupVals(i) = str2double(tok{1});

    % Load data
    S = load(fpath);
    if ~isfield(S, varName)
        error('Variable "%s" not found in %s', varName, fname);
    end
    v = S.(varName)(:);   % column

    n = sum(~isnan(v));
    sigma = std(v, 0, 'omitnan');
    
    fileMeans(i) = mean(v, 'omitnan');
    fileStds(i)  = sigma / sqrt(n);   % <-- SEM as "error"
end

% ---------- Plot: 12 points with error bars ----------

% X-axis = group values (10, 30, 50, 70).
% This will give multiple points at the same x (replicates).
x = xLabels(:);   % ensure column
y = fileMeans;
e = fileStds;

% Optional: small jitter so points at same x are slightly separated horizontally
jitterAmp = 0.4;
xJitter = x + (rand(size(x)) - 0.5) * 2 * jitterAmp;

% Colors per group using your zfcolors
colors = zfcolors();          % make sure zfcolors.m is on the path
uGroups = unique(groupVals);  % unique conditions, e.g. [10 30 50 70]

% Map group value -> color index
colorMap = containers.Map('KeyType','double','ValueType','any');
for k = 1:numel(uGroups)
    colorMap(uGroups(k)) = colors(k,:);   % first 4 colors of zfcolors
end

fig = figure('Color','w');
ax  = axes('Parent',fig); hold(ax,'on'); box(ax,'on');

% Plot each point with its own color based on group
for i = 1:nFiles
    c = colorMap(groupVals(i));
    errorbar(ax, xJitter(i), y(i), e(i), ...
        'o', ...
        'MarkerFaceColor', c, ...
        'MarkerEdgeColor', c, ...
        'MarkerSize', 12, ...
        'Color', c, ...
        'LineWidth', 2, ...
        'CapSize', 12);
end

% ---- Axes styling ----
xlabel(ax, '\Phi_{MF}', 'FontName','Times New Roman','FontSize',40);
ylabel(ax, 'v_{defect}(\mum/h)', 'FontName','Times New Roman','FontSize',40);
xlim([0 65]);
ylim([35 55]);

set(ax, 'FontName','Times New Roman','FontSize',35, ...
        'LineWidth', 4, 'TickDir','out', 'Box','on');

% Use the actual group values as xticks
set(ax, 'XTick', uGroups, 'XTickLabel', arrayfun(@num2str, uGroups, 'UniformOutput', false));

% Optional legend: one entry per group
legLabels = arrayfun(@(g) sprintf('%d', g), uGroups, 'UniformOutput', false);
for k = 1:numel(uGroups)
    hDummy(k) = plot(ax, NaN, NaN, 'o', ...
        'MarkerFaceColor', colors(k,:), ...
        'MarkerEdgeColor', 'k', ...
        'Color', colors(k,:), ...
        'LineWidth', 2);
end
%legend(hDummy, legLabels, 'Location','best', 'FontSize',20, 'Title','Condition');

% Save figure
outPNG = fullfile(resultDir, 'defect_velocity_mean.tif');
exportgraphics(fig, outPNG, 'Resolution', 300);
fprintf('Saved figure to: %s\n', outPNG);