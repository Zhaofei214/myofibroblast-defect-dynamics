%% velocity_hdf_vs_mf_500.m
% This repository contains MATLAB scripts for analyzing cell motility from
% linked nucleus tracking data and comparing the velocity distributions of
% HDF and myofibroblasts MF.
%
% The analysis reads a CSV file, which is captured from trackmate in imageJ, containing linked cell tracks and cell-type
% assignments, computes frame-to-frame cell velocities based on nucleus
% displacement, and compares HDF and MF motility using:
%
%   - Histogram distributions
%   - Boxplots
%   - Kernel density estimation curves
%   - Two-sample t-test statistics
%
% Required input:
%   data/sample/velocity_hdf_vs_mf/linked_tracks_with_celltype.csv
%
% Output:
%   Results/velocity_hdf_vs_mf/
%       velocity_histogram_hdf_vs_mf.tif
%       velocity_boxplot_hdf_vs_mf.tif
%       velocity_density_hdf_vs_mf.tif
%       velocity_HDF_vs_MF.csv
%       velocity_statistics.txt
%
% Copyright (c) 2025 Zhaofei Zheng
% Released under the MIT License. See LICENSE in the repository root.

close all; clear; clc;

% -------------------- Repo-safe paths (do not use pwd) --------------------
repoRoot   = fileparts(fileparts(mfilename("fullpath"))); % scripts/ → repo root
resultPath = fullfile(repoRoot, 'Results', 'velocity_hdf_vs_mf_500');
if ~exist(resultPath, 'dir'); mkdir(resultPath); end

% -------------------- Input paths --------------------
runName  = 'velocity_hdf_vs_mf_500';
dataDir  = fullfile(repoRoot, 'data', 'sample', runName);
inputCsv = fullfile(dataDir, 'linked_tracks_with_celltype.csv');

% -------------------- Output paths --------------------
outputCsv        = fullfile(resultPath, 'velocity_HDF_vs_MF.csv');
statsTxt         = fullfile(resultPath, 'velocity_statistics.txt');
histFigFile      = fullfile(resultPath, 'velocity_histogram_hdf_vs_mf.tif');
boxFigFile       = fullfile(resultPath, 'velocity_boxplot_hdf_vs_mf.tif');
densityFigFile   = fullfile(resultPath, 'velocity_density_hdf_vs_mf.tif');

% -------------------- User parameters --------------------
pixelSize_um = 2.344;   % micron/pixel
dt_hr        = 0.5;     % time between frames (hours)

% -------------------- Load data --------------------
T = readtable(inputCsv);

requiredCols = {'TRACK_ID','Linked_X_px','Linked_Y_px','Linked_Frame_1based','CellType'};
assert(all(ismember(requiredCols, T.Properties.VariableNames)), ...
    'Missing required columns in input CSV.');

% -------------------- Sort by track and time --------------------
T = sortrows(T, {'TRACK_ID','Linked_Frame_1based'});

% -------------------- Compute velocities --------------------
velocities = [];
cellTypes  = [];

uniqueTracks = unique(T.TRACK_ID);

for i = 1:numel(uniqueTracks)
    trackID = uniqueTracks(i);
    idx = T.TRACK_ID == trackID;
    track = T(idx, :);

    if height(track) < 2
        continue;
    end

    x = track.Linked_X_px;
    y = track.Linked_Y_px;

    dx = diff(x);
    dy = diff(y);

    disp_um = sqrt(dx.^2 + dy.^2) * pixelSize_um;
    v = disp_um / dt_hr;   % um/hr

    % Assign track cell type using first row
    type = string(track.CellType(1));

    velocities = [velocities; v];
    cellTypes  = [cellTypes; repmat(type, numel(v), 1)];
end

% -------------------- Keep only valid cell types --------------------
validIdx = (cellTypes == "HDF") | (cellTypes == "MF");
velocities = velocities(validIdx);
cellTypes  = cellTypes(validIdx);

v_HDF = velocities(cellTypes == "HDF");
v_MF  = velocities(cellTypes == "MF");

% -------------------- Save per-step velocity table --------------------
outTable = table(velocities, cellTypes, 'VariableNames', {'Velocity_um_per_hr','CellType'});
writetable(outTable, outputCsv);

% -------------------- Statistics --------------------
meanHDF = mean(v_HDF, 'omitnan');
meanMF  = mean(v_MF,  'omitnan');

[~, p, ci, stats] = ttest2(v_HDF, v_MF);

if     p < 0.001
    stars = '***';
elseif p < 0.01
    stars = '**';
elseif p < 0.05
    stars = '*';
else
    stars = 'n.s.';
end

fid = fopen(statsTxt, 'w');
fprintf(fid, 'Velocity comparison: HDF vs MF\n\n');
fprintf(fid, 'Number of HDF velocity samples: %d\n', numel(v_HDF));
fprintf(fid, 'Number of MF velocity samples : %d\n\n', numel(v_MF));
fprintf(fid, 'Mean velocity of HDF: %.4f um/hr\n', meanHDF);
fprintf(fid, 'Mean velocity of MF : %.4f um/hr\n\n', meanMF);
fprintf(fid, 'Two-sample t-test results:\n');
fprintf(fid, 't-statistic = %.6f\n', stats.tstat);
fprintf(fid, 'degrees of freedom = %.6f\n', stats.df);
fprintf(fid, 'p-value = %.6g\n', p);
fprintf(fid, '95%% CI = [%.6f, %.6f]\n', ci(1), ci(2));
fprintf(fid, 'Significance = %s\n', stars);
fclose(fid);

fprintf('\nMean velocity:\n');
fprintf('HDF: %.2f um/hr\n', meanHDF);
fprintf('MF : %.2f um/hr\n', meanMF);
fprintf('p-value = %.4g\n', p);
fprintf('Significance: %s\n', stars);

% -------------------- Plot 1: Histogram comparison --------------------
fig1 = figure('Color', 'w', 'Visible', 'off');
ax1  = axes('Parent', fig1); hold(ax1, 'on'); box(ax1, 'on');

histogram(ax1, v_HDF, ...
    'Normalization', 'probability', ...
    'BinWidth', 2, ...
    'FaceColor', [1 0 1], ...
    'EdgeColor', 'k');

histogram(ax1, v_MF, ...
    'Normalization', 'probability', ...
    'BinWidth', 2, ...
    'FaceColor', [0 0.6 0], ...
    'EdgeColor', 'k');

xlabel(ax1, 'Velocity (\mum/hr)', 'FontSize', 40, 'FontName', 'Times New Roman');
ylabel(ax1, 'Probability',        'FontSize', 40, 'FontName', 'Times New Roman');

legend(ax1, {'HDF','MF'}, ...
    'FontSize', 28, ...
    'FontName', 'Times New Roman', ...
    'Location', 'best', ...
    'Box', 'off');

set(ax1, 'FontName', 'Times New Roman', ...
         'FontSize', 35, ...
         'LineWidth', 4, ...
         'TickDir', 'out', ...
         'Box', 'on');

exportgraphics(fig1, histFigFile, 'Resolution', 300);

% -------------------- Plot 2: Boxplot comparison --------------------
fig2 = figure('Color', 'w', 'Visible', 'off');
ax2  = axes('Parent', fig2); hold(ax2, 'on'); box(ax2, 'on');

group = [repmat("HDF", numel(v_HDF), 1); repmat("MF", numel(v_MF), 1)];
vals  = [v_HDF; v_MF];

boxplot(ax2, vals, group, 'Symbol', '');

ylabel(ax2, 'Velocity (\mum/hr)', 'FontSize', 40, 'FontName', 'Times New Roman');

set(ax2, 'FontName', 'Times New Roman', ...
         'FontSize', 35, ...
         'LineWidth', 4, ...
         'TickDir', 'out', ...
         'Box', 'on');

h = findobj(ax2, 'Tag', 'Box');

colors = [ ...
    0 0.6 0; ...
    1 0 1];

for j = 1:numel(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors(j,:), ...
        'FaceAlpha', 0.5, ...
        'EdgeColor', colors(j,:));
end

exportgraphics(fig2, boxFigFile, 'Resolution', 300);

% -------------------- Plot 3: Density comparison --------------------
fig3 = figure('Color', 'w', 'Visible', 'off');
ax3  = axes('Parent', fig3); hold(ax3, 'on'); box(ax3, 'on');

[f1, xi1] = ksdensity(v_HDF);
[f2, xi2] = ksdensity(v_MF);

plot(ax3, xi1, f1, 'Color', [1 0 1],   'LineWidth', 4);
plot(ax3, xi2, f2, 'Color', [0 0.6 0], 'LineWidth', 4);

xlabel(ax3, 'Velocity (\mum/hr)', 'FontSize', 40, 'FontName', 'Times New Roman');
ylabel(ax3, 'Probability',        'FontSize', 40, 'FontName', 'Times New Roman');

legend(ax3, {'HDF','MF'}, ...
    'FontSize', 28, ...
    'FontName', 'Times New Roman', ...
    'Location', 'best', ...
    'Box', 'off');

xlim(ax3, [0 35]);

set(ax3, 'FontName', 'Times New Roman', ...
         'FontSize', 35, ...
         'LineWidth', 4, ...
         'TickDir', 'out', ...
         'Box', 'on');

exportgraphics(fig3, densityFigFile, 'Resolution', 300);