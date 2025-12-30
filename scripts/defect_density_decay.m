%% defect_density_decay.m
% This repository contains MATLAB scripts for analyzing and fitting the temporal decay
% of +1/2 topological defect density in dense cell monolayers.
% 
% The analysis reads preprocessed CSV files containing defect counts and time indices,
% normalizes defect density by the physical imaging area, and fits the decay using an
% exponential model of the form:
% 
%     ρ(t) = A exp(-t / B) + C
%
% Outputs:
%   Results/defectdensitydecay.tif
%
% The script generates:
% - A publication-quality plot of defect density versus time with fitted curves
% - Text files summarizing fitted parameters, confidence intervals, and goodness-of-fit
%
% Copyright (c) 2025 Zhaofei Zheng
% Released under the MIT License. See LICENSE in the repository root.


clear; clc;

% -------------------- Repo-safe paths (do not use pwd) --------------------
repoRoot = fileparts(fileparts(mfilename("fullpath"))); % scripts/ → repo root
resultPath = fullfile(repoRoot, 'Results');
if ~exist(resultPath, 'dir'); mkdir(resultPath); end

% --------- Input folders (EDIT HERE) ----------
runName = 'defect_density_decay';
dataDir = fullfile(repoRoot, 'data/sample', runName);
folderPaths = { ...
    fullfile(dataDir, 'Processed_Results_0215_750celldensity_10_03'), ...
    fullfile(dataDir, 'Processed_Results_0215_750celldensity_50_03') ...
};

% Optional labels for legend
datasetLabels = { ...
    '10% MF', ...
    '50% MF' ...
};
markerList = {'o', 's', '^', 'd'};  

% Colors for each dataset (N x 3)
cols = zfcolors(); 



% --------- Constants ----------
AreaPixel = 800*800;
Area      = AreaPixel * 2.344 * 2.344 * 1e-6;   % um^2
frameRate = 0.5;           % h/frame

model       = fittype('A*exp(-x/B) + C', 'independent','x','dependent','y');
startPoints = [2, 7, 1.5];

% --------- Figure for density + fit ----------
fig = figure('Color','w');
ax  = axes('Parent',fig); hold(ax,'on'); box(ax,'on');

hFit  = gobjects(numel(folderPaths),1);
hData = gobjects(numel(folderPaths),1);

for k = 1:numel(folderPaths)

    folderPath = folderPaths{k};
    inputFile  = fullfile(folderPath, 'processed_results.csv');
    T = readtable(inputFile);

    % ----- Data -----
    x = T{:,1} * frameRate;     % time (h)
    y = T{:,2} / Area;          % +1/2 defect density (1/mm^2)
    orderParameter = T{:,5};

    % ----- Fit -----
    [fitResult, gof] = fit(x, y, model, 'StartPoint', startPoints);
    ci = confint(fitResult, 0.95);    % 95% CI

    xFit = linspace(min(x), max(x), 400);
    yFit = feval(fitResult, xFit);

    % ----- Plot on shared axes -----
    % shadow points
    hData(k) = scatter(ax, x, y, 160, ...
        'Marker', markerList{2*k-1}, ...
        'MarkerFaceColor', cols(2*k-1,:), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.25);

    % fit line
    hFit(k) = plot(ax, xFit, yFit, '-', ...
        'Color', cols(2*k-1,:), ...
        'LineWidth', 6);

    % ----- Save fit stats for this dataset -----
    outTXT = fullfile(resultPath, sprintf('fitting_results_%d.txt', k));
    fid = fopen(outTXT,'w');
    fprintf(fid, 'Exponential Fit: y = A*exp(-x/B) + C\n\n');

    pnames  = coeffnames(fitResult);
    pvalues = coeffvalues(fitResult);
    fprintf(fid, 'Fitting Parameters (95%% CI):\n');
    for i = 1:numel(pnames)
        fprintf(fid, '%s = %.6g (%.6g – %.6g)\n', ...
            pnames{i}, pvalues(i), ci(1,i), ci(2,i));
    end
    fprintf(fid, '\nGoodness of Fit:\n');
    gofFields = fieldnames(gof);
    for i = 1:numel(gofFields)
        fprintf(fid, '%s = %.6g\n', gofFields{i}, gof.(gofFields{i}));
    end
    fclose(fid);

end

% ----- Styling for shared density+fit figure -----
xlabel(ax, 't (h)', 'FontName','Times New Roman', 'FontSize',40);
ylabel(ax, '\rho (1/mm^2)', 'FontName','Times New Roman', 'FontSize',40);
set(ax, 'FontName','Times New Roman', 'FontSize',40, ...
        'LineWidth', 5, ...
        'TickDir','out', 'Box','on');
set(ax, 'LooseInset', max(get(ax,'TightInset'), 0.02));

% Legend (using fit lines or data points)
legend(ax, hFit, datasetLabels, ...
    'Location','northeast', ...
    'FontName','Times New Roman', ...
    'FontSize',30, ...
    'Box','off');

% Save the figure (PNG, 300 dpi)
resultFilename = fullfile(resultPath, 'defectdensitydecay.tif');
exportgraphics(gcf, resultFilename, 'Resolution', 300);