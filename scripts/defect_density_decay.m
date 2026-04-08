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


close all; clear;

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

%% Draw plot curves of decay constant and platue value
% Load the Excel file
data = readtable(fullfile(dataDir,'Data for scatter plot.xlsx'));

% Extract variables
concentration = data.Concentration;

AreaPixel = [800*800, 800*800];
Area      = AreaPixel * 2.344 * 2.344 * 1e-6;   % um^2


% For B
B          = data.B;
B_err_low  = B - data.B_low;
B_err_high = data.B_high - B;

% For C
C          = data.C;
C_err_low  = C - data.C_low;
C_err_high = data.C_high - C;

% --------- Sort everything by concentration (ascending) ---------
[concentration, sortIdx] = sort(concentration);  % sort concentrations
B          = B(sortIdx);
B_err_low  = B_err_low(sortIdx);
B_err_high = B_err_high(sortIdx);

C          = C(sortIdx);
C_err_low  = C_err_low(sortIdx);
C_err_high = C_err_high(sortIdx);

% ---------------------------------------------------------------

% Load your custom colormap
addpath('/Users/zhaofei/Documents/MATLAB');
cols = zfcolors();   % at least 8 rows

n    = numel(concentration);
gSz  = 3;                          % group size (3 per color)
nGrp = ceil(n / gSz);

% --- Main figure ---
fig = figure('Color','w', 'Visible','off');
ax  = axes('Parent',fig); hold(ax,'on'); box(ax,'on'); grid(ax,'off');


% === Main plot: B with error bars ===
for i = 1:n
    gi  = ceil(i / gSz);           % group index
    clr = cols(gi, :);             % pick color from zfcolors

    errorbar(concentration(i), B(i), B_err_low(i), B_err_high(i), ...
        'o', ...
        'MarkerFaceColor', clr, ...
        'MarkerEdgeColor', clr, ...
        'Color', clr, ...
        'LineStyle', 'none', ...
        'CapSize', 12, ...
        'LineWidth', 2.5, ...
        'MarkerSize', 12);
end

xlabel('\Phi_{MF}', ...
    'FontName', 'Times New Roman', 'FontSize', 18);
ylabel('$\tau_D$ (h)', ...
    'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 18);
ylim([3 25])
xlim([0 70])
pbaspect([4 3 1])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 40);
set(ax, 'FontName','Times New Roman', 'FontSize', 40, ...
        'LineWidth', 5, 'TickDir','out', 'Box','on');

% === Inset plot: C with error bars ===
inset_pos  = [0.26, 0.65, 0.25, 0.25];  % [x, y, width, height]
inset_axes = axes('Position', inset_pos); hold(inset_axes, 'on');

for i = 1:n
    gi  = ceil(i / gSz);
    clr = cols(gi, :);

    errorbar(concentration(i), C(i), C_err_low(i), C_err_high(i), ...
        'o', ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k', ...
        'Color', 'k', ...
        'LineStyle', 'none', ...
        'CapSize', 10, ...
        'LineWidth', 1, ...
        'MarkerSize', 10);
end

xlabel('\Phi_{MF}', 'FontSize', 6, 'FontName', 'Times New Roman');
ylabel('$\rho_\infty$', 'Interpreter', 'latex', ...
    'FontSize', 8, 'FontName', 'Times New Roman');
xlim([0 70]);
ylim([1 6]);
pbaspect([5 4 1])

set(inset_axes, 'FontName', 'Times New Roman', 'FontSize', 25);
grid(inset_axes, 'off'); box(inset_axes, 'on');


% --- Plot figures seperately ---
%% --- Figure 1: B with error bars ---
fig1 = figure('Color','w');
ax1  = axes('Parent',fig1); hold(ax1,'on'); box(ax1,'on');

for i = 1:n
    gi  = ceil(i / gSz);
    clr = cols(gi, :);

    h = errorbar(concentration(i), B(i), B_err_low(i), B_err_high(i), ...
        'o', ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k', ...
        'Color', [0.6 0.6 0.6], ...   % gray error bars
        'LineStyle', 'none', ...
        'CapSize', 12, ...
        'LineWidth', 2.5, ...
        'MarkerSize', 12);

h.MarkerEdgeColor = 'k';  % restore marker edge to black
end

xlabel('\Phi_{MF}', 'FontName', 'Times New Roman', 'FontSize', 40);
ylabel('$\tau_D$ (h)', 'Interpreter','latex', ...
    'FontName','Times New Roman','FontSize',40);
ylim([3 20]); xlim([0 70]);
pbaspect([4 3 1])
set(ax1, 'FontName','Times New Roman','FontSize',35, ...
         'LineWidth',4,'TickDir','out','Box','on');

%% --- Figure 2: C with error bars ---
fig2 = figure('Color','w');
ax2  = axes('Parent',fig2); hold(ax2,'on'); box(ax2,'on');

for i = 1:n
    gi  = ceil(i / gSz);
    clr = cols(gi, :);

    h = errorbar(concentration(i), C(i), C_err_low(i), C_err_high(i), ...
        'o', ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k', ...
        'Color', [0.6 0.6 0.6], ...   % gray error bars
        'LineStyle', 'none', ...
        'CapSize', 12, ...
        'LineWidth', 2.5, ...
        'MarkerSize', 12);

h.MarkerEdgeColor = 'k';  % restore marker edge to black
end

xlabel('\Phi_{MF}', 'FontName','Times New Roman','FontSize',40);
ylabel('\rho_{\infty} (1/mm^2)', 'FontName','Times New Roman','FontSize',40);
xlim([0 70]); ylim([1 6]);
pbaspect([4 3 1])
set(ax2, 'FontName','Times New Roman','FontSize',35, ...
         'LineWidth',4,'TickDir','out','Box','on');