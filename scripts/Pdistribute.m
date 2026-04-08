%% Pdistribution_averaged_with_error.m
% Compute and plot the averaged orientation distribution P(theta)
% with shaded error band from replicate director files.

clear; close all;

% -------------------- Repo-safe paths (do not use pwd) --------------------
repoRoot = fileparts(fileparts(mfilename("fullpath"))); % scripts/ → repo root

% --------- User-adjustable dataset name (relative to Data/) --------------
runName = 'Pdistribution';

% Input directory containing director_data_*.mat
dataDir = fullfile(repoRoot, 'data/sample', runName);

% Output directory for final figures
resultDir = fullfile(repoRoot, 'Results');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

% List of file pairs to average
directorList = {
    {'director_data_10_01.mat', 'director_data_10_02.mat'}, ...
    {'director_data_30_01.mat', 'director_data_30_02.mat'}, ...
    {'director_data_50_01.mat', 'director_data_50_02.mat'}, ...
    {'director_data_70_01.mat', 'director_data_70_02.mat'}
};

legendNames = {'10% MF', '30% MF', '50% MF', '70% MF'};
markerList  = {'o', 's', '^', 'd'};

% Histogram bin setup
numBins    = 50;
edges      = linspace(-pi, pi, numBins + 1);
binCenters = edges(1:end-1) + diff(edges)/2;
binWidth   = edges(2) - edges(1);

% === Use zfcolors ===
addpath('/Users/zhaofei/Documents/MATLAB');
cols = zfcolors();
colorOrder = cols(1:numel(directorList), :);

% Create figure
fig = figure('Color','w');
ax  = axes('Parent', fig);
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'off');

hPlot = gobjects(1, numel(directorList));

for k = 1:numel(directorList)

    nRep = numel(directorList{k});
    pdfAll = zeros(nRep, numBins);   % store PDF from each replicate

    for j = 1:nRep
        filePath = fullfile(dataDir, directorList{k}{j});
        data = load(filePath);

        cosTheta = data.directors(:,:,1);
        sinTheta = data.directors(:,:,2);

        theta = atan2(sinTheta, cosTheta);
        thetaVec = theta(:);

        % Remove NaN if needed
        thetaVec = thetaVec(~isnan(thetaVec));

        % Mean subtraction using circular mean
        meanTheta = atan2(mean(sin(thetaVec)), mean(cos(thetaVec)));
        deltaTheta = mod(thetaVec - meanTheta + pi, 2*pi) - pi;

        % Histogram -> PDF
        counts = histcounts(deltaTheta, edges);
        totalSamples = numel(deltaTheta);
        pdfVals = counts ./ (totalSamples * binWidth);

        pdfAll(j, :) = pdfVals;
    end

    % Mean and error across replicates
    pdfAvg = mean(pdfAll, 1);
    pdfStd = std(pdfAll, 0, 1);
    pdfSem = pdfStd ./ sqrt(nRep);

    % Choose error type for shadow
    errVals = pdfSem;   % <-- use SEM
    % errVals = pdfStd; % <-- use STD instead if preferred

    % Shaded error band
    xFill = [binCenters, fliplr(binCenters)];
    yFill = [pdfAvg - errVals, fliplr(pdfAvg + errVals)];

    fill(xFill, yFill, colorOrder(k,:), ...
        'FaceAlpha', 0.20, ...
        'EdgeColor', 'none');

    % Mean curve
    hPlot(k) = plot(binCenters, pdfAvg, ...
        '-', ...
        'Marker', markerList{k}, ...
        'LineWidth', 3, ...
        'MarkerSize', 8, ...
        'Color', colorOrder(k,:), ...
        'MarkerFaceColor', colorOrder(k,:), ...
        'MarkerEdgeColor', colorOrder(k,:));
    
    
end

% Axis formatting
xlabel('\delta\theta (radians)', 'FontName', 'Times New Roman', 'FontSize', 40);
ylabel('P(\delta\theta)', 'FontName', 'Times New Roman', 'FontSize', 40);

xlim([-5*pi/8, 5*pi/8]);
ylim([0 1.7]);
xticks([-pi/2, -pi/4, 0, pi/4, pi/2]);
xticklabels({'-\pi/2','-\pi/4','0','\pi/4','\pi/2'});

set(ax, 'FontName', 'Times New Roman', ...
        'FontSize', 40, ...
        'LineWidth', 5, ...
        'TickDir', 'out', ...
        'Box', 'on');

legend(hPlot, legendNames, ...
       'Location', 'northeast', ...
       'FontName', 'Times New Roman', ...
       'FontSize', 32, ...
       'Box', 'off');

exportgraphics(fig, fullfile(resultDir, 'Pdistribution_with_error.tif'), 'Resolution', 300);