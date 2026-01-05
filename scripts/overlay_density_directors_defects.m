%% overlay_density_directors_defects.m
% Overlay density maps (G/R/combined) with director field and defect markers,
% and export publication-ready TIFFs.
%
% Expected input directory structure (relative to repo root):
%   data/sample/overlay_density_directors_defects/
%       defectData.mat
%       segmentation/density_maps_frame*.mat
%       defectimages/director_data_*.mat
%
% Outputs:
%   Results/overlay_green_directors_frame*.tif
%   Results/overlay_red_directors_frame*.tif
%   Results/overlay_combined_directors_frame*.tif
%
% Notes:
% - Repo-safe paths via mfilename("fullpath")
% - Export uses FIGURE handle (not axes) to avoid exportgraphics invalid-figure warning
%
% Copyright (c) 2025 Zhaofei Zheng
% Released under the MIT License. See LICENSE in the repository root.

close all;
clear;

% -------------------- Repo-safe paths --------------------
repoRoot = fileparts(fileparts(mfilename("fullpath"))); % scripts/ -> repo root

% -------------------- User inputs --------------------
frameNum = 80; % example frame

runName   = 'overlay_density_directors_defects';
runFolder = fullfile('data', 'sample', runName);
baseDir   = fullfile(repoRoot, runFolder);


% Output folder
resultPath = fullfile(repoRoot, 'Results');
if ~exist(resultPath, 'dir')
    mkdir(resultPath);
end

% -------------------- File paths --------------------
densityMatPath = fullfile(baseDir, sprintf('density_maps_frame%d.mat', frameNum));
directorsPath  = fullfile(baseDir, sprintf('director_data_%d.mat', frameNum));
defectMatPath  = fullfile(baseDir, 'defectData.mat');

% -------------------- Visual params --------------------
densityMin = 1e-3;    % cells/um^2
densityMax = 20e-3;   % cells/um^2

% Director overlay params
gridSpacing = 60;     % px
lineLength  = 30;     % px
lineWidth   = 2;      % px
lineColor   = 'k';

% Defect params
markerSize = 8;
edgeWidth  = 1;

% -------------------- Load density maps --------------------
assert(exist(densityMatPath,'file')==2, 'Missing file: %s', densityMatPath);
S = load(densityMatPath);

if ~isfield(S,'G_cells_per_um2') || ~isfield(S,'R_cells_per_um2')
    error('density_maps_frame%d.mat missing required variables G_cells_per_um2 / R_cells_per_um2.', frameNum);
end
G_cells_per_um2 = S.G_cells_per_um2;
R_cells_per_um2 = S.R_cells_per_um2;

[H, W] = size(G_cells_per_um2);

% -------------------- Load directors --------------------
assert(exist(directorsPath,'file')==2, 'Missing file: %s', directorsPath);
Ddir = load(directorsPath);
if ~isfield(Ddir,'directors')
    error('director_data_%d.mat does not contain variable "directors".', frameNum);
end
directors = Ddir.directors;

% -------------------- Load defect data --------------------
assert(exist(defectMatPath,'file')==2, 'Missing file: %s', defectMatPath);
Ddef = load(defectMatPath);

requiredDefFields = {'defectX','defectY','defectCharge'};
for k = 1:numel(requiredDefFields)
    if ~isfield(Ddef, requiredDefFields{k})
        error('defectData.mat missing required variable "%s".', requiredDefFields{k});
    end
end

x = Ddef.defectX(:);
y = Ddef.defectY(:);
q = Ddef.defectCharge(:);

% -------------------- Colormaps --------------------
n = 256;
cmG = [linspace(1,0,n)', linspace(1,1,n)', linspace(1,0,n)'];    % white->green
cmR = [linspace(1,1,n)', linspace(1,0.5,n)', linspace(1,1,n)'];  % white->pink

% -------------------- Output names --------------------
outNameG = fullfile(resultPath, sprintf('overlay_green_directors_frame%d.tif', frameNum));
outNameR = fullfile(resultPath, sprintf('overlay_red_directors_frame%d.tif', frameNum));
outNameC = fullfile(resultPath, sprintf('overlay_combined_directors_frame%d.tif', frameNum));

% -------------------- Combined RGB background --------------------
G_norm = mat2gray(G_cells_per_um2);
R_norm = mat2gray(R_cells_per_um2);

combinedRGB = zeros(H, W, 3);
combinedRGB(:,:,2) = 1 - R_norm;   % green
combinedRGB(:,:,1) = 1 - G_norm;   % red
combinedRGB(:,:,3) = 1 - G_norm;   % blue -> pink-ish

%% -------------------- Draw & export: Green --------------------
[f, ax] = newCleanFig();
imagesc(ax, G_cells_per_um2, [densityMin densityMax]);
colormap(ax, cmG);

drawDirectorfield(directors, 'Axes', ax, ...
    'gridSpacing', gridSpacing, 'lineLength', lineLength, ...
    'lineWidth', lineWidth, 'color', lineColor, 'scaleTo', [H, W]);

plotDefects(ax, x, y, q, markerSize, edgeWidth);

drawnow; % ensure render is complete before export
exportgraphics(f, outNameG, 'Resolution', 300);
close(f);

%% -------------------- Draw & export: Red --------------------
[f, ax] = newCleanFig();
imagesc(ax, R_cells_per_um2, [densityMin densityMax]);
colormap(ax, cmR);

drawDirectorfield(directors, 'Axes', ax, ...
    'gridSpacing', gridSpacing, 'lineLength', lineLength, ...
    'lineWidth', lineWidth, 'color', lineColor, 'scaleTo', [H, W]);

plotDefects(ax, x, y, q, markerSize, edgeWidth);

drawnow;
exportgraphics(f, outNameR, 'Resolution', 300);
close(f);

%% -------------------- Draw & export: Combined --------------------
[f, ax] = newCleanFig();
imshow(combinedRGB, 'Parent', ax);

drawDirectorfield(directors, 'Axes', ax, ...
    'gridSpacing', gridSpacing, 'lineLength', lineLength, ...
    'lineWidth', lineWidth, 'color', lineColor, 'scaleTo', [H, W]);

plotDefects(ax, x, y, q, 3*markerSize, edgeWidth);

drawnow;
exportgraphics(f, outNameC, 'Resolution', 300);
close(f);

%% -------------------- Helper functions --------------------
function [f, ax] = newCleanFig()
f = figure('Color','w', 'Visible','off');
f.InvertHardcopy = 'off';
ax = axes('Parent', f);
axis(ax,'image'); axis(ax,'ij'); axis(ax,'off');
hold(ax,'on');
end

function plotDefects(ax, x, y, q, markerSize, edgeWidth)
idxPlus  = q > 0;
idxMinus = q < 0;

plot(ax, x(idxPlus), y(idxPlus), 'o', ...
    'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor','k', ...
    'MarkerSize',markerSize, 'LineWidth',edgeWidth);

plot(ax, x(idxMinus), y(idxMinus), 'o', ...
    'MarkerFaceColor',[0 0 1], 'MarkerEdgeColor','k', ...
    'MarkerSize',markerSize, 'LineWidth',edgeWidth);
end

function h = drawDirectorfield(directors, varargin)
p = inputParser;
addParameter(p,'Axes',[],@(x) isempty(x) || ishghandle(x,'axes'));
addParameter(p,'gridSpacing',60);
addParameter(p,'lineLength',45);
addParameter(p,'lineWidth',2);
addParameter(p,'color','y');
addParameter(p,'scaleTo',[],@(x) isempty(x) || (isnumeric(x) && numel(x)==2));
parse(p,varargin{:});

ax = p.Results.Axes; if isempty(ax); ax = gca; end

theta = atan2(directors(:,:,2), directors(:,:,1));

if ~isempty(p.Results.scaleTo)
    tgt_rows = p.Results.scaleTo(1);
    tgt_cols = p.Results.scaleTo(2);
    theta_resamp = imresize(theta, [tgt_rows, tgt_cols], 'bilinear');
else
    [tgt_rows, tgt_cols] = size(theta);
    theta_resamp = theta;
end

[X, Y] = meshgrid(1:p.Results.gridSpacing:tgt_cols, ...
                  1:p.Results.gridSpacing:tgt_rows);
theta_samp = theta_resamp(1:p.Results.gridSpacing:end, ...
                          1:p.Results.gridSpacing:end);

h = gobjects(numel(X),1); k = 0;
for i = 1:size(X,1)
    for j = 1:size(X,2)
        ang = theta_samp(i,j);
        if isnan(ang), continue; end
        dx = p.Results.lineLength * cos(ang);
        dy = p.Results.lineLength * sin(ang);
        k = k + 1;
        h(k) = line(ax, [X(i,j)-dx/2, X(i,j)+dx/2], ...
                         [Y(i,j)-dy/2, Y(i,j)+dy/2], ...
                         'Color', p.Results.color, ...
                         'LineWidth', p.Results.lineWidth);
    end
end
h = h(1:k);
end