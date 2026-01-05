%===========================================
% Velocity magnitude + directors + defects
%===========================================
% Description:
%   This script loads PIV velocity (u_filtered, v_filtered), computes velocity
%   magnitude (converted to µm/h), resizes the velocity map to full-resolution
%   image coordinates, and overlays:
%     (1) director field (from directors(:,:,1:2) = [cosθ, sinθ])
%     (2) defect positions and charges (from defectData.mat)
%   The final figure is exported as a high-resolution TIFF (300 dpi).
%
% Inputs (relative to repo root):
%   Data/<datasetName>/PIVlab.mat
%   Data/<datasetName>/director_data_<frameNum>.mat
%   Data/<datasetName>/defectData.mat
%
% Outputs:
%   Results/<datasetName>/velocity_magnitude_<datasetName>.mat
%   Results/<datasetName>/velocity_overlay_frame<frameNum>.tif
%
% License:
%   MIT License (see LICENSE file in the repository root).
%   Copyright (c) 2025 Zhaofei Zheng
%
%===========================================

close all; clear;

% -------------------- Repo-safe paths (do not use pwd) --------------------
repoRoot = fileparts(fileparts(mfilename("fullpath"))); % scripts/ → repo root
dataRoot = fullfile(repoRoot, 'data');
resRoot  = fullfile(repoRoot, 'Results');

runName = 'overlay_velocity_director';

% ---------- Inputs ----------
frameNum = 60;
datasetName = fullfile('sample', runName);
fileTag = regexprep(datasetName, '[\\/]', '_');   % for filenames only

baseDir = fullfile(dataRoot, datasetName);
pivFile       = fullfile(baseDir, 'PIVlab.mat');                     % has u_filtered, v_filtered
directorsPath = fullfile(baseDir, sprintf('director_data_%d.mat', frameNum));
defectMatPath = fullfile(baseDir, 'defectData.mat');                 % has defectX, defectY, defectCharge

% Optional: add dependency folder if it exists inside the repo
cthetaDir = fullfile(repoRoot, 'Ctheta');
if exist(cthetaDir, 'dir')
    addpath(cthetaDir);
end

% ---------- Outputs ----------
outDir = fullfile(resRoot, datasetName);
if ~exist(outDir, 'dir'), mkdir(outDir); end

saveVelMat = fullfile(outDir, sprintf('velocity_magnitude_%s.mat', fileTag));
savePath   = fullfile(outDir, sprintf('velocity_overlay_%s_frame%d.tif', fileTag, frameNum));

% value is calculated later based on the image size
% lineLength  = gridSpacing;      % px line length, same for this one
lineWidth   = 2;       % line thickness
lineColor   = 'w';     % director line color (e.g., 'y' for yellow)
% Defect params
markerSize = 20;       % circle size
edgeWidth  = 0.1;      % circle edge thickness (unused here, kept as-is)

% ---------- Conversion factors (from your code 1) ----------
px_per_mm   = 2.344;     % example scaling in your pipeline (not used directly here)
h_per_frame = 2.5;       % [h/frame]
convFactor  = (px_per_mm) * (1/h_per_frame);  % μm/h per (px/frame)

% -------------------------------
% Load PIV velocity and compute magnitude
% -------------------------------
P = load(pivFile);
u = P.u_filtered;   % size ~ [Ny_piv, Nx_piv]
v = P.v_filtered;

% Here u{1}means extract frame 60
velocity_magnitude = sqrt(u{1}.^2 + v{1}.^2) * convFactor;
meanVelocity = mean(velocity_magnitude, 'all');
stdVelocity  = std(velocity_magnitude, 0, 'all');

fprintf('The mean velocity is %.2f with std %.2f\n', meanVelocity, stdVelocity);
save(saveVelMat, 'velocity_magnitude');

% -------------------------------
% Load director field and defects (full-res pixel coords)
% -------------------------------
Ddir = load(directorsPath);
if ~isfield(Ddir,'directors')
    error('director_data_%d.mat must contain variable "directors".', frameNum);
end
directors = Ddir.directors;        % size ~ [H, W, 2] where (:,:,1)=cosθ, (:,:,2)=sinθ
[H, W, ~] = size(directors);

% Now get the grid size
gridSpacing = max(1, round(H / size(u{1},1)));
lineLength  = gridSpacing;

Ddef = load(defectMatPath);
f = Ddef.defectFrames(:);
index = (f == frameNum);
x = Ddef.defectX(:);               % pixel coords in [1..W]
x = x(index);
y = Ddef.defectY(:);               % pixel coords in [1..H]
y = y(index);
q = Ddef.defectCharge(:);
q = q(index);

% -------------------------------
% Resize velocity map to match full-res image size
% -------------------------------
vel_full = imresize(velocity_magnitude, [H, W], 'nearest');

% -------------------------------
% Plot: velocity map + directors + defects
% -------------------------------
fig = figure('Color','w', ...
             'Units','pixels', ...
             'Position',[100 100 1200 1200]);   % [left bottom width height]
ax  = axes('Parent', fig);
ax.Toolbar.Visible = 'off'; % suppress "axes toolbar" export warning

imagesc(ax, vel_full);
axis(ax, 'image');
axis(ax, 'ij');                % top-left origin to match your coords
colormap(ax, parula);

% Colorbar (horizontal, keep units)
c = colorbar(ax, 'southoutside');
c.Label.String      = 'Velocity (\mum/h)';
c.Label.Interpreter = 'tex';
c.Label.FontName    = 'Times New Roman';
c.Label.FontSize    = 50;
% Tick label font
c.FontSize          = 50;   % colorbar tick numbers
c.FontName          = 'Times New Roman';
clim([0, 8]);                 % adjust if needed

% Remove all x/y labels and numbers
set(ax, 'XTick', [], 'YTick', [], 'XColor','none','YColor','none', ...
        'Box','on', 'LineWidth', 2);

hold(ax, 'on');

% Overlay director field (resampled to full-res target = [H, W])
drawDirectorfield(directors, ...
    'Axes', ax, ...
    'gridSpacing', gridSpacing, ...
    'lineLength',  lineLength, ...
    'lineWidth',   lineWidth, ...
    'color',       lineColor, ...
    'scaleTo',     [H, W]);

% Overlay defects
idxPlus  = q > 0;   % +1/2 → red
idxMinus = q < 0;   % -1/2 → blue
plot(ax, x(idxPlus),  y(idxPlus),  'o', 'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor','k', ...
    'MarkerSize',markerSize);
plot(ax, x(idxMinus), y(idxMinus), 'o', 'MarkerFaceColor',[0 0 1], 'MarkerEdgeColor','k', ...
    'MarkerSize',markerSize);

% Make sure figure is saved at on-screen size, high resolution
set(fig, 'PaperPositionMode', 'auto');

% Save as 300 dpi TIFF
print(fig, savePath, '-dtiff', '-r300');


% --- helper ---
function h = drawDirectorfield(directors, varargin)
% directors: Ny x Nx x 2 ((:,:,1)=cosθ, (:,:,2)=sinθ)

p = inputParser;
addParameter(p,'Axes',[],@(x) isempty(x) || ishghandle(x,'axes'));
addParameter(p,'gridSpacing',60);
addParameter(p,'lineLength',45);
addParameter(p,'lineWidth',2);
addParameter(p,'color','w');
addParameter(p,'scaleTo',[],@(x) isempty(x) || (isnumeric(x) && numel(x)==2));
parse(p,varargin{:});
ax = p.Results.Axes; if isempty(ax); ax = gca; end

theta = atan2(directors(:,:,2), directors(:,:,1));   % radians

% Resample theta to background size (if provided)
if ~isempty(p.Results.scaleTo)
    tgt_rows = p.Results.scaleTo(1);
    tgt_cols = p.Results.scaleTo(2);
    theta_resamp = imresize(theta, [tgt_rows, tgt_cols], 'bilinear');
else
    [tgt_rows, tgt_cols] = size(theta);
    theta_resamp = theta;
end

% Sample grid (gridSpacing must be integer for indexing)
step = max(1, round(p.Results.gridSpacing));

[X, Y] = meshgrid(1:step:tgt_cols, 1:step:tgt_rows);
theta_samp = theta_resamp(1:step:end, 1:step:end);

X = X + 0.5 * step;
Y = Y + 0.5 * step;
hold(ax,'on');
axis(ax,'ij'); axis(ax,'image');
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