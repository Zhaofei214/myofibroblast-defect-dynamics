%===========================================================
% Theta colormap (NO director overlay) with physical units
%===========================================================
% Description:
%   This script loads precomputed director fields (stored as a variable named
%   "directors" in .mat files), computes the nematic orientation angle
%   theta = mod(atan2(ny, nx), pi), and exports theta colormap images using
%   physical coordinates.
%
% Inputs:
%   - Data/sample/angle_color_plot/*.mat files containing variable: directors (Ny x Nx x 2)
%
% Outputs:
%   - Results/theta_maps/theta_map_XX.png (300 dpi)
%
% License:
%   MIT License (see LICENSE file in the repository root).
%   Copyright (c) 2025 Zhaofei Zheng
%
% Notes:
%   - Paths are resolved relative to the repository root (no pwd dependency).
%   - Optional cropping code is included but commented out.
%   - Further image processing was performed in ImageJ/Fiji to generate cell masks using thresholding.
%===========================================================

clear; clc; close all;

% -------------------- Repo-safe paths (do not use pwd) --------------------
repoRoot = fileparts(fileparts(mfilename("fullpath"))); % scripts/ → repo root

% Input dataset folder (relative to repo root)
folderName = fullfile('data/sample', 'angle_color_plot');
dataDir    = fullfile(repoRoot, folderName);
                 % input .mat files live here
outDir   = fullfile(repoRoot, 'Results', 'theta_maps'); % outputs go here
if ~exist(outDir, 'dir'), mkdir(outDir); end

% -------- Inputs --------
fileList = { ...
    'director_data_10.mat', ...
    'director_data_30.mat', ...
    'director_data_50.mat', ...
    'director_data_70.mat' ...
};
labels = {'10% MF','30% MF','50% MF','70% MF'};

% -------- Calibration --------
pixel_size_mm = 2.344;  % size per pixel in millimeters
% Desired crop size
cropSize = 2000;  % This is for crop the full image to 2000 by 2000 t0 remove the edges

for f = 1:numel(fileList)
    matPath = fullfile(dataDir, fileList{f});
    if ~isfile(matPath), warning('Missing file: %s', matPath); continue; end
    S = load(matPath);
    if ~isfield(S,'directors'), warning('No "directors" in %s', matPath); continue; end

    directors = S.directors;
    theta = atan2(directors(:,:,2), directors(:,:,1));
    theta = mod(theta, pi);

    [Ny, Nx] = size(theta);

    % % Compute start and end indices (centered)
    % yStart = floor(Ny/2 - cropSize/2) + 1;
    % yEnd   = yStart + cropSize - 1;
    %
    % xStart = floor(Nx/2 - cropSize/2) + 1;
    % xEnd   = xStart + cropSize - 1;
    %
    % % Crop out central region
    % theta = theta(yStart:yEnd, xStart:xEnd);
    % [Ny, Nx] = size(theta);

    % ---- Physical coordinates ----
    x_mm = (0:Nx-1) * pixel_size_mm;
    y_mm = (0:Ny-1) * pixel_size_mm;

    % ---- Plot ----
    fig = figure('Color','w','Position',[100 100 900 750]);
    ax = axes('Parent', fig); hold(ax,'on');
    set(ax,'LineWidth',3);   % makes all 4 sides of the box thicker

    imagesc(x_mm, y_mm, theta, 'Parent', ax);
    axis(ax,'equal','tight'); set(ax,'YDir','normal');
    colormap(ax, hsv); caxis(ax, [0 pi]);
    ax.Toolbar.Visible = 'off';

    % ---- Colorbar on the right side ----
    % c = colorbar(ax, 'eastoutside');  % moved from southoutside → eastoutside
    % c.Ticks = [0 pi/4 pi/2 3*pi/4 pi];
    % c.TickLabels = {'0','\pi/4','\pi/2','3\pi/4','\pi'};
    % c.FontName = 'Times New Roman'; c.FontSize = 40;
    % c.Label.String = '\theta (rad)';
    % c.Label.FontName = 'Times New Roman'; c.Label.FontSize = 40;

    % ---- Axes styling ----
    % ---- Remove all labels, ticks, and axis lines ----
    axis off;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLabel = [];
    ax.YLabel = [];
    ax.Title = [];

    title(ax, strrep(labels{f}, '%','%'), 'FontName','Times New Roman','FontSize',20);

    outPath = fullfile(outDir, sprintf('theta_map_%02d.png', f));
    exportgraphics(fig, outPath, 'Resolution', 300);
    close(fig);
end