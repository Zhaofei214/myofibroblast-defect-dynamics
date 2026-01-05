%% velocity_spatial_correlation_two_conditions.m
% Compute and compare isotropic spatial velocity correlation C_v(r) from two datasets
% (e.g., 750 vs 500 cells/mm^2) using precomputed PIVlab outputs.
%
% This script:
% - Loads PIV velocity fields from two folders (each containing PIVlab.mat)
% - Computes mean speed (um/h) for each condition
% - Computes isotropic spatial velocity correlation C_v(r) via direct lag-binning
% - Generates a publication-quality plot and saves it to Results/
%
% Expected repo structure (relative to repo root):
%   data/sample/velocity_spatial_correlation/
%       750_50/PIVlab.mat
%       500_50/PIVlab.mat
%
% Outputs:
%   Results/velocity_correlation.png
%
% Notes:
% - This script uses repo-safe paths and does NOT rely on pwd or absolute paths.
% - If zfcolors is not available on the MATLAB path, it falls back to MATLAB's "lines".
%
% Copyright (c) 2025 Zhaofei Zheng
% Released under the MIT License. See LICENSE in the repository root.

clear; close all

%% -------------------- Repo-safe paths (do not use pwd) --------------------
repoRoot   = fileparts(fileparts(mfilename("fullpath"))); % scripts/ -> repo root
resultPath = fullfile(repoRoot, 'Results');
if ~exist(resultPath, 'dir'); mkdir(resultPath); end

%% -------------------- Input folders (EDIT HERE) --------------------
runName = 'velocity_spatial_correlation';
dataDir = fullfile(repoRoot, 'data', 'sample', runName);

% Two conditions (each must contain PIVlab.mat)
baseDir1 = fullfile(dataDir, '750_50');
baseDir2 = fullfile(dataDir, '500_50');

%% -------------------- Parameters --------------------
um_per_px = 2.344;          % micron conversion
% velocity vectors may be on a coarse grid (e.g., pixel size effectively scaled)
time_per_h = 0.5 * 5;       % 5 frames
dr   = 1;                   % radial bin size (px)
rmax = 400;                 % max radius (px)

%% -------------------- Load colormap (repo-friendly) --------------------
% Option A (recommended): put zfcolors.m into repo (e.g., scripts/utils/) and add that folder:
% addpath(fullfile(repoRoot, 'scripts', 'utils'));
% Then zfcolors will be found automatically.
%
% Here we try to call zfcolors if it exists; otherwise fall back.
if exist('zfcolors', 'file') == 2
    cols = zfcolors;
else
    cols = lines(7); % fallback
end

%% ===== Load first dataset =====
PIVPath1 = fullfile(baseDir1, 'PIVlab.mat');
assert(exist(PIVPath1,'file')==2, 'Missing file: %s', PIVPath1);

S1 = load(PIVPath1, 'u_filtered', 'v_filtered');
Vx1 = S1.u_filtered{1};
Vy1 = S1.v_filtered{1};

% ===== Mean velocity magnitude =====
speed1 = sqrt(Vx1.^2 + Vy1.^2) * um_per_px / time_per_h;  % um/h
meanV1 = mean(speed1(:), 'omitnan');

[r1, Cr1, ~] = velocitySpatialCorrelation_direct(Vx1, Vy1, dr, rmax);

%% ===== Load second dataset =====
PIVPath2 = fullfile(baseDir2, 'PIVlab.mat');
assert(exist(PIVPath2,'file')==2, 'Missing file: %s', PIVPath2);

S2 = load(PIVPath2, 'u_filtered', 'v_filtered');
Vx2 = S2.u_filtered{1};
Vy2 = S2.v_filtered{1};

% ===== Mean velocity magnitude =====
speed2 = sqrt(Vx2.^2 + Vy2.^2) * um_per_px / time_per_h;  % um/h
meanV2 = mean(speed2(:), 'omitnan');

[r2, Cr2, ~] = velocitySpatialCorrelation_direct(Vx2, Vy2, dr, rmax);

%% ===== Plot results =====
figure('Color','w','Units','inches','Position',[1 1 8 8]);
hold on;

Nx = size(Vx1,1);
px_per_vec = 800 / Nx;

plot(r1 * um_per_px * px_per_vec, Cr1, ...
    'LineWidth', 3, ...
    'Color', cols(3,:), ...
    'LineStyle', '-', ...
    'Marker', '^', ...
    'MarkerSize', 14, ...
    'MarkerFaceColor', cols(3,:));

plot(r2 * um_per_px * px_per_vec, Cr2, ...
    'LineWidth', 3, ...
    'Color', cols(3,:), ...
    'LineStyle', '--', ...
    'Marker', '*', ...
    'MarkerSize', 14, ...
    'MarkerFaceColor', cols(3,:));

yline(0, '--k', 'LineWidth', 2);
xlim([0 max(r1) * um_per_px]);

xlabel('r (\mum)', 'FontName','Times New Roman', 'FontSize',40);
ylabel('C_v(r)', 'FontName','Times New Roman', 'FontSize',40);
legend({'750 cells/mm^2', '500 cells/mm^2'}, 'Location','northeast', ...
       'FontName','Times New Roman', 'FontSize',30, 'Box','off');

xticks([100 500]);
xticklabels({'100','500'});

axis square;
set(gca, 'FontName','Times New Roman', 'FontSize',36, ...
         'LineWidth',2, 'TickDir','out', 'Box','on');

hold off;

outFig = fullfile(resultPath, 'velocity_correlation.png');
exportgraphics(gcf, outFig, 'Resolution', 300);
fprintf('Saved: %s\n', outFig);

fprintf('Mean |V_750| = %.1f um/h, Mean |V_500| = %.1f um/h\n', meanV1, meanV2);

%% ========================================================================
function [r, Cr, Npairs] = velocitySpatialCorrelation_direct(Vx, Vy, dr, rmax)
% velocitySpatialCorrelation_direct
% Compute isotropic spatial velocity correlation C(r) from definition,
% without FFT, using integer lags and radial binning.
%
% Inputs
%   Vx, Vy : velocity components (Ny x Nx)
%   dr     : radial bin width in pixels (default 1)
%   rmax   : max radius in pixels (default = floor(min(Nx,Ny)/2))
%
% Outputs
%   r      : bin centers (pixels)
%   Cr     : correlation C(r) normalized so C(0)=1
%   Npairs : number of dot-product samples contributing to each bin

    if nargin < 3 || isempty(dr),   dr = 1; end
    [Ny, Nx] = size(Vx);
    if nargin < 4 || isempty(rmax), rmax = floor(min(Nx,Ny)/2); end

    % ---- 2) Precompute denominator <|v|^2> for normalization ----
    v2 = Vx.^2 + Vy.^2;
    denom = mean(v2(:), 'omitnan');       % equals C(0) by definition
    if denom == 0
        r = (0:dr:rmax).';
        Cr = zeros(size(r));
        Npairs = zeros(size(r));
        return
    end

    % ---- 3) Set up radial bins and accumulators ----
    edges = 0:dr:(rmax+dr);               % right-open bins [r, r+dr)
    r = edges(1:end-1) + dr/2;            % bin centers
    S = zeros(size(r));                   % sum of <v(x)·v(x+lag)>
    N = zeros(size(r));                   % count of samples

    % ---- 4) Loop over integer lags (dx,dy) and accumulate radially ----
    rmax2 = rmax^2;
    for dy = -rmax:rmax
        for dx = -rmax:rmax
            if dx*dx + dy*dy > rmax2
                continue
            end

            % Overlapping submatrices for this lag
            if dy >= 0
                y1 = 1:(Ny-dy);   y2 = (1+dy):Ny;
            else
                y1 = (1-dy):Ny;   y2 = 1:(Ny+dy);
            end
            if dx >= 0
                x1 = 1:(Nx-dx);   x2 = (1+dx):Nx;
            else
                x1 = (1-dx):Nx;   x2 = 1:(Nx+dx);
            end

            % Dot products v(x)·v(x+lag) for all overlapping pixels
            dp = Vx(y1,x1).*Vx(y2,x2) + Vy(y1,x1).*Vy(y2,x2);
            dp_norm = dp / denom; % normalization inside the loop

            % Radial bin index for this lag magnitude
            rr = sqrt(double(dx*dx + dy*dy));
            b = find(rr >= edges(1:end-1) & rr < edges(2:end), 1, 'first');
            if isempty(b), continue; end

            % Accumulate sum and counts (use mean over the overlap)
            S(b) = S(b) + mean(dp_norm(:), 'omitnan');
            N(b) = N(b) + 1;
        end
    end

    % ---- 5) Average over all lags that fell into each radial bin ----
    Cr = S ./ max(N,1);
    Npairs = N;

    % Ensure r=0 bin is exactly 1 (numerically)
    if ~isempty(Cr) && ~isnan(Cr(1)) && ~isinf(Cr(1))
        Cr(1) = 1;
    end

    % Optionally trim empty bins at the tail
    lastNonZero = find(Npairs>0, 1, 'last');
    if ~isempty(lastNonZero)
        r      = r(1:lastNonZero);
        Cr     = Cr(1:lastNonZero);
        Npairs = Npairs(1:lastNonZero);
    end
end