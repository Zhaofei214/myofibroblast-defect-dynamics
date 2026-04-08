function run_all()
%RUN_ALL  Entry point to run the full defect-dynamics analysis pipeline.
%
% Description:
%   Initializes the repository environment and (optionally) executes key
%   analysis scripts in a standard order. This file is intended to be the
%   single command users run after cloning the repo.
%
% Expected repo structure:
%   repo/
%     scripts/        (this file + analysis scripts)
%     Data/
%       raw/          (raw .czi inputs)
%       ...           (processed/intermediate outputs)
%     Results/        (final figures/tables)
%
% Usage:
%   run_all
%
% License:
%   MIT License (see LICENSE file in the repository root).
%   Copyright (c) 2026 Zhaofei Zheng
close all; clear; clc;


% -------------------- Repo-safe paths (do not use pwd) --------------------
thisFile   = mfilename("fullpath");
scriptsDir = fileparts(thisFile);
repoRoot   = fileparts(scriptsDir);

% -------------------- Initialize environment --------------------
if exist('startup.m','file') == 2
    startup;  % your startup.m should add paths and validate dependencies
else
    % Minimal fallback: add scripts and common dependency folders
    addpath(scriptsDir);
    if exist(fullfile(repoRoot,'Ctheta'),'dir'); addpath(fullfile(repoRoot,'Ctheta')); end
    if exist(fullfile(repoRoot,'external'),'dir'); addpath(genpath(fullfile(repoRoot,'external'))); end
end

fprintf("Repository initialized.\n");
fprintf("Repo root: %s\n", repoRoot);

% -------------------- Ensure standard folders exist --------------------
dataRawDir = fullfile(repoRoot, 'Data', 'raw');
resultsDir = fullfile(repoRoot, 'Results');

if ~exist(dataRawDir, 'dir'); mkdir(dataRawDir); end
if ~exist(resultsDir, 'dir'); mkdir(resultsDir); end

fprintf("Place raw .czi files in: %s\n", dataRawDir);

% -------------------- Optional: run pipeline scripts --------------------
% Toggle these to true/false depending on what you want to run by default.
RUN_CZI_PROCESSING    = false;  % e.g., generates processed_results.csv, director_data_*.mat, etc.
RUN_DEFECT_FITTING    = false;  % defect density decay fitting/plot
RUN_COLOR_THETA_MAPS  = false;  % theta colormap export
RUN_VELOCITY_OVERLAY  = false;  % velocity + directors + defects overlay
RUN_P_DISTRIBUTION    = false;  % P(theta) / angular distribution analysis
RUN_DEFECT_VELOCITY   = false;  % defect velocity analysis
RUN_COUNT_DEFECTTYPE  = false;  % count defects by type/charge
RUN_CELL_RATIO        = false;  % red/green area ratio check
RUN_OVERLAY           = false;  % overlay of cell density, director field and defect
RUN_V_CORRELATION     = false;  % velocity correlation based on PIV output
RUN_DEFECT_DOMINANCE  = false;  % classify defect type based on cell dominance (G/R)
RUN_VELOCITY_HDF_MF   = false;  % compare the HDF and MD cell velocity at 500 cell density

% Build an ordered pipeline list (keeps numbering consistent even if toggles change)
steps = {};

if RUN_CZI_PROCESSING
    steps(end+1,:) = {"CZI processing / defect counting", "cziprocessing_defectcounting"}; %#ok<AGROW>
end
if RUN_DEFECT_FITTING
    steps(end+1,:) = {"Defect density decay fitting", "defect_density_decay"}; %#ok<AGROW>
end
if RUN_COLOR_THETA_MAPS
    steps(end+1,:) = {"Theta colormap export", "color_angle_plot"}; %#ok<AGROW>
end
if RUN_VELOCITY_OVERLAY
    steps(end+1,:) = {"Velocity/director/defect overlay", "overlay_velocity_director"}; %#ok<AGROW>
end
if RUN_P_DISTRIBUTION
    steps(end+1,:) = {"P(theta) distribution analysis", "Pdistribution"}; %#ok<AGROW>
end
if RUN_DEFECT_VELOCITY
    steps(end+1,:) = {"Defect velocity analysis", "defectvelocity"}; %#ok<AGROW>
end
if RUN_COUNT_DEFECTTYPE
    steps(end+1,:) = {"Count defect types", "countdefecttype"}; %#ok<AGROW>
end
if RUN_CELL_RATIO
    steps(end+1,:) = {"Cell ratio verification", "cell_ratio_verification"}; %#ok<AGROW>
end
if RUN_OVERLAY
    steps(end+1,:) = {"Overlay density + directors + defects", "overlay_density_directors_defects"}; %#ok<AGROW>
end
if RUN_V_CORRELATION
    steps(end+1,:) = {"Velocity spatial correlation (PIV)", "velocity_spatial_correlation"}; %#ok<AGROW>
end
if RUN_DEFECT_DOMINANCE
    steps(end+1,:) = {"Defect dominance classification (G/R)", "defect_color_dominance"}; %#ok<AGROW>
end
if RUN_VELOCITY_HDF_MF
    steps(end+1,:) = {"HDF MF velocity comparasion at 500 cell density", "velocity_hdf_vs_mf_500"}; %#ok<AGROW>
end

% -------------------- Run pipeline --------------------
if isempty(steps)
    fprintf("\nNo steps enabled. Set RUN_* toggles to true in run_all.m.\n");
    return;
end

nSteps = size(steps, 1);
for i = 1:nSteps
    fprintf("\n[%d/%d] %s...\n", i, nSteps, steps{i,1});
    run_script_if_exists(steps{i,2});
end

fprintf("\nDone. Outputs are saved under: %s\n", resultsDir);

end


% -------------------- helper --------------------
function run_script_if_exists(scriptName)
%RUN_SCRIPT_IF_EXISTS Run a script by name if it exists on the MATLAB path.
% Runs in base workspace to prevent scripts that call "clear" from breaking run_all.

scriptFile = [char(scriptName) '.m'];

if exist(scriptFile, 'file') ~= 2
    warning('Script not found on path: %s (skipping)', scriptFile);
    return;
end

fprintf("  -> Running %s\n", scriptFile);

try
    evalin('base', sprintf("run('%s')", scriptName));
catch ME
    warning("Failed running %s: %s", scriptFile, ME.message);
end
end