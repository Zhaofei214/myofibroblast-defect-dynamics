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
%   Copyright (c) 2025 Zhaofei Zheng

% -------------------- Repo-safe paths (do not use pwd) --------------------
thisFile = mfilename("fullpath");
scriptsDir = fileparts(thisFile);
repoRoot = fileparts(scriptsDir);

% -------------------- Initialize environment --------------------
if exist('startup.m','file') == 2
    startup;  % your startup.m should add paths and validate dependencies
else
    % Minimal fallback: add scripts and common dependency folders
    addpath(scriptsDir);
    if exist(fullfile(repoRoot,'Ctheta'),'dir'); addpath(fullfile(repoRoot,'Ctheta')); end
end

fprintf("Repository initialized.\n");
fprintf("Repo root: %s\n", repoRoot);

% -------------------- Ensure standard folders exist --------------------
dataRawDir = fullfile(repoRoot, 'Data', 'raw');
resultsDir = fullfile(repoRoot, 'Results');

if ~exist(dataRawDir, 'dir'), mkdir(dataRawDir); end
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

fprintf("Place raw .czi files in: %s\n", dataRawDir);


% -------------------- Optional: run pipeline scripts --------------------
% Toggle these to true/false depending on what you want to run by default.
RUN_CZI_PROCESSING    = true;   % e.g., generates processed_results.csv, director_data_*.mat, etc.
RUN_COUNT_DEFECTTYPE  = true;   % count defects by type/charge
RUN_DEFECT_VELOCITY   = true;   % defect velocity analysis
RUN_P_DISTRIBUTION    = true;   % P(distribution) / angular distribution analysis
RUN_DEFECT_FITTING    = true;   % defect density decay fitting/plot
RUN_COLOR_THETA_MAPS  = true;   % theta colormap export
RUN_VELOCITY_OVERLAY  = true;   % velocity + directors + defects overlay
RUN_CELL_RATIO        = true;   % red/green area ratio check


% --- CZI processing / defect counting ---
if RUN_CZI_PROCESSING
    fprintf("\n[1/8] Running CZI processing / defect counting...\n");
    run_script_if_exists('cziprocessing_defectcounting');
end

% --- Defect density decay fit/plot ---
if RUN_DEFECT_FITTING
    fprintf("\n[2/8] Running defect density decay fitting...\n");
    run_script_if_exists('defect_density_decay');
end

% --- Theta maps ---
if RUN_COLOR_THETA_MAPS
    fprintf("\n[3/8] Exporting theta maps...\n");
    run_script_if_exists('color_angle_plot');  % rename to your actual script name
end

% --- Velocity overlay ---
if RUN_VELOCITY_OVERLAY
    fprintf("\n[4/8] Generating velocity/director/defect overlays...\n");
    run_script_if_exists('overlay_velocity_director'); % rename to your actual script name
end

% --- P-distribution analysis ---
if RUN_P_DISTRIBUTION
    fprintf("\n[5/8] Running P-distribution analysis...\n");
    run_script_if_exists('Pdistribution');
end

% --- Defect velocity analysis ---
if RUN_DEFECT_VELOCITY
    fprintf("\n[6/8] Running defect velocity analysis...\n");
    run_script_if_exists('defectvelocity');
end

% --- Count defect type (e.g., +1/2 vs -1/2, integer vs half-integer, etc.) ---
if RUN_COUNT_DEFECTTYPE
    fprintf("\n[7/8] Running count defect type...\n");
    run_script_if_exists('countdefecttype');
end

% --- Cell ratio verification ---
if RUN_CELL_RATIO
    fprintf("\n[8/8] Running cell ratio verification...\n");
    run_script_if_exists('cell_ratio_verification'); % rename to your actual script name
end

fprintf("\nDone. Outputs (if enabled) are saved under: %s\n", resultsDir);

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
    % Use literal text here so it still works even if something was cleared in base
    warning("Failed running %s: %s", scriptFile, ME.message);
end
end