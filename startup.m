function startup()
%STARTUP  Initialize paths for the myofibroblast defect dynamics repository
%
% This function adds all required source and script directories to the
% MATLAB path. It should be run once per MATLAB session after cloning
% the repository.
%
% Usage:
%   startup
%
% Author: Zhaofei Zheng
% License: MIT (see LICENSE file)

% Determine repository root (location of this file)
repoRoot = fileparts(mfilename("fullpath"));

% Define key subdirectories
srcDir     = fullfile(repoRoot, "src");
scriptsDir = fullfile(repoRoot, "scripts");

% Add paths 
if exist(srcDir, "dir")
    addpath(genpath(srcDir));
else
    warning("src directory not found: %s", srcDir);
end

if exist(scriptsDir, "dir")
    addpath(genpath(scriptsDir));
else
    warning("scripts directory not found: %s", scriptsDir);
end

fprintf("Repository paths added:\n");
fprintf("  %s\n", srcDir);
fprintf("  %s\n", scriptsDir);

end