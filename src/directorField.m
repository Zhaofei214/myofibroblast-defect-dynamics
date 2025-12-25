function [defect_points_m, defect_points_p, orderParameter, directors, vars] = directorField(img, outputDir, vars)
%--------------------------------------------------------------------------
% directorField
%
% Purpose
%   Compute the local director field from an input image, detect +1/2 and -1/2
%   topological defects, estimate a global order parameter, and update defect-
%   associated metadata (e.g., orientations) for downstream analysis.
%
% Description
%   This function is designed for the CZI processing pipeline used in the PNAS
%   manuscript "Myofibroblasts slow down defect recombination dynamics in mixed
%   cell monolayers." Given a (typically merged) fluorescence image, it:
%     1) Computes a director field (via visualizeCellOrientation.m)
%     2) Detects defects (via get_all_defects.m)
%     3) Computes defect orientations / bookkeeping (via defectOrient1.m)
%   Processing parameters are written to a text file in the output directory.
%
% Inputs
%   img        (numeric 2D array) Input image (grayscale). Recommended range [0,1]
%              but any numeric matrix is accepted.
%   outputDir  (char/string) Path to an output directory where parameters are logged.
%   vars       (struct) Accumulator struct used across frames. Expected fields:
%              defectFrames, defectCharge, defectX, defectY, defectOrientation
%              (fields may be empty at initialization).
%
% Outputs
%   defect_points_m (Nx2 double) [x,y] positions of -1/2 defects (pixels)
%   defect_points_p (Mx2 double) [x,y] positions of +1/2 defects (pixels)
%   orderParameter  (double) Scalar order parameter (if computed), otherwise NaN
%   directors       (struct/array) Director field representation returned by
%                   visualizeCellOrientation.m (format defined by your pipeline)
%   vars            (struct) Updated accumulator struct
%
% Example
%   img = im2double(imread('example.png'));
%   vars = struct('defectFrames',[],'defectCharge',[],'defectX',[],'defectY',[],'defectOrientation',[]);
%   [dm, dp, S, directors, vars] = directorField(img, "output", vars);
%
% Dependencies
%   - visualizeCellOrientation.m  (your repo)  : computes directors (and optionally orderParameter)
%   - get_all_defects.m           (your repo)  : detects defect positions
%   - defectOrient1.m             (your repo)  : computes defect orientations / updates vars
%   - MATLAB toolboxes commonly used:
%       Image Processing Toolbox (recommended)
%
% Citation
%   If you use or adapt this code, please cite the associated manuscript:
%   Z. Zheng et al., “Myofibroblasts slow down defect recombination dynamics
%   in mixed cell monolayers,” PNAS (in submission).
%
% Author
%   Zhaofei Zheng
%
% License
%   MIT License. See the LICENSE file in the repository root.
%--------------------------------------------------------------------------

% -----------------------------
% Input validation
% -----------------------------
if nargin < 3
    error("directorField requires three inputs: img, outputDir, vars.");
end

validateattributes(img, {'numeric','logical'}, {'2d','nonempty'}, mfilename, 'img', 1);

if ~(ischar(outputDir) || isstring(outputDir))
    error("outputDir must be a char or string.");
end
outputDir = char(outputDir);

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

if ~isstruct(vars)
    error("vars must be a struct accumulator.");
end

% -----------------------------
% Processing parameters
% -----------------------------
params = struct();
params.gridSpacing    = 60;   % Grid spacing for visualization / defect sampling (pixels)
params.lineLength     = 45;   % Length of orientation lines (pixels)
params.lineWidth      = 2;    % Line width for visualization
params.smoothingSigma = 20;   % Gaussian smoothing sigma
params.lineColor      = 'y';  % Line color for visualization
params.defthres       = 0.1;  % Defect detection threshold

% -----------------------------
% Log parameters (reviewer-friendly)
% -----------------------------
try
    paramFile = fullfile(outputDir, 'parameters.txt');
    fid = fopen(paramFile, 'w');
    if fid ~= -1
        fprintf(fid, "directorField parameters\n");
        fprintf(fid, "gridSpacing: %g\n", params.gridSpacing);
        fprintf(fid, "lineLength: %g\n", params.lineLength);
        fprintf(fid, "lineWidth: %g\n", params.lineWidth);
        fprintf(fid, "smoothingSigma: %g\n", params.smoothingSigma);
        fprintf(fid, "lineColor: %s\n", params.lineColor);
        fprintf(fid, "defthres: %g\n", params.defthres);
        fclose(fid);
    end
catch
    % Non-fatal: if file logging fails, continue.
end

% -----------------------------
% 1) Compute directors (+ orderParameter if available)
% -----------------------------
orderParameter = NaN;
directors = [];

if exist('visualizeCellOrientation', 'file') ~= 2
    error("Missing dependency: visualizeCellOrientation.m is not on the MATLAB path.");
end

% Try to call with two outputs if supported
try
    if nargout('visualizeCellOrientation') >= 2
        [directors, orderParameter] = visualizeCellOrientation(img, params);
    else
        directors = visualizeCellOrientation(img, params);
    end
catch ME
    error("visualizeCellOrientation failed: %s", ME.message);
end

% -----------------------------
% 2) Detect defects
% -----------------------------
defect_points_p = zeros(0,2);
defect_points_m = zeros(0,2);

if exist('get_all_defects', 'file') ~= 2
    error("Missing dependency: get_all_defects.m is not on the MATLAB path.");
end

try
    % Some implementations may optionally return vars as a 3rd output.
    if nargout('get_all_defects') >= 3
        [defect_points_p, defect_points_m, vars] = get_all_defects(directors, params.defthres, params.gridSpacing, vars);
    else
        [defect_points_p, defect_points_m] = get_all_defects(directors, params.defthres, params.gridSpacing);
    end
catch ME
    error("get_all_defects failed: %s", ME.message);
end

% Ensure outputs are Nx2
if ~isempty(defect_points_p)
    defect_points_p = double(defect_points_p);
    if size(defect_points_p,2) ~= 2
        error("defect_points_p must be an Nx2 array of [x,y] positions.");
    end
end
if ~isempty(defect_points_m)
    defect_points_m = double(defect_points_m);
    if size(defect_points_m,2) ~= 2
        error("defect_points_m must be an Nx2 array of [x,y] positions.");
    end
end

% -----------------------------
% 3) Defect orientation bookkeeping (optional but recommended)
% -----------------------------
if exist('defectOrient1', 'file') == 2
    try
        if ~isempty(defect_points_p)
            xP = defect_points_p(:,1); yP = defect_points_p(:,2);
        else
            xP = []; yP = [];
        end
        if ~isempty(defect_points_m)
            xM = defect_points_m(:,1); yM = defect_points_m(:,2);
        else
            xM = []; yM = [];
        end

        vars = defectOrient1(directors, xP, yP, xM, yM, params.gridSpacing, vars);
    catch ME
        % Non-fatal: keep defects even if orientation step fails
        warning("defectOrient1 failed (continuing): %s", ME.message);
    end
else
    warning("defectOrient1.m not found on path. Skipping defect orientation calculation.");
end

end