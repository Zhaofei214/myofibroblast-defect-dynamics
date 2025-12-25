function [defect_points_p, defect_points_m, vars] = get_all_defects(directors, defthres, gridSpacing, vars)
%--------------------------------------------------------------------------
% get_all_defects
%
% Purpose
%   Identify +1/2 and -1/2 topological defects from a sampled nematic
%   director field using winding number analysis.
%
% Description
%   This function analyzes a coarse-grained director field and detects
%   topological defects by evaluating angular changes around closed loops
%   on the director grid. Regions exhibiting large angular discontinuities
%   are classified as defect cores. Each detected defect is assigned a
%   topological charge of +1/2 or -1/2.
%
%   The function is designed to be called iteratively across time frames
%   and can update a cumulative accumulator structure (`vars`) that stores
%   defect positions and charges for downstream analysis.
%
% Inputs
%   directors (struct)
%       Director field structure with fields:
%         directors.X            x-coordinates of grid points
%         directors.Y            y-coordinates of grid points
%         directors.orientation  director angles (radians)
%
%   defthres (double)
%       Threshold for detecting sharp angular changes between neighboring
%       grid points. Typical value: 0.1–0.2 radians.
%
%   gridSpacing (double)
%       Spacing of the director grid in pixels. Used for spatial scaling.
%
%   vars (struct, optional)
%       Accumulator structure used across frames. Expected fields:
%         defectCharge, defectX, defectY
%       If empty or omitted, a new structure is created.
%
% Outputs
%   defect_points_p (Np x 2 double)
%       Pixel coordinates [x,y] of +1/2 defects.
%
%   defect_points_m (Nm x 2 double)
%       Pixel coordinates [x,y] of -1/2 defects.
%
%   vars (struct)
%       Updated accumulator structure containing appended defect metadata.
%
% Example
%   [dp, dm, vars] = get_all_defects(directors, 0.1, 60, vars);
%
% Dependencies
%   - MATLAB base functions
%   - Assumes nematic symmetry: theta and theta + pi are equivalent
%
% Citation
%   If you use or adapt this code, please cite:
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
    error("get_all_defects requires at least directors, defthres, gridSpacing.");
end

if nargin < 4 || isempty(vars)
    vars = struct('defectCharge', [], 'defectX', [], 'defectY', []);
end

requiredFields = {'X','Y','orientation'};
for f = requiredFields
    if ~isfield(directors, f{1})
        error("directors.%s is missing.", f{1});
    end
end

theta = directors.orientation;
X = directors.X;
Y = directors.Y;

validateattributes(theta, {'numeric'}, {'2d','nonempty'}, mfilename, 'directors.orientation');
validateattributes(defthres, {'numeric'}, {'scalar','positive'}, mfilename, 'defthres');
validateattributes(gridSpacing, {'numeric'}, {'scalar','positive'}, mfilename, 'gridSpacing');

% -----------------------------
% Initialize outputs
% -----------------------------
defect_points_p = zeros(0,2);
defect_points_m = zeros(0,2);

% Grid size
[nRows, nCols] = size(theta);

% -----------------------------
% Winding number detection
% -----------------------------
for i = 1:nRows-1
    for j = 1:nCols-1
        % Extract orientation around a plaquette
        loopTheta = [ ...
            theta(i,j), ...
            theta(i,j+1), ...
            theta(i+1,j+1), ...
            theta(i+1,j), ...
            theta(i,j) ];

        if any(isnan(loopTheta))
            continue;
        end

        % Compute angle differences with nematic symmetry
        dtheta = diff(loopTheta);
        dtheta = mod(dtheta + pi/2, pi) - pi/2;

        % Sum of angular changes (winding number)
        winding = sum(dtheta);

        % Classify defects
        if abs(winding - pi) < defthres
            % +1/2 defect
            x = mean([X(i,j), X(i,j+1), X(i+1,j+1), X(i+1,j)]);
            y = mean([Y(i,j), Y(i,j+1), Y(i+1,j+1), Y(i+1,j)]);
            defect_points_p(end+1,:) = [x,y];
            vars.defectCharge(end+1,1) = +0.5;
            vars.defectX(end+1,1) = x;
            vars.defectY(end+1,1) = y;

        elseif abs(winding + pi) < defthres
            % -1/2 defect
            x = mean([X(i,j), X(i,j+1), X(i+1,j+1), X(i+1,j)]);
            y = mean([Y(i,j), Y(i,j+1), Y(i+1,j+1), Y(i+1,j)]);
            defect_points_m(end+1,:) = [x,y];
            vars.defectCharge(end+1,1) = -0.5;
            vars.defectX(end+1,1) = x;
            vars.defectY(end+1,1) = y;
        end
    end
end

end