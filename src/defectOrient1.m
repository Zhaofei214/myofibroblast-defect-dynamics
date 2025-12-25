function orientation = defectOrient1(directorField, defectPos, radius)
%--------------------------------------------------------------------------
% defectOrient1
%
% Purpose
%   Compute the average local director orientation surrounding a topological
%   defect within a specified neighborhood.
%
% Description
%   This function extracts the director orientations in a circular region
%   around a detected defect and computes a representative mean orientation.
%   The averaging accounts for nematic symmetry (θ ≡ θ + π).
%
%   The resulting orientation can be used to characterize defect polarity,
%   defect–cell alignment, or defect-induced anisotropy in active nematic
%   systems.
%
% Inputs
%   directorField (struct)
%       Director field structure with fields:
%         directorField.X            x-coordinates of grid points
%         directorField.Y            y-coordinates of grid points
%         directorField.orientation  director angles (radians)
%
%   defectPos (1x2 double)
%       Defect position specified as [x, y] in pixel coordinates.
%
%   radius (double)
%       Radius (in pixels) defining the neighborhood around the defect over
%       which orientations are averaged.
%
% Outputs
%   orientation (double)
%       Mean director orientation (radians) around the defect, accounting
%       for nematic symmetry.
%
% Example
%   ori = defectOrient1(directors, [x0, y0], 100);
%
% Notes
%   - NaN values in the director field are ignored.
%   - The nematic mean is computed using ⟨exp(2iθ)⟩.
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
requiredFields = {'X','Y','orientation'};
for f = requiredFields
    if ~isfield(directorField, f{1})
        error("directorField.%s is missing.", f{1});
    end
end

validateattributes(defectPos, {'numeric'}, {'numel',2}, ...
    mfilename, 'defectPos');

validateattributes(radius, {'numeric'}, {'scalar','positive'}, ...
    mfilename, 'radius');

X = directorField.X;
Y = directorField.Y;
theta = directorField.orientation;

% -----------------------------
% Distance mask
% -----------------------------
dist = sqrt((X - defectPos(1)).^2 + (Y - defectPos(2)).^2);
mask = dist <= radius;

localTheta = theta(mask);
localTheta = localTheta(~isnan(localTheta));

if isempty(localTheta)
    orientation = NaN;
    return;
end

% -----------------------------
% Nematic averaging
% -----------------------------
orientation = 0.5 * angle(mean(exp(1i * 2 * localTheta)));

end