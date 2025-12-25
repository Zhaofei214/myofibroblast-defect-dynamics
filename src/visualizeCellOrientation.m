function [directors, S] = visualizeCellOrientation(inputImage, params)
%--------------------------------------------------------------------------
% visualizeCellOrientation
%
% Purpose
%   Compute and visualize the local director (orientation) field of elongated
%   cells from a grayscale microscopy image using a structure-tensor approach.
%
% Description
%   This function estimates the local cell orientation by computing spatial
%   gradients of the input image, constructing and smoothing the structure
%   tensor, and extracting the principal orientation angle at each pixel.
%   The orientation field is sampled on a coarse grid and visualized as
%   line segments (directors) overlaid on the input image.
%
%   A scalar order parameter S is also computed to quantify the degree of
%   global nematic alignment in the field.
%
% Inputs
%   inputImage (numeric 2D array)
%       Grayscale image representing cell morphology or cytoskeletal signal.
%       Recommended to be normalized to [0, 1], but not strictly required.
%
%   params (struct)
%       Structure containing visualization and smoothing parameters:
%         params.gridSpacing    (double) Grid spacing for director sampling (pixels)
%         params.lineLength     (double) Length of director lines (pixels)
%         params.lineWidth      (double) Line width for visualization
%         params.smoothingSigma (double) Gaussian smoothing sigma (pixels)
%         params.lineColor      (char)   Line color (e.g., 'y')
%
% Outputs
%   directors (struct)
%       Structure containing the director field with fields:
%         directors.X           x-coordinates of grid points
%         directors.Y           y-coordinates of grid points
%         directors.orientation Local director angle (radians)
%
%   S (double)
%       Global nematic order parameter computed from the director field.
%
% Example
%   params.gridSpacing = 60;
%   params.lineLength = 45;
%   params.lineWidth = 2;
%   params.smoothingSigma = 20;
%   params.lineColor = 'y';
%
%   [directors, S] = visualizeCellOrientation(img, params);
%
% Dependencies
%   - Image Processing Toolbox (imgradientxy, imgaussfilt)
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
validateattributes(inputImage, {'numeric','logical'}, {'2d','nonempty'}, ...
    mfilename, 'inputImage', 1);

requiredFields = {'gridSpacing','lineLength','lineWidth','smoothingSigma','lineColor'};
for f = requiredFields
    if ~isfield(params, f{1})
        error("Missing params.%s", f{1});
    end
end

inputImage = double(inputImage);

% -----------------------------
% Structure tensor computation
% -----------------------------
[Gx, Gy] = imgradientxy(inputImage, "sobel");

Axx = Gx.^2;
Ayy = Gy.^2;
Axy = Gx .* Gy;

% Smooth tensor components
Axx = imgaussfilt(Axx, params.smoothingSigma);
Ayy = imgaussfilt(Ayy, params.smoothingSigma);
Axy = imgaussfilt(Axy, params.smoothingSigma);

% Orientation angle (nematic director)
orientation = 0.5 * atan2(2*Axy, Axx - Ayy);

% -----------------------------
% Grid sampling
% -----------------------------
[rows, cols] = size(inputImage);
x = params.gridSpacing:params.gridSpacing:cols;
y = params.gridSpacing:params.gridSpacing:rows;
[X, Y] = meshgrid(x, y);

orientationInterp = interp2(orientation, X, Y, 'linear');

% -----------------------------
% Visualization
% -----------------------------
imshow(inputImage, 'Colormap', gray, 'DisplayRange', [0 0.5]);
axis ij;
hold on;

for i = 1:size(X,1)
    for j = 1:size(X,2)
        angle = orientationInterp(i,j);
        if isnan(angle)
            continue;
        end

        dx = params.lineLength * cos(angle);
        dy = params.lineLength * sin(angle);

        line([X(i,j)-dx, X(i,j)+dx], ...
             [Y(i,j)-dy, Y(i,j)+dy], ...
             'Color', params.lineColor, ...
             'LineWidth', params.lineWidth);
    end
end

hold off;

% -----------------------------
% Output structure
% -----------------------------
directors = struct();
directors.X = X;
directors.Y = Y;
directors.orientation = orientationInterp;

% -----------------------------
% Global nematic order parameter
% -----------------------------
angles = orientationInterp(:);
angles = angles(~isnan(angles));

if isempty(angles)
    S = NaN;
else
    S = abs(mean(exp(1i * 2 * angles)));
end

end