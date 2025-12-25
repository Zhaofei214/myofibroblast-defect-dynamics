function cmap = zfcolors(n)
%--------------------------------------------------------------------------
% zfcolors
%
% Purpose
%   Generate a custom perceptually ordered colormap used throughout the
%   analysis and visualization of nematic cell orientation and defect
%   dynamics.
%
% Description
%   This function returns a smoothly interpolated colormap designed for
%   visualizing orientation fields, density maps, and defect-related
%   quantities in cell monolayers. The colormap is defined by a set of
%   anchor RGB values and interpolated to the desired number of entries.
%
%   If no input is provided, the colormap defaults to the current figure's
%   colormap length.
%
% Inputs
%   n (integer, optional)
%       Number of color levels in the output colormap.
%       Default: size(get(gcf,'colormap'),1)
%
% Outputs
%   cmap (n x 3 double)
%       Colormap array with RGB values scaled to [0, 1].
%
% Example
%   colormap(zfcolors(256));
%   colorbar;
%
%   imagesc(data);
%   colormap(zfcolors);
%
% Notes
%   - This colormap is optimized for smooth visual transitions and is
%     suitable for both qualitative and quantitative visualization.
%   - Color values are interpolated linearly in RGB space.
%
% Citation
%   If you use or adapt this colormap, please cite:
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
% Default number of colors
% -----------------------------
if nargin < 1 || isempty(n)
    n = size(get(gcf,'colormap'), 1);
end

validateattributes(n, {'numeric'}, {'scalar','integer','positive'}, ...
    mfilename, 'n');

% -----------------------------
% Anchor color definition
% -----------------------------
% Original color anchors (RGB in [0,255])
baseColors = [ ...
     68,  1, 84;
     72, 35,116;
     64, 67,135;
     52, 94,141;
     41,120,142;
     32,144,140;
     34,167,132;
     68,190,112;
    121,209, 81;
    189,223, 38;
    253,231, 37 ];

% Normalize to [0,1]
baseColors = baseColors / 255;

% -----------------------------
% Interpolation
% -----------------------------
x = linspace(0, 1, size(baseColors,1));
xq = linspace(0, 1, n);

cmap = interp1(x, baseColors, xq, 'linear');

% Safety: ensure bounds
cmap = max(min(cmap,1),0);

end