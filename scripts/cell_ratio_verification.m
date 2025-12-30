%===========================================
% Cell ratio verification (magenta / green area)
%===========================================
% Description:
%   This script estimates the relative cell ratio between two fluorescent
%   channels (magenta and green) by thresholding each channel and computing
%   the relative area occupied by each cell population.
%
%   The ratio is reported as:
%       greenArea / (greenArea + magentaArea)
%
% Inputs (relative to repo root):
%   Data/<datasetName>/*.tif or *.png
%
% Outputs:
%   - Printed cell area ratios in the MATLAB command window
%
% License:
%   MIT License (see LICENSE file in the repository root).
%   Copyright (c) 2025 Zhaofei Zheng
%
% Notes:
%   - Thresholding uses Otsu’s method (graythresh).
%   - Small objects are removed using bwareaopen.
%===========================================

clear; clc;

% -------------------- Repo-safe paths (do not use pwd) --------------------
repoRoot = fileparts(fileparts(mfilename("fullpath"))); % scripts/ → repo root
dataRoot = fullfile(repoRoot, 'Data');

% ---------- Inputs ----------
datasetName = fullfile('sample', 'cell_ratio_verification');
inputDir = fullfile(dataRoot, datasetName);

fileNames = { ...
    '10%image.tif', ...
    '30%image.tif', ...
    '50%image.tif', ...
    '70%image.tif', ...
};

% ---------- Analysis ----------
for i = 1:length(fileNames)

    imgPath = fullfile(inputDir, fileNames{i});
    if ~isfile(imgPath)
        warning('Missing image: %s', imgPath);
        continue;
    end

    img = imread(imgPath);   % Read RGB image

    magentaChannel   = img(:,:,1);
    greenChannel = img(:,:,2);

    % ---- Thresholding (Otsu) ----
    magentaBW   = imbinarize(magentaChannel,   graythresh(magentaChannel));
    greenBW = imbinarize(greenChannel, graythresh(greenChannel));

    % ---- Remove small noise ----
    magentaBW   = bwareaopen(magentaBW,   50);
    greenBW = bwareaopen(greenBW, 50);

    % ---- Compute areas ----
    magentaArea   = sum(magentaBW(:));
    greenArea = sum(greenBW(:));

    relativeDensity = greenArea / (greenArea + magentaArea);

    fprintf('%s: Cell area ratio = %.1f%%\n', ...
        fileNames{i}, relativeDensity * 100);
end