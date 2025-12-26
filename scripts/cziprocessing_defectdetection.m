%% cziprocessing_defectcounting.m
% Pipeline script: read a multi-channel CZI time series, compute director field
% and defects for each frame, and export summary CSV + per-frame outputs.
%
% Expected repo structure:
%   external/bfmatlab/   (Bio-Formats MATLAB toolbox; not tracked unless licensed)
%   data/raw/            (raw .czi files; gitignored)
%   output/              (generated outputs; gitignored)
%
% Author: Zhaofei Zheng
% License: MIT (see LICENSE in repo root)

close all; clear; clc;

%% -------------------- Paths (portable) --------------------
repoRoot = fileparts(fileparts(mfilename("fullpath"))); % scripts/ -> repo root

% Bio-Formats
bfPath = fullfile(repoRoot, "external", "bfmatlab");
assert(exist(bfPath, "dir") == 7, ...
    "Bio-Formats not found. Put bfmatlab in: %s", bfPath);
addpath(bfPath);

% Input / Output directories
dataRawDir = fullfile(repoRoot, "data", "raw");
outRootDir = fullfile(repoRoot, "output");

if ~exist(outRootDir, "dir"); mkdir(outRootDir); end
if ~exist(dataRawDir, "dir")
    mkdir(dataRawDir);
    error("Created %s. Please place your .czi file there and rerun.", dataRawDir);
end

%% -------------------- User settings --------------------
% Input file (place under data/raw/)
cziName = "live_cell_imaging_every30min_10x_4tiles_750_50_03.czi";
cziFile = fullfile(dataRawDir, cziName);
assert(exist(cziFile, "file") == 2, "CZI file not found: %s", cziFile);

% Channel indices for bfGetPlane indexing:
% NOTE: Bio-Formats uses 0-based channel index in reader.getIndex().
MagentaChannel   = 0;   % adjust if needed
greenChannel = 1;   % adjust if needed (set to [] if not used)

% Output naming
runTag = "Processed_Results_750celldensity_50";
outputDir  = fullfile(outRootDir, runTag);
outputDir2 = fullfile(outputDir, "defectimages");
if ~exist(outputDir,  "dir"); mkdir(outputDir);  end
if ~exist(outputDir2, "dir"); mkdir(outputDir2); end

% Background correction parameter (pixels)
bgSigma = 80;

% Contrast stretch limits
stretchLim = [0.005 0.995];

%% -------------------- Init accumulator --------------------
vars = struct('defectFrames', [], ...
              'defectCharge', [], ...
              'defectX', [], ...
              'defectY', [], ...
              'defectOrientation', []);

save(fullfile(outputDir, "defectData.mat"), "-struct", "vars");

%% -------------------- Read CZI metadata --------------------
reader = bfGetReader(char(cziFile));

numSeries = reader.getSeriesCount();
reader.setSeries(0);

numChannels = reader.getSizeC();
numFrames   = reader.getSizeT();

fprintf("Bio-Formats reader loaded.\n");
fprintf("File: %s\nSeries: %d | Channels: %d | Frames: %d\n", cziName, numSeries, numChannels, numFrames);

%% -------------------- CSV output --------------------
csvFile = fullfile(outputDir, "processed_results.csv");
fid = fopen(csvFile, "w");
assert(fid ~= -1, "Cannot open CSV for writing: %s", csvFile);

fprintf(fid, "Frame,defectP,defectM,defectTotal,orderParameter\n");

%% -------------------- Analysis loop --------------------
for t = 1:numFrames
    % Read channels for this frame (Z=0, T=t-1)
    imgR = read_and_preprocess_plane(reader, MagentaChannel,   t, bgSigma, stretchLim);
    if ~isempty(greenChannel)
        imgG = read_and_preprocess_plane(reader, greenChannel, t, bgSigma, stretchLim);
        mergedImg = (imgR + imgG) / 2;
    else
        mergedImg = imgR;
    end

    % Create invisible figure for any overlay drawings done inside directorField
    figHandle = figure("Name","defectVisualization", ...
        "Position",[400 100 800 600], "Visible","off");
    axis tight;

    % Core analysis
    [defect_points_m, defect_points_p, orderParameter, directors, vars] = ...
        directorField(mergedImg, outputDir, vars);

    % Counts
    defectP = size(defect_points_p, 1);
    defectM = size(defect_points_m, 1);
    defectTotal = defectP + defectM;

    % Track frames for each defect (append one entry per defect)
    vars.defectFrames = [vars.defectFrames; t * ones(defectTotal, 1)];

    % Write CSV row
    fprintf(fid, "%d,%d,%d,%d,%.6f\n", t, defectP, defectM, defectTotal, orderParameter);

    % Save processed overlay image
    imgFilename = fullfile(outputDir2, sprintf("Processed_Merged_RG_Frame%d.png", t));
    saveas(figHandle, imgFilename);
    close(figHandle);

    % Save director field for this frame
    directorFieldName = fullfile(outputDir2, sprintf("director_data_%d.mat", t));
    save(directorFieldName, "directors");

    if mod(t, 10) == 0 || t == numFrames
        fprintf("Processed frame %d / %d\n", t, numFrames);
    end
end

%% -------------------- Finalize --------------------
save(fullfile(outputDir, "defectData.mat"), "-struct", "vars");
fclose(fid);

fprintf("Processing complete!\nCSV: %s\nMAT: %s\n", csvFile, fullfile(outputDir,"defectData.mat"));


%% ==================== Helper functions ====================

function imgOut = read_and_preprocess_plane(reader, ch, t, bgSigma, stretchLim)
% Read a single plane using Bio-Formats and apply background correction + contrast.
%
% Inputs:
%   ch: 0-based channel index (as used by reader.getIndex(0,ch,t-1))
%   t : 1-based frame index (MATLAB loop counter)
%
% Output:
%   imgOut: double image in [0,1]

    planeIndex = reader.getIndex(0, ch, t-1) + 1; % bfGetPlane is 1-based
    raw = bfGetPlane(reader, planeIndex);

    corrected = correct_gradient(raw, bgSigma);
    corrected = mat2gray(double(corrected));

    imgOut = imadjust(corrected, stretchlim(corrected, stretchLim), []);
end

function corrected = correct_gradient(img, sigma)
% Estimate a smooth background with Gaussian filtering and normalize it out.
    img = double(img);
    background = imgaussfilt(img, sigma);
    corrected = img ./ background;
    corrected = mat2gray(corrected);
end