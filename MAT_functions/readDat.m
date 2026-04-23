function curves = readDat(Settings)
% Marco Nabacino PoliMi
%READDAT read .dat output from TRS.   
%   CURVES = READDAT(SETTINGS) reads the data according to the specified
%   settings and stores the result in CURVES.
%   
%
%   Input:
%   SETTINGS is a structure that must have the following fields:
%   - inputFile: the path to the .dat file to read
%   - nChannels: the number of channels of the acquisition board
%   - nCurvesTot: the number of curves in the file
%
%   Additionally, the following fields can be specified. If not, default
%   values are used:
%   - dataType: the encoding of the curve data (default 'uint16')
%   - headerSize: the number of elements in the header (default 764)
%   - headerType: the encoding of the header (default 'uint8')
%   - subheaderSize: the number of elements in each subheader (default 204)
%   - subheaderType: the encoding of each subheader (default 'uint8')
%
%
%   Output:
%   CURVES is a nChannels-by-nCurves matrix whose columns are the
%   time-resolved curves stored in the file. 

%% Read settings
nChannels = Settings.nChannels;
nCurvesTot = Settings.nCurvesTot;
inputFile = Settings.inputFile;

if (isfield(Settings, 'headerSize'))
    headerSize = Settings.headerSize;
else
    headerSize = 764;
end
if (isfield(Settings, 'headerType'))
    headerType = Settings.headerType;
else
    headerType = 'uint8';
end
if (isfield(Settings, 'subheaderSize'))
    subheaderSize = Settings.subheaderSize;
else
    subheaderSize = 204;
end
if (isfield(Settings, 'subheaderType'))
    subheaderType = Settings.subheaderType;
else
    subheaderType = 'uint8';
end
if (isfield(Settings, 'dataType'))
    dataType = Settings.dataType;
else
    dataType = 'uint16';
end

%% Read all curves.
% Initialize data matrix.
curves = zeros(nChannels, nCurvesTot);

% Open data file.
fidIn = fopen(inputFile, 'r');

% Read header.
fread(fidIn, headerSize, headerType);

% Read curves. Every curve is preceded by a subheader.
for curveIdx = 1:nCurvesTot
    fread(fidIn, subheaderSize, subheaderType);
    curves(:, curveIdx) = fread(fidIn, nChannels, dataType);
end
fclose(fidIn);
end