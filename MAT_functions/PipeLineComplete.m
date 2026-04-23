%% HOMER3 CW-fNIRS pipeline using native Homer3 Block Average
% This script:
% 1) Loads a SNIRF file using Homer3
% 2) Converts intensity to OD
% 3) Optionally performs motion artifact detection
% 4) Optionally performs wavelet correction
% 5) Optionally performs spline correction
% 6) Optionally performs bandpass filtering
% 7) Converts OD to concentration (HbO / HbR / HbT)
% 8) Computes native Homer3 Block Average using hmrR_BlockAvg
% 9) Plots single-channel and all-channel Block Average
% 10) Builds HRF heatmap and HRF 3D surface
% 11) Plots probe geometry and an HRF metric on the probe
%


clear; close all; clc

%% ================= PATHS =================
HOMER3_DIR = 'C:\Users\catea\Downloads\Homer3-master';
addpath(genpath(HOMER3_DIR));

%% ================= INPUT =================
snirfFile = 'C:\Users\catea\OneDrive - Politecnico di Milano\Corsi\METHODS FOR THE FUNCTIONAL EVALUATION OF THE CENTRAL NERVOUS SYSTEM\2024-2025\DATI\CW_Alessandro Lia\2025-04-02 Fabio Negretti\2025-04-02_000\2025-04-02_001.snirf';

keepStims = ["1","2"];
% Only these stimulus labels will be kept if they exist in the file.

%% ================= USER OPTIONS =================
exampleCh   = 1;   % Example S-D pair index used for single-channel plots
hbTypeIdx   = 1;   % 1 = HbO, 2 = HbR, 3 = HbT
condIdxPlot = 1;   % Condition used for all-channel BlockAvg, heatmap, 3D plots

plotAllBlockAvgChannels = true;
channelsPerFigure = 9;

showProbe3D            = true;
showProbe3D_HRF        = true;
probeMetricHbTypeIdx   = 1;      % 1 = HbO, 2 = HbR, 3 = HbT
probeMetricTimeWindow  = [4 8];  % Time window (s) for probe HRF metric
probeUseWavelengthIdx  = 1;      % Wavelength index used to map channels on probe
showChannelLinks       = true;

%% ================= DENOISING OPTIONS =================
useMotionArtifactDetect    = true;
useMotionArtifactByChannel = false;
useWaveletMotionCorrection = true;
useSplineMotionCorrection  = true;
useBandpassFilter          = true;

% ---- Motion artifact detection parameters ----
tMotion     = 0.5;
tMask       = 1.0;
STDEVthresh = 50;
AMPthresh   = 5;

% ---- Wavelet parameters ----
iqr_wavelet    = 1.5;
turnon_wavelet = 1;

% ---- Spline parameters ----
pSpline = 0.99;

% ---- Bandpass parameters ----
hpf = 0.01;
lpf = 0.2;

% ---- MBLL parameters ----
dpf = [6 6];

%% ================= BLOCK AVERAGE PARAMETERS =================
trangeBA = [-2 20];

%% ================= SAVE =================
outMat = 'Homer3_results_nativeBlockAvg.mat';

%% ================= LOAD SNIRF =================
assert(exist(snirfFile,'file')==2, 'SNIRF not found: %s', snirfFile);
% Check that the file exists before trying to load it.
% Input:
%   snirfFile -> full path to the .snirf file
% Output:
%   none
% If the file is missing, execution stops here.

snirf   = SnirfClass(snirfFile);
% Homer3 class constructor.
% Input:
%   snirfFile -> path to a SNIRF file
% Output:
%   snirf -> a SNIRF container object with fields such as:
%            - data
%            - probe
%            - stim
%            - aux (if present)
%
% This is the main object that represents the recording.

dataObj = snirf.data;
% Extract the SNIRF.data container from the file.
% Input:
%   snirf -> SNIRF container
% Output:
%   dataObj -> Homer3 data container
%
% This contains the actual recorded time series.

probe   = snirf.probe;
% Extract the probe definition.
% Input:
%   snirf
% Output:
%   probe -> Homer3 ProbeClass object
%
% The probe stores:
%   - source positions
%   - detector positions
%   - wavelengths
%   - landmarks
%   - coordinate system


% Extract the first data object robustly
% Why needed:
% Homer3 object storage can vary depending on file/version.
d0 = get_first_data_element(dataObj);

% Extract time and numeric data
t  = get_time_from_data(d0);

% Custom helper.
% Input:
%   d0 -> one Homer3 data object
% Output:
%   Y0 -> numeric matrix [time x measurements]
%
% The helper tries:
%   - GetDataTimeSeries()
%   - .dataTimeSeries
%   - GetDataTimeSeries('reshape')
Y0 = get_data_matrix(d0);

nMeas = size(Y0, 2);

% Extract measurement list
% Typical columns in raw OD/intensity measurement list:
%   col 1 -> source index
%   col 2 -> detector index
%   col 3 -> data type / placeholder
%   col 4 -> wavelength index
ml0 = get_meas_list_safe(d0);

% Build active measurement lists in Homer3 format:
% [source detector active wavelength_or_type]
if size(ml0,2) >= 4
    mlActMan  = { [ml0(:,1:2), ones(size(ml0,1),1), ml0(:,4)] };
    mlActAuto = { [ml0(:,1:2), ones(size(ml0,1),1), ml0(:,4)] };
elseif size(ml0,2) >= 3
    mlActMan  = { [ml0(:,1:2), ones(size(ml0,1),1), ml0(:,3)] };
    mlActAuto = { [ml0(:,1:2), ones(size(ml0,1),1), ml0(:,3)] };
else
    error('Unexpected measurement list format in raw data.');
end

fprintf('Loaded: %s\n', snirfFile);
fprintf('Time points: %d | nMeas: %d\n', numel(t), nMeas);
fprintf('Wavelengths: %s\n', mat2str(get_probe_wavelengths(probe)));
fprintf('MeasList size: %d x %d\n', size(ml0,1), size(ml0,2));

%% ================= KEEP STIMS =================
hasStim = false;
try
    hasStim = ~isempty(snirf.stim);
catch
end

if hasStim
    stimNamesAll = get_stim_names(snirf.stim);
    % Custom helper.
% Input:
%   snirf.stim -> SNIRF stimulus container
% Output:
%   stimNamesAll -> string array with all condition labels
    fprintf('Stim names found in file:\n');
    disp(stimNamesAll)

    keepMask = ismember(lower(strtrim(stimNamesAll)), lower(strtrim(keepStims)));
% Input:
%   stimNamesAll -> names found in the file
%   keepStims    -> names you want to keep
% Output:
%   keepMask     -> logical mask

    if any(keepMask)
        snirf.stim = snirf.stim(keepMask);
        fprintf('Kept stims: %s\n', strjoin(get_stim_names(snirf.stim), ', '));
    else
        warning('Requested stimulus labels were not found. Using all available stimuli.');
        fprintf('Available stimuli: %s\n', strjoin(stimNamesAll, ', '));
    end
else
    warning('No stimuli found in SNIRF.');
end

condNames = get_stim_names(snirf.stim);
fprintf('Final number of conditions: %d\n', numel(condNames));
disp(condNames)

if isempty(condNames)
    error('No valid conditions available after stimulus selection.');
end

if condIdxPlot > numel(condNames)
    warning('condIdxPlot=%d is out of range. Resetting to 1.', condIdxPlot);
    condIdxPlot = 2;
end

%% ================= PIPELINE =================
% -------- 1) Intensity -> Optical Density --------
dodRawObj = hmrR_Intensity2OD(dataObj);

% -------- 2) Motion artifact detection (optional) --------
tInc   = [];
tIncCh = [];

if useMotionArtifactDetect
    try
        tInc = hmrR_MotionArtifact(dodRawObj, mlActMan, mlActAuto, ...
                                   tMotion, tMask, STDEVthresh, AMPthresh);
%tInc = inclusion mask over time:  1 = keep, 0 = artifact region
% Needed: Spline correction often needs a temporal mask telling it where the bad segments are.
    catch ME
        warning('Global motion artifact detection failed: %s', ME.message);
        tInc = [];
    end
end

if useMotionArtifactByChannel
    try
        tIncCh = hmrR_MotionArtifactByChannel(dodRawObj, mlActMan, mlActAuto, ...
                                              tMotion, tMask, STDEVthresh, AMPthresh);
    catch ME
        warning('Channel-wise motion artifact detection failed: %s', ME.message);
        tIncCh = [];
    end
end

% -------- 3) Wavelet correction (optional) --------
if useWaveletMotionCorrection
    dodWaveObj = hmrR_MotionCorrectWavelet(dodRawObj, mlActMan, mlActAuto, ...
                                           iqr_wavelet, turnon_wavelet);
else
    dodWaveObj = dodRawObj;
end

% -------- 4) Spline correction (optional) --------

if useSplineMotionCorrection
    try
        % If no motion mask exists, create a full inclusion mask
        if isempty(tInc)
            tInc_use = {ones(length(t),1)};
        else
            % Make sure tInc is a cell, because some Homer3 versions expect { ... }
            if iscell(tInc)
                tInc_use = tInc;
            else
                tInc_use = {tInc};
            end
        end

        dodSplineObj = hmrR_MotionCorrectSpline(dodWaveObj, t, tInc_use, pSpline);

    catch ME
        warning('Spline correction failed: %s. Using wavelet output instead.', ME.message);
        dodSplineObj = dodWaveObj;
    end
else
    dodSplineObj = dodWaveObj;
end
% -------- 5) Bandpass filtering (optional) --------
if useBandpassFilter
    dodFiltObj = hmrR_BandpassFilt(dodSplineObj, hpf, lpf);
else
    dodFiltObj = dodSplineObj;
end

% -------- 6) OD -> concentration --------
wl = get_probe_wavelengths(probe);
assert(numel(dpf)==numel(wl), ...
    'dpf length (%d) must match number of wavelengths (%d).', numel(dpf), numel(wl));

dcObj = hmrR_OD2Conc(dodFiltObj, probe, dpf);

%% ================= EXTRACT OD DATA =================
dodRaw    = get_data_matrix(get_first_data_element(dodRawObj));
dodWave   = get_data_matrix(get_first_data_element(dodWaveObj));
dodSpline = get_data_matrix(get_first_data_element(dodSplineObj));
dodFilt   = get_data_matrix(get_first_data_element(dodFiltObj));
ml        = get_meas_list_safe(get_first_data_element(dodRawObj));

nCh = size(dodRaw,2);
assert(exampleCh >= 1 && exampleCh <= max(1, estimate_num_sd_pairs(ml)), ...
    'exampleCh out of range.');

fprintf('OD raw size      : %d x %d\n', size(dodRaw,1), size(dodRaw,2));
fprintf('OD wavelet size  : %d x %d\n', size(dodWave,1), size(dodWave,2));
fprintf('OD spline size   : %d x %d\n', size(dodSpline,1), size(dodSpline,2));
fprintf('OD filtered size : %d x %d\n', size(dodFilt,1), size(dodFilt,2));

%% ================= OD PLOT: CHOSEN CHANNEL =================
% The example OD channel here refers to the raw measurement column index.
figure('Name', sprintf('OD example measurement %d', exampleCh), 'Color', 'w');
plot(t, dodRaw(:,exampleCh),    'k-', 'LineWidth', 1); hold on;
plot(t, dodWave(:,exampleCh),   'b-', 'LineWidth', 1.1);
plot(t, dodSpline(:,exampleCh), 'm-', 'LineWidth', 1.1);
plot(t, dodFilt(:,exampleCh),   'r-', 'LineWidth', 1.3);
grid on
xlabel('Time (s)');
ylabel('OD');
title(sprintf('OD example measurement | %s', make_channel_title(ml, exampleCh)), 'Interpreter', 'none');
legend({'Raw OD','After wavelet','After spline','After bandpass'}, 'Location', 'best');

%% ================= OD PLOT: RAW VS FINAL =================
figure('Name', sprintf('OD raw vs filtered - measurement %d', exampleCh), 'Color', 'w');
plot(t, dodRaw(:,exampleCh),  'k-', 'LineWidth', 1); hold on;
plot(t, dodFilt(:,exampleCh), 'r-', 'LineWidth', 1.3);
grid on
xlabel('Time (s)');
ylabel('OD');
title(sprintf('OD raw vs final filtered | %s', make_channel_title(ml, exampleCh)), 'Interpreter', 'none');
legend({'Raw OD','Final filtered OD'}, 'Location', 'best');

%% ================= OD HEATMAPS =================
figure('Name','Raw OD heatmap','Color','w');
imagesc(t, 1:nCh, dodRaw');
axis xy
xlabel('Time (s)');
ylabel('Measurement');
title('Raw OD - all measurements');
colorbar

figure('Name','Wavelet-corrected OD heatmap','Color','w');
imagesc(t, 1:nCh, dodWave');
axis xy
xlabel('Time (s)');
ylabel('Measurement');
title('Wavelet-corrected OD - all measurements');
colorbar

figure('Name','Spline-corrected OD heatmap','Color','w');
imagesc(t, 1:nCh, dodSpline');
axis xy
xlabel('Time (s)');
ylabel('Measurement');
title('Spline-corrected OD - all measurements');
colorbar

figure('Name','Final filtered OD heatmap','Color','w');
imagesc(t, 1:nCh, dodFilt');
axis xy
xlabel('Time (s)');
ylabel('Measurement');
title('Final filtered OD - all measurements');
colorbar

%% ================= NATIVE HOMER3 BLOCK AVERAGE =================
% Each condition is processed separately using the native Homer3 stim container.
% The output is a Homer3 data object. We then extract:
% - numeric matrix
% - time vector
% - measurement list
% and reconstruct HbO/HbR/HbT channel-by-channel.

HRF_by_condition  = cell(numel(condNames),1);
tHRF_by_condition = cell(numel(condNames),1);
mlBA_by_condition = cell(numel(condNames),1);

figure('Name','BlockAvg single channel HbO/HbR','Color','w');
nCond = numel(condNames);

for c = 1:numel(condNames)

    stim_single = snirf.stim(c);

    [dataAvgBA, dataStdBA, nTrialsBA, dataSum2BA, yTrialsBA] = ...
        hmrR_BlockAvg(dcObj, stim_single, trangeBA);

    yBA = get_data_matrix(get_first_data_element(dataAvgBA));
    tHRF_BA = get_time_from_data(get_first_data_element(dataAvgBA));
    mlBA = get_meas_list_safe(get_first_data_element(dataAvgBA));

    HRF_by_condition{c}  = yBA;
    tHRF_by_condition{c} = tHRF_BA(:);
    mlBA_by_condition{c} = mlBA;

    subplot(nCond,1,c); hold on; grid on

    % Reconstruct HbO/HbR/HbT for the chosen S-D pair index
    idxHbO = find_hb_channels(mlBA, exampleCh, 'HbO');
    idxHbR = find_hb_channels(mlBA, exampleCh, 'HbR');
    idxHbT = find_hb_channels(mlBA, exampleCh, 'HbT');

    if ~isempty(idxHbO)
        plot(tHRF_BA, yBA(:,idxHbO(1)), 'r-', 'LineWidth', 2);
    end
    if ~isempty(idxHbR)
        plot(tHRF_BA, yBA(:,idxHbR(1)), 'b-', 'LineWidth', 2);
    end
    if ~isempty(idxHbT)
        plot(tHRF_BA, yBA(:,idxHbT(1)), 'k--', 'LineWidth', 1.5);
        legend({'HbO','HbR','HbT'}, 'Location', 'best');
    else
        legend({'HbO','HbR'}, 'Location', 'best');
    end

    xlabel('Time (s)');
    ylabel('\DeltaHb (\muM)');
    title(sprintf('BlockAvg | %s | channel %d', condNames(c), exampleCh), ...
        'Interpreter','none');
end

%% ================= BLOCKAVG PLOT FOR ALL CHANNELS =================
if plotAllBlockAvgChannels
    yBA_plot  = HRF_by_condition{condIdxPlot};
    tHRF_plot = tHRF_by_condition{condIdxPlot};
    mlBA_plot = mlBA_by_condition{condIdxPlot};

    if isempty(yBA_plot)
        warning('Selected condition has no BlockAvg output.');
    else
        nPairs = estimate_num_sd_pairs(mlBA_plot);

        nFig = ceil(nPairs / channelsPerFigure);
        nRow = ceil(sqrt(channelsPerFigure));
        nCol = ceil(channelsPerFigure / nRow);

        for f = 1:nFig
            figure('Name', sprintf('BlockAvg channels %d', f), 'Color', 'w');

            chStart = (f-1)*channelsPerFigure + 1;
            chEnd   = min(f*channelsPerFigure, nPairs);

            for ch = chStart:chEnd
                subplot(nRow, nCol, ch-chStart+1)
                hold on; grid on

                idxHbO = find_hb_channels(mlBA_plot, ch, 'HbO');
                idxHbR = find_hb_channels(mlBA_plot, ch, 'HbR');
                idxHbT = find_hb_channels(mlBA_plot, ch, 'HbT');

                if ~isempty(idxHbO)
                    plot(tHRF_plot, yBA_plot(:,idxHbO(1)), 'r-', 'LineWidth', 1.5);
                end
                if ~isempty(idxHbR)
                    plot(tHRF_plot, yBA_plot(:,idxHbR(1)), 'b-', 'LineWidth', 1.5);
                end
                if ~isempty(idxHbT)
                    plot(tHRF_plot, yBA_plot(:,idxHbT(1)), 'k--', 'LineWidth', 1.0);
                end

                xlabel('Time (s)');
                ylabel('\DeltaHb');
                title(sprintf('Ch %d', ch));
            end

            sgtitle(sprintf('BlockAvg per channel | %s', condNames(condIdxPlot)), ...
                'Interpreter', 'none');
        end
    end
end

%% ================= HRF HEATMAP + HRF 3D =================
yBA_plot  = HRF_by_condition{condIdxPlot};
tHRF_plot = tHRF_by_condition{condIdxPlot};
mlBA_plot = mlBA_by_condition{condIdxPlot};

if isempty(yBA_plot)
    warning('HRF for selected condition is empty. Skipping HRF heatmap and 3D plot.');
else
    hbName = hb_idx_to_name(hbTypeIdx);

    nPairs = estimate_num_sd_pairs(mlBA_plot);
    HRF_all = nan(length(tHRF_plot), nPairs);

    for ch = 1:nPairs
        idxHb = find_hb_channels(mlBA_plot, ch, hbName);
        if ~isempty(idxHb)
            HRF_all(:,ch) = yBA_plot(:,idxHb(1));
        end
    end

    figure('Name','HRF heatmap','Color','w');
    imagesc(tHRF_plot, 1:nPairs, HRF_all');
    axis xy
    xlabel('Time (s)');
    ylabel('Channel');
    title(sprintf('HRF heatmap | %s | %s', condNames(condIdxPlot), hbName), ...
        'Interpreter','none');
    colorbar

    [T_hrf, Ch_hrf] = meshgrid(tHRF_plot, 1:nPairs);

    figure('Name','HRF 3D surface','Color','w');
    surf(T_hrf, Ch_hrf, HRF_all', 'EdgeColor', 'none');
    xlabel('Time (s)');
    ylabel('Channel');
    zlabel('\DeltaHb (\muM)');
    title(sprintf('HRF 3D | %s | %s', condNames(condIdxPlot), hbName), ...
        'Interpreter','none');
    colorbar
    view(45,30)
    grid on
end

%% ================= PROBE GEOMETRY 2D/3D =================
if showProbe3D
    srcPos  = get_probe_srcPos(probe);
    detPos  = get_probe_detPos(probe);
    landPos = get_probe_landmarkPos(probe);

    if isempty(srcPos) || isempty(detPos)
        warning('Probe plot skipped: source/detector coordinates not available.');
    else
        is3D = (size(srcPos,2) >= 3) && (size(detPos,2) >= 3);

        if is3D
            figure('Name','Probe geometry 3D','Color','w'); hold on

            if ~isempty(landPos) && size(landPos,2) >= 3
                scatter3(landPos(:,1), landPos(:,2), landPos(:,3), 30, [0.7 0.7 0.7], 'filled');
            end

            scatter3(srcPos(:,1), srcPos(:,2), srcPos(:,3), 90, 'r', 'filled');
            scatter3(detPos(:,1), detPos(:,2), detPos(:,3), 90, 'b', 'filled');

            for i = 1:size(srcPos,1)
                text(srcPos(i,1), srcPos(i,2), srcPos(i,3), sprintf(' S%d', i), ...
                    'Color', 'r', 'FontSize', 9);
            end

            for i = 1:size(detPos,1)
                text(detPos(i,1), detPos(i,2), detPos(i,3), sprintf(' D%d', i), ...
                    'Color', 'b', 'FontSize', 9);
            end

            if showChannelLinks && ~isempty(ml)
                for ch = 1:size(ml,1)
                    sIdx = ml(ch,1);
                    dIdx = ml(ch,2);

                    if sIdx <= size(srcPos,1) && dIdx <= size(detPos,1)
                        pS = srcPos(sIdx,:);
                        pD = detPos(dIdx,:);
                        plot3([pS(1) pD(1)], [pS(2) pD(2)], [pS(3) pD(3)], ...
                            '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
                    end
                end
            end

            axis equal
            grid on
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title('Probe geometry 3D');
            view(3)

        else
            figure('Name','Probe geometry 2D','Color','w'); hold on

            if ~isempty(landPos) && size(landPos,2) >= 2
                scatter(landPos(:,1), landPos(:,2), 30, [0.7 0.7 0.7], 'filled');
            end

            scatter(srcPos(:,1), srcPos(:,2), 90, 'r', 'filled');
            scatter(detPos(:,1), detPos(:,2), 90, 'b', 'filled');

            for i = 1:size(srcPos,1)
                text(srcPos(i,1), srcPos(i,2), sprintf(' S%d', i), ...
                    'Color', 'r', 'FontSize', 9);
            end

            for i = 1:size(detPos,1)
                text(detPos(i,1), detPos(i,2), sprintf(' D%d', i), ...
                    'Color', 'b', 'FontSize', 9);
            end

            if showChannelLinks && ~isempty(ml)
                for ch = 1:size(ml,1)
                    sIdx = ml(ch,1);
                    dIdx = ml(ch,2);

                    if sIdx <= size(srcPos,1) && dIdx <= size(detPos,1)
                        pS = srcPos(sIdx,:);
                        pD = detPos(dIdx,:);
                        plot([pS(1) pD(1)], [pS(2) pD(2)], ...
                            '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
                    end
                end
            end

            axis equal
            grid on
            xlabel('X'); ylabel('Y');
            title('Probe geometry 2D');
        end
    end
end

%% ================= PROBE HRF MAP 2D/3D =================
if showProbe3D_HRF
    srcPos  = get_probe_srcPos(probe);
    detPos  = get_probe_detPos(probe);
    landPos = get_probe_landmarkPos(probe);

    yBA_probe  = HRF_by_condition{condIdxPlot};
    tHRF_probe = tHRF_by_condition{condIdxPlot};
    mlBA_probe = mlBA_by_condition{condIdxPlot};

    if isempty(srcPos) || isempty(detPos)
        warning('Probe HRF map skipped: source/detector coordinates not available.');
    elseif isempty(yBA_probe)
        warning('Probe HRF map skipped: HRF not available for selected condition.');
    else
        hbName = hb_idx_to_name(probeMetricHbTypeIdx);

        tMask = tHRF_probe >= probeMetricTimeWindow(1) & tHRF_probe <= probeMetricTimeWindow(2);
        if ~any(tMask)
            warning('No HRF samples in selected time window.');
        else
            nPairs = estimate_num_sd_pairs(mlBA_probe);
            metricAllChannels = nan(1,nPairs);

            for ch = 1:nPairs
                idxHb = find_hb_channels(mlBA_probe, ch, hbName);
                if ~isempty(idxHb)
                    metricAllChannels(ch) = mean(yBA_probe(tMask, idxHb(1)), 1, 'omitnan');
                end
            end

            idxWl = get_meas_idx_for_wavelength(ml, probeUseWavelengthIdx);
            pairMap = unique(ml(:,1:2), 'rows', 'stable');

            % Use all available S-D pairs from BlockAvg output.
            idxUsePairs = 1:nPairs;
            metricUse = metricAllChannels(idxUsePairs);

            [midPos, srcDetPairs] = compute_pair_midpoints(pairMap, srcPos, detPos, idxUsePairs);

            is3D = (size(srcPos,2) >= 3) && (size(detPos,2) >= 3);

            if is3D
                figure('Name','Probe HRF map 3D','Color','w'); hold on

                if ~isempty(landPos) && size(landPos,2) >= 3
                    scatter3(landPos(:,1), landPos(:,2), landPos(:,3), 30, [0.7 0.7 0.7], 'filled');
                end

                scatter3(srcPos(:,1), srcPos(:,2), srcPos(:,3), 70, 'r', 'filled');
                scatter3(detPos(:,1), detPos(:,2), detPos(:,3), 70, 'b', 'filled');

                if showChannelLinks
                    for i = 1:size(srcDetPairs,1)
                        sIdx = srcDetPairs(i,1);
                        dIdx = srcDetPairs(i,2);
                        pS = srcPos(sIdx,:);
                        pD = detPos(dIdx,:);
                        plot3([pS(1) pD(1)], [pS(2) pD(2)], [pS(3) pD(3)], ...
                            '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
                    end
                end

                scatter3(midPos(:,1), midPos(:,2), midPos(:,3), 120, metricUse(:), 'filled');

                axis equal
                grid on
                xlabel('X'); ylabel('Y'); zlabel('Z');
                title(sprintf('Probe HRF map 3D | %s | %s | mean %.1f-%.1f s', ...
                    condNames(condIdxPlot), hbName, ...
                    probeMetricTimeWindow(1), probeMetricTimeWindow(2)), ...
                    'Interpreter', 'none');
                colorbar
                view(3)

            else
                figure('Name','Probe HRF map 2D','Color','w'); hold on

                if ~isempty(landPos) && size(landPos,2) >= 2
                    scatter(landPos(:,1), landPos(:,2), 30, [0.7 0.7 0.7], 'filled');
                end

                scatter(srcPos(:,1), srcPos(:,2), 70, 'r', 'filled');
                scatter(detPos(:,1), detPos(:,2), 70, 'b', 'filled');

                if showChannelLinks
                    for i = 1:size(srcDetPairs,1)
                        sIdx = srcDetPairs(i,1);
                        dIdx = srcDetPairs(i,2);
                        pS = srcPos(sIdx,:);
                        pD = detPos(dIdx,:);
                        plot([pS(1) pD(1)], [pS(2) pD(2)], ...
                            '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
                    end
                end

                scatter(midPos(:,1), midPos(:,2), 120, metricUse(:), 'filled');

                axis equal
                grid on
                xlabel('X'); ylabel('Y');
                title(sprintf('Probe HRF map 2D | %s | %s | mean %.1f-%.1f s', ...
                    condNames(condIdxPlot), hbName, ...
                    probeMetricTimeWindow(1), probeMetricTimeWindow(2)), ...
                    'Interpreter', 'none');
                colorbar
            end
        end
    end
end

%% ================= SAVE =================
save(outMat, ...
    'snirfFile','keepStims', ...
    'useMotionArtifactDetect','useMotionArtifactByChannel', ...
    'useWaveletMotionCorrection','useSplineMotionCorrection','useBandpassFilter', ...
    'tMotion','tMask','STDEVthresh','AMPthresh', ...
    'iqr_wavelet','pSpline','hpf','lpf','dpf', ...
    'exampleCh','hbTypeIdx','condIdxPlot', ...
    'dodRaw','dodWave','dodSpline','dodFilt','ml', ...
    'HRF_by_condition','tHRF_by_condition','mlBA_by_condition','probe');

fprintf('Saved: %s\n', outMat);

%% ===================== LOCAL HELPERS =====================

function d0 = get_first_data_element(dataObj)
% Return the first data element robustly.
try
    d0 = dataObj(1);
    return
catch
end

try
    if iscell(dataObj)
        d0 = dataObj{1};
        return
    end
catch
end

d0 = dataObj;
end

function Y = get_data_matrix(dataObj)
% Extract numeric data matrix [time x measurements].
try
    Y = dataObj.GetDataTimeSeries();
    return
catch
end

try
    Y = dataObj.dataTimeSeries;
    return
catch
end

try
    Y = dataObj.GetDataTimeSeries('reshape');
    return
catch
end

error('Unable to extract data matrix from dataObj.');
end

function t = get_time_from_data(dataObj)
% Extract time vector robustly.
try
    t = dataObj.GetTime();
    t = t(:);
    return
catch
end

try
    t = dataObj.time;
    t = t(:);
    return
catch
end

try
    t = dataObj.t;
    t = t(:);
    return
catch
end

error('Unable to extract time vector from dataObj.');
end

function ml = get_meas_list_safe(dataObj)
% Extract measurement list robustly.
try
    ml = dataObj.GetMeasList('reshape');
    return
catch
end

try
    ml = dataObj.GetMeasList();
    return
catch
end

try
    ml = dataObj.measList;
    return
catch
end

error('Unable to extract measurement list from dataObj.');
end

function wl = get_probe_wavelengths(probe)
% Extract wavelength vector from probe.
try
    wl = probe.GetWavelengths();
catch
    try
        wl = probe.wavelengths;
    catch
        try
            wl = probe.lambda;
        catch
            wl = [];
        end
    end
end
wl = wl(:)';
end

function names = get_stim_names(stimArr)
% Return stimulus names as a string array.
if isempty(stimArr)
    names = string.empty(0,1);
    return
end

names = strings(numel(stimArr),1);
for i = 1:numel(stimArr)
    names(i) = get_stim_name_one(stimArr(i));
end
end

function name = get_stim_name_one(stim)
% Extract one stimulus name robustly.
try
    name = string(stim.name);
    return
catch
end

try
    name = string(stim.GetName());
    return
catch
end

name = "stim";
end

function data = get_stim_data_one(stim)
% Extract stimulus data, usually [onset duration amplitude].
try
    data = stim.data;
    return
catch
end

try
    data = stim.GetData();
    return
catch
end

data = [];
end

function ttl = make_channel_title(ml, ch)
% Build a readable title for one measurement.
if isempty(ml)
    ttl = sprintf('Ch %d', ch);
    return
end

if size(ml,2) >= 4
    ttl = sprintf('Ch %d | S%d-D%d | wlIdx %d', ch, ml(ch,1), ml(ch,2), ml(ch,4));
elseif size(ml,2) >= 2
    ttl = sprintf('Ch %d | S%d-D%d', ch, ml(ch,1), ml(ch,2));
else
    ttl = sprintf('Ch %d', ch);
end
end

function srcPos = get_probe_srcPos(probe)
% Extract source positions, preferring 3D then 2D.
srcPos = [];

try
    srcPos = probe.sourcePos3D;
    if ~isempty(srcPos), return; end
catch
end

try
    srcPos = probe.GetSrcPos();
    if ~isempty(srcPos), return; end
catch
end

try
    srcPos = probe.sourcePos2D;
    if ~isempty(srcPos), return; end
catch
end
end

function detPos = get_probe_detPos(probe)
% Extract detector positions, preferring 3D then 2D.
detPos = [];

try
    detPos = probe.detectorPos3D;
    if ~isempty(detPos), return; end
catch
end

try
    detPos = probe.GetDetPos();
    if ~isempty(detPos), return; end
catch
end

try
    detPos = probe.detectorPos2D;
    if ~isempty(detPos), return; end
catch
end
end

function landPos = get_probe_landmarkPos(probe)
% Extract landmark positions, preferring 3D then 2D.
landPos = [];

try
    landPos = probe.landmarkPos3D;
    if ~isempty(landPos), return; end
catch
end

try
    landPos = probe.landmarkPos2D;
    if ~isempty(landPos), return; end
catch
end
end

function idx = get_meas_idx_for_wavelength(ml, wlIdxWanted)
% Return indices of raw OD measurements belonging to a given wavelength index.
idx = [];
if isempty(ml)
    return
end

if size(ml,2) >= 4
    idx = find(ml(:,4) == wlIdxWanted);
end
end

function idx = find_hb_channels(ml, chPairIdx, hbName)
% Find columns corresponding to one S-D pair and one Hb type.
%
% This helper assumes that BlockAvg concentration output is organized so that
% each unique source-detector pair appears up to 3 times in order:
% HbO, HbR, HbT.
%
% This is a practical reconstruction step for Homer3 native BlockAvg output.

idx = [];
if isempty(ml)
    return
end

pairs = unique(ml(:,1:2), 'rows', 'stable');

if chPairIdx > size(pairs,1)
    return
end

targetPair = pairs(chPairIdx,:);
pairIdx = find(ml(:,1)==targetPair(1) & ml(:,2)==targetPair(2));

switch hbName
    case 'HbO'
        if numel(pairIdx) >= 1, idx = pairIdx(1); end
    case 'HbR'
        if numel(pairIdx) >= 2, idx = pairIdx(2); end
    case 'HbT'
        if numel(pairIdx) >= 3, idx = pairIdx(3); end
end
end

function nPairs = estimate_num_sd_pairs(ml)
% Estimate the number of unique source-detector pairs.
if isempty(ml)
    nPairs = 0;
    return
end
pairs = unique(ml(:,1:2), 'rows', 'stable');
nPairs = size(pairs,1);
end

function hbName = hb_idx_to_name(hbTypeIdx)
% Map Hb index to Hb name.
switch hbTypeIdx
    case 1
        hbName = 'HbO';
    case 2
        hbName = 'HbR';
    case 3
        hbName = 'HbT';
    otherwise
        hbName = 'HbO';
end
end

function [midPos, srcDetPairs] = compute_pair_midpoints(pairMap, srcPos, detPos, idxPairs)
% Compute midpoints for source-detector pairs.
n = numel(idxPairs);
midPos = nan(n,3);
srcDetPairs = nan(n,2);

for i = 1:n
    pairIdx = idxPairs(i);

    if pairIdx > size(pairMap,1)
        continue
    end

    sIdx = pairMap(pairIdx,1);
    dIdx = pairMap(pairIdx,2);

    if sIdx <= size(srcPos,1) && dIdx <= size(detPos,1)
        pS = srcPos(sIdx,:);
        pD = detPos(dIdx,:);

        if size(pS,2) == 2
            pS(3) = 0;
        end
        if size(pD,2) == 2
            pD(3) = 0;
        end

        midPos(i,:) = (pS + pD)/2;
        srcDetPairs(i,:) = [sIdx dIdx];
    end
end
end