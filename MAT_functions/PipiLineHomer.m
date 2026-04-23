%% HOMER3 CW-fNIRS pipeline (single script) + BlockAvg + GLM
clear; close all; clc

%% -------- PATHS (metti SOLO Homer3; evita Homer2 per conflitti) --------
% HOMER3_DIR = 'C:\Users\catea\Downloads\Homer3-master';  % <-- ok
HOMER3_DIR = 'C:\Users\catea\Downloads\Homer3-master'; 
% addpath(genpath('C:\Users\Public\homer3\'))
addpath(genpath(HOMER3_DIR));

%% -------- INPUT --------
snirfFile = 'C:\Users\catea\OneDrive - Politecnico di Milano\Corsi\METHODS FOR THE FUNCTIONAL EVALUATION OF THE CENTRAL NERVOUS SYSTEM\2024-2025\CW_Augusto\2025-04-02_001\2025-04-02_001.snirf';

keepStims  = ["back1","back2","control"];


%% -------- PREPROC PARAMS --------
iqr_wavelet = 1.5;
turnon_wavelet = 1;

hpf = 0.01;
lpf = 0.2;

ppf = [6 6];     % deve matchare #wavelengths

%% -------- GLM PARAMS --------
trange = [-2 20];
glmSolveMethod = 1;    % OLS
idxBasis = 2;          % gamma canonical
paramsBasis = [0.1 3.0 10.0  1.8 3.0 10.0];
rhoSD_ssThresh = 15;   % mm
flagSSmethod = 0;
driftOrder = 3;
flagMotionCorrect = 0;

outMat = 'Homer3_results.mat';

%% ================= LOAD SNIRF =================
assert(exist(snirfFile,'file')==2, 'SNIRF not found: %s', snirfFile);

sn = SnirfClass(snirfFile);
methods(sn.stim)
snirf = SnirfClass(snirfFile);   % invece di SnirfLoad
dcObj = snirf.data;              % oppure snirf.data(1) se array
dataObj = snirf.data;      % OGGETTO (non matrice)
probe   = snirf.probe;

t = get_time_from_data(dataObj);

% nMeas: numero misure direttamente dai dati (robusto con la tua versione)
if numel(dataObj) > 1
    d0 = dataObj(1);
else
    d0 = dataObj;
end
nMeas = size(d0.GetDataTimeSeries(), 2);

mlActMan  = {ones(nMeas,1)};
mlActAuto = {ones(nMeas,1)};

fprintf('Loaded: %s\n', snirfFile);
fprintf('Time points: %d | nMeas: %d\n', numel(t), nMeas);
fprintf('Wavelengths: %s\n', mat2str(get_probe_wavelengths(probe)));

%% ================= KEEP STIMS =================
if isfield(snirf,'stim') && ~isempty(snirf.stim)
    stimNames = get_stim_names(snirf.stim);
    snirf.stim = snirf.stim(ismember(stimNames, keepStims));
    fprintf('Kept stims: %s\n', strjoin(get_stim_names(snirf.stim), ', '));
else
    warning('No stim found in SNIRF.');
end

% stim matrix s (time x nCond)
s = stim2sMatrix(snirf.stim, t);
condNames = get_stim_names(snirf.stim);

%% ================= PIPELINE =================
% 1) Intensity -> OD
dodObj = hmrR_Intensity2OD(dataObj);

% 2) Motion correction wavelet (tua signature)
dodObj = hmrR_MotionCorrectWavelet(dodObj, mlActMan, mlActAuto, iqr_wavelet, turnon_wavelet);

% 3) Bandpass
dodObj = hmrR_BandpassFilt(dodObj, hpf, lpf);

% 4) OD -> Hb (MBLL)
wl = get_probe_wavelengths(probe);
assert(numel(ppf)==numel(wl), 'ppf length (%d) must match #wavelengths (%d).', numel(ppf), numel(wl));
dcObj = hmrR_OD2Conc(dodObj, probe, ppf);

%% ================= BLOCK AVERAGE =================
trangeBA = [-2 20];
chIdx = 1;

figure; hold on; grid on

for c = 1:numel(condNames)

    s_single = s(:,c);

    [dataAvgBA, ~, ~, tHRF_BA] = ...
        BlockAvg_NoSnirf(dcObj, s_single, trangeBA);

    yBA = dataAvgBA(1).GetDataTimeSeries();

    plot(tHRF_BA, yBA(:,1,chIdx), 'LineWidth',2)

end

legend(condNames,'Interpreter','none')
xlabel('Time (s)')
ylabel('\DeltaHbO (\muM)')
title('BlockAvg HbO - all conditions')
%% ================= GLM =================

% ----- REQUIRED CELL FORMAT FOR HOMER3 -----
% Forza il measlist completo nel data object
ml_full = snirf.data(1).GetMeasList('reshape');
dcObj(1).SetMeasList(ml_full);
ml = dcObj(1).GetMeasList('reshape');   % NON SrcDetPairs

nMeas = size(ml,1);

mlActAuto = {ones(nMeas,1)};
tIncAuto  = {ones(length(t),1)};
tIncAuto  = {ones(length(t),1)};
rcMap     = {};
c_vector  = 0;
Aaux      = [];

[dcAvg, dcAvgStd, nTrialsGLM, dcNew, dcResid, dcSum2, beta, R, hmrstats] = ...
    hmrR_GLM(dcObj, snirf.stim, probe, ...
             mlActAuto, Aaux, tIncAuto, rcMap, ...
             trange, glmSolveMethod, idxBasis, paramsBasis, ...
             rhoSD_ssThresh, flagSSmethod, driftOrder, c_vector);

%% ================= SAVE =================
save(outMat, ...
    'snirfFile','keepStims','iqr_wavelet','hpf','lpf','ppf','trange', ...
    'condNames','nMeas', ...
    'dcAvg','dcAvgStd','beta','R','probe');

fprintf('Saved: %s\n', outMat);

%% ===================== LOCAL HELPERS =====================
function t = get_time_from_data(dataObj)
try
    t = dataObj.GetTime();
catch
    t = dataObj.time;
end
t = t(:);
end

function wl = get_probe_wavelengths(probe)
try
    wl = probe.GetWavelengths();
catch
    try
        wl = probe.wavelengths;
    catch
        wl = probe.lambda;
    end
end
wl = wl(:)';
end

function names = get_stim_names(stimArr)
if isempty(stimArr)
    names = string.empty(0,1); return
end
names = strings(numel(stimArr),1);
for i = 1:numel(stimArr)
    names(i) = get_stim_name_one(stimArr(i));
end
end

function name = get_stim_name_one(stim)
try
    name = string(stim.name);
catch
    try
        name = string(stim.GetName());
    catch
        name = "stim";
    end
end
end

function data = get_stim_data_one(stim)
try
    data = stim.data;   % [onset duration amplitude]
catch
    data = stim.GetData();
end
end

function s = stim2sMatrix(stimArr, t)
if isempty(stimArr)
    s = zeros(numel(t),0); return
end
nCond = numel(stimArr);
s = zeros(numel(t), nCond);
for c = 1:nCond
    stimData = get_stim_data_one(stimArr(c));
    if isempty(stimData), continue; end
    onsets = stimData(:,1);
    for k = 1:numel(onsets)
        [~, it] = min(abs(t - onsets(k)));
        s(it,c) = 1;
    end
end
end