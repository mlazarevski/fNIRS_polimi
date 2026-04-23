clear all, close all, clc 
addpath(genpath('C:\Users\catea\OneDrive - Politecnico di Milano\Corsi\METHODS FOR THE FUNCTIONAL EVALUATION OF THE CENTRAL NERVOUS SYSTEM\2025-2026\Codici\functions'));
addpath(genpath('C:\Users\catea\OneDrive - Politecnico di Milano\Corsi\METHODS FOR THE FUNCTIONAL EVALUATION OF THE CENTRAL NERVOUS SYSTEM\2025-2026\Codici\nirs-toolbox-master'))
HOMER3_DIR = 'C:\Users\catea\OneDrive - Politecnico di Milano\Corsi\METHODS FOR THE FUNCTIONAL EVALUATION OF THE CENTRAL NERVOUS SYSTEM\2024-2025\homer2';
addpath(genpath(HOMER3_DIR));
%% data import

raw = nirs.io.loadDirectory('nirsdataC', {'group' 'subject'});   
demotbl = nirs.createDemographicsTable( raw )

figure;
raw(1).draw;
title('intensity data');
drawnow;

figure; 
raw(1).probe.draw;
drawnow;

%% stimuli correction

% stim_channel1 --> 1-BACK          = stim1
% stim_channel2 --> 2-BACK          = stim2
% stim_channel4 --> PUSH FOR CENTER = stim0
 
job = nirs.modules.DiscardStims();                   
job.listOfStims = {'stim_aux1' 'stim_aux2' 'stim_channel128', 'stim_channel8'};
raw = job.run( raw );

job = nirs.modules.RenameStims();              
job.listOfChanges = {'stim_channel1' 'stim1'; 'stim_channel2' 'stim2'; 'stim_channel4' 'stim0'};
raw = job.run( raw );

raw = nirs.design.change_stimulus_duration(raw, {'stim1'}, 28);
raw = nirs.design.change_stimulus_duration(raw, {'stim2'}, 28);
raw = nirs.design.change_stimulus_duration(raw, {'stim0'}, 28);

figure;
raw(1).draw;
title('intensity data');
drawnow;

%% prune channels

tbl = nirs.util.scalp_coupling_index( raw(1) );

%% pre-processing

job = nirs.modules.Resample();
raw = job.run( raw );

job = nirs.modules.OpticalDensity(); 
dod = job.run( raw );

job = nirs.modules.TDDR(); 
dod = job.run( dod );

job = nirs.modules.RunMatlabCode();            
job.FunctionHandle = @(x) BPF_NIRSToolbox(x, 0.01, 0.5);
job.rungroup = 1;
dod = job.run( dod );

job = nirs.modules.RunMatlabCode();             
job.FunctionHandle = @(x) PCA_NIRSToolbox(x, 1);
dod = job.run( dod );

figure;
dod(1).draw;
title('optical density data');
drawnow;

%% Modified Beer-Lambert law

job = nirs.modules.BeerLambertLaw();       
job.PPF = DPFcorrection(dod(1).probe.types, 55)';  
hb = job.run( dod );

figure;
hb(1).draw;
title('hb data');
drawnow;

%% GLM

job = nirs.modules.LabelShortSeperation();     
hb = job.run( hb );

% GLM function
%         Barker, Jeffrey W., Ardalan Aarabi, and Theodore J. Huppert."Autoregressive model based algorithm fo
%         r correcting motion and serially correlated errors in fNIRS." Biomedical optics express 4.8 (2013): 
%         1366-1379.
job = nirs.modules.GLM();                
job.verbose = 1;
job.AddShortSepRegressors = 1;  

jb = nirs.design.basis.Canonical();
jb.peakTime = 6;
job.basis('default') = jb;

SubjStats = job.run( hb );

%% GLM results

SubjStats.draw('tstat', [], 'p < 0.05');
drawnow;


%% GLM contrast
close all;

C = {'stim1-stim0', 'stim2-stim1', 'stim2-stim0'}               
cc = SubjStats.ttest( C );                  
      
cc.draw('tstat', [], 'p < 0.05');
drawnow;

%% summary results according to ROI
close all;

load ROI_cognitive_configuration.mat;

[ ccTBL,ccFIGURE ] = nirs.util.roiAverage(cc, ROI, ROIname);
disp( ccTBL );                          % 


tmp = ccFIGURE.draw('tstat');
% tmp = ccFIGURE.draw('beta');

figure(1);
% set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(gca, 'XTickLabelRotation', 45);
title('HbO contrast');
drawnow;

figure(2);
% set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(gca, 'XTickLabelRotation', 45);
drawnow;