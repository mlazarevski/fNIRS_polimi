function [ nirsData ] = correctMotorImagineryStimuli( nirsData )

% This function corrects stimuli timing and name of data from motor
% imaginery task through NIRSToolbox.
% Moreover, it add the post-stimulus baseline to each left or right
% tasks.
% 
% It is suggested to run this function on raw data classes before any other
% analyses. 
% 
% 
%       % % % % % % % % % % % %
%       % STIMULI DESCRIPTION %
%       % % % % % % % % % % % %
% 
% stim_channel1 --> LEFT finger tapping task's onset
% stim_channel2 --> RIGHT finger tapping task's onset
%
% stim_channel4 --> LEFT motor imaginery task's onset
% stim_channel8 --> RIGHT motor imaginery task's onset
%

%% stimuli

% remove auxillary measurements as stimuli
job = nirs.modules.DiscardStims();
job.listOfStims = {'stim_aux1' 'stim_aux2' 'stim_aux3' ...
    'stim_aux4' 'stim_aux5' 'stim_aux6'};
nirsData = job.run( nirsData );

% rename task and change task duration
job = nirs.modules.RenameStims();
job.listOfChanges = {'stim_channel1' 'left'; ...
    'stim_channel2' 'right'; ...
    'stim_channel4' 'leftImagine'; ...
    'stim_channel8' 'rightImagine'};
nirsData = job.run( nirsData );

nirsData = nirs.design.change_stimulus_duration(nirsData, {'left'}, 10);
nirsData = nirs.design.change_stimulus_duration(nirsData, {'right'}, 10);
nirsData = nirs.design.change_stimulus_duration(nirsData, {'leftImagine'}, 10);
nirsData = nirs.design.change_stimulus_duration(nirsData, {'rightImagine'}, 10);

%% stimuli baseline 

% add post-task resting as stimulus
for ind=1:length(nirsData)
    tt = nirsData( ind ).time;
    stim1 = nirsData( ind ).stimulus('left');
    stim2 = nirsData( ind ).stimulus('right');
    stim3 = nirsData( ind ).stimulus('leftImagine');
    stim4 = nirsData( ind ).stimulus('rightImagine');

    [ ~,idx1 ] = findVectorClosestValue(tt, stim1.onset + 10);
    [ ~,idx2 ] = findVectorClosestValue(tt, stim2.onset + 10);
    [ ~,idx3 ] = findVectorClosestValue(tt, stim3.onset + 10);
    [ ~,idx4 ] = findVectorClosestValue(tt, stim4.onset + 10);

    stim1.onset = tt( idx1 );
    stim2.onset = tt( idx2 );
    stim3.onset = tt( idx3 );
    stim4.onset = tt( idx4 );
    
    stim1.name = 'leftBaseline';
    stim2.name = 'rightBaseline';
    stim3.name = 'leftImagineBaseline';
    stim4.name = 'rightImagineBaseline';
    
    nirsData( ind ).stimulus('leftBaseline') = stim1;
    nirsData( ind ).stimulus('rightBaseline') = stim2;
    nirsData( ind ).stimulus('leftImagineBaseline') = stim3;
    nirsData( ind ).stimulus('rightImagineBaseline') = stim4;
end

nirsData = nirs.design.change_stimulus_duration(nirsData, {'leftBaseline'}, 15);
nirsData = nirs.design.change_stimulus_duration(nirsData, {'rightBaseline'}, 15);
nirsData = nirs.design.change_stimulus_duration(nirsData, {'leftImagineBaseline'}, 15);
nirsData = nirs.design.change_stimulus_duration(nirsData, {'rightImagineBaseline'}, 15);

end