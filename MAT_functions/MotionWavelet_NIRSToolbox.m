function [ dod ] = MotionWavelet_NIRSToolbox( dod,iqr,strNIRSToolbox )

% NOTE:
% NIRSToolbox needs to be removed from current path when
% "hmrMotionCorrectWavelet" is used (solve for inconsistency).
%
% NOTE:
% depending on the dataset, this operation takes a lot of time. A waitbar
% object is displayed.
%

if not( exist('strNIRSToolbox', 'var') ) | isempty( strNIRSToolbox )
    load strNIRSToolbox.mat;
end

f = waitbar(0, sprintf('MotionWavelet:  %u / %u', 0, length(dod)));

for ind=1:length(dod)
    val = dod(ind).data;
    SD = nirs.util.probe2sd( dod(ind).probe );

    waitbar(ind/length(dod), f, sprintf('MotionWavelet:  %u / %u', ind, length(dod)));
    
    rmpath( genpath(strNIRSToolbox) )
    val = hmrMotionCorrectWavelet(val, SD, iqr);
    addpath( genpath(strNIRSToolbox) )

    dod(ind).data = val;
end

close(f);

end

