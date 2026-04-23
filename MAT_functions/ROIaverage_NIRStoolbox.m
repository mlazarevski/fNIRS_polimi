function [ nirsDataAvg ] = ROIaverage_NIRStoolbox( nirsData,ROI,ROIname )

% 23/07/2024
% This fcn performs ROI averaging on fNIRS signals (i.e., either raw, dod or 
% hb data classes) using NIRSToolbox. 
% However, I noticed that ROI averaged signals according to toolbox functions 
% are divided by the number of channels in the ROI itself, hence providing a 
% different scaling compared to perform this operation analytically. 
% Therefore, this function corrects for this scaling on ROI averaged
% signals. 
% 
% This function is also valid for block averaged data. 
% 
% This function requires that input variable has at least two wavelengths
% (i.e., raw and dod data) or two chromophores (i.e., hb data). Otherwise, 
% computation is stopped. 
%
% This function is not valid for nirs.core.ChannelStats data class or any
% other statistical data classes. 
% 

%% 

% compute ROI averaging
nirsDataAvg = nirs.util.roiAverage(nirsData, ROI, ROIname);

% check number of wavelengths / chromophores
if length( unique(nirsData(1).probe.link.type) ) < 2
    fprintf('\n\n\t\t*************\n');
    fprintf('\t\t****Error****\n');
    fprintf('\t\t*************\n\n');
    fprintf('The number of wavelengths or chromophores is not sufficient. Computation is stopped.\n\n');
    return;
end

% correct scaling according to number of channels per ROI
for ind=1:length(nirsDataAvg)
    fprintf('\n\nComputing ROI average: %u / %u\n\n', ind, length(nirsDataAvg));
    val = nirsDataAvg( ind ).data;
    
    count = 1;
    for indr=1:length(ROI)
        nch = height( ROI{indr} );
        
        % scale wavelength / chromophore #1
        val( :,count ) = nch .* val( :,count ); 
        
        % scale wavelength / chromophore #2
        val( :,count+1 ) = nch .* val( :,count+1 );
        
        % update index
        count = count + 2;
    end
    
    % update data
    nirsDataAvg( ind ).data = val;
end

end