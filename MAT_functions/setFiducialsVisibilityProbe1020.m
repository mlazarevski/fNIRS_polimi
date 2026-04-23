function [ nirsData ] = setFiducialsVisibilityProbe1020( nirsData,fidVisibility )

% This function sets fiducials visibility of nirs.core.Probe1020 data
% class. 
% The input (nirsData) can be one of the following data classes:
% - nirs.core.Data
% - nirs.core.ChannelStats
% - nirs.core.Probe1020
%
% Fiducials visibility (fidVisibility) can be set as:
% - true    : display 
% - false   : no display
% 

%%

if isa(nirsData, 'nirs.core.Data') | isa(nirsData, 'nirs.core.ChannelStats')
    
    for ind=1:length(nirsData)
        pp = nirsData(ind).probe;
        mesh = pp.getmesh;
        mesh(1).fiducials.Draw(:) = fidVisibility;
        pp = pp.register_mesh2probe(mesh);
        nirsData(ind).probe = pp;
    end
    
elseif isa(nirsData, 'nirs.core.Probe1020')
   
    mesh = nirsData.getmesh;
    mesh(1).fiducials.Draw(:) = fidVisibility;
    nirsData = nirsData.register_mesh2probe(mesh);
    
end

end
