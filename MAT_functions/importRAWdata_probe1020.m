function [ raw ] = importRAWdata_probe1020( folderData,folderHierarchy,probeVisualizationMethod )

% This function imports .nirs data throught nirs.io.loadDirectory method
% and computes a subject-specific nirs.core.Probe1020 data class from
% digpts.txt file.        
%
% Therefore, this function works under the hypothesis that each subject has 
% only one corresponding digpts.txt file within the subject-specific folder. 
% Morethan one run per subject is allowed. 
% An example of data directory structure is provided:
%     folderData
%         sbj001 
%             sbj001_data_run01.nirs
%             sbj001_data_run02.nirs
%             digpts.txt
%         sbj002
%             sbj002_data.nirs
%             digpts.txt
%         sbj003
%             sbj003_data.nirs
%             digpts.txt
%         ...
%
%
% INPUT:
% - folderData: pathanme of .nirs data folder 
% - folderHierarchy: cell array describing the folder hierarchy of
%               folderData; see nirs.io.loadDirectory for further details
% - probeVisualizationMethod: string defining probe visualization; see
%               "defaultdrawfcn" property in nirs.core.Probe1020 for
%               further details; default method is '10-20'
%
% OUTPUT:
% - raw: nirs.core.Data class containing raw data in folderData
%       
% NOTE:
% the definition of subject-specific nirs.core.Probe1020 data classes is
% based on the registration of digpts.txt data over the Colin27 mesh
% available in Brain AnalyzIR toolbox.
% 

%% import

% select .nirs data folder 
if not( exist('folderData', 'var') ) | isempty( folderData ) 
    folderData = uigetdir(pwd, 'Select .nirs data folder: ');
end
    
% set raw data folder hierarchy
if not( exist('folderHierarchy', 'var') ) | isempty( folderHierarchy ) 
    prompt = {'Set folder hierarchy (whitespace-separated entries): '};
    defAns = {'group subject'};

    folderHierarchy = inputdlg(prompt, '', [1 60], defAns);
end
        
% set probe visualization method
if not( exist('probeVisualizationMethod', 'var') ) | isempty( probeVisualizationMethod )
    probeVisualizationMethod = '10-20';
end

% retrieve .nirs data
pathnameNIRS = dir( fullfile(folderData, '**', '*.nirs') );

% retrieve digpts.txt data
% NOTE: each subject-specific folder must contain at least one digpts.txt
% file 
pathnameDIGPTS = pathnameNIRS;
for ind=1:length(pathnameNIRS)
    pathnameDIGPTS( ind ).name = 'digpts.txt';
end

%% compute

% import raw data class
raw = nirs.io.loadDirectory(folderData, folderHierarchy);

for ind=1:length(pathnameNIRS)
    fprintf('\n\nComputing probe1020: %u / %u\n\n', ind, length(pathnameNIRS));
    
    % register digpts.txt on Colin27 mesh in NIRS toolbox
    % NOTE: a new digpts_NIRSToolbox.txt file is created within the
    % corresponding folder
    % NOTE: registration is based on affine transform between Nz, LPA, RPA
    % and Cz points
    [ srcpos,detpos,fidpos,fidpos_label ] = registerDIGPTStoMesh_NIRSToolbox(pathnameDIGPTS(ind).folder, ...
        pathnameDIGPTS(ind).name);
    if isempty( srcpos )
        return;
    end
    
    % substitute probe to raw data class 
    SD = nirs.util.probe2sd( raw(ind).probe );

    SD = rmfield(SD, 'MeasListAct');
    SD = rmfield(SD, 'Description');

    SD.SrcPos = srcpos;
    SD.DetPos = detpos;

    SD.SrcPos3D = srcpos;
    SD.DetPos3D = detpos;

    SD.Landmarks.pos = fidpos;
    SD.Landmarks.labels = fidpos_label;
    SD.Landmarks2D.pos = fidpos;
    SD.Landmarks2D.labels = fidpos_label;
    SD.Landmarks3D.pos = fidpos;
    SD.Landmarks3D.labels = fidpos_label;

    probe = nirs.util.sd2probe( SD );
    probe.defaultdrawfcn = probeVisualizationMethod;
    
    raw(ind).probe = probe;

end

end




function [ srcpos,detpos,fidpos,fidpos_label ] = registerDIGPTStoMesh_NIRSToolbox( digptsDIR,digptsNAME )

% internal fcn 
% 

%% 

% import digpts.txt data
pathnameDIGPTS = fullfile(digptsDIR, digptsNAME);

try
    [ srcpos,detpos,fidpos,fidpos_label ] = openDigptsTXT( pathnameDIGPTS );
    fidpos_label = strrep(fidpos_label, ':', '');

catch
    fprintf('\n\nMissing digpts.txt file at %s.\nComputation is stopped.\n\n', digptsDIR);
    srcpos = [];
    detpos = [];
    fidpos = [];
    fidpos_label = [];
    return;

end

% register digpts.txt on Colin27 mesh in NIRS toolbox
% NOTE: registration is based on affine transform between Nz, LPA, RPA
% and Cz points
mesh = nirs.registration.Colin27.mesh;
fidlst = mesh(1).fiducials;

if not( isempty(find( strcmp(fidpos_label, {'Nz'}))) ) & ...
        not( isempty(find( strcmp(fidpos_label, {'LPA'}))) ) & ...
        not( isempty(find( strcmp(fidpos_label, {'RPA'}))) ) & ...
        not( isempty(find( strcmp(fidpos_label, {'Cz'}))) )

    tmp = fidlst.Name;
    tmp = strrep(tmp, 'spmnas', 'Nz');
    tmp = strrep(tmp, 'spmlpa', 'LPA');
    tmp = strrep(tmp, 'spmrpa', 'RPA');
    fidlst.Name = tmp;
    
elseif not( isempty(find( strcmp(fidpos_label, {'nz'}))) ) & ...
        not( isempty(find( strcmp(fidpos_label, {'al'}))) ) & ...
        not( isempty(find( strcmp(fidpos_label, {'ar'}))) ) & ...
        not( isempty(find( strcmp(fidpos_label, {'cz'}))) )
    
    tmp = fidlst.Name;
    tmp = strrep(tmp, 'spmnas', 'nz');
    tmp = strrep(tmp, 'spmlpa', 'al');
    tmp = strrep(tmp, 'spmrpa', 'ar');
    tmp = strrep(tmp, 'Cz', 'cz');
    tmp = strrep(tmp, 'Iz', 'iz');
    fidlst.Name = tmp;
    
end

fidpos_mesh = [];
for ind=1:length(fidpos_label)
    pos = find( strcmp(fidlst.Name, fidpos_label{ind}) );
    fidpos_mesh = [ fidpos_mesh; fidlst.X(pos) fidlst.Y(pos) fidlst.Z(pos) ];
end

% %         subplot(121)
% %         plot3DdataMatrix(fidpos_mesh, 'g')
% %         hold on
% %         plot3DdataMatrix(fidpos, 'm')
% %         plot3DdataMatrix(srcpos, 'r')
% %         plot3DdataMatrix(detpos, 'b')
        
[ T,fidpos ] = computeAffineTransformFromPointcloud(fidpos, fidpos_mesh);
srcpos = xform_apply(srcpos, T);
detpos = xform_apply(detpos, T);

% %         subplot(122)
% %         plot3DdataMatrix(fidpos_mesh, 'g')
% %         hold on
% %         plot3DdataMatrix(fidpos, 'm')
% %         plot3DdataMatrix(srcpos, 'r')
% %         plot3DdataMatrix(detpos, 'b')
      
end