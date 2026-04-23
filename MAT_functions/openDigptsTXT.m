function [ srcpos,detpos,fidpos,fidpos_label ] = openDigptsTXT( pathDigpts )

% [ self-explanatory ]
%

%%

% import configuration
if nargin < 1
    
    [ filename,pathname ] = uigetfile('digpts*.txt', 'import configuration');
    pathDigpts = fullfile(pathname, filename);

end

fid = fopen( pathDigpts );
c = textscan(fid, '%s %f %f %f');
fclose( fid );
optLabel = c{1};
xyzIN = [ c{2} c{3} c{4} ];

% source - detector - refpts 
poss = strfind(optLabel, 's');
posd = strfind(optLabel, 'd');
poss = find( cell2mat( cellfun(@(x) not(isempty(x)), poss, 'UniformOutput', false) ) );
posd = find( cell2mat( cellfun(@(x) not(isempty(x)), posd, 'UniformOutput', false) ) );
posr = find( not( ismember([1:length(optLabel)]', [poss; posd]) ) );

fidpos = xyzIN( posr,: );
srcpos = xyzIN( poss,: );
detpos = xyzIN( posd,: );

fidpos_label = optLabel( posr );

end
