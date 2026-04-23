function [ hb ] = MBLL_NIRSToolbox_Homer2( dod,ppf )

% This function applies MBLL to a nirs.core.Data class dataset using Homer
% functions. 
% Output "hb" is a nirs.core.data class dataset as well.
% 
% If no ppf value is specified, default value is set to [0.1 0.1] (i.e.,
% default value of nirs.modules.BeerLambertLaw method). 
%
% NOTE: this fcn is valid only for dual-wavelength data. 
% 

%%

% set ppf value
if not( exist('ppf', 'var') ) | isempty( ppf )
    ppf = [ 0.1 0.1 ];
end
    
% compute output probe 
job = nirs.modules.BeerLambertLaw();
tmp = job.run( dod(1) );
pp = tmp.probe;

%%

hb = dod;

for ind=1:length(dod)
    % compute homer-compatible probe and data 
    SD = nirs.util.probe2sd( dod( ind ).probe );
    val = dod( ind ).data;
    
    % apply homer MBLL 
    dc = hmrOD2Conc( val, SD, ppf );

    % arrange output data matrix 
    dcnew = [];
    
    for indc=1:size(dc,3)
        dcnew = [ dcnew dc(:,1,indc) dc(:,2,indc) ];
    end
    
    % update 
    hb( ind ).data = dcnew;
    hb( ind ).probe = pp;
end

end
