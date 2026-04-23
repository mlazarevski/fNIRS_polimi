function [ DPF ] = DPFcorrection( wl,age )

% This script computes the age- and wavelength-specific DPF according to 
% [Scholkmann et al., 2013].
%
% 23/07/2021
% This fcn is able to cope with multiple wavelengths at a time.
% 

%%

% NOTE: [Scholkmann et al., 2013] reports that interpolated age- and
% wavelength-dependent DPF values are:
% - referred to frontal cortex --> in our case, we just need an indicative
%    value, hence it is fine
% - not completely reliable for subjects whos age is above 50. Therefore,
%    as a simplification, the DPF subjects who are older than 50yo are 
%    automatically clipped to 50yo-case itself
if age > 50
    age = 50;
end

alpha = 223.3;
beta = 0.05624;
gamma = 0.8493;
delta = -5.723*( 10^-7 );
epsilon = 0.001245;
zeta = -0.9025;

DPF = alpha + beta*( age^gamma ) + delta*( wl.^3 ) + epsilon*( wl.^2 ) + zeta.*wl; 

end

