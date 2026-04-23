% BASED UPON PREVIOUS BANDPASS FUNCTION
%
% y2 = hmrHighPassFiltConc( y, fs, hpf, lpf )
% Perform a high pass filtering by additional separation of [HbO] and [HbR]
% concentration types
%
% INPUT:
% y - data to filter #time points x #channels of data
% fs - sample frequency (Hz). If length(fs)>1 then this is assumed to be a time
%      vector from which fs is estimated
% hpf - high pass filter frequency (Hz)
%       Typical value is 0 to 0.02.
%
% OUTPUT:
% y2 - filtered data

function [y2] = hmrHighPassFiltConc( y, fs, SD, hpf )

ytemp = y;
y2temp = zeros(size(y));

% retrieve [HbO] and [HbR] data
sdList = SD.MeasList;
iO = find( sdList(:,4) == 1 );
iR = find( sdList(:,4) == 2 );

% convert t to fs
% assume fs is a time vector if length>1
if length(fs)>1
    fs = 1/(fs(2)-fs(1));
end

% % %
% [HbO]
% % % 

y = ytemp(:,iO);

% high pass filter
FilterType = 1;
FilterOrder = 5;
if FilterType==1 | FilterType==5
    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,hpf,'high');
elseif FilterType==4
%    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,hpf,'high',Filter_Rp,Filter_Rs);
else
%    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,hpf,'high',Filter_Rp);
end

if FilterType~=5
    y2=filtfilt(fb,fa,y); 
else
    y2 = y;
end

y2temp(:,iO) = y2;

% % %
% [HbR]
% % %

y = ytemp(:,iR);

% high pass filter
FilterType = 1;
FilterOrder = 5;
if FilterType==1 | FilterType==5
    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,hpf,'high');
elseif FilterType==4
%    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,hpf,'high',Filter_Rp,Filter_Rs);
else
%    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,hpf,'high',Filter_Rp);
end

if FilterType~=5
    y2=filtfilt(fb,fa,y); 
else
    y2 = y;
end

y2temp(:,iR) = y2;
clear y2
y2 = y2temp;

end