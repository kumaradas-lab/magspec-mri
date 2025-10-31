function w = Hamming(L, sflag)
%% Create Hamming filter window
%
%   w = Hamming(L, sflag)
%
% INPUT:
%   L       number of samples
%   sflag   'symmetric' or 'periodic' (default: 'symmetric')
%
%
% ------------------------------------------------------------------------
% (C) Copyright 2018 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------

if nargin == 1
    sflag = 'symmetric';
end   

L = double(L);
if strcmp(sflag, 'periodic'), L = L+1; end;
wp = (0:L-1).';
a = 25/46;
b = 21/46;
w=  + a ...
    - b .* cos(2.*pi.*wp./(L-1));
if strcmp(sflag, 'periodic'), w = w(1:end-1); end;
if numel(w) == 1, w = 1; end

end
