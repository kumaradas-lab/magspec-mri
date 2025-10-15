function w = Blackman(L, sflag)
%% Create Blackman filter window
%
%   w = Blackman(L, sflag)
%
% INPUT:
%   L       number of samples
%   sflag   'symmetric' or 'periodic' (default: 'symmetric')
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2019-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if nargin == 1
  sflag = 'symmetric';
end

L = double(L);
if strcmp(sflag, 'periodic')
  L = L+1;
end
wp = (0:L-1).';
a0 = 7938/18608;
a1 = 9240/18608;
a2 = 1430/18608;
w = ...
  + a0 ...
  - a1 .* cos(2.*pi.*wp./(L-1)) ...
  + a2 .* cos(4.*pi.*wp./(L-1));

if strcmp(sflag, 'periodic')
  w = w(1:end-1);
end
if isscalar(w)
  w = 1;
end

end
