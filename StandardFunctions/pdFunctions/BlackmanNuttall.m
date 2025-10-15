function w = BlackmanNuttall(L, sflag)
%% Create Blackman-Nuttall filter window
%
%   w = BlackmanNuttall(L, sflag)
%
% The Blackman-Nuttall filter window is -- apart from the four almost identical
% coefficients -- identical to the Blackman-Harris filter windows. That
% illustrates the necessary accuracy for the implementatino of these
% coefficients when it comes to this class of window functions.
%
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
a0 = 0.3635819;
a1 = 0.4891775;
a2 = 0.1365995;
a3 = 0.0106411;
w = ...
  + a0 ...
  - a1 .* cos(2.*pi.*wp./(L-1)) ...
  + a2 .* cos(4.*pi.*wp./(L-1)) ...
  - a3 .* cos(6.*pi.*wp./(L-1));


if strcmp(sflag, 'periodic')
  w = w(1:end-1);
end
if isscalar(w)
  w = 1;
end

end
