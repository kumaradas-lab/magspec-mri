function w = Hamming(L, sflag)
%% Create Hamming filter window
%
%   w = Hamming(L, sflag)
%
% INPUT:
%   L       Sample points (in the interval [-1;1]) if `sflag` is set to
%           'sampled'; number of samples otherwise.
%   sflag   'symmetric', 'periodic', or 'sampled'. (default: 'symmetric')
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if nargin == 1
  sflag = 'symmetric';
end

% parameters for Hamming filter
a = 25/46;
b = 21/46;

if strcmp(sflag, 'sampled')
  % L is between -1 and 1
  w = + a ...
      - b .* cos(pi .* (L+1));
else
  % L is number of equally spaced samples
  L = double(L);
  if strcmp(sflag, 'periodic')
    L = L+1;
  end
  wp = (0:L-1).';
  w = + a ...
      - b .* cos(2.*pi.*wp./(L-1));
  if strcmp(sflag, 'periodic')
    w = w(1:end-1);
  end
end

if isscalar(w)
  w = 1;
end

end
