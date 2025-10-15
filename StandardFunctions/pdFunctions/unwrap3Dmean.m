function unwrapped = unwrap3D(toUnwrap, tol)
%% Unwrap 3-dimensional angle data (in radians)
%
%   unwrapped = unwrap3D(toUnwrap, tol)
%
% INPUT:
%   toUnwrap    Phase angles in radians
%   tol         Jump tolerance passed to Matlab's "unwrap"
%
% OUTPUT:
%   unwrapped   Unwrapped data where multiples of 2pi are added where jumps are
%               larger than pi (see: unwrap). The data is unwrapped starting
%               with the largest dimension. Unwrapping along subsequent
%               dimensions is done such that the average jump of unwrapped
%               dimensions is minimal. This works best for continuous volumes.
%
% ------------------------------------------------------------------------
% (C) Copyright 2011-2017 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------

%% input check
if nargin < 2
  tol = [];
end

[~, sortOrder] = sort(size(toUnwrap));

if numel(sortOrder) ~= 3
  error('unwrap3d: Input argument "toUnwrap" must have 3 non-singular dimensions.');
end

%% unwrap data
% first unwrap largest dimension
uw1 = unwrap(toUnwrap, tol, sortOrder(3));
% unwrap second largest dimension
uw2a = mean(uw1, sortOrder(3), 'omitnan');
uw2 = bsxfun(@minus, uw1, uw2a-unwrap(uw2a, tol, sortOrder(2)));
% unwrap shortest dimension
uw3a = mean(mean(uw2, sortOrder(3), 'omitnan'), sortOrder(2), 'omitnan');
unwrapped = bsxfun(@minus, uw2, uw3a-unwrap(uw3a, tol, sortOrder(1)));

end
