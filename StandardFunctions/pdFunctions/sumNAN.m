function s = sumNAN(A, dim)
% do a sum over a matrix even if it is including nan.

warning('PD:sumNAN:deprecated', ...
  ['The function "sumNAN" is deprecated and will be removed in a future version of OpenMatlab.\n', ...
  'Use "sum(..., ''omitnan'')" instead.']');

[s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9),s(10)]=size(A);
if nargin == 1
  dim = find(s>1, 1, 'first');
  if isempty(dim)
    dim = 1;
  end
end
A(isnan(A)) = 0;
s = sum(A, dim);


%% -----------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
