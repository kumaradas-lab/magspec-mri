function y = sinc(x, a)
%% sinc function
%
%   y = sinc(x, a)
%
% INPUT:
%   x     numeric array with x values
%   a     scalar defining the location of the first zero-crossing (default: 1)
%
% OUTPUT:
%   y     numeric array with the values of
%           sin(pi*x/a) ./ (pi*x/a)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if nargin < 2, a = 1; end

iZeros = find(abs(x) <= eps(a));
x(iZeros) = 1;  % avoid division by zero
y = sin(pi*x/a) ./ (pi*x/a);
y(iZeros) = 1;  % limit to 0 is 1
y(abs(y) <= eps(1)) = 0;  % numeric deviations

end
