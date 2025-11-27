function ComplexDataAnglewMean = get_MeanPhaseWeighted(ComplexData, dim)
%% Calculate mean phase of complex signal weighted with the signal amplitude
%
%   ComplexDataAnglewMean = get_MeanPhaseWeighted(ComplexData, dim)
%
%
% INPUT:
%
%   ComplexData
%       Array with arbitrary dimensions with complex valued data.
%
%   dim
%       Dimension along which the average weighted phase should be calculated.
%       (Default: the first non-singleton dimension of ComplexData)
%
%
% OUTPUT:
%
%   ComplexDataAnglewMean
%       Array with the same dimensions as the input ComplexData (apart from the
%       dimension dim) which contains the weighted average phase in radians
%       along dimension dim. The magnitude of the input ComplexData is used as
%       the weight.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if nargin < 2
  dim = find(size(ComplexData)>1, 1, 'first');
end


%% calculate weighted average phase along dim
ComplexDataw = abs(ComplexData);  % magnitude used as weight

% phase of input (unwrapped along dim)
ComplexDataAngle = unwrap(angle(ComplexData), [], dim);

% weighted average
ComplexDataAnglew = ComplexDataAngle .* ComplexDataw;
ComplexDataAnglewMean = sum(ComplexDataAnglew, dim) ./ sum(ComplexDataw, dim);

ComplexDataAnglewMean(all(ComplexDataw==0, dim)) = 0;  % avoid 0/0

end
