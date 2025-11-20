function output = logspaceMatrix(start, stop, numSteps, OneStepValueStop)
%% Create matrix with logarithmically spaced (differing) columns
%
%   output = logspaceMatrix(start, stop, numSteps, OneStepValueStop)
%
% Each input argument must be either a (row) vector of the same size or a
% scalar.
%
% INPUT:
%
%   start
%       Row vector or scalar with the start values for each column in the result
%       matrix.
%       (Default: 1)
%
%   stop
%       Row vector or scalar with the stop values for each colums in the result
%       matrix.
%       (Default: 100)
%
%   numSteps
%       Row vector or scalar with the number of steps for each column in the
%       result matrix. If the number of steps is different for different
%       columns, the remaining elements are filled with NaN values.
%       If the number of steps is set to NaN, it is set to
%       floor(abs(log10(stop)-log10(start)+1)) for that column (i.e.,
%       approximately one value per order of magnitude).
%       (Default: NaN)
%
%   OneStepValueStop
%       Logical value that decides which limit is returned if the number of
%       steps for a columns is 1. If true, the stop value is used. If false, the
%       start value is used.
%       (Default: true)
%
%
% OUTPUT:
%
%   output
%       Matrix with the size maximum number of elements for each column times
%       the maximum number of elements in any of the input arguments. Each
%       column contains logarithmically spaced elements according to the input
%       arguments.
%
%
% EXAMPLE:
%
% >> output = logspaceMatrix([1,1,5,1,1], [100,0.001,60,2,2], [NaN,NaN,3,1,1], [0,0,0,0,1])
% output =
%     1.0000    1.0000    5.0000    1.0000    2.0000
%    10.0000    0.1000   17.3205       NaN       NaN
%   100.0000    0.0100   60.0000       NaN       NaN
%        NaN    0.0010       NaN       NaN       NaN
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if nargin < 1 || isempty(start)
  start = 1;
end
if nargin < 2 || isempty(stop)
  stop = 100;
end
if nargin < 3 || isempty(numSteps)
  numSteps = NaN;  % i.e., number of steps is floor(abs(log10(stop)-log10(start)+1))
end
if nargin < 4 || isempty(OneStepValueStop)
  % if numSteps is 1:
  % true => stop value used for output
  % false => start value used for output
  OneStepValueStop = true;
end


%% preprocess input (to have row vectors of consistent sizes)
% This block will fail if arguments have inconsistent sizes
start = log10(reshape(start, 1, []));
stop = log10(reshape(stop, 1, []));
numSteps = reshape(numSteps, 1, []);

range = stop - start;
stepSize = range ./ (numSteps-1);

start = start + zeros(size(stepSize));
stop = stop + zeros(size(stepSize));
numSteps = numSteps + zeros(size(stepSize));

% handle values for columns for which numSteps was set to NaN
numSteps(isnan(stepSize)) = floor(abs(stop(isnan(stepSize)) - start(isnan(stepSize))) + 1);
stepSize(isnan(stepSize)) = sign(range(isnan(stepSize)) ./ (numSteps(isnan(stepSize)) - 1));


%% create matrix with logspace columns
output = repmat((1:max(numSteps)).', 1, numel(stop));
output(bsxfun(@gt, output, numSteps)) = NaN;
output = output .* stepSize + start - stepSize;
OneStepValueStop = logical(OneStepValueStop) | false(1, size(output, 2));

if any(numSteps)
  output(1, OneStepValueStop & (numSteps==1))=  stop(1, OneStepValueStop & (numSteps==1));
  output(1,~OneStepValueStop & (numSteps==1))= start(1,~OneStepValueStop & (numSteps==1));
  output = 10.^output;
end


end
