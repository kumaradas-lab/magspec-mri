function [phaseDiffWeightedMean, phaseDiffWeightedStd] = get_MeanPhaseDiffWeighted(complexData, dim, nanflag)
%% Weighted mean phase difference
%
%   [phaseDiffWeightedMean, phaseDiffWeightedStd] = get_MeanPhaseDiffWeighted(complexData, dim, nanflag)
%
% Calculate the mean "phaseDiffWeightedMean" and standard error
% "phaseDiffWeightedStd" of the phase difference weighted by the squared
% absolute amplitude of "ComplexData" in direction "dim".
% If nanflag is 'includenan', the mean and standard error evaluate to NaN if the
% input includes NaN values in the specified dimension (default). If it is set to
% 'omitnan', NaN values are ignored.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2018 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if nargin < 2, dim = find(size(complexData)>1, 1, 'first'); end
if nargin < 3, nanflag = 'includenan'; end

nsize = cumsum(size(complexData)>=1);
nsize(1) = nsize(dim);
nsize(dim) = 1;
complexData = permute(complexData, nsize);
weights = min(min(abs(complexData(1:end-1,:,:,:,:,:)), abs(complexData(2:end,:,:,:,:,:))), ...
  abs((complexData(1:end-1,:,:,:,:,:) + complexData(2:end,:,:,:,:,:))/2)).^2;
complexData = permute(complexData, nsize); % FIXME: Should this be ipermute?
weights = permute(weights, nsize); % FIXME: Should this be ipermute?
if (dim==1) && isvector(complexData) && all(isfinite(complexData(:)))
  PhaseDiff = diff(unwrapVector(angle(complexData), pi), 1, dim);
else
  PhaseDiff = diff(unwrap(angle(complexData), [], dim), 1, dim);
end
if strcmpi(nanflag, 'omitnan')
  PhaseDiff(isnan(PhaseDiff)) = 0;
  weights(isnan(weights)) = 0;
end
PhaseDiffWeighted = PhaseDiff .* weights;
phaseDiffWeightedMean = sum(PhaseDiffWeighted,dim) ./ sum(weights, dim);

if nargout > 1
  phaseDiffWeightedStd = sqrt(sum(bsxfun(@minus, PhaseDiff, phaseDiffWeightedMean).^2 .* weights) ./ ...
                              sum(weights, dim));
end

end
