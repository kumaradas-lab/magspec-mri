function [phaseDiffWeightedMean, phaseDiffWeightedStd] = get_MeanPhaseDiffWeighted(complexData, dim, nanflag)
%% Weighted mean phase difference
%
%   [phaseDiffWeightedMean, phaseDiffWeightedStd] = ...
%     get_MeanPhaseDiffWeighted(complexData, dim, nanflag)
%
% Calculate the mean slope "phaseDiffWeightedMean" and standard error
% "phaseDiffWeightedStd" of the phase weighted by the squared absolute amplitude
% of "ComplexData" in direction "dim".
% This is done by a weighted linear regression of the phase of complexData.
%
% If nanflag is 'includenan', the mean and standard error evaluate to NaN if the
% input includes NaN values in the specified dimension (default). If it is set
% to 'omitnan', (trailing) NaN values are ignored.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% FIXME: Check if "non-trailing" NaN values are handled correctly.

if nargin < 2, dim = find(size(complexData)>1, 1, 'first'); end
if nargin < 3, nanflag = 'includenan'; end

nsize = cumsum(size(complexData)>=1);
nsize(1) = nsize(dim);
nsize(dim) = 1;
complexData = permute(complexData, nsize);
if (dim==1) && isvector(complexData) && all(isfinite(complexData(:)))
  phaseUnwrap = unwrapVector(angle(complexData), pi);
else
  phaseUnwrap = unwrap(angle(complexData), [], 1);
end
absDataSq = abs(complexData).^2;

[nSamples, nAQs] = size(complexData);

permSize = size(phaseUnwrap);
permSize(1) = 1;
phaseDiffWeightedMean = zeros(permSize);
if nargout > 1
  phaseDiffWeightedStd = zeros(permSize);
end

omitNaN = strcmpi(nanflag, 'omitnan');

x = [ones(1,nSamples)*nSamples; 1:nSamples].';
% FIXME: Can we avoid this loop?
for iData = 1:nAQs
  if omitNaN
    iUsed = ~isnan(phaseUnwrap(:,iData));
    numUsed = find(iUsed, 1, 'last');
    iUsed = iUsed(1:numUsed);
    x = [ones(1,numUsed)*numUsed; 1:numUsed].';
  else
    iUsed = true(nSamples,1);
  end

  if ~any(iUsed)
    % all amplitudes are NaN
    phaseDiffWeightedMean(1,iData) = NaN;
  else
    % use signal amplitude as weight in fit
    % b = x(iUsed,:) \ phaseUnwrap(iUsed,iData);  % unweighted fit
    weights = absDataSq(iUsed,iData) / max(absDataSq(iUsed,iData));
    if all(isnan(weights))
      % occurs when all absDataSq are zero (0/0)
      phaseDiffWeightedMean(1,iData) = NaN;
      b = [1; 0];  % dummy value to avoid error if standard error is requested
    else
      b = bsxfun(@times, weights, x(iUsed,:)) \ (weights .* phaseUnwrap(iUsed,iData));
      phaseDiffWeightedMean(1,iData) = b(2);  % only take slope, offset doesn't matter
    end
  end

  if nargout > 1
    if ~any(iUsed)
      % all amplitudes are NaN
      phaseDiffWeightedStd(1,iData) = NaN;
    else
      % See: https://de.wikipedia.org/wiki/Standardfehler_der_Regression
      % sigma = sqrt(sum(((x(iUsed,:)*b-phaseUnwrap(iUsed,iData)).^2)) / ...
      %         (nSamples-2));
      % together with:
      % https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Statistical_properties
      sigma = sqrt(1/(sum(iUsed)-2) * ...
        sum(((phaseUnwrap(iUsed,iData) - x(iUsed,:)*b) .* (weights*sum(iUsed)/sum(weights))).^2));
      phaseDiffWeightedStd(1,iData) = ...
        sigma / sqrt(mean((x(iUsed,2) - mean(x(iUsed,2))).^2));
    end
  end
end

phaseDiffWeightedMean = ipermute(phaseDiffWeightedMean, nsize);

if nargout > 1
  phaseDiffWeightedStd = ipermute(phaseDiffWeightedStd, nsize);
end

end
