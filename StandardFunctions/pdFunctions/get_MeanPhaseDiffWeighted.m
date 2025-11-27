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
% to 'omitnan', (trailing!) NaN values are ignored.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if nargin < 2, dim = find(size(complexData)>1, 1, 'first'); end
if nargin < 3, nanflag = 'includenan'; end

omitNaN = strcmpi(nanflag, 'omitnan');

% FIXME: Check if "non-trailing" NaN values are handled correctly.


%% prepare for linear fit
% permute DIM to first dimension
nsize = cumsum(size(complexData)>=1);
nsize(1) = nsize(dim);
nsize(dim) = 1;
complexData = permute(complexData, nsize);
[nSamples, nAQs] = size(complexData);

% un-wrap phase
if (dim==1) && isvector(complexData) && all(isfinite(complexData(:)))
  phaseUnwrap = unwrapVector(angle(complexData), pi);
else
  phaseUnwrap = unwrap(angle(complexData), [], 1);
end

% normalize data
meanAbsData = mean(abs(complexData(:)), 'omitnan');
absDataSq = abs(complexData / meanAbsData).^2;


%% perform linear fit
% initialize vectors for results
permSize = size(phaseUnwrap);
permSize(1) = 1;
phaseDiffWeightedMean = zeros(permSize);
if nargout > 1
  phaseDiffWeightedStd = zeros(permSize);
end

% prepare equation for linear fit
x = [ones(1,nSamples)*nSamples; 1:nSamples].' / nSamples;

if ~omitNaN
  % values are the same for each acquisition window
  numUsed = nSamples;
  iUsed = true(nSamples, 1);
end

% FIXME: Can we avoid this loop?
for iData = 1:nAQs
  if omitNaN
    % update equation for linear fit for each acquisition window
    iUsed = ~isnan(phaseUnwrap(:,iData));
    numUsed = find(iUsed, 1, 'last');
    iUsed = iUsed(1:numUsed);
    x = [ones(1,numUsed)*numUsed; 1:numUsed].' / numUsed;
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
      phaseDiffWeightedMean(1,iData) = 0;
      b = [1; 0];  % dummy value to avoid error if standard error is requested
    else
      % weighted least-squares-fit with linear equation
      b = bsxfun(@times, weights, x(iUsed,:)) \ (weights .* phaseUnwrap(iUsed,iData));
      phaseDiffWeightedMean(1,iData) = b(2) / numUsed;  % only take slope, offset doesn't matter
    end
  end

  if nargout > 1
    % calculate standard error
    if ~any(iUsed)
      % all amplitudes are NaN
      phaseDiffWeightedStd(1,iData) = NaN;
    else
      % See: https://de.wikipedia.org/wiki/Standardfehler_der_Regression
      % sigma = sqrt(sum(((x(iUsed,:)*b-phaseUnwrap(iUsed,iData)).^2)) / ...
      %         (nSamples-2));
      % together with:
      % https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Statistical_properties
      sigma = sqrt(...
        1/(numUsed-2) ...
        * sum(((phaseUnwrap(iUsed,iData) - x(iUsed,:)*b) ...
               .* (weights/sum(weights)*numUsed)).^2));
      phaseDiffWeightedStd(1,iData) = sigma ...
        / sqrt(sum((x(iUsed,2) - mean(x(iUsed,2))).^2)) ...
        / numUsed;
    end
  end
end


%% reverse (potential) initial permutation of the input
phaseDiffWeightedMean = ipermute(phaseDiffWeightedMean, nsize);

if nargout > 1
  phaseDiffWeightedStd = ipermute(phaseDiffWeightedStd, nsize);
end


end
