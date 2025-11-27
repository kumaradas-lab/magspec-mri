function window = LorentzToGauss(windowSize, LB, Wmax, fSample, doPlot)
%% Window function that transforms Lorentz frequency profile to a normal distribution
%
%       window = LorentzToGauss(windowSize, LB, Wmax, fSample, doPlot)
%
% Applying the window function that is returned by this function to an
% exponentially decaying (time-domain) function leads to it corresponding to
% having a normal (Gaussian) distribution in the frequency-domain (instead of a
% Lorentz distribution).
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if nargin < 1 || isempty(windowSize)
  windowSize = 2^10;
end
if nargin < 2 || isempty(LB)
  LB = -(windowSize/fSample) * 2;
end
if nargin < 3 || isempty(Wmax)
  Wmax=1/3;
end
if nargin < 4 || isempty(fSample)
  fSample = windowSize;
end
if nargin < 5
  doPlot = 0;
end


%% create window function
alfa = -LB*pi/2/Wmax/(windowSize/fSample);
k = (0:windowSize-1).';
window = exp(-pi*LB*k./fSample) .* exp(-alfa.*(k./fSample).^2);


%% optionally plot window function
if doPlot
  hf = figure(234);
  clf(hf);
  hax = axes(hf);
  plot(hax, k./fSample, window);
  xlabel(hax, 'time / s');
  ylabel(hax, 'window amplitude');
end

end
