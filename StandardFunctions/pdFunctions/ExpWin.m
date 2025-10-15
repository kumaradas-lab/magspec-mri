function window = ExpWin(windowSize, LW, fSample, doPlot)
%% Windows function for exponential decay
%
%   window = ExpWin(windowSize, LW, fSample, doPlot)
%
% INPUT:
%   windowSize      number of elements for "window"
%   LW              linewidth (FWHM) of FFT of window
%   fSample         sample frequency
%   doPlot          show figure with window function
%
% OUTPUT:
%   window          column vector with window function

if nargin < 1 || isempty(windowSize), windowSize = 2^10; end
if nargin < 3 || isempty(fSample),    fSample = windowSize; end
if nargin < 2 || isempty(LW),         LW = 2*windowSize/fSample; end
if nargin < 4 || isempty(doPlot),     doPlot = false; end

k = (0:windowSize-1).' / fSample;
window = exp(-pi*LW*k);

if doPlot
    figure(111);
    plot(k, window);
    xlabel('time in s');
    ylabel('window amplitude');
    ylim([0 max(ylim)]);
end

end
