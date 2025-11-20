function [sigma, lamda, kappa] = get_DWGslk(epsilon, s, doPlot)
%% Calculate form factors for gradient shapes in diffusion encoding
%
%   [sigma, lamda, kappa] = get_DWGslk(epsilon, s, doPlot)
%
% Calculate gradient pulse form factor parameters sigma, lambda, and kappa
% according to equations 9, 10, and 11 in D. Sinnaeve, Concepts Magn Reson Part
% A 40A: 39–65, 2012 [1].
%
%
% INPUT:
%
%   epsilon
%     Time for gradient shape scaled to interval from t = [0...1].
%
%   s
%     Amplitude for gradient shape at sampling points epsilon.
%
%   doPlot
%     Plot shape and (intermediate) results.
%
%
% OUTPUT:
%
%   sigma, lambda, kappa
%     Gradient pulse form factors for calculation of diffusion encoding
%     amplitude.
%
%
% [1]: https://onlinelibrary.wiley.com/doi/epdf/10.1002/cmr.a.21223
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if nargin < 3, doPlot = 0; end

% example input (trapezoid)
% epsilon = [-Seq.DWG.Duration/2, -Seq.DWG.Duration/2+Seq.DWG.tRamp, ...
%   Seq.DWG.Duration/2-Seq.DWG.tRamp, Seq.DWG.Duration/2]./Seq.DWG.Duration + 0.5;
% s = [0, 1, 1, 0];


%% Integrate pulse shape function

% Linear distribution from 0 to 1, with 1000001 samples
% Additionally add the input sample points "epsilon". ("unique" also sorts)
epsilonN = unique([linspace(epsilon(1),epsilon(end),1000001), epsilon], 'sorted');
% re-sample input amplitude "s" to finer grid
sN = interp1(epsilon, s, epsilonN, 'linear');
sN = abs(sN);  % FIXME: Is this correct?

% numerically integrate

% equation 9
% FIXME: Do we need to make sure that epsilon(end) == 1?
S = cumsum([0, diff(epsilon).*(s(1:end-1)+s(2:end))./2], 2);
sigma = S(end);

SN = cumsum([0, diff(epsilonN).*(sN(1:end-1)+sN(2:end))./2], 2);  % equation 12
lamda = 1./sigma .* cumsum([0,diff(epsilonN).*(SN(1:end-1)+SN(2:end))./2], 2);  % equation 10
kappa = 1./sigma.^2 .* cumsum([0,diff(epsilonN).*((SN(1:end-1).^2)+(SN(2:end).^2))./2], 2);  % equation 11

if doPlot
  S2N = cumsum([0, diff(epsilonN).*(SN(1:end-1).^2+SN(2:end).^2)./2], 2);  % ???
  b = cumsum([0, diff(epsilonN).*(S2N(1:end-1)+S2N(2:end))./2], 2);  % ???
  hf = figure(11);
  clf(hf);
  hax = axes(hf);
  plot(hax, epsilon, s, ...  % pulse shape
    epsilonN, sN, ...  % interpolated pulse shape
    epsilonN, SN, ...  % integrated pulse shape
    epsilonN, lamda, ...
    epsilonN, kappa, ...
    epsilonN, S2N, ...
    epsilonN, b);
  legend(hax, {'s', 'sN', 'SN', 'lamda', 'kappa', 'S2N', 'b'});
end

lamda = lamda(end);
kappa = kappa(end);

end
