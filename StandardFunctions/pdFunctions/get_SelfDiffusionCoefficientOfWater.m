function SelfDiffusionCoefficient = get_SelfDiffusionCoefficientOfWater(Temperature_C)
%% Get approximate self diffusion coefficient of water at temperature
%
%       SelfDiffusionCoefficient = get_SelfDiffusionCoefficientOfWater(Temperature_C)
%
%
% INPUT:
%
%   Temperature_C
%       Temperature in degrees C for which the approximate self diffusion
%       coefficient is calculated.
%
% OUTPUT:
%
%   SelfDiffusionCoefficient
%       Approximate self diffusion coefficent in m^2/s of water at the given
%       temperature. The Speedy-Angell power law (formula (2) in [1]) is used.
%       This formula is in good agreement with measured values above
%       approximately 250 K (approx. -23 degrees C).
%
% [1]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9697084/
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

D0 = 1.635e-8;  % D0 = 1.635e-8 +/- 2.242e-11 (m^2/s)
Ts = 215.05;  % Ts = 215.05 +/- 1.20 K
g = 2.063;  % gamma = 2.063 +/-0.051
% https://pmc.ncbi.nlm.nih.gov/articles/PMC9697084/ Speedy-Angell power law (2)
SelfDiffusionCoefficient = D0 .* (((Temperature_C+273.15) ./ Ts)-1).^g;  

end


% Table 1: Diffusion coefficients D, diaphragm cell; R. Mills (1973) [1]
% Temperature          D in um2/ms
%   1 deg C  274.15 K  1.149
%   4 deg C  277.15 K  1.276
%   5 deg C  278.15 K  1.313
%  15 deg C  288.15 K  1.777
%  25 deg C  298.15 K  2.299
%  35 deg C  308.15 K  2.919
%  45 deg C  318.15 K  3.575
% Table 2: Diffusion coefficients D, diaphragm cell; A.J. Easteal et al. (1989) [2]
% Temperature          D in um^2/ms
%   0 deg C  273.15 K  1.130
%  10 deg C  283.15 K  1.536
%  20 deg C  293.15 K  2.022
%  25 deg C  298.15 K  2.296
%  30 deg C  303.15 K  2.590
%  40 deg C  313.15 K  3.240
%  50 deg C  323.15 K  3.968
%  60 deg C  333.15 K  4.772
%  70 deg C  343.15 K  5.646
%  80 deg C  353.15 K  6.582
%  90 deg C  363.15 K  7.578
% 100 deg C  373.15 K  8.623
% Table 3: Diffusion Coefficients D, NMR; K.R. Harris and L.A. Woolf (1980) [3]
% Temperature          D in um2/ms
%   4 deg C  277.15 K  1.27
%  10 deg C  283.15 K  1.55
%  25 deg C  298.15 K  2.30
%  45 deg C  318.15 K  3.55
%  60 deg C  333.15 K  4.70
% Table 4: Diffusion coefficients D, NMR; P.S. Tofts et al. (2000) [4]
% Temperature          D in um2/ms
%  15 deg C  288.15 K  1.756
%  20 deg C  293.15 K  2.023
%  25 deg C  298.15 K  2.317
%  30 deg C  303.15 K  2.616
% Table 5: Diffusion coefficients D, NMR; Holz et al. (2000) [5]
% Temperature          D in um2/ms
%  15 deg C  288.15 K  1.766
%  20 deg C  293.15 K  2.025
%  25 deg C  298.15 K  2.299
%  30 deg C  303.15 K  2.597
%  35 deg C  308.15 K  2.895
%  40 deg C  313.15 K  3.222
%  45 deg C  318.15 K  3.601
%  50 deg C  323.15 K  3.983
%  56 deg C  329.15 K  4.444
% Table 6: Diffusion coefficients D, NMR; K.T. Gillen et al. (1972) [6]
% Temperature            D in um2/ms
% -30.65 deg C  242.5 K  0.187
% -28.75 deg C  244.4 K  0.219
% -26.85 deg C  246.3 K  0.263
% -24.95 deg C  248.2 K  0.321
% -23.15 deg C  250.0 K  0.341
% -21.35 deg C  251.8 K  0.395
% -19.15 deg C  254.0 K  0.438
% -17.35 deg C  255.8 K  0.477
% -14.45 deg C  258.7 K  0.553
% -11.65 deg C  261.5 K  0.633
%  -9.45 deg C  263.7 K  0.700
%   0.35 deg C  273.5 K  1.05
%  12.15 deg C  285.3 K  1.58
%  25.05 deg C  298.2 K  2.23
% Table 7: Diffusion Coefficients D, NMR; various studies
% Temperature           D/(um2/ms)     Ref.
% 18.4 deg C  291.55 K  1.98 +/- 0.14  [7]
% 18.5 deg C  291.65 K  1.95 +/- 0.02  [8]
% 18.5 deg C  291.65 K  2.03 +/- 0.01  [8]
% 25.0 deg C  298.15 K  2.36 +/- 0.04  [8]
% 25.0 deg C  298.15 K  2.23 +/- 0.06  [9]
