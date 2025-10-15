%% Properties of coil system of portable TWMPI.

if ~exist('iDevice', 'var'), iDevice = 1; end

% Most of these settings are in LoadSystem_Specific (device #231) now.
% FIXME: Can this file be deleted?

%% gradient coil system
if all(HW.Grad(iDevice).PaCurrentControlled)
  % for DC-600
  HW.Grad(iDevice).LoadIin2Amp(1:4) = 1;  % amplitude in amp-units/Ampere

  % adjust when efficiency is known
  HW.Grad(iDevice).AmpUnit(1:4) = {'A'};
  HW.Grad(iDevice).AmpUnitScale(1:4) = 1;
else
  % for audio amplifier TA2400 MK-X
  HW.Grad(iDevice).LoadUin2Amp(1:4) = 1;  % amplitude in amp-units/Volts

  % adjust when efficiency is known
  HW.Grad(iDevice).AmpUnit(1:4) = {'V'};
  HW.Grad(iDevice).AmpUnitScale(1:4) = 1;
end
% FIXME: Measure resistance of coils
HW.Grad(iDevice).LoadRin(1:4) = [2, 2, 2, 2];  % resistance of coils - approx. 1.5 - 2.5 Ohms
% FIXME: Measure inductance of coils
% HW.Grad(iDevice).Inductance = [];

% model parameters for heat in coil groups
HW.Grad(iDevice).CoilThermalGroup = [1, 1, 1, 1];  % assignment of gradient coils to thermal groups
HW.Grad(iDevice).CoilPowerDissipation = 2; % Watt (maximum approx. 5A^2*2Ohm = 50 Watt for gradient coils)
HW.Grad(iDevice).CoilTemperatur = 20;  % °C
HW.Grad(iDevice).CoilMaxTemperature = 60;  % °C
mass = (0.4e-3/2)^2*pi * (0.12 * 20 * 4 * (2 + 3) * 2) * 8.96e3;  % mass of copper in coils in kg
heatcap_copper = 385; % heat capacity of copper in J/kg/K
% thermal capacity for each coil group in J/K
HW.Grad(iDevice).CoilThermalCapacity = mass*heatcap_copper;  % imaging gradients

% model parameters for fuse blow estimation
HW.Grad(iDevice).CoilMaxDcCurrent = [1, 1, 1, 1] * 100 * 0.9;  % A
HW.Grad(iDevice).CoilCurrentSquareTime = [0.9, 0.9, 0.9, 0.9] * 100 * 0.9;  % A^2*sec


% simple gradient to DAC channel relation
HW.Grad(iDevice).x = 1;
HW.Grad(iDevice).y = 2;
HW.Grad(iDevice).z = 3;
HW.Grad(iDevice).B = 4;
% simple sign
HW.Grad(iDevice).xyzBDir = [1, 1, 1, 1];  % x y z B0 sign
% gradient names
HW.Grad(iDevice).Name = {'1:2x3', '2:3x2', '3:4y3', '4:5y2'};
% HW.Grad(iDevice).Name = {'x3', 'x2', 'y3', 'y2'};
% HW.Grad(iDevice).Name = {'x_cos', 'x_sin', 'y_cos', 'y_sin'};

shimIdx = sum([HW.Grad(1:(iDevice-1)).n]) + (1:HW.Grad(iDevice).n);
HW.MagnetShim(shimIdx) = zeros(1, HW.Grad(iDevice).n);
