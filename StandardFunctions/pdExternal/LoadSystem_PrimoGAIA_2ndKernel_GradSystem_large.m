%% Settings for small gradient system for the PrimoGaia measurement project
%
% ------------------------------------------------------------------------------
% (C) Copyright 2022-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

if ~exist('iDevice', 'var')
  iDevice = 1;
end

if numel(HW.Grad) >= iDevice
  %% By default first drive-l for gradient system

  % 4-channel gradient system (S. Lotter, p. 56)
  % S1 --> channel 2 (B0 direction)
  % S2 --> channel 1 (tube direction)
  % S3 --> channel 4 (unused)
  % S4 --> channel 3 (perpendicular to tube and B0 direction)

  % S2 --> channel 1 (tube direction)
  % S1 --> channel 2 (B0 direction)
  % S4 --> channel 3 (perpendicular to tube and B0 direction)
  % S3 --> channel 4 (unused)


  % Channel names and units for sequence plots
  HW.Grad(iDevice).Name = {'Grad x', 'Grad y', 'Grad z', 'Grad 4'};
  HW.Grad(iDevice).AmpUnit = {'mT/m', 'mT/m', 'mT/m', 'mT/m'};
  HW.Grad(iDevice).AmpUnitScale = [1e-3, 1e-3, 1e-3, 1e-3];

  maxGrad_Dev1 = 4;

  if 1
    %% large gradient system

    % correlation channels to encoding directions (to match coordinate system of B0 cage)
    HW.Grad(iDevice).x = 1;
    HW.Grad(iDevice).y = 3;
    HW.Grad(iDevice).z = 4;
    HW.Grad(iDevice).B = 2;  % 4th gradient channel

    % coil efficiency
    % Measured in (T/m)/A (Tesla per Meter per Ampere);
    % or more generally: [amplitude unit]/[amplifier input unit]
    % signs without any filters
    % 4th channel not calibrated!
     HW.Grad(iDevice).LoadIin2Amp(1:3) = [...
       0.000166, ...
       0.000166*1.5, ...
       0.000166/2*54/60];  % x y z

    % signs for right-handed coordinate system
    HW.Grad(iDevice).xyzBDir(1:4) = [1, 1, 1, 1];  % x y z prepol - sign (for inversing direction)

    % resistance of gradient coils (including cabling)
    HW.Grad(iDevice).LoadRin(1:4) = [7.9, 8.3, 7.9, 6.9];  % only gradients (without cabling)
    HW.Grad(iDevice).LoadRin(1:4) = [8.4, 8.3+0.5, 8.3, 6.9+0.5];  % gradients + cable

  end


  % gradients used for auto shimming routine
  HW.Grad(iDevice).ShimGradients = [1, 1, 1, 0];

  % extent (and location) of the image volume in meters
  HW.Grad(iDevice).ImageVol = [-0.06, 0.06, -0.04, 0.04, -0.04, 0.04]*1.5;  % FIXME: Anpassen
  HW.Grad(iDevice).ImageVolOffset = [0, 0, 0];  % offset of tube

  % FIXME: Alles anpassen
%   HW.Grad(iDevice).tRamp = 20e-3;  % minimum ramp time (for imaging gradients) in s

  % die folgenden drei funktionieren für Gradienten um die 50 µT/m
  HW.Grad(iDevice).tRamp = 2e-3;  % minimum ramp time (for imaging gradients) in s
    HW.Grad(iDevice).tEC = 6e-3+2e-3;  % eddy current time in s (i.e. time until gradient is stable)
%   HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [0.0055-0.8e-3, 0.0015+1.4e-3+1e-3+0.75e-3, 0.0015+2.4e-3+2e-3*2-4.5e-3+0.5e-3];  % Time delay of grad amp
%   HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [0.0055-0.8e-3-0.5e-3, 0.0015+1.4e-3+1e-3, 0.0015+2.4e-3+2e-3*2-4.5e-3-1.2e-3/2];  % Time delay of grad amp
%   HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [0.0055-0.8e-3-0.5e-3-1e-3, 0.0015+1.4e-3+1e-3, 0.0015+2.4e-3+2e-3*2-4.5e-3-1.2e-3/2];  % Time delay of grad amp
%   HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [0.0032, 0.0039, 0.0028];  % Time delay of grad amp
%   HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [0.00675, 0.0065, 0.0065];  % 1 mm Aluroehre Time delay of grad amp
%   HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [0.0075, 0.009, 0.01]; % schwer einzustellen da moeglicher fSample Fehler
%   HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [0.0075-1.5e-3, 0.006-2e-3, 0.009-4e-3]; % schwer einzustellen da moeglicher fSample Fehler
%   HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [0.015, 0.015,  0.015];  % Time delay of grad amp mit großer 3 mm Aluröhre

  % gradient time delay with 3mm aluminum shield
  HW.Grad(iDevice).SystemTimeDelay(1:4) = [10.103, 0.050, 8.931, 9.427]*1e-3;  % gradient time delay in s - 2024-01-11 11:16:02

%   HW.Grad(iDevice).tRamp = 5e-3;  % minimum ramp time (for imaging gradients) in s
%   HW.Grad(iDevice).tEC = 5e-3;  % eddy current time in s (i.e. time until gradient is stable)
%   HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [0.0055-0.8e-3, 0.0015+1.4e-3+1e-3, 0.0015+2.4e-3+2e-3*2-4.5e-3];  % Time delay of grad amp

  % maximum total power (coil heating, model parameters)
  % FIXME: Adapt to properties of used coils.
  HW.Grad(iDevice).CoilPowerDissipation(1) = 5;  % power dissipation at maximum temperature for each coil group in Watt
  HW.Grad(iDevice).CoilTemperatur(1) = 20;  % resting temperature for each coil group in degrees C
  HW.Grad(iDevice).CoilMaxTemperature(1) = 60;  % maximum temperature for each coil group in degrees C
  HW.Grad(iDevice).CoilThermalCapacity(1) = 0.06 * 0.5 * 0.5 * 3 * 0.8 * 1000;  % thermal capacity for each coil group in J/K

  % properties for "fuse blow" model
  % There is actually no fuse to protect gradient system and other component.
  % The "weakest" component is probably the inductance in the gradient filter.
  % Setting a fuse current of 1.5 A is probably reasonable to protect that
  % component (Toni, 30-Sep-2024).
  HW.Grad(iDevice).CoilMaxDcCurrent(1:maxGrad_Dev1) = 1.5 * 0.9;  % maximum (nominal) DC current of fuse in Ampere
  HW.Grad(iDevice).CoilCurrentSquareTime(1:maxGrad_Dev1) = 0.9 * 0.9;  % maximum "accumulated heat" in A^2*sec


  HW.Grad(iDevice).MaxAmpSlice = 50e-6;
else
  warning('PD:B0Cage:DeviceNotConnected', ...
    'The selected device %d for the gradient system is not connected.', iDevice);
end
