%% Settings for B0 field cage for PrimoGAIA 2nd kernel
%
% By default a second console is used for the B0 cage.
%
% The coordinate system is:
%   X: from left to right (along gradient tube)
%   Y: from front to back
%   Z: from bottom to top
%
% ------------------------------------------------------------------------------
% (C) Copyright 2022-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------


if ~exist('iDevice', 'var')
  iDevice = 2;
end

if numel(HW.Grad) >= iDevice
  %% By default second drive-l for B0 field cage

  LoadGradAmp_DC600_SN_39;

  % For simplicity, use a simple correlation:
  HW.Grad(iDevice).x = 1;  % corresponds to Grad(5) in pulse programs
  HW.Grad(iDevice).y = 3;  % corresponds to Grad(6) in pulse programs
  HW.Grad(iDevice).z = 4;  % corresponds to Grad(7) in pulse programs
  HW.Grad(iDevice).B = 2;  % corresponds to Grad(8) in pulse programs
  % HW.Grad(2).xyzBDir = [-1, -1, -1, 1];  % x y z prepol - sign (for inversing direction)
  % exchange y-z coordinates (9th March 2021)
  HW.Grad(iDevice).xyzBDir = [1, 1, 1, 1];  % x y z prepol - sign (for inversing direction)

  % Channel names and units for sequence plots
  HW.Grad(iDevice).Name = {'B0 x', 'B0 y', 'B0 z', 'unused'};
  HW.Grad(iDevice).AmpUnit(1:3) = {[char(181) 'T']};
  HW.Grad(iDevice).AmpUnit{4} = 'A';
  HW.Grad(iDevice).AmpUnitScale = [1e-6, 1e-6, 1e-6, 1];

  HW.Grad(iDevice).HoldShim = 1;

  %% coil efficiency
  % Measured in T/A (Tesla per Ampere);
  % or more generally: [amp unit]/[controller input unit]

  % Alternatively, amplitude in Ampere:
  % HW.Grad(iDevice).LoadIin2Amp = [1, 1, 1, 1];  % x y z XXX?

%   HW.Grad(iDevice).LoadIin2Amp = [5.096746300113425e-05, 6.026636720829278e-05, 5.844179958134018e-05, 1];  % B0 cage (x y z XXX?) calibrated in T/A with MR on 24/02/2022
%   HW.Grad(iDevice).LoadIin2Amp = [5.188719483334751e-05, 5.265215425384763e-05, 5.777788441818475e-05, 1];  % B0 cage (x y z XXX?) calibrated in T/A with MR on 01/04/2022
%   HW.Grad(iDevice).LoadIin2Amp = [5.197075475924804e-05, 5.265180502984102e-05, 5.080516422217215e-05, 1];  % B0 cage (x y z XXX?) calibrated in T/A with MR on 01/04/2022
  HW.Grad(iDevice).LoadIin2Amp = [5.197075475924804e-05, 5.080516422217215e-05, 5.265180502984102e-05, 1];  % B0 cage (x y z XXX?) calibrated in T/A with MR on 01/04/2022

  %% resistance of coils in Ohm
  HW.Grad(iDevice).LoadRin = [14, 14, 14, 14];  % noch nicht gemessen!

  % gradients used for auto shimming routine
  HW.Grad(iDevice).ShimGradients(1:4) = 0;

  % extent (and location) of the image volume in meters
  HW.Grad(iDevice).ImageVol = HW.Grad(1).ImageVol;  % FIXME: Would it make sense if we set something different?
  HW.Grad(iDevice).ImageVolOffset = HW.Grad(1).ImageVolOffset;  % FIXME: Would it make sense if we set something different?

  % FIXME: Alles anpassen
  HW.Grad(iDevice).tRamp = 2e-3;  % minimum ramp time (for B0 cage fields) in s
  HW.Grad(iDevice).SystemTimeDelay = [163e-6, 217e-6, 215e-6, 215e-6]-0e-6;  % time delay of grad amp in s (FIXME: Can this be measured?)
  HW.Grad(iDevice).tEC = max(HW.Grad(iDevice).SystemTimeDelay(1:3))*0.5;  % eddy current time in s (i.e. time until gradient is stable)

  % maximum total power (coil heating, model parameters)
  % FIXME: Adapt to properties of used coils.
  % 4th channel is currently unused. Set some values anyway.
  HW.Grad(iDevice).CoilThermalGroup = 1:4;
  HW.Grad(iDevice).CoilPowerDissipation(1:4) = 220;  % Watt
  HW.Grad(iDevice).CoilTemperatur(1:4) = 20;  % temperature of coil at start of measurement in degrees C
  HW.Grad(iDevice).CoilMaxTemperature(1:4) = 60;  % maximum temperature (in model) before damage occurs in degrees C
  % thermal capacity in J/K;  FIXME: Adjust values?
  specificHeatCapacity = 500; % J/kg/K (estimated for mix of copper wires and aluminum frame)
  volumeCoil = 2*4*1*0.015^2;  % volume of coil in m^3 (2 coils with 4 sides of ~1 m and a cross-section of ~1.5 cm^2)
  mass = volumeCoil*5e3;  % in kg (estimated density for mix of copper wires and aluminum frame)
  HW.Grad(iDevice).CoilThermalCapacity(1:4) = specificHeatCapacity*mass;

  % properties for "fuse blow" model
  HW.Grad(iDevice).CoilMaxDcCurrent = [8, 8.0, 8.0, 8.0] * 0.9;  % maximum (nominal) DC current of fuse in Ampere
  HW.Grad(iDevice).CoilCurrentSquareTime = [0.9, 0.9, 0.9, 0.9] * 0.9;  % maximum "accumulated heat" in A^2*sec

  fprintf('B0 field via cage in z-direction %3.1f %cT (in %s).\n', HW.B0*1e6, char(181), mfilename());
  HW.Grad(iDevice).AmpOffsetExtra(3) = -HW.B0;
%   HW.Grad(iDevice).AmpOffsetExtra(2) = -HW.B0;
%   HW.Grad(iDevice).AmpOffsetExtra(1) = -HW.B0;
%   HW.Grad(iDevice).AmpOffsetExtra(1) = 50e-6;
  % HW.Grad(2).AmpOffsetExtra(1:3) = -cross(HW.Grad(2).AmpOffset(1:3), [1,0,0]) / norm(HW.Grad(2).AmpOffset(1:3),2) * HW.B0;
%   HW.Grad(2).AmpOffsetExtra(1:3) = [9.988, -42.897, -21.812]*1e-6;


%   HW.Grad(iDevice).AmpOffsetExtra(1) = 0e-6; % Polarität falsch
%   HW.Grad(iDevice).AmpOffsetExtra(2) = 0e-6;
%  HW.Grad(iDevice).AmpOffsetExtra(3) = 50e-6;


else
  warning('PD:B0Cage:DeviceNotConnected', ...
    'The selected device %d for the B0 cage is not connected.', iDevice);
end
