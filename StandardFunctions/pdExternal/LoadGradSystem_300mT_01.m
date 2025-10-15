% Parameters for gradient system in 300 mT magnet

if ~exist('iDevice', 'var'), iDevice = 1; end

% assignment gradient direction to DAC channel
HW.Grad(iDevice).x = 1;
HW.Grad(iDevice).y = 2;
HW.Grad(iDevice).z = 3;
HW.Grad(iDevice).B = 4;

% gradient sign
HW.Grad(iDevice).xyzBDir = [1, 1, 1, 1];  % x y z y2

% gradient names and units for pulse program plot
HW.Grad(iDevice).Name = {'Grad x', 'Grad y', 'Grad z', 'Grad y2'};
HW.Grad(iDevice).AmpUnit = {'mT/m', 'mT/m', 'mT/m', 'mT/m'};
HW.Grad(iDevice).AmpUnitScale = [1e-3, 1e-3, 1e-3, 1e-3];

if ~exist('combinedChannels', 'var')
  combinedChannels = false;  % set to true if y-gradient is connected to 2 channels (2 and 4)
end
if combinedChannels
  HW.Grad(iDevice).CombineCurrentOutputs = [2; 4];  % channels, each column will be treated as a combined output
else
  HW.Grad(iDevice).CombineCurrentOutputs = [];
end

% resistance of gradients
% HW.Grad(iDevice).LoadRin = [5.8, 3.2, 5.9, 3.2];  % x y z y2
% HW.Grad(iDevice).LoadRin = [6.2, 3.9, 6.4, 3.2];  % x y z y2 am 9.8.2017 nur mit kurzen 1.5mm Kabeln
HW.Grad(iDevice).LoadRin = [6.2, 3.9, 4.7, 3.9];  % x y z y2 am 9.8.2017 nur mit kurzen 1.5mm Kabeln

% gradient efficiency in (T/(m*A)) Tesla per Meter per Ampere
% HW.Grad(iDevice).LoadIin2Amp = [0.058448, 0.031694, 0.045762, 0.031694];  % x y z B0  (T/(m*A)) am 18.05.2018 mit Calibrate_GradientSystem (y nur am Kanal 2)
HW.Grad(iDevice).LoadIin2Amp = [0.061952, 0.063381/2, 0.071987, 0.065668/2];  % x y z B0  (T/(m*A)) am 30.05.2018 mit Calibrate_GradientSystem (y an Kanal 2&4)

% Independently of what was set before, element 2 and 4 must always have the
% same values for this setup! (Y-connector or just mirror)
HW.Grad.LoadRin(4) = HW.Grad.LoadRin(2);
HW.Grad.LoadIin2Amp(4) = HW.Grad.LoadIin2Amp(2);

if combinedChannels
  HW.Grad(iDevice).LoadIin2Amp(HW.Grad.CombineCurrentOutputs(1,:)) = 2*HW.Grad(iDevice).LoadIin2Amp(HW.Grad(iDevice).CombineCurrentOutputs(1,:));
  HW.Grad(iDevice).LoadIin2Amp(HW.Grad.CombineCurrentOutputs(2,:)) = HW.Grad(iDevice).LoadIin2Amp(HW.Grad(iDevice).CombineCurrentOutputs(1,:));
  HW.Grad(iDevice).LoadRin(HW.Grad.CombineCurrentOutputs(2,:)) = HW.Grad(iDevice).LoadRin(HW.Grad(iDevice).CombineCurrentOutputs(1,:));
end

% gradients used for shimming
HW.Grad(iDevice).ShimGradients = [1 1 1 0];

% maximum image volume
% HW.Grad(iDevice).ImageVol = [-0.04, 0.04, -0.06, 0.06, -0.04, 0.04];  % [xmin xmax ymin ymax zmin zmax] % "default"-38mm coil
HW.Grad(iDevice).ImageVol = [-0.025, 0.025, -0.025, 0.025, -0.025, 0.025];  % [xmin xmax ymin ymax zmin zmax] % 38mm coil H-free
% HW.Grad(iDevice).ImageVol = [-0.04, 0.04, -0.03, 0.03, -0.04, 0.04];  % [xmin xmax ymin ymax zmin zmax] % 38mm coil H-free
HW.Grad(iDevice).ImageVolOffset = [0, 0, 0];  % offset of tube

% Time delay of gradient amplifier
% HW.Grad(iDevice).SystemTimeDelay = [154e-6,160e-6,190e-6,202e-6]-0e-6;
% HW.Grad(iDevice).SystemTimeDelay = [280e-6,220e-6,320e-6,202e-6]-0e-6;  % am 9.8.2017 mit GradAmp
% HW.Grad(iDevice).SystemTimeDelay = [240e-6,220e-6,320e-6,202e-6]-0e-6;  % am 10.8.2017 mit GradAmp
% HW.Grad(iDevice).SystemTimeDelay = [0.000241344, 0.000277104, 0.000348752, 0.000277104];  % Time delay of grad amp, am 18.05.2018 mit Calibrate_GradDelay (y nur am Kanal 2)
% HW.Grad(iDevice).SystemTimeDelay(1:3) = [0.00023888, 0.000388304, 0.000258664];  % Time delay of grad amp, am 30.05.2018 mit Calibrate_GradDelay (y an Kanal 2&4)
HW.Grad(iDevice).SystemTimeDelay(1:3) = [0.000176902, 0.000371194, 0.000185846];  % 13.10.2021 mit 10 mm Spule

% HW.Grad(iDevice).SystemTimeDelay = [1020e-6, 730e-6, 750e-6, 202e-6]-0e-6;

% time until gradient is stable after ramp (eddy current time)
HW.Grad(iDevice).tEC = max(HW.Grad.SystemTimeDelay(1:3))*.5;
% HW.Grad(iDevice).tEC = max(250e-6)*.5;

HW.Grad(iDevice).tRamp = 250e-6;  % minimum ramp time in s

% model parameters for heat in coil groups
HW.Grad(iDevice).CoilPowerDissipation(:) = 5; % Watt
HW.Grad(iDevice).CoilTemperatur(:) = 20; % °C
HW.Grad(iDevice).CoilMaxTemperature(:) = 60; % °C
HW.Grad(iDevice).CoilThermalCapacity(:) = 0.06*0.5*0.5*3*0.8*1000; % J/K
% model parameters for fuse blow estimation
HW.Grad(iDevice).CoilMaxDcCurrent = [0.75, 6.5, 0.75, 6.5] * 0.9; % A
HW.Grad(iDevice).CoilCurrentSquareTime = [0.9, 0.9, 0.9, 0.9] * 0.9; % A^2*sec

% HW.Grad(iDevice).MaxAmpSlice = 0.002;
% HW.Grad(iDevice).MaxAmp = 0.02;
