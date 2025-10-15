%% Parameters for gradient system

if ~exist('iDevice', 'var'), iDevice = 1; end

% assignment gradient direction to DAC channel
HW.Grad(iDevice).x = 1;
HW.Grad(iDevice).y = 2;
HW.Grad(iDevice).z = 3;
HW.Grad(iDevice).B = 4;

% gradient sign
HW.Grad(iDevice).xyzBDir(1:4) = [1, 1, 1, 1];  % x y z y2

% gradient names and units for pulse program plot
HW.Grad(iDevice).Name(1:4) = {'Grad x', 'Grad y', 'Grad z', 'Grad y2'};
HW.Grad(iDevice).AmpUnit(1:4) = {'mT/m', 'mT/m', 'mT/m', 'mT/m'};
HW.Grad(iDevice).AmpUnitScale(1:4) = [1e-3, 1e-3, 1e-3, 1e-3];
HW.Grad(iDevice).LengthUnit = 'mm';
HW.Grad(iDevice).LengthUnitScale = 1e-3;

if ~exist('combinedChannels', 'var')
  combinedChannels = false;  % set to true if y-gradient is connected to 2 channels (2 and 4)
end
if combinedChannels
  HW.Grad(iDevice).CombineCurrentOutputs = [2; 4];  % channels, each column will be treated as a combined output
else
  HW.Grad(iDevice).CombineCurrentOutputs = [];
end

% resistance of gradients
HW.Grad(iDevice).LoadRin(1:4) = [6.2, 3.9, 4.7, 3.9];  % x y z y2

% inductance of gradient in Henry
HW.Grad(iDevice).Inductance(1:4) = [0, 0, 0, 0];  % x y z y2

% gradient efficiency in (T/(m*A)) Tesla per Meter per Ampere
HW.Grad(iDevice).LoadIin2Amp(1:4) = [0.061952, 0.063381/2, 0.071987, 0.063381/2];  % x y z B0  (T/(m*A))

if combinedChannels
  % Independently of what was set before, element 2 and 4 must always have the
  % same values for this setup! (Y-connector or just mirror)
  HW.Grad.LoadRin(4) = HW.Grad.LoadRin(2);
  HW.Grad.LoadIin2Amp(4) = HW.Grad.LoadIin2Amp(2);

  HW.Grad(iDevice).LoadIin2Amp(HW.Grad.CombineCurrentOutputs(1,:)) = 2*HW.Grad(iDevice).LoadIin2Amp(HW.Grad(iDevice).CombineCurrentOutputs(1,:));
  HW.Grad(iDevice).LoadIin2Amp(HW.Grad.CombineCurrentOutputs(2,:)) = HW.Grad(iDevice).LoadIin2Amp(HW.Grad(iDevice).CombineCurrentOutputs(1,:));
  HW.Grad(iDevice).LoadRin(HW.Grad.CombineCurrentOutputs(2,:)) = HW.Grad(iDevice).LoadRin(HW.Grad(iDevice).CombineCurrentOutputs(1,:));
end

% maximum amplitude in T/(m*A)
HW.Grad(iDevice).MaxAmp(1:4) = Inf(1, 4);  % x y z y2

% gradients used for shimming
HW.Grad(iDevice).ShimGradients(1:4) = [1 1 1 0];

% maximum image volume
HW.Grad(iDevice).ImageVol = [-0.025, 0.025, -0.025, 0.025, -0.025, 0.025];  % [xmin xmax ymin ymax zmin zmax]
HW.Grad(iDevice).ImageVolOffset = [0, 0, 0];  % [x y z] offset of tube in m

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

%% Might also be set in gradient amplifier configuration file

% time delay of gradient amplifier and gradient system
HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [50e-6, 50e-6, 50e-6];  % gradient output channel ordering
% HW.Grad(iDevice).MaxAmpSlice = 0.002;
