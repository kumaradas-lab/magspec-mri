%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 10;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [0.0002897053, 0.0005113268, -0.001826162, -0.001214733];  % offset current in A, 07-Nov-2016 15:07:18

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3339845, 0.3331593, 0.3335963, 0.3334057] - HW.Grad(iDevice).PaOffsetI(1:4)) ./ 1;  % amplification in A/V, 07-Nov-2016 15:08:35

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 18e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s
HW.Grad(iDevice).SystemTimeDelay(1:3) = [35.328, 50.680, 47.446]*1e-6;  % time delay of gradient amplifier in s
HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok

% load MRE Setting
HW.Grad(iDevice).LoadRin(4) = 100;  % resistance parallel to Piezo in Ohm
% set correct sequence plot names
HW.Grad(iDevice).Name(4) = {'Piezo Current'};
HW.Grad(iDevice).AmpUnit(4) = {'A'};
HW.Grad(iDevice).AmpUnitScale(4) = 1;
% set correct gradient efficiency
HW.Grad(iDevice).LoadIin2Amp(4) = 1;
