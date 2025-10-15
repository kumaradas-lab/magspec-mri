%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 44;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

% correlation of DAC channels to gradient orientations
HW.Grad(iDevice).x = 1;
HW.Grad(iDevice).y = 2;
HW.Grad(iDevice).z = 3;
HW.Grad(iDevice).B = 4; 

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [0.0007978848, -0.001510718, 0.001026152, -0.0003211803];  % offset current in A  14-Feb-2022 14:49:11

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3347015, 0.3320718, 0.3341042, 0.3338977] - HW.Grad(iDevice).PaOffsetI) ./ 1;  % amplification in A/V  14-Feb-2022 14:49:11

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tRamp = 500e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s
HW.Grad(iDevice).tEC = 200e-6;                                 % eddy current time in s

HW.Grad(iDevice).SystemTimeDelay(1:3) = [160e-6, 140e-6, 160e-6];  % time delay of gradient amplifier in s
HW.Grad(iDevice).MaxAmpSlice = 0.05;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
