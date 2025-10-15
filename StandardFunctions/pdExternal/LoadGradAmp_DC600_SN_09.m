%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 9;                               % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [1.178049e-05, -0.00272666, -0.001894867, 0.004445753];  % offset current in A

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3356205, 0.3322447, 0.3330938, 0.3380657] - HW.Grad(iDevice).PaOffsetI(1:4)) ./ 1;  % amplification in A/V

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 18e-6*5;                             % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s
% HW.Grad(iDevice).SystemTimeDelay(1:3) = [49.8e-6, 57e-6, 58.28e-6];  % time delay of gradient amplifier in s
HW.Grad(iDevice).SystemTimeDelay(1:3) = [49.8e-6, 57e-6, 58.28e-6];  % time delay of gradient amplifier in s
HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
