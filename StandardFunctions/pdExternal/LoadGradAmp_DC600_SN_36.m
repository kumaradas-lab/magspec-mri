%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 36;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [-0.001173692, -0.005148218, -0.003356433, -0.0002601287];  % offset current in A, 16-Mar-2021 12:59:20

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3310783, 0.3277, 0.3296105, 0.3329217] - HW.Grad(iDevice).PaOffsetI(1:4)) ./ 1;  % amplification in A/V, 16-Mar-2021 12:59:20

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s

HW.Grad.SystemTimeDelay(1:3) = [2.2148e-05  2.69016e-05  2.40048e-05];  % time delay of gradient amplifier in s
HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
