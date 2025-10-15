%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 38;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [0.000418361, -0.002221, -0.00314363, -0.005902653];  % offset current in A, 16-Mar-2021 14:21:08

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3327872, 0.3303737, 0.3289998, 0.3271452] - HW.Grad(iDevice).PaOffsetI(1:4)) ./ 1;  % amplification in A/V, 16-Mar-2021 14:21:08

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s

HW.Grad(iDevice).SystemTimeDelay(1:3) = [2.37896e-05, 3.25888e-05, 2.58152e-05];  % Time delay of grad amp
HW.Grad(iDevice).MaxAmpSlice=0.1;                             % Max Grad Amp ?

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
