%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 8;                               % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [0.002781935, 0.0007123647, -0.004950025, -0.002564188];  % offset current in A, 11-Nov-2016 16:11:09

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3359287, 0.3340262, 0.3302193, 0.3319968] - HW.Grad(iDevice).PaOffsetI(1:4)) ./ 1;  % amplification in A/V, 11-Nov-2016 16:12:48

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 42e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s

HW.Grad(iDevice).SystemTimeDelay(1:3) = [1.7648e-05, 1.9408e-05, 1.756e-05];  % time delay of gradient amplifier in s
HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
