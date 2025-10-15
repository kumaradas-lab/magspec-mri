%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 22;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [-0.002464537, 0.003317852, -0.005837413, -0.001497292];  % offset current in A, 05-Nov-2018 12:37:21

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.33094, 0.335337, 0.3272828, 0.3297582] - HW.Grad(iDevice).PaOffsetI(1:4)) ./ 1;  % amplification in A/V, 05-Nov-2018 12:37:21

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s

switch HW.UserName
  case 'magnet_01_probe_1H' % Magnet 164
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [2.12384e-05, 2.54496e-05, 2.54552e-05];  % time delay of gradient amplifier in s
  case 'sampleHeater' % Magnet 170
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [7.492e-05, 0.000105309, 0.000118237];  % time delay of gradient amplifier in s
  otherwise % Magnet 164
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [2.12384e-05, 2.54496e-05, 2.54552e-05];  % time delay of gradient amplifier in s
end

HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
