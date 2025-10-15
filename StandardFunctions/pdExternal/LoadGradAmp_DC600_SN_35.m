%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 35;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [0.001402637, 0.005817095, -0.003808113, -0.002240312];  % offset current in A, 15-Mar-2021 16:13:37

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.333951, 0.3382715, 0.3286292, 0.3303183] - HW.Grad(iDevice).PaOffsetI(1:4)) ./ 1;  % amplification in A/V, 15-Mar-2021 16:13:37

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s

switch HW.UserName
  case 'magnet_01_probe_1H'
    % Magnet 1: 1H Wechselprobenkopf
    % Magnet 165
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [1.68496e-05, 1.9252e-05, 1.6668e-05];  % time delay of gradient amplifier in s
  case 'magnet_01_probe_31P'
    % Magnet 1: 31P Wechselprobenkopf
    % Magnet 165
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [1.68496e-05, 1.9252e-05, 1.6668e-05];  % time delay of gradient amplifier in s
  case 'magnet_02'
    % Magnet 2: 10mm Magnet KEIN Wechselprobenkopf
    % Magnet 166
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [2.27648e-05, 2.64584e-05, 2.35296e-05];  % time delay of gradient amplifier in s
end

HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
