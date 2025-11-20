%% Calibration values for gradient amplifier DC-600
% original calibration

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 41;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
% HW.Grad(iDevice).PaOffsetI(1:4) = [0.0015263,0.002125053,0.000573302,-0.0005810385]; % Offset Current 11-Feb-2022 10:01:41
HW.Grad(iDevice).PaOffsetI(1:4) = [-0.0002590115, 0.001552958, -0.001099738, -0.001513675];  % offset current in A  24-Oct-2022 10:44:39

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3326805, 0.3349972, 0.332308, 0.3322325] - HW.Grad(iDevice).PaOffsetI) ./ 1;  % amplification in A/V  24-Oct-2022 10:44:39
% HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3346045, 0.335718, 0.3340415, 0.3333802] - HW.Grad(iDevice).PaOffsetI) ./ 1;  % Input Voltage to output Current ratio% Input Voltage to output Current ratio 11-Feb-2022 10:01:41

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s

HW.Grad(iDevice).SystemTimeDelay(1:3) = [41.034, 55.630, 47.734]*1e-6;  % time delay of gradient amplifier in s
HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok

%% Calibration values for gradient amplifier DC-600
% RESISTOR SWITCHING BOX calibration
%
% if ~exist('iDevice', 'var'), iDevice = 1; end
%
% HW.Grad(iDevice).ExtGradSN = 41;                              % serial number
%
% HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
% HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier
%
% HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 0];     % if current controlled, set to 1; if voltage controlled, set to 0
%
% HW.Grad(iDevice).PaRin(1:3) = [24e3, 24e3, 24e3];       % INA137 input impedance
%
% HW.Grad(iDevice).PaOffsetU(1:3) = [0, 0, 0];               % offset voltage in V
% HW.Grad(iDevice).PaOffsetI(1:3) = [-0.0002590115, 0.001552958, -0.001099738];  % offset current in A  24-Oct-2022 10:44:39
%
% HW.Grad(iDevice).PaUin2PaIout(1:3) = ([0.3326805, 0.3349972, 0.332308] - HW.Grad(iDevice).PaOffsetI(1:3))./ 1;  % amplification in A/V  24-Oct-2022 10:44:39
%
% HW.Grad(iDevice).PaPmaxInt(1:3) = [100, 100, 100];       % maximum internal power dissipation in W
%
% HW.Grad(iDevice).PaRout(1:3) = [15000, 15000, 15000];  % output impedance in Ohm
%
% HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
% HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s
%
% % 13C
% % HW.Grad(iDevice).SystemTimeDelay(1:3) = [41.614, 54.732, 49.752]*1e-6;  % time delay of gradient amplifier in s
% % 15N
% % HW.Grad(iDevice).SystemTimeDelay(1:3) = [42.672, 54.432, 45.558]*1e-6;  % time delay of gradient amplifier in s
% % 31P
% HW.Grad(iDevice).SystemTimeDelay(1:3) = [42.820, 52.664, 42.680]*1e-6;  % time delay of gradient amplifier in s
%
% HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m
%
% HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
% HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
