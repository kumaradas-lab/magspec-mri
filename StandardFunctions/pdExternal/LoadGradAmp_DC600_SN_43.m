%% Calibration values for gradient amplifier DC-600
% don't touch Ch4 for resistor switching box

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 43;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 0];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:3) = [24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:3) = [0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:3) = [-0.0003164087, 0.005041358, 0.00740063];  % offset current in A  14-Feb-2022 13:51:36

HW.Grad(iDevice).PaUin2PaIout(1:3) = ([0.3336597, 0.3388468, 0.3400932] - HW.Grad(iDevice).PaOffsetI(1:3)) ./ 1;  % amplification in A/V  14-Feb-2022 13:51:36

HW.Grad(iDevice).PaPmaxInt(1:3) = [100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:3) = [15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s
% User specific gradient delays are set in LoadSystem_Specific
% HW.Grad(iDevice).SystemTimeDelay(1:3) = [3.43072e-05, 3.87304e-05, 4.58992e-05];  % time delay of gradient amplifier in s
HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok

%% Calibration values for gradient amplifier DC-600
% original calibration
%
% if ~exist('iDevice', 'var'), iDevice = 1; end
%
% HW.Grad(iDevice).ExtGradSN = 43;                              % serial number
%
% HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
% HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier
%
% HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0
%
% HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance
%
% HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
% HW.Grad(iDevice).PaOffsetI(1:4) = [-0.0003164087, 0.005041358, 0.00740063, -0.001645523];  % offset current in A  14-Feb-2022 13:51:36
%
% HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3336597, 0.3388468, 0.3400932, 0.331954] - HW.Grad(iDevice).PaOffsetI) ./ 1;  % amplification in A/V  14-Feb-2022 13:51:36
%
% HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W
%
% HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm
%
% HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
% HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s
% % TODO: Run Grad Delay with DC-600 calibration!
% HW.Grad(iDevice).SystemTimeDelay(1:3) = [50e-6, 50e-6, 50e-6];  % time delay of gradient amplifier in s
% HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m
%
% HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
% HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
