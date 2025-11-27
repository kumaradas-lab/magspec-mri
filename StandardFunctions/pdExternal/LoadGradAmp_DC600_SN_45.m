%% Calibration values for gradient amplifier DC-600
% don't touch Ch4 for resistor switching box

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 45;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 0];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:3) = [24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:3) = [0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:3) = [0.001642142, 0.0004985007, 0.002987378];  % offset current in A  15-Feb-2022 08:13:42

HW.Grad(iDevice).PaUin2PaIout(1:3) = ([0.3349747, 0.3342055, 0.3361173] - HW.Grad(iDevice).PaOffsetI(1:3)) ./ 1;  % amplification in A/V  15-Feb-2022 08:13:42

HW.Grad(iDevice).PaPmaxInt(1:3) = [100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:3) = [15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s
% User specific gradient delays are set in LoadSystem_Specific
HW.Grad(iDevice).SystemTimeDelay(1:3) = [77.472, 103.170, 120.346]*1e-6;  % time delay of gradient amplifier in s
HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok

%% Calibration values for gradient amplifier DC-600
% original calibration
% 
% if ~exist('iDevice', 'var'), iDevice = 1; end
% 
% HW.Grad(iDevice).ExtGradSN = 45;                              % serial number
% 
% HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
% HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier
% 
% HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0
% 
% HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance
% 
% HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
% HW.Grad(iDevice).PaOffsetI(1:4) = [0.001642142, 0.0004985007, 0.002987378, -0.0007462725];  % offset current in A  15-Feb-2022 08:13:42
% 
% HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3349747, 0.3342055, 0.3361173, 0.3319795] - HW.Grad(iDevice).PaOffsetI) ./ 1;  % amplification in A/V  15-Feb-2022 08:13:42
% 
% HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W
% 
% HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm
% 
% HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
% HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s
% % TODO: Run Grad Delay with DC-600 calibration!
% HW.Grad(iDevice).SystemTimeDelay(1:3) = [50e-6, 50e-6, 50e-6];  % time delay of gradient amplifier in s
% 
% % magnet_188 - Konsole 35: 
% HW.Grad(iDevice).SystemTimeDelay(1:3) = [7.4316e-05   0.00010474  0.000124528]; % Time delay of grad amp
% 
% HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m
% 
% HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
% HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
