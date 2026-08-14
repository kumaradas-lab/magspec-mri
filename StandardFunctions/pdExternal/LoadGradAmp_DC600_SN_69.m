%% Calibration values for gradient amplifier DC-600 V2
% CH4 5 Ohm, max. output current CH4 0.5A
% output port: RJ-45

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 69;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [0.0009001385, 0.004552188, 0.00517718, -1.264858e-05];  % offset current in A  20-Oct-2025 09:33:30

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3344253, 0.337384, 0.3382685, 0.03333673] - HW.Grad(iDevice).PaOffsetI) ./ 1;  % amplification in A/V  20-Oct-2025 09:33:30

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

% long eddy currents for 32mm bore magnet 233
HW.Grad(iDevice).tRamp = 100e-6;                              % minimum ramp time in s
HW.Grad(iDevice).tEC = 500e-6;                                % eddy current time in s
% time delay of gradient amplifier is in LoadCoil_MultiCore_12
HW.Grad(iDevice).MaxAmpSlice = 50e-3;                         % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok

% fuses: T3.0A in magnet; T1.0A in DC-600
HW.Grad(iDevice).CoilMaxDcCurrent(:) = 1.0 * 0.9;  % ampere rating of fuse in A
HW.Grad(iDevice).CoilCurrentSquareTime(:) = 1.98;  % nominal melting I^2t of fuse in A^2*sec
