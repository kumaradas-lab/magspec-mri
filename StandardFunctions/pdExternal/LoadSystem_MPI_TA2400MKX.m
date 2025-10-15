%% Calibration values for audio amplifiers TA2400 MK-X

if ~exist('iDevice', 'var'), iDevice = 1; end


%% gradient audio amplifier TA2400 MK-X
HW.Grad(iDevice).ExtGradSN = 0;  % serial number

HW.Grad(iDevice).PowerDown = 0;  % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;  % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [0, 0, 0, 0];  % if current controlled, set to 1; if voltage controlled, set to 0
HW.Grad(iDevice).AmpCurrentDependent(1:4) = [0, 0, 0, 0];  % voltage dependent output

HW.Grad(iDevice).PaRin(1:4) = [20e3, 20e3, 20e3, 20e3];  % input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0.01, 0, 0];  % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [0, 0, 0, 0];  % offset current in A

HW.Grad(iDevice).PaUin2PaUout(1:4) = (91.2 - HW.Grad(iDevice).PaOffsetU) ./ 1;  % amplification in A/V  30-Aug-2022 13:05:54

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];  % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = 0.1;  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;  % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;  % eddy current time in s
HW.Grad(iDevice).SystemTimeDelay(1:4) = [0e-6, 0e-6, 0e-6, 0e-6];  % time delay of gradient amplifier in s
% HW.Grad(iDevice).MaxAmpSlice = 0.1;  % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 0;  % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 0;  % gradient and temperature of DC-600 ok

% HW.Grad(iDevice).PaUoutMax(1:4) = 5*HW.Grad(iDevice).PaUin2PaUout;  % 5 Volts
% maximum at amplifier input -- doesn't work for voltage controlled amplifiers
