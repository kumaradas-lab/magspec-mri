%% Calibration values for audio amplifiers Amp1k2

if ~exist('iDevice', 'var'), iDevice = 1; end


%% gradient audio amplifier AmpN1k2 SN 003 (A:CH3 & B:CH4) and SN 004 (A:CH1 & B:CH2)
HW.Grad(iDevice).ExtGradSN = 003004;  % serial number
HW.Grad(iDevice).ExtGradType = 'N1k2';  % type of external gradient amplifier

HW.Grad(iDevice).PowerDown = 0;  % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;  % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [0, 0, 0, 0];  % if current controlled, set to 1; if voltage controlled, set to 0
% HW.Grad(iDevice).AmpCurrentDependent(1:4) = [0, 0, 0, 0];  % voltage dependent output - set in LoadSystem_Specific

HW.Grad(iDevice).PaRin(1:4) = [20e3, 20e3, 20e3, 20e3];  % input impedance (data sheet: 94 kOhm)

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0.01, 0.01, 0.01];  % offset voltage in V - 25-06-2024 after 60 min pre-heat
HW.Grad(iDevice).PaOffsetI(1:4) = [0, 0, 0, 0];  % offset current in A

HW.Grad(iDevice).PaUin2PaUout(1:4) = (24.55 - HW.Grad(iDevice).PaOffsetU(1:4)) ./ 1;  % amplification in V/V  (data sheet: 27.6-28 dB)

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];  % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = 0;  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;  % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;  % eddy current time in s
HW.Grad(iDevice).SystemTimeDelay(1:4) = [0e-6, 0e-6, 0e-6, 0e-6];  % time delay of gradient amplifier in s
% HW.Grad(iDevice).MaxAmpSlice = 0.1;  % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;  % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 2;  % gradient and temperature of DC-600 ok

% HW.Grad(iDevice).PaUoutMax(1:4) = 5*HW.Grad(iDevice).PaUin2PaUout;  % 5 Volts
% maximum at amplifier input -- doesn't work for voltage controlled amplifiers
