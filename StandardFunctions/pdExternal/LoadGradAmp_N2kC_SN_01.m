%% Calibration values for audio amplifiers Amp N2k

if ~exist('iDevice', 'var'), iDevice = 1; end

%% gradient audio amplifier Amp N2k Compact SN 01
HW.Grad(iDevice).ExtGradSN = 01;  % serial number
HW.Grad(iDevice).ExtGradType = 'N2kC';  % type of external gradient amplifier

HW.Grad(iDevice).PowerDown = 0;  % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;  % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(3) = [0];  % if current controlled, set to 1; if voltage controlled, set to 0
% HW.Grad(iDevice).AmpCurrentDependent(3) = [0];  % voltage dependent output - set in LoadSystem_Specific

HW.Grad(iDevice).PaRin(3) = [20e3];  % input impedance (data sheet: 94 kOhm)

HW.Grad(iDevice).PaOffsetU(3) = [0.0];  % offset voltage in V - 25-06-2024 after 60 min pre-heat
HW.Grad(iDevice).PaOffsetI(3) = [0];  % offset current in A

HW.Grad(iDevice).PaUin2PaUout(3) = 10;

HW.Grad(iDevice).PaPmaxInt(3) = [100];  % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(3) = 0;  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;  % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;  % eddy current time in s
HW.Grad(iDevice).SystemTimeDelay(3) = [50e-6];  % time delay of gradient amplifier in s

HW.Grad(iDevice).MaxAmpSlice = 0.1;  % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;  % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 2;  % gradient and temperature of DC-600 ok

% HW.Grad(iDevice).PaUoutMax(3) = 5*HW.Grad(iDevice).PaUin2PaUout;  % 5 Volts
% maximum at amplifier input -- doesn't work for voltage controlled amplifiers
