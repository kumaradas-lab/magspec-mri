%% Combined gradient amplifiers DC-600 SN_23 and SN_24
%
% amplifier SN_24: output 1&3 for x and output 2&4 for y
% amplifier SN_23: output 1&3 for z and output 2&4 for slice
%
% Y-adapter for connection to Ext. Grad:
% DAC channel 1 to SN_24 input 1&3
% DAC channel 2 to SN_24 input 2&4
% DAC channel 3 to SN_23 input 1&3
% DAC channel 4 to SN_23 input 2&4
%
% Take care to turn on both(!) DC-600. Power on detection is not reliable.

HW.Grad.ExtGradSN = 23024;
HW.Grad.ExtGradType = 'DC600';

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3]/2;                  % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
% 22: HW.Grad.PaOffsetI=[-0.002464537,0.003317852,-0.005837413,-0.001497292]; % Offset Current 05-Nov-2018 12:37:21
% 23: HW.Grad.PaOffsetI=[-0.001519135,-0.00123852,-0.00452556,0.0007421785]; % Offset Current 12-Nov-2018 13:18:49
HW.Grad.PaOffsetI=[sum([-0.002464537,-0.005837413]),sum([0.003317852,-0.001497292]), ... % for 22 (ch1&3 for x and ch2&4 for y)
  sum([-0.001519135,-0.00452556]),sum([-0.00123852,0.0007421785])]; % for 23 (ch1&3 for z and ch2&4 for slice)

% Input Voltage to output Current ratio
% 22: HW.Grad.PaUin2PaIout=([0.33094,0.335337,0.3272828,0.3297582]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio 05-Nov-2018 12:37:21
% 23: HW.Grad.PaUin2PaIout=([0.3319863,0.331572,0.3288497,0.3324137]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio 12-Nov-2018 13:18:49
HW.Grad.PaUin2PaIout=([sum([0.33094,0.3272828]),sum([0.335337,0.3297582]), ... % for 22 (ch1&3 for x and ch2&4 for y)
  sum([0.3319863,0.3288497]),sum([0.331572,0.3324137])]-HW.Grad.PaOffsetI)./1; % for 23 (ch1&3 for z and ch2&4 for slice)


HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Gradient and temperature ok of DC600

HW.Grad.PaPmaxInt=[100,100,100,100]*2;                  % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Output impedance

%% Gradient system

% ACHTUNG: Kopie in LoadSystem_Specific!
% Gradient efficiency must be reset here because it might have been changed when
% switching to the external power supply.
HW.Grad.LoadIin2Amp = [0.0303227, 0.0342713, 0.0334524, 0.102177];  % 2020-07-16T15:52:40, x y z B0  (T/(m*A)) Tesla per Meter per Ampere

% Zuordnung Gradienten Richtung und DAC Channel
HW.Grad.x = 1;
HW.Grad.y = 2;
HW.Grad.z = 3;
HW.Grad.B = 4; % slice gradient

% Zuordnung Gradienten Namen für Plot
HW.Grad.Name = {'Grad x','Grad y','Grad z','Grad slice'};
HW.Grad.AmpUnit = {'mT/m', 'mT/m', 'mT/m', 'mT/m'};
HW.Grad.AmpUnitScale = [1e-3, 1e-3, 1e-3, 1e-3];

HW.Grad.tRamp = 200e-6;  % minimum ramp time in s
HW.Grad.tEC = 100e-6;  % eddy current time in s
HW.Grad.MaxAmpSlice = 0.040;


% parameters for heat development model
HW.Grad.CoilThermalGroup = [1, 1, 1, 2];  % assignment of gradient coils to thermal groups
% Group 1: imaging gradient system
HW.Grad.CoilPowerDissipation(1) = 5;  % power dissipation at maximum temperature in Watt
HW.Grad.CoilTemperatur(1) = 20;  % resting temperature for each coil group in degrees C
HW.Grad.CoilMaxTemperature(1) = 60;  % maximum temperature for each coil group in degrees C
% Group 2: slice gradient
HW.Grad.CoilTemperatur(2) = 40;  % resting temperature for each coil group in degrees C
HW.Grad.CoilMaxTemperature(2) = 60;  % maximum temperature for each coil group in degrees C
% mass of copper in coils in kg
mass_x = 200e-3;
mass_y = 450e-3;
mass_z = 250e-3;
% mass_slice = 1.6;
heatcap_copper = 385; % heat capacity of copper in J/kg/K
% thermal capacity for each coil group in J/K
HW.Grad.CoilThermalCapacity(1) = (mass_x+mass_y+mass_z)*heatcap_copper;  % imaging gradients
% HW.Grad.CoilThermalCapacity(2) = mass_slice*heatcap_copper;  % slice gradient
% 1000 W * 10 s / 2.5 K (measured at gradient surface)
% Estimate that temperature of coils might rise twice as fast.
HW.Grad.CoilThermalCapacity(2) = 2000;  % slice gradient
% parameters for fuse blow model
HW.Grad.CoilMaxDcCurrent(1:3) = [2, 2, 2];  % maximum (nominal) DC current of fuse in Ampere
HW.Grad.CoilCurrentSquareTime(1:3) = [0.9, 0.9, 0.9] * 0.95;  % time-lag fuse parameter (maximum "accumulated heat") in A^2*sec

HW.Grad.PaUoutMax = [60, 60, 60, 64];
HW.Grad.PaUoutMin = -HW.Grad.PaUoutMax;

HW.FindPulseDuration.tPulse90Estimated = []; % use current PaUout2Amplitude as estimate
HW.FindPulseDuration.AQSlice.alfa  = 0.5*pi;
HW.FindPulseDuration.AQSlice.phi   = 0.0*pi;
HW.FindPulseDuration.AQSlice.theta = 0.0*pi;
HW.FindPulseDuration.AQSlice.HzPerPixMin = 600;

% Settings for slice gradient
HW.Grad.Slice.channel = 4; % channel that is used for the slice gradient
HW.Grad.Slice.tRamp = 250e-3; % ramp time for the slice gradient in s
HW.Grad.Slice.tEC = 30e-3; % time for eddy currents to decay in s
