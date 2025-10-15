%% Combined gradient amplifiers DC-600 SN #30 (rack) and SN #19 (extern)
%
% amplifier SN #30: output 1&3 for x and output 2&4 for y
% amplifier SN #19: output 1&3 for z and output 2&4 for B0(?)
%
% Y-adapter for connection to "Ext. Grad" connector at console:
% DAC channel 1 to SN #30 input 1&3
% DAC channel 2 to SN #30 input 2&4
% DAC channel 3 to SN #19 input 1&3
% DAC channel 4 to SN #19 input 2&4
%
% Take care to turn on both(!) DC-600. Power-on detection is not reliable.

HW.Grad.ExtGradSN = 30019;
HW.Grad.ExtGradType = 'DC600';

HW.Grad.PowerDown = 0;
HW.Grad.PaEnable = 1;

HW.Grad.PaCurrentControlled = [1, 1, 1, 1];  % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad.PaRin = [24e3, 24e3, 24e3, 24e3]/2;  % INA137 input impedance in Ohm  (divided by 2 for parallel amplifiers)

HW.Grad.PaOffsetU = [0, 0, 0, 0];  % offset voltage in Volt
% #30: HW.Grad.PaOffsetI = [-0.000342823, 0.00331074, -0.0002937403, -0.002939028];  % offset current in Ampere - 14-May-2020 14:44:27
% #19: HW.Grad.PaOffsetI = [-0.004124942, -0.002204735, 0.0006412072, 0.0004767803];  % offset current in Ampere - 25-Oct-2018 16:24:38
HW.Grad.PaOffsetI = ...
  [sum([-0.000342823,-0.0002937403]), sum([0.00331074,-0.002939028]), ...  % for #30 (ch1&3 for x and ch2&4 for y)
   sum([-0.004124942,0.0006412072]), sum([-0.002204735,0.0004767803])];   % for #19 (ch1&3 for z and ch2&4 for B0)

% amplification factor from input voltage to output current
% #30: HW.Grad.PaUin2PaIout = ([0.3317603, 0.3349532, 0.3327892, 0.329336] - HW.Grad.PaOffsetI)./1;  % Input Voltage to output Current ratio 14-May-2020 14:44:27
% #19: HW.Grad.PaUin2PaIout = ([0.3292542, 0.3304975, 0.3332017, 0.3328445] - HW.Grad.PaOffsetI)./1;  % Input Voltage to output Current ratio 25-Oct-2018 16:24:38
HW.Grad.PaUin2PaIout = ...
  ([sum([0.3317603,0.3327892]), sum([0.3349532,0.3293360]), ...  % for #30 (ch1&3 for x and ch2&4 for y)
    sum([0.3292542,0.3332017]), sum([0.3304975,0.3328445])] ...  % for #19 (ch1&3 for z and ch2&4 for B0)
   - HW.Grad.PaOffsetI)./1;

HW.Grad.PaPmaxInt = [100, 100, 100, 100]*2;  % maximum internal power dissipation (times 2 for parallel amplifiers)

HW.Grad.PaRout = [15000, 15000, 15000, 15000];  % output impedance in Ohm

% FIXME: Adapt to match connected gradient system
HW.Grad.tRamp = 50e-6;  % minimum ramp time in s
HW.Grad.tEC = 50e-6;  % setting time in s (for eddy currents)
HW.Grad.SystemTimeDelay(1:3) = [50e-6, 50e-6, 50e-6];  % [x y z] time delay of gradient amplifier in s
HW.Grad.MaxAmpSlice = 0.1;  % maximum gradient amplitude in T/m for slice selective pulses

HW.Grad.Status1 = 1;  % power supply ok of DC600
HW.Grad.Status2 = 1;  % gradient and temperature ok of DC600
