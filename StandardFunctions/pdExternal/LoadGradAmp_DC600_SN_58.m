%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 58;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [0.0007468187, -0.001478997, -0.0004166522, -0.0003359438];  % offset current in A  26-Jan-2024 09:28:58

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3339077, 0.3317162, 0.3341623, 0.3330982] - HW.Grad(iDevice).PaOffsetI) ./ 1;  % amplification in A/V  26-Jan-2024 09:28:58

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s

switch HW.UserName
  case {'probe_10mm_H1', 'teach'}
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [35.134, 41.592, 41.674]*1e-6;  % time delay of gradient amplifier in s
  case 'probe_15mm_H1'
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [60.158, 90.810, 101.916]*1e-6;  % time delay of gradient amplifier in s
  otherwise
    error('PD:LoadGradAmp:UnknownUser', ...
      'No settings for user "%s" specified in "%s"', HW.UserName, mfilename());
end

HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
