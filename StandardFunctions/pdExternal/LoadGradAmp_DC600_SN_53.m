%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 53;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [0.001653053, -1.873647e-05, 0.00178375, 0.00171609];  % offset current in A  31-Aug-2022 09:46:56

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3347352, 0.3334587, 0.335439, 0.335004] - HW.Grad(iDevice).PaOffsetI) ./ 1;  % amplification in A/V  31-Aug-2022 09:46:56

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s

switch HW.UserName
  case {'probe_H1', 'teach'}
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [36.486, 49.410, 49.832]*1e-6;  % time delay of gradient amplifier in s
  case 'probe_H1C13'
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [40.720, 54.550, 50.652]*1e-6;  % time delay of gradient amplifier in s
  case 'probe_H1N15'
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [40.378, 54.232, 51.268]*1e-6;  % time delay of gradient amplifier in s
  otherwise
    error('PD:LoadGradAmp:UnknownUser', ...
      'No settings for user "%s" specified in "%s"', HW.UserName, mfilename());
end

HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
