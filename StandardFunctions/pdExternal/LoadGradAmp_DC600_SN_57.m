%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 57;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [-0.001044605, -0.003801982, -0.001567125, 0.0002391333];  % offset current in A  26-Jan-2024 08:51:53

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3320125, 0.3298305, 0.3321587, 0.3335088] - HW.Grad(iDevice).PaOffsetI) ./ 1;  % amplification in A/V  26-Jan-2024 08:51:53

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s

switch HW.UserName
  case {'probe_10mm_H1', 'teach'}
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [32.770, 42.710, 39.710]*1e-6;  % time delay of gradient amplifier in s
  case 'probe_15mm_H1'
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [61.128, 89.308, 102.806]*1e-6;  % time delay of gradient amplifier in s
  case 'probe_MRE'
    HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [42.826, 40.330, 32.022]*1e-6;  % time delay of gradient amplifier in s
    HW.Grad(iDevice).tRamp = 42e-6;  % minimum ramp time in s
  otherwise
    error('PD:LoadGradAmp:UnknownUser', ...
      'No settings for user "%s" specified in "%s"', HW.UserName, mfilename());
end

HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok
