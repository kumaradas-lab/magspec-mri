%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 40;                              % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [-0.003558705, -0.0009927628, 0.0014262, -0.001224293];  % offset current in A  11-Feb-2022 08:54:35

HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.3297928, 0.3322115, 0.334648, 0.3324415] - HW.Grad(iDevice).PaOffsetI) ./ 1;  % amplification in A/V  11-Feb-2022 08:54:35

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s

switch HW.UserName
  case '10mm_single'
    % Magnet # 183
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [3.1112e-05, 4.3672e-05, 4.0624e-05];  % Time delay of grad amp
  case '15mm_single'
    % Magnet # 180
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [7.6704e-05, 0.000101776, 0.0001192];  % Time delay of grad amp
  case '10mm_both'
    % both connected - active: Magnet # 183
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [3.3792e-05, 4.9016e-05, 4.5424e-05];  % Time delay of grad amp
  case '15mm_both'
    % both connected - active: Magnet # 180
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [8.0728e-05, 0.000105616, 0.000122448];  % Time delay of grad amp
%   case '14.8mm_light_single'
%     % Magnet # 180
%     HW.Grad(iDevice).SystemTimeDelay(1:3) = [60.792, 86.778, 94.342]*1e-6;  % time delay of gradient amplifier in s
  case '15mm_light_single'
    % Magnet # 180
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [61.838, 89.026, 104.196]*1e-6;  % DC-600 time delay of gradient amplifier in s
  case {'15mm_light_single_b1', '15mm_light_single_b1_DC600', 'teach'}
    % - kopiert von '15mm_light_single'
    % Magnet # 180
    HW.Grad(iDevice).SystemTimeDelay(1:3) = [61.838, 89.026, 104.196]*1e-6;  % DC-600 time delay of gradient amplifier in s
  otherwise
    error('PD:LoadGradAmp:UnknownUser', ...
      'No settings for user "%s" specified in "%s"', HW.UserName, mfilename());
end

HW.Grad(iDevice).MaxAmpSlice = 0.1;  % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;  % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;  % gradient and temperature of DC-600 ok
