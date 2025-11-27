%% Calibration values for gradient amplifier

% all vector properties are in gradient output channel ordering

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 9999;                            % serial number
HW.Grad(iDevice).ExtGradType = 'XYZ';                         % type of external gradient amplifier

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % input impedance in Ohm

% current controlled
HW.Grad(iDevice).PaOffsetI(1:4) = [0, 0, 0, 0];               % offset current in A
HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.33, 0.33, 0.33, 0.33] - HW.Grad(iDevice).PaOffsetI(1:4)) ./ 1;  % amplification in A/V

% voltage controlled
HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
% HW.Grad(iDevice).PaUin2PaUout(1:4) = ([2, 2, 2, 2] - HW.Grad(iDevice).PaOffsetU) ./ 1;  % amplification in V/V

HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 100, 100];       % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 50e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s

% check if Status pins of "Ext. Grad" connector are high
% Status1 pin at high level indicates that the power supply is ok.
% Status2 pin at high level indicates that there is no general error.
% The values of the following properties can be set to 0, 1, or 2.
% If set to 0, the corresponding status pin is not checked.
% If set to 1, the corresponding status pin is checked immediately before
% starting the measurement and while polling the measurement data.
% If set to 2, the corresponding status pin is checked while polling the
% measurement data (but not immediately before starting the measurement).
HW.Grad(iDevice).Status1 = 0;
HW.Grad(iDevice).Status2 = 0;


%% Might also be set in gradient system configuration file

% time delay of gradient amplifier and gradient system
HW.Grad(iDevice).SystemTimeDelay(1:4) = [50e-6, 50e-6, 50e-6, 50e-6];  % in s
HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m
