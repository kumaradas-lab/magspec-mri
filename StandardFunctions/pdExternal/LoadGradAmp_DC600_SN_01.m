%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 1;                               % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 0];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, -0.0034];         % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [0.0041773-450e-6, 0.0033587, 0.0054375, 0];  % offset current in A

HW.Grad(iDevice).PaUin2PaUout(1:4) = [0, 0, 0, 4];            % Voltage gain in V/V
HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.33827*996.41934e-003, 0.33788, 0.33792, 0] - HW.Grad(iDevice).PaOffsetI(1:4)) ./ 1;  % amplification in A/V
HW.Grad(iDevice).PaPmaxInt(1:4) = [100, 100, 90, 60];         % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 1.1];    % output impedance in Ohm

HW.Grad(iDevice).tRamp = 18e-6*4;                             % minimum ramp time in s
HW.Grad(iDevice).tEC = 36e-6*2;                               % eddy current time in s
% HW.Grad(iDevice).SystemTimeDelay(1:4) = [4.8208e-05, 5.6208e-05, 5.452e-05, 1.95e-05];  % time delay of gradient amplifier in s
HW.Grad(iDevice).SystemTimeDelay(1:3) = 00e-6;                % time delay of gradient amplifier in s
% HW.Grad(iDevice).SystemTimeDelay(1:3) =[4.1904e-05, 6.1216e-05, 6.1536e-05];  % time delay of gradient amplifier in s
% HW.Grad(iDevice).SystemTimeDelay(1:3) =[9.6288e-05, 0.000112736, 0.00011904];  % time delay of gradient amplifier in s
HW.Grad(iDevice).SystemTimeDelay(1:3) = [4.6296e-05, 5.6568e-05, 5.9056e-05];  % time delay of gradient amplifier in s
HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok

% properties for "fuse blow" model
HW.Grad(iDevice).CoilMaxDcCurrent(1:4) = [0.75,0.75,0.75,0.75];  % current for which fuse starts to meld in A
HW.Grad(iDevice).CoilCurrentSquareTime(1:4) = [0.9,0.9,0.9,0.9];  % fuse characteristics (I2t rating) in A^2*sec
