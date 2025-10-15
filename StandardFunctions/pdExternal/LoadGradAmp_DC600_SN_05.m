%% Calibration values for gradient amplifier DC-600

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 5;                               % serial number
HW.Grad(iDevice).ExtGradType = 'DC600';

HW.Grad(iDevice).PowerDown = 0;                               % power down amplifier after some time (sleep)
HW.Grad(iDevice).PaEnable = 1;                                % un-mute the amplifier

HW.Grad(iDevice).PaCurrentControlled(1:4) = [1, 1, 1, 1];     % if current controlled, set to 1; if voltage controlled, set to 0

HW.Grad(iDevice).PaRin(1:4) = [24e3, 24e3, 24e3, 24e3];       % INA137 input impedance

HW.Grad(iDevice).PaOffsetU(1:4) = [0, 0, 0, 0];               % offset voltage in V
HW.Grad(iDevice).PaOffsetI(1:4) = [-0.0021933, 0.003974, 0.00016173, 0.0040016];  % offset current in A

HW.Grad(iDevice).PaUin2PaUout(1:4) = [0, 0, 0, 0];            % Voltage gain in V/V
HW.Grad(iDevice).PaUin2PaIout(1:4) = ([0.33108,0.33731,0.33135,0.33866] - HW.Grad(iDevice).PaOffsetI(1:4)) ./ 1;  % amplification in A/V

HW.Grad(iDevice).PaPmaxInt(1:4) = [80, 80, 80, 80];           % maximum internal power dissipation in W

HW.Grad(iDevice).PaRout(1:4) = [15000, 15000, 15000, 15000];  % output impedance in Ohm

HW.Grad(iDevice).tRamp = 18e-6;                               % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                 % eddy current time in s
HW.Grad(iDevice).SystemTimeDelay(HW.Grad.xyzB(1:3)) = [51.614, 48.260, 37.014]*1e-6;  % time delay of gradient amplifier in s
HW.Grad(iDevice).MaxAmpSlice = 0.1;                           % maximum gradient amplitude for slice selection in T/m

HW.Grad(iDevice).Status1 = 1;                                 % power supply of DC-600 ok
HW.Grad(iDevice).Status2 = 1;                                 % gradient and temperature of DC-600 ok

% Piezo
HW.Grad.LoadRin(4) = 100;  % resistance parallel to Piezo in Ohm
% name for Plot
HW.Grad(iDevice).Name(4) = {'Piezo Current'};
HW.Grad(iDevice).AmpUnit(4) = {'A'};
HW.Grad(iDevice).AmpUnitScale(4) = 1;
% efficiency
HW.Grad(iDevice).LoadIin2Amp(4) = 1;  % 1 A per Ampere
