if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 26;

HW.Grad(iDevice).PowerDown = 0;
HW.Grad(iDevice).PaEnable = 1;

% Grad Arrays
HW.Grad(iDevice).PaCurrentControlled = [1, 1, 1, 1];             % If current controlled set 1, if voltage controlled set 0

HW.Grad(iDevice).PaRin = [24e3, 24e3, 24e3, 24e3];               % INA137 input impedance

HW.Grad(iDevice).PaOffsetU = [0, 0, 0, 0];                       % Offset voltage
HW.Grad(iDevice).PaOffsetI = [0.001548337, 0.001442667, 0.0008353942, -0.00323112];  % Offset Current 13-May-2019 13:12:48

HW.Grad(iDevice).PaUin2PaIout = ([0.333947, 0.3341918, 0.3335547, 0.3284165] - HW.Grad(iDevice).PaOffsetI)./1;  % Input Voltage to output Current ratio% Input Voltage to output Current ratio 13-May-2019 13:12:48

HW.Grad(iDevice).PaPmaxInt = [100, 100, 100, 100];               % Maximum internal power dissipation

HW.Grad(iDevice).PaRout = [15000, 15000, 15000, 15000];          % Ouput impedance

HW.Grad(iDevice).tRamp = 50e-6;                                  % minimum ramp time in s
HW.Grad(iDevice).tEC = 50e-6;                                    % Setting time;
HW.Grad(iDevice).SystemTimeDelay(1:3) = [50e-6, 50e-6, 50e-6];   % Time delay of grad amp
HW.Grad(iDevice).MaxAmpSlice = 0.1;                              % Max Grad Amp ?

HW.Grad(iDevice).Status1 = 1;                                    % Power supply ok of DC600
HW.Grad(iDevice).Status2 = 1;                                    % Gradient and temperature ok of DC600
