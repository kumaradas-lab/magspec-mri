%% gradient amplifier settings for Aspect system
% used with MMRT #172

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.Grad(iDevice).ExtGradSN = 10172;  % dummy serial number

HW.Grad(iDevice).PowerDown = 0;
HW.Grad(iDevice).PaEnable = 1;

HW.Grad(iDevice).PaCurrentControlled = [1, 1, 1, 1];  % If current controlled set 1, if Voltage controlled set 0
HW.Grad(iDevice).PaRin = [24e3, 24e3, 24e3, 24e3];  % INA137 input impedance

HW.Grad(iDevice).PaOffsetU = [0, 0, 0, 0];  % offset voltage in Volts
HW.Grad(iDevice).PaOffsetI = [0, 0, 0, 0];  % offset current in Amperes

HW.Grad(iDevice).PaUin2PaUout = [0, 0, 0, 0];  % voltage gain (for voltage controlled mode)
% input voltage to output current ratio (for current controlled mode)
HW.Grad(iDevice).PaUin2PaIout = ([0.1,0.1,0.1,0.1] - HW.Grad(iDevice).PaOffsetI)./1;  
HW.Grad(iDevice).PaPmaxInt = [1000, 1000, 1000, 1000];  % maximum internal power dissipation in Watts

HW.Grad(iDevice).PaRout = [15000, 15000, 15000, 15000];  % output impedance in Ohms

HW.Grad(iDevice).tRamp = 200e-6;  % minimum ramp time in seconds (used, e.g., in imaging sequences)
HW.Grad(iDevice).tEC = 200e-6;   % settling time of eddy currents in seconds

HW.Grad(iDevice).SystemTimeDelay(1:4) = [52e-06 50e-06 31e-06 50e-06];  % time delay of grad amp in seconds

% HW.Grad(iDevice).Status1 = 1;  % power supply ok of DC600
% HW.Grad(iDevice).Status2 = 1;  % gradient and temperature ok of DC600

% properties for "fuse blow" model
HW.Grad(iDevice).CoilMaxDcCurrent = [1, 1, 1, 1] ./ 10;  % current for which fuse starts to meld in A
HW.Grad(iDevice).CoilCurrentSquareTime = [1, 1, 1, 1];  % fuse characteristics (I2t rating) in A^2*sec


%% gradient system settings
% (Also set here because this amplifier is "coupled" with that gradient system.)

% extent (and location) of the image volume in meters
HW.Grad(iDevice).ImageVol = [-0.05, 0.05, -0.05, 0.05, -0.05, 0.05];  % FOV in m
HW.Grad(iDevice).ImageVolOffset = [0, 0, 0];  % offset of sample in m

% correlation of output channel to gradient direction
HW.Grad(iDevice).x = 1;
HW.Grad(iDevice).y = 2;
HW.Grad(iDevice).z = 3;
HW.Grad(iDevice).B = 4;  % B0 / unused
HW.Grad(iDevice).xyzBDir = [1, 1, 1, 1];  % x y z B0 current sign

% to be determined!
HW.Grad(iDevice).LoadRin = [1, 1, 1, 1];  % x y z B0 resistance of gradients

% gradient efficiency (to be determined!)
HW.Grad(iDevice).LoadIin2Amp = [0.3, 0.3, 0.3, 1];  % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
% HW.MagnetShim = [0,0,0,0];  % x y z in T/m and B0 in T

% reasonable parameters for auto shim
HW.FindShim.ShimStep = 0.2e-3*ones(1,4);
HW.FindShim.ShimStart = [0.000155185, -0.000442921, -0.000677680, 0];
