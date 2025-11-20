%% calibration properties and settings for the probe TWMPI #4

HW.MPI.probeType = 'TWMPI';
HW.MPI.probeSN = 4;

if ~exist('iDevice', 'var'), iDevice = 1; end

% active gradiometer
HW.TX(iDevice).PaUout2Amplitude(1) = 1;  % output in Volts

% MPI pulse program settings
% encoding periods
HW.MPI.kf1 = 1;  % number of periods during AQ for rotational frequency (f1, "envelope" Grad)
HW.MPI.kf2 = 13;  % number of periods during AQ for travel frequency (f2, "fast" Grad)
HW.MPI.kf3 = 189;  % number of periods during AQ for swipe frequency (f3, TX)
% phase correction (device dependent)
HW.MPI.phaseCH14 = deg2rad(35);  % 25-06-2024 - additional phase correction in rad for channels 1-4 (gradient channels)
HW.MPI.phaseCH5 = deg2rad(80);  % 30-07-2024 - additional phase correction in rad for channel 5 (rf channel)
% other encoding parameters
HW.MPI.swipeRadius = 8.6e-3;  % maximum displacement of the FFL from the center by the swipe frequency in m

% limits for sequence parameters
HW.MPI.ampGrad_max = [ 3, 3, 3, 3 ];  % [T/m]
HW.MPI.ampSwipe_max = 10;  % [T/m]
HW.MPI.kf1_min = 1;  % [periods/AQ]
HW.MPI.kf1_max = 25;  % [periods/AQ]
HW.MPI.kf2_min = 1;  % [periods/AQ]
HW.MPI.kf2_max = 25;  % [periods/AQ]
HW.MPI.kf3_min = 189;  % [periods/AQ]
HW.MPI.kf3_max = 189;  % [periods/AQ]

HW.MPI.tSampleGrad = 24e-6;

% gradient coils (channels 1-4)
HW.Grad(iDevice).Name = {'1:DHC3y', '2:DHC2y', '3:DHC3x', '4:DHC2x'};
HW.Grad(iDevice).AmpUnit(1:4) = {'T/m'};
HW.Grad(iDevice).AmpUnitScale(1:4) = 1;
% FIXME: Measure resistance of coils
HW.Grad(iDevice).LoadRin(1:4) = [2.4, 1.6, 2.7, 1.8];  % resistance of coils - approx. 1.5 - 2.5 Ohms
ws = warning('off', 'PD:Grad:PaCurrentControlledChanged');  % temporarily disable warning
HW.Grad(iDevice).AmpCurrentDependent(1:4) = 0;  % voltage dependent output
HW.Grad(iDevice).LoadUin2Amp = [0.071, 0.105, 0.064/1.3, 0.093/1.3]*5/6.7;  % 25-05-2024 - gradient coil efficiencies in T/m/V
warning(ws);
% FIXME: Set reasonable values for the following parameters
% model parameters for fuse blow estimation
HW.Grad(iDevice).CoilMaxDcCurrent = [100, 100, 100, 100];  % maximum (nominal) DC current of fuse in Ampere
HW.Grad(iDevice).CoilCurrentSquareTime = [100, 100, 100, 100];  % time-lag fuse parameter (maximum "accumulated heat") in A^2*sec
% parameters for heat development model
HW.Grad(iDevice).CoilThermalGroup = [1, 2, 3, 4];  % assignment of gradient coils to thermal groups
HW.Grad(iDevice).CoilTemperatur(1:4) = 22;  % resting temperature for each coil group in degrees C
HW.Grad(iDevice).CoilMaxTemperature(1:4) = 60;  % maximum temperature for each coil group in degrees C
massCopper = [23.9, 16.0, 26.9, 17.9]*1e-3;  % mass of gradient coils from simulation
specificHeatCapacity = 385;  % J/kg/K (copper wire)
HW.Grad(iDevice).CoilThermalCapacity(1:4) = massCopper*specificHeatCapacity;  % thermal capacity of gradient coils in J/K
HW.Grad(iDevice).CoilPowerDissipation(1:4) = 0.1;  % power dissipation at maximum temperature in Watt
% set shim to 0 on all channels
shimIdx = sum([HW.Grad(1:(iDevice-1)).n]) + (1:HW.Grad(iDevice).n);
HW.MagnetShim(shimIdx) = zeros(1, HW.Grad(iDevice).n);

% rf coil (channel 5)
PaIout2Amplitude = 1/6.7*2*HW.MPI.swipeRadius;  % coil efficiency in T/A (as simulated)
% factor added on 04-Jul-2024 to account for the match difference of TX2 and amplifier input
TX_PaRin = 5e3;  % impedance of rf amplifier at input in Ohm
TX_UoutRout = 50;  % impedance of TX2 port of MMRT in Ohm
matchFactor = 2 * TX_PaRin/(TX_PaRin+TX_UoutRout);  % factor compared to 50 Ohm termination
voltageFor5Ampere = 4*matchFactor;  % 25-06-2024 - PaUout to reach 5 Ampere through the coil at resonant frequency (kf3) (measured with clamp-on ammeter)
HW.TX(iDevice).PaUout2Amplitude(2) = 5/voltageFor5Ampere * PaIout2Amplitude;  % rf coil efficiency in T/V
HW.TX(iDevice).AmplitudeName = '5:sol';
HW.TX(iDevice).AmplitudeUnit = 'mT';
HW.TX(iDevice).AmplitudeUnitScale = 1e-3;
% parameters for heat development model
% FIXME: Set reasonable parameters
HW.TX(iDevice).CoilTemperature = 22;  % resting temperature of coil in degrees C
HW.TX(iDevice).CoilMaxTemperature = 60;  % maximum temperature of coil before damage in degrees C
massCopper = 47.8e-3;  % mass of solenoid from simulation
specificHeatCapacity = 385;  % J/kg/K (copper wire)
HW.TX(iDevice).CoilThermalCapacity = massCopper*specificHeatCapacity;  % thermal capacity of coil in J/K
HW.TX(iDevice).CoilPowerDissipation = 1;  % thermal power dissipation in W at maximum temperature
HW.TX(iDevice).LoadRin(2) = voltageFor5Ampere/5;  % effective resistance at resonant frequency (kf3) in Ohm
