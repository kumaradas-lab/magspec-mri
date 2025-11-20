%% calibration properties and settings for the probe TWMPI #1
% portable TWMPI device (DC-600)

HW.MPI.probeType = 'TWMPI';
HW.MPI.probeSN = 1;

if ~exist('iDevice', 'var'), iDevice = 1; end

% active gradiometer
HW.TX(iDevice).PaUout2Amplitude(1) = 1;  % output in Volts

% MPI pulse program settings
% encoding periods
HW.MPI.kf1 = 1;  % number of periods during AQ for rotational frequency (f1, "envelope" Grad)
HW.MPI.kf2 = 11;  % number of periods during AQ for travel frequency (f2, "fast" Grad)
HW.MPI.kf3 = 111;  % number of periods during AQ for swipe frequency (f3, TX)
% phase correction (device dependent)
HW.MPI.phaseCH14 = deg2rad(35);  % 10-07-2024 - additional phase correction in rad for channels 1-4 (gradient channels)
HW.MPI.phaseCH5 = deg2rad(135);  % 10-07-2024 - additional phase correction in rad for channel 5 (rf channel)
% other encoding parameters
HW.MPI.swipeRadius = 8.6e-3;  % maximum displacement of the FFL from the center by the swipe frequency in m

% limits for sequence parameters
HW.MPI.ampGrad_max = [ 1, 1, 1, 1 ];  % [T/m]
HW.MPI.ampSwipe_max = 1;  % [T/m]
HW.MPI.kf1_min = 1;  % [periods/AQ]
HW.MPI.kf1_max = 200;  % [periods/AQ]
HW.MPI.kf2_min = 1;  % [periods/AQ]
HW.MPI.kf2_max = 200;  % [periods/AQ]
HW.MPI.kf3_min = 1;  % [periods/AQ]
HW.MPI.kf3_max = 200;  % [periods/AQ]

HW.MPI.tSampleGrad = 24e-6;

% gradient coils (channels 1-4)
HW.Grad(iDevice).Name = {'1:DHC3y', '2:DHC2y', '3:DHC3x', '4:DHC2x'};
HW.Grad(iDevice).AmpUnit(1:4) = {'T/m'};
HW.Grad(iDevice).AmpUnitScale(1:4) = 1;


% Calibration is only necessary if DC-600 amplification is different
% with sinusoidal output.
I_1Tm = [5.9, 5.9, 5.9*1.25, 5.9*1.25];  % [A/(T/m)] I_p from simulation (1 strand!!!)
isCalibratedGrad = true;
ampSet = 0.2;  % amplitude in T/m set for calibration measurement
if ~isCalibratedGrad
  calFactor = 1;
else
  I_meas(1) = 2.35/2;  % [A] I_p measured with current probe in CH1
  I_meas(2) = 2.36/2;  % [A] I_p measured with current probe in CH2
  I_meas(3) = 2.95/2;  % [A] I_p measured with current probe in CH3
  I_meas(4) = 2.94/2;  % [A] I_p measured with current probe in CH4
  calFactor = I_1Tm*ampSet ./ I_meas;
end
num_strands = 1;  % number of strands in series
HW.Grad(iDevice).LoadIin2Amp(1:4) = num_strands .* calFactor ./ I_1Tm;  % gradient coil efficiencies in T/m/A from simulation

% FIXME: Measure resistance of coils
HW.Grad(iDevice).LoadRin(1:4) = [4.7, 3.2, 5.2, 3.5];  % resistance of coils in Ohm from simulation
% FIXME: Measure inductance of coils
% HW.Grad(iDevice).Inductance = [];

% FIXME: Set reasonable values for the following parameters
% model parameters for fuse blow estimation
HW.Grad(iDevice).CoilMaxDcCurrent = [1, 1, 1, 1] * 100 * 0.9;  % maximum (nominal) DC current of fuse in Ampere
HW.Grad(iDevice).CoilCurrentSquareTime = [0.9, 0.9, 0.9, 0.9] * 100 * 0.9;  % time-lag fuse parameter (maximum "accumulated heat") in A^2*sec
% parameters for heat development model
HW.Grad(iDevice).CoilThermalGroup = [1, 1, 1, 1];  % assignment of gradient coils to thermal groups
HW.Grad(iDevice).CoilTemperatur(1) = 22;  % resting temperature for each coil group in degrees C
HW.Grad(iDevice).CoilMaxTemperature(1) = 60;  % maximum temperature for each coil group in degrees C
massCopper = sum([23.9, 16.0, 26.9, 17.9]*1e-3);  % mass of gradient coils in kg - copied from TWMPI #3
specificHeatCapacity = 385;  % J/kg/K (copper wire)
HW.Grad(iDevice).CoilThermalCapacity(1) = massCopper*specificHeatCapacity;  % thermal capacity of gradient coils in J/K
HW.Grad(iDevice).CoilPowerDissipation(1) = 2;  % power dissipation at maximum temperature in Watt

% set shim to 0 on all channels
shimIdx = sum([HW.Grad(1:(iDevice-1)).n]) + (1:HW.Grad(iDevice).n);
HW.MagnetShim(shimIdx) = zeros(1, HW.Grad(iDevice).n);


% rf coil (channel 5)

% Calibration is only necessary if DC-600 amplification is different
% with sinusoidal output or factor of differential converter is not
% exactly 2.
num_layers = 1;  % number of (identical) layers of coil
I5_1Tm = 5.9/num_layers;  % [A/(T/m)] I_p from simulation : target current
% FIXME: measure the resistance
LoadRin = 2.1;  % resistance of CH5 coil in Ohm from simulation
factorDifferential = 2;  % factor for differential converter
PaUin2PaIout = mean([0.3326805, 0.332308]);  % amplifier gain in A/V - copied from LoadGradAmp_DC600_SN_41 CH1&3
dualChannel = 2;  % factor for dual channel at DC-600 output
isCalibratedCH5 = true;
ampSet = 0.6*1.25;  % amplitude in T/m set for calibration measurement - 12-07-2024: factor 1.25 added after calibration to increase rf amplitude for swipe
if ~isCalibratedCH5
  calFactor = 1;
else
  I5_meas = 4.33;  % [A] I_p measured with current probe
  calFactor = I5_meas / ( ampSet * I5_1Tm );
end
Uout2PaIout = calFactor * factorDifferential * PaUin2PaIout * dualChannel;  % amplification in A/V
PaIout2Amplitude = 1/I5_1Tm * HW.MPI.swipeRadius;  % coil efficiency in T/A (as simulated)
% Use resistance as simulated
HW.TX(iDevice).LoadRin(2) = LoadRin;  % resistance of CH5 coil at frequency in Ohm
HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).LoadRin(2) * Uout2PaIout;  % amplification in V/V
HW.TX(iDevice).PaUout2Amplitude(2) = PaIout2Amplitude / HW.TX(iDevice).LoadRin(2);  % rf coil efficiency in T/V

HW.TX(iDevice).AmplitudeName = '5:sol';
HW.TX(iDevice).AmplitudeUnit = 'mT';
HW.TX(iDevice).AmplitudeUnitScale = 1e-3;
