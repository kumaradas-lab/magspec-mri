%% calibration properties and settings for the probe iMPI #1

HW.MPI.probeType = 'iMPI';
HW.MPI.probeSN = 1;

if ~exist('iDevice', 'var'), iDevice = 1; end

% active gradiometer
HW.TX(iDevice).PaUout2Amplitude(1) = 1;  % output in Volts

% MPI pulse program settings
% encoding periods

HW.MPI.kf1 = 0;  % number of periods during AQ for rotational frequency (f1, "envelope" Grad)
HW.MPI.kf2 = 3;  % number of periods during AQ for travel frequency (f2, "fast" Grad)
HW.MPI.kf3 = 132; % number of periods during AQ for swipe frequency (f3, TX)
% phase correction (device dependent)
HW.MPI.phaseCH14 = deg2rad(0);  % TODO - additional phase correction in rad for channels 1-4 (gradient channels)
HW.MPI.phaseCH5 = deg2rad(20);  % TODO - additional phase correction in rad for channel 5 (rf channel)
% other encoding parameters
HW.MPI.swipeRadius = 130e-3;  % maximum displacement of the FFL from the center by the swipe frequency in m

% limits for sequence parameters
HW.MPI.ampGrad_max = [ 0.25, 0.25, 0.25, 0.25 ];  % [T/m]
HW.MPI.ampSwipe_max = 0.25;  % [T/m]
HW.MPI.kf1_min = 0;  % [periods/AQ]
HW.MPI.kf1_max = 0;  % [periods/AQ]
HW.MPI.kf2_min = 1;  % [periods/AQ]
HW.MPI.kf2_max = 9;  % [periods/AQ]
HW.MPI.kf3_min = 100;  % [periods/AQ]
HW.MPI.kf3_max = 140;  % [periods/AQ]

% device specific settings
HW.MPI.tRep = 0.07;
HW.MPI.tAQ = 0.05;
HW.MPI.tSampleGrad = 72e-6;
HW.MPI.tRampUpGrad = 10e-3;
HW.MPI.tRampDownGrad = 10e-3;
HW.MPI.tRampUpHoldTX = 30e-3;
HW.MPI.RX_Uin_max = 20e-3;
HW.MPI.nPhaseTX = 1;

% gradient coils (channels 1-2)
HW.Grad(iDevice).Name = {'1:DHC1', '2:DHC2', '3:-', '4:-'};
HW.Grad(iDevice).AmpUnit(1:4) = {'T/m'};
HW.Grad(iDevice).AmpUnitScale(1:4) = 1;

% for calibration: Measure resistance of coils at frequency
I_1Tm = [500, 500, 500, 500];  % [A/(T/m)] I_p from simulation
LoadIin2Amp = 1./I_1Tm;  % gradient coil efficiencies in T/m/A from simulation
LoadRin = [0.5, 0.5, 0.5, 0.5];  % resistance of coils in Ohm from simulation
amp_set = 0.01*ones(1, 4);  % amplitude in T/m set for calibration measurement
isCalibratedGrad = true;  % set to true when the gradient coils are calibrated
if ~isCalibratedGrad
  % Use resistance as simulated
  HW.Grad(iDevice).LoadRin(1:4) = LoadRin;
else
  % measure peak current in Ampere with current probe when amplitude
  % is set to amp_set
  % calibrated on TODO (date):
  I_meas(1) = 36.5/2;  % I_p @ 0.01 T/m from measurement CH1
  I_meas(2) = 36.7/2;  % I_p @ 0.01 T/m from measurement CH2
  I_meas(3) = 1/2;  % I_p @ 0.01 T/m from measurement CH3
  I_meas(4) = 1/2;  % I_p @ 0.01 T/m from measurement CH4

  HW.Grad(iDevice).LoadRin(1:4) = LoadRin .* I_meas ./ amp_set .* LoadIin2Amp;  % resistance of coils in Ohm at frequency
  % HW.Grad(iDevice).LoadRin(1:4) = [1.8, 1.8, 1.8, 1.8];
  % sanity check:
  % The new LoadRin should be only slightly larger than the LoadRin
  % from the simulation.
end

ws = warning('off', 'PD:Grad:AmpCurrentDependentChanged');  % temporarily disable warning
HW.Grad(iDevice).AmpCurrentDependent(1:4) = 0;  % voltage dependent output
HW.Grad(iDevice).LoadUin2Amp = HW.Grad(iDevice).LoadRin .* LoadIin2Amp;  % gradient coil efficiencies in T/m/V
warning(ws);
%--------------------------------------------------------------------

% FIXME: Set reasonable values for the following parameters
% model parameters for fuse blow estimation
HW.Grad(iDevice).CoilMaxDcCurrent = [1000, 1000, 1000, 1000];  % maximum (nominal) DC current of fuse in Ampere
HW.Grad(iDevice).CoilCurrentSquareTime = [1000, 1000, 1000, 1000];  % time-lag fuse parameter (maximum "accumulated heat") in A^2*sec
% parameters for heat development model
HW.Grad(iDevice).CoilThermalGroup = [1, 2, 3, 4];  % assignment of gradient coils to thermal groups
HW.Grad(iDevice).CoilTemperatur(1:4) = 220;  % resting temperature for each coil group in degrees C
HW.Grad(iDevice).CoilMaxTemperature(1:4) = 600;  % maximum temperature for each coil group in degrees C
massCopper = [3000, 3000, 3000, 3000]*1e-3;  % mass of gradient coils from simulation
specificHeatCapacity = 385;  % J/kg/K (copper wire)
HW.Grad(iDevice).CoilThermalCapacity(1:4) = massCopper*specificHeatCapacity;  % thermal capacity of gradient coils in J/K
HW.Grad(iDevice).CoilPowerDissipation(1:4) = 0.1;  % power dissipation at maximum temperature in Watt
% set shim to 0 on all channels
shimIdx = sum([HW.Grad(1:(iDevice-1)).n]) + (1:HW.Grad(iDevice).n);
HW.MagnetShim(shimIdx) = zeros(1, HW.Grad(iDevice).n);


% rf coil (channel 5)
% factor added on 04-Jul-2024 to account for the match difference of TX2 and amplifier input
TX_PaRin = 0.6;  % impedance of rf amplifier at input in Ohm
TX_UoutRout = 50;  % impedance of TX2 port of MMRT in Ohm
matchFactor = 2 * TX_PaRin/(TX_PaRin+TX_UoutRout);  % factor compared to 50 Ohm termination

% for calibration
num_layers = 1;  % number of (identical) layers of coil
I5_1Tm = 400/num_layers;  % [A/(T/m)] I_p from simulation : target current - value from simulation
PaIout2Amplitude = 1/I5_1Tm*HW.MPI.swipeRadius;  % coil efficiency in T/A (as simulated)
isCalibratedCH5 = true;  % set to true when the CH5 coil is calibrated
amp_set = 0.01;  % amplitude in T/m set for calibration measurement
if ~isCalibratedCH5
  % step 1: set coil efficiency as simulated
  HW.TX(iDevice).PaUout2Amplitude(2) = 1 / matchFactor * PaIout2Amplitude;  % rf coil efficiency in T/V
else
  % step 2: measure I5_1Tm_meas
  % calibrated on 05-07-2024:
  I5_1Tm_meas = 1.0/2;  % [A/(T/m)] I_p measured with current probe
  HW.TX(iDevice).PaUout2Amplitude(2) = I5_1Tm_meas/I5_1Tm / amp_set / matchFactor * PaIout2Amplitude;  % rf coil efficiency in T/V
end

HW.TX(iDevice).AmplitudeName = '3:sol';
HW.TX(iDevice).AmplitudeUnit = 'mT';
HW.TX(iDevice).AmplitudeUnitScale = 1e-3;

% parameters for heat development model
% FIXME: Set reasonable parameters
HW.TX(iDevice).CoilTemperature = 220;  % resting temperature of coil in degrees C
HW.TX(iDevice).CoilMaxTemperature = 600;  % maximum temperature of coil before damage in degrees C
massCopper = 470.8e-3;  % mass of solenoid from simulation
specificHeatCapacity = 385;  % J/kg/K (copper wire)
HW.TX(iDevice).CoilThermalCapacity = massCopper*specificHeatCapacity;  % thermal capacity of coil in J/K
HW.TX(iDevice).CoilPowerDissipation = 1;  % thermal power dissipation in W at maximum temperature
HW.TX(iDevice).LoadRin(2) = 1.5;  % effective resistance at resonant frequency (kf3) in Ohm
