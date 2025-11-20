%% calibration properties and settings for the probe TWMPI #5 (large bore)

HW.MPI.probeType = 'TWMPI';
HW.MPI.probeSN = 5;

if ~exist('iDevice', 'var'), iDevice = 1; end

% active gradiometer
HW.TX(iDevice).PaUout2Amplitude(1) = 1;  % output in Volts

% MPI pulse program settings
% encoding periods
HW.MPI.kf1 = 1;  % number of periods during AQ for rotational frequency (f1, "envelope" Grad)
HW.MPI.kf2 = 13;  % number of periods during AQ for travel frequency (f2, "fast" Grad)
HW.MPI.kf3 = 175;  % number of periods during AQ for swipe frequency (f3, TX)
% phase correction (device dependent)
HW.MPI.phaseCH14 = deg2rad(20);  % 16-01-2025 - additional phase correction in rad for channels 1-4 (gradient channels)
HW.MPI.phaseCH5 = deg2rad(-42);  % 16-01-2025 - additional phase correction in rad for channel 5 (rf channel)
% other encoding parameters
HW.MPI.swipeRadius = 19.0e-3;  % maximum displacement of the FFL from the center by the swipe frequency in m

HW.MPI.tSampleGrad = 24e-6;

HW.MPI.RX_Uin_max = 50e-3;

% limits for sequence parameters
HW.MPI.ampGrad_max = [ 1, 1, 1, 1 ];  % [T/m]
HW.MPI.ampSwipe_max = 1.3;  % [T/m]
HW.MPI.kf1_min = 0;  % [periods/AQ]
HW.MPI.kf1_max = 25;  % [periods/AQ]
HW.MPI.kf2_min = 1;  % [periods/AQ]
HW.MPI.kf2_max = 25;  % [periods/AQ]
HW.MPI.kf3_min = 1;%175;  % [periods/AQ]
HW.MPI.kf3_max = 400;%175;  % [periods/AQ]

% gradient coils (channels 1-4)
HW.Grad(iDevice).Name = {'1:DHC3y', '2:DHC2y', '3:DHC3x', '4:DHC2x'};
HW.Grad(iDevice).AmpUnit(1:4) = {'T/m'};
HW.Grad(iDevice).AmpUnitScale(1:4) = 1;

% if 1  % remove this block when calibrating next time
% % FIXME: Measure resistance of coils
% HW.Grad(iDevice).LoadRin(1:4) = [1.3, 0.7, 1.5, 0.9];  % resistance of coils - approx. 1.5 - 2.5 Ohms
% ws = warning('off', 'PD:Grad:PaCurrentControlledChanged');  % temporarily disable warning
% HW.Grad(iDevice).AmpCurrentDependent(1:4) = 0;  % voltage dependent output
%
% % for calibration
% % step 1: set the output to a specific value
% % !! -> comment out the latter ...LoadUin2Amp
% % HW.Grad(iDevice).LoadUin2Amp = [0.05, 0.05, 0.05, 0.05];  % attention: 20V
% % step 2: define from simulation the target current
% %         measure with current probe the amplitude
% I1_1Tm = 23.3;     % [A] I_p from simu
% I2_1Tm = 23.3;     % [A] I_p from simu
% I3_1Tm = 23.3*1.26; % [A] I_p from simu
% I4_1Tm = 23.3*1.26; % [A] I_p from simu
% I1_1Tm_meas = 21.46/2; % I_p from measurement
% I2_1Tm_meas = 30.80/2; % I_p from measurement
% I3_1Tm_meas = 19.38/2; % I_p from measurement
% I4_1Tm_meas = 28.26/2; % I_p from measurement
% % step 3: use the calibrated value
% HW.Grad(iDevice).LoadUin2Amp = [0.05*I1_1Tm_meas/I1_1Tm, 0.05*I2_1Tm_meas/I2_1Tm, 0.05*I3_1Tm_meas/I3_1Tm, 0.05*I4_1Tm_meas/I4_1Tm];  % gradient coil efficiencies in T/m/V
%
% warning(ws);
% end

% ---------- use this block next time for calibration ----------------
% for calibration: Measure resistance of coils at frequency
I_1Tm = [23.3, 23.3, 23.3*1.26, 23.3*1.26];  % [A/(T/m)] I_p from simulation
LoadIin2Amp = 1./I_1Tm;  % gradient coil efficiencies in T/m/A from simulation
LoadRin = [1.3, 0.7, 1.5, 0.9];  % resistance of coils in Ohm from simulation
amp_set = 0.1*ones(1, 4);  % amplitude in T/m set for calibration measurement
isCalibratedGrad = true;  % set to true when the gradient coils are calibrated
if ~isCalibratedGrad
  % Use resistance as simulated
  HW.Grad(iDevice).LoadRin(1:4) = LoadRin;
else
  % measure peak current in Ampere with current probe when amplitude
  % is set to amp_set
  % calibrated on TODO (date):
  I_meas(1) = 4.6/2;  % I_p @ 1 T/m from measurement CH1
  I_meas(2) = 11.3/2;  % I_p @ 1 T/m from measurement CH2
  I_meas(3) = 4.9/2;  % I_p @ 1 T/m from measurement CH3
  I_meas(4) = 10.9/2;  % I_p @ 1 T/m from measurement CH4

  HW.Grad(iDevice).LoadRin(1:4) = LoadRin .* I_meas ./ amp_set .* LoadIin2Amp;  % resistance of coils in Ohm at frequency
  % sanity check:
  % The new LoadRin should be only slightly larger than the LoadRin
  % from the simulation.
end

ws = warning('off', 'PD:Grad:AmpCurrentDependentChanged');  % temporarily disable warning
HW.Grad(iDevice).AmpCurrentDependent(1:4) = 0;  % voltage dependent output
HW.Grad(iDevice).LoadUin2Amp = HW.Grad(iDevice).LoadRin .* LoadIin2Amp;  % gradient coil efficiencies in T/m/V
warning(ws);
% --------------------------------------------------------------------

% FIXME: Set reasonable values for the following parameters
% model parameters for fuse blow estimation
HW.Grad(iDevice).CoilMaxDcCurrent = [100, 100, 100, 100];  % maximum (nominal) DC current of fuse in Ampere
HW.Grad(iDevice).CoilCurrentSquareTime = [100, 100, 100, 100];  % time-lag fuse parameter (maximum "accumulated heat") in A^2*sec
% parameters for heat development model
HW.Grad(iDevice).CoilThermalGroup = [1, 2, 3, 4];  % assignment of gradient coils to thermal groups
HW.Grad(iDevice).CoilTemperatur(1:4) = 22;  % resting temperature for each coil group in degrees C
HW.Grad(iDevice).CoilMaxTemperature(1:4) = 60;  % maximum temperature for each coil group in degrees C
massCopper = [149.5, 99.7, 165.9, 110.6]*1e-3;  % mass of gradient coils from simulation
specificHeatCapacity = 385;  % J/kg/K (copper wire)
HW.Grad(iDevice).CoilThermalCapacity(1:4) = massCopper*specificHeatCapacity;  % thermal capacity of gradient coils in J/K
HW.Grad(iDevice).CoilPowerDissipation(1:4) = 0.1;  % power dissipation at maximum temperature in Watt

% set shim to 0 on all channels
shimIdx = sum([HW.Grad(1:(iDevice-1)).n]) + (1:HW.Grad(iDevice).n);
HW.MagnetShim(shimIdx) = zeros(1, HW.Grad(iDevice).n);


% rf coil (channel 5)
% factor added on 04-Jul-2024 to account for the match difference of TX2 and amplifier input
TX_PaRin = 5e3;  % impedance of rf amplifier at input in Ohm
TX_UoutRout = 50;  % impedance of TX2 port of MMRT in Ohm
matchFactor = 2 * TX_PaRin/(TX_PaRin+TX_UoutRout);  % factor compared to 50 Ohm termination

% for calibration
num_layers = 2;  % number of (identical) layers of coil
I5_1Tm = 23.3/num_layers;  % [A/(T/m)] I_p from simulation : target current - value from simulation
PaIout2Amplitude = 1/I5_1Tm*HW.MPI.swipeRadius;  % coil efficiency in T/A (as simulated)
isCalibratedCH5 = true;  % set to true when the CH5 coil is calibrated
amp_set = 0.1;  % amplitude in T/m set for calibration measurement
if ~isCalibratedCH5
  % step 1: set coil efficiency as simulated
  HW.TX(iDevice).PaUout2Amplitude(2) = 1 / matchFactor * PaIout2Amplitude;  % rf coil efficiency in T/V
else
  % step 2: measure I5_1Tm_meas
  % calibrated on 08-08-2024:
  I5_1Tm_meas = 8.2/2;  % [A] I_p measured with current probe
  HW.TX(iDevice).PaUout2Amplitude(2) = I5_1Tm_meas/I5_1Tm/amp_set / matchFactor * PaIout2Amplitude;  % rf coil efficiency in T/V
end

HW.TX(iDevice).AmplitudeName = '5:sol';
HW.TX(iDevice).AmplitudeUnit = 'mT';
HW.TX(iDevice).AmplitudeUnitScale = 1e-3;

% parameters for heat development model
% FIXME: Set reasonable parameters
HW.TX(iDevice).CoilTemperature = 22;  % resting temperature of coil in degrees C
HW.TX(iDevice).CoilMaxTemperature = 60;  % maximum temperature of coil before damage in degrees C
massCopper = 229.0e-3;  % mass of solenoid from simulation
specificHeatCapacity = 385;  % J/kg/K (copper wire)
HW.TX(iDevice).CoilThermalCapacity = massCopper*specificHeatCapacity;  % thermal capacity of coil in J/K
HW.TX(iDevice).CoilPowerDissipation = 2;  % thermal power dissipation in W at maximum temperature
HW.TX(iDevice).LoadRin(2) = 2.0;  % effective resistance at resonant frequency (kf3) in Ohm
