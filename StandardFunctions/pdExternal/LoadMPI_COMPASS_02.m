%% calibration properties and settings for the probe COMPASS #2

HW.MPI.probeType = 'COMPASS';
HW.MPI.probeSN = 2;

if ~exist('iDevice', 'var'), iDevice = 1; end

% active gradiometer
HW.TX(iDevice).PaUout2Amplitude(1) = 1;  % output in Volts
HW.TX(iDevice).PaUout2Amplitude(2) = 1;  % output in Volts

HW.MPI.identifier = 'COMPASS_V2';

% limits for sequence parameters
HW.MPI.ACFrequency_min = 1000;  % minimum frequency in Hz
HW.MPI.ACFrequency_max = 50000;  % maximum frequency in Hz
HW.MPI.DC_max = 25e-3;  % maximum DC amplitude in T
HW.MPI.AC_max = 25e-3;  % maximum AC amplitude in T
HW.MPI.DutyCycle_max = 0.75;

% set shim to 0 on all channels
shimIdx = sum([HW.Grad(1:(iDevice-1)).n]) + (1:HW.Grad(iDevice).n);
HW.MagnetShim(shimIdx) = zeros(1, HW.Grad(iDevice).n);

% calibration
% step 0: set coil efficiency, e.g. from simulation
HW.MPI.COMPASS_CoilEfficiency = 1.25e-3;  % efficiency of the TX coil in T/A

% step 1: run DC-only sequence and measure voltage from preview and
% current from current-probe to calculate LoadRIn = 0.652 V/(1.4291 A/2)
HW.TX(iDevice).LoadRin(2) = 3.65/(8.15/2);  % resistance of COMPASS TX coil in Ohm

% step 2: run AC-only sequence with frequency range 1kHz to
% 25kHz, e.g. 5 steps, AC, DC steps = 1, no DC-value to set up inductance!
% -> all amplitudes should show the same value
HW.TX(iDevice).Inductance(2) = 61e-6;  % inductance of COMPASS TX coil in Henry

% step 3: set isCalibratedCH5=false
%         set ACmax and DCmac to I_peak_exp=5e-3 [T]
%         perform measurement and get I_peak values of AC
%         and DC max ! @10kHz
% step 4: set isCalibratedCH5=true
isCalibratedCH5 = true;
if ~isCalibratedCH5
  HW.MPI.COMPASS_ACcal = 1;  % calibration factor
  HW.MPI.COMPASS_DCcal = 1;  % calibration factor
else
  I_peak_exp = 5e-3 / HW.MPI.COMPASS_CoilEfficiency;  % [A]
  I_peak_ACmeas = 8.3/2;  % [A]
  I_peak_DCmeas = 8.3/2;  % [A]
  HW.MPI.COMPASS_ACcal = 1 * I_peak_exp / I_peak_ACmeas;  % calibration factor
  HW.MPI.COMPASS_DCcal = 1 * I_peak_exp / I_peak_DCmeas;  % calibration factor
end

% step 4: get the current for specific frequencies @5mT
TXACcal=[1,   1000, 7.98/2];
TXACcal=[2,   2000, 7.95/2];
TXACcal=[3,   5000, 7.93/2];
TXACcal=[4,   7500, 7.96/2];
TXACcal=[5,  10000, 8.00/2];
TXACcal=[6,  15000, 8.06/2];
TXACcal=[7,  20000, 8.07/2];
TXACcal=[8,  25000, 8.00/2];
TXACcal=[9,  30000, 7.93/2];
TXACcal=[10, 40000, 7.65/2];

HW.TX(iDevice).AmplitudeName = '5:sol';
HW.TX(iDevice).AmplitudeUnit = 'V';
HW.TX(iDevice).AmplitudeUnitScale = 1;

% parameters for heat development model
% FIXME: Set reasonable parameters
HW.TX(iDevice).CoilTemperature = 22;  % resting temperature of coil in degrees C
HW.TX(iDevice).CoilMaxTemperature = 60;  % maximum temperature of coil before damage in degrees C
massCopper = 26.7e-3;  % mass of solenoid from simulation
specificHeatCapacity = 385;  % J/kg/K (copper wire)
HW.TX(iDevice).CoilThermalCapacity = massCopper*specificHeatCapacity;  % thermal capacity of coil in J/K
HW.TX(iDevice).CoilPowerDissipation = 60;  % thermal power dissipation in W at maximum temperature
% HW.TX(iDevice).Rout = HW.TX(iDevice).LoadRin(2);  % effective resistance at resonant frequency (kf3) in Ohm
