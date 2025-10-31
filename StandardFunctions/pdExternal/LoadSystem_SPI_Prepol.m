%% System settings for SPI prepolarization amplifier

iDevice = 1;  % device number with the SPI connection

HW.Grad(iDevice).B = [4, 5];

%% properties for plot
HW.Grad(iDevice).Name{5} = 'B0 Prepol';
HW.Grad(iDevice).AmpUnit{5} = 'mT';
HW.Grad(iDevice).AmpUnitScale(5) = 1e-3;


%% DAC properties
HW.Grad(iDevice).DacBits(5) = 18;  % number of DAC bits

% dummy values (not needed for current controlled amplifiers?)
HW.Grad(iDevice).Dac2ExtGradUout(5) = 10/2^HW.Grad(1).DacBits(5);  % DAC output efficiency without load
HW.Grad(iDevice).Dac2ExtGradUout500Ohm(5) = 10/2^HW.Grad(1).DacBits(5);  % DAC output efficiency with 500 Ohm load

HW.Grad(iDevice).ExtGradOffsetDAC(5) = 0;  % offset in DAC for 0V output  % needs to be measured!

HW.DigitalIO(iDevice).InvertChannelOut = sum(2.^([1,3,4,5]-1));


%% amplifier properties
warning('off', 'PD:Grad:PaCurrentControlledChanged');
HW.Grad(iDevice).PaCurrentControlled(5) = 1;  % current controlled amplifier
HW.Grad(iDevice).AmpCurrentDependent(5) = 1;  % current dependent amplifier output

HW.Grad(iDevice).LoadRin(5) = 0.053;  % resistance of coil in Ohm, email by G. Ferrante 08-Jul-21

HW.Grad(iDevice).PaRout(5) = 15000;  % output impedance in Ohm
HW.Grad(iDevice).PaPmaxInt(5) = 100;  % maximum internal power dissipation (ever used???)

HW.Grad(iDevice).PaUin2PaIout(5) = 200/2.5;  % amplification in A/V  % needs to be measured!
HW.Grad(iDevice).PaUoutMax(5) = 20;  % maximum voltage of amplifier driver in V  % from Gianni Ferrante during Skype call on 22/07/2021

HW.Grad(iDevice).PaOffsetU(5) = 0; % amplifier offset in V  (set to 0 for current controlled amplifier)
HW.Grad(iDevice).PaOffsetI(5) = 0; % amplifier offset in A  % needs to be measured!


%% coil properties
HW.Grad(iDevice).LoadIin2Amp(5) = 20e-3/200;  % coil efficiency in T/A  % needs to be measured!
HW.Grad(iDevice).Inductance(5) = 15e-3;  % inductance of coil in H, email by G. Ferrante 08-Jul-21

% properties for coil overheat model
iCoilGroup = 3;
HW.Grad(iDevice).CoilThermalGroup(5) = iCoilGroup;
HW.Grad(1).CoilPowerDissipation(iCoilGroup) = 2000*sqrt(0.5);  % power dissipation at maximum temperature in Watt (50% duty cycle at 2 kW)
HW.Grad(1).CoilTemperatur(iCoilGroup) = 20;  % resting temperature in degrees C
% thermal capacity in J/K;  FIXME: Adjust values?
heatCapacityAlu = 0.9e3; % J/kg/K
massAlu = 1.2e3; % kg
HW.Grad(1).CoilThermalCapacity(iCoilGroup) = heatCapacityAlu*massAlu;
% maximum temperature in degrees C
% 2 kW for 10s with a 50% duty cycle should be the limit. Estimate the
% corresponding steady-state coil temperature and set it as the limit:
HW.Grad(1).CoilMaxTemperature(iCoilGroup) = HW.Grad(1).CoilTemperatur(iCoilGroup) + ...
  2000 * 10 / HW.Grad(1).CoilThermalCapacity(iCoilGroup) / sqrt(0.5);
% properties for "fuse blow" model
% maximum (nominal) DC current of fuse in Ampere, 50% duty cycle at 200 A
HW.Grad(iDevice).CoilMaxDcCurrent(5) = 200*sqrt(0.5);
% maximum "accumulated heat" in A^2*sec, max. 10 sec at 200 A
HW.Grad(iDevice).CoilCurrentSquareTime(5) = (200^2 - HW.Grad(iDevice).CoilMaxDcCurrent(5)^2) * 10 * 1.1;  % add 10% margin


%% pulse program defaults
HW.DefSeqValues.Prepol.useGradSPI = true;  % use SPI signal to switch prepolarization coil
HW.DefSeqValues.Prepol.useGrad = false;  % use signal on gradient channel to switch prepolarization coil
HW.DefSeqValues.Prepol.useDigitalIO = false;  % use digital output signal to switch prepolarization coil

HW.DefSeqValues.Prepol.tPrepol = 3;  % pre-polarization time in seconds
HW.DefSeqValues.Prepol.SPI.Channel = 5;  % gradient channel for SPI (should be 5)
HW.DefSeqValues.Prepol.SPI.Device = iDevice;  % device number (1 or 2)
HW.DefSeqValues.Prepol.SPI.tRampUp = 100e-3;  % ramp up time in seconds  % FIXME: Adjust value
HW.DefSeqValues.Prepol.SPI.nRampUp = 28;  % number of segments in ramp up to reduce load on amplifier
HW.DefSeqValues.Prepol.SPI.ampPrepol = 10e-3;  % amplitude for prepolarization in T  % FIXME: Adjust value
HW.DefSeqValues.Prepol.SPI.tRampDownNonAdiabatic = 100e-3;  % non-adiabatic ramp down time in seconds  % FIXME: Adjust value
HW.DefSeqValues.Prepol.SPI.nRampDownNonAdiabatic = 28;  % number of segments in non-adiabatic ramp down to reduce load on amplifier
HW.DefSeqValues.Prepol.SPI.ampNonAdiabatic = 200e-6;  % amplitude at end of non-adiabatic ramp in T  % FIXME: Adjust value
HW.DefSeqValues.Prepol.SPI.maxNegU = -10;  % "maximum" negative voltage of amplifier for ramp down in V  % FIXME: Adjust value
HW.DefSeqValues.Prepol.SPI.tRampDownAdiabatic = 20e-3;  % adiabatic ramp down time in seconds  % FIXME: Adjust value
% FIXME: The postset times of the two IO signals must(!) be different.
HW.DefSeqValues.Prepol.SPI.TransformerOnOutput = 2^0;  % digital output to turn transformer on
HW.DefSeqValues.Prepol.SPI.TransformerOnTOffset = 0.15;  % offset time for transformer on signal in seconds  % FIXME: Adjust value
HW.DefSeqValues.Prepol.SPI.TransformerOnTPostset = 0.15;  % postset time (time between pulse end and transformer off signal) in seconds  % FIXME: Adjust value
HW.DefSeqValues.Prepol.SPI.ClampOnOutput = 2^1;  % digital output to clamp the prepolarization coil
HW.DefSeqValues.Prepol.SPI.ClampOffTOffset = 0.12;  % offset time for clamp off signal in seconds  % FIXME: Adjust value
HW.DefSeqValues.Prepol.SPI.ClampOffTPostset = 10e-3;  % postset time (time between pulse end and clamp on signal) in seconds  % FIXME: Adjust value
HW.DefSeqValues.Prepol.tEnd = -(HW.DefSeqValues.Prepol.SPI.tRampDownNonAdiabatic/2 + ...
  HW.DefSeqValues.Prepol.SPI.tRampDownAdiabatic + 50e-3);


HW.Function_Before_Measurement = @check_SPI_Prepol_Amplifier;
HW.Function_While_Measurement = @check_SPI_Prepol_Amplifier;
HW.DefSeqValues.Prepol.SPI.ErrorInputDevice = iDevice;  % device for amplifier error signal
HW.DefSeqValues.Prepol.SPI.ErrorInputChannel = 6;  % digital input channel for amplifier error signal
