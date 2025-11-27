%% System settings for the PrimoGaia measurement project
%
% Set a sensible Larmor frequency before calling this script (from the
% LoadMySystem file of that system).
%
% ------------------------------------------------------------------------------
% (C) Copyright 2021-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------


%% Sequence settings

HW.FindFrequencySweep.maxTime = Inf;  % Don't run automatic frequency sweeps

HW.FindShim.AQEcho = 0.7;
HW.FindShim.fSample = 150e3;

HW.Grad(1).PaEnable = 1;

HW.FindPulseDuration.AQSlice.nRead = 16;
HW.FindPulseDuration.AQSlice.sizeRead = 0.05;
HW.FindPulseDuration.AQSlice.alfa = 0;
HW.FindPulseDuration.AQSlice.phi = 0.5 * pi;
HW.FindPulseDuration.AQSlice.theta = 0;


%% TX settings
HW.TX(1).ChannelDef = 2;
HW.TX(1).Def.Amplitude(HW.TX(1).ChannelDef) = 0.99 * HW.TX(1).Max.Amplitude(HW.TX(1).ChannelDef);  % send with maximum voltage by default
% HW.TX(1).Def.Uout(HW.TX(1).ChannelDef) = HW.TX(1).Max.Uout(HW.TX(1).ChannelDef);  % send with maximum voltage by default
HW.TX(1).Def.Uout(HW.TX(1).ChannelDef) = 3.8;  % send with 2 volts by default
HW.TX(1).Max.Uout(HW.TX(1).ChannelDef) = 3.8;  % send with 2 volts by default

% mit Thomann Behringer A800 Audio Verstärker
HW.TX(1).Uout2PaUout(2) = 20; % x2 wegen fehlendem 50 Ohm am Eingang
HW.TX(1).Def.Uout(HW.TX(1).ChannelDef) = 0.25;  %
HW.TX(1).Def.Uout(HW.TX(1).ChannelDef) = 0.5;  %
% HW.TX(1).Def.Uout(HW.TX(1).ChannelDef) = 0.5;  %
% HW.TX(1).Def.Uout(HW.TX(1).ChannelDef) = 0.25;  %
% HW.TX(1).Max.Uout(HW.TX(1).ChannelDef) = 4 / 2;  %
HW.TX(1).Max.Uout(HW.TX(1).ChannelDef) = 1;  %
HW.TX(1).Def.PaUout(HW.TX(1).ChannelDef) = 20;         % whole body TX coil
HW.TX(1).Def.PaUout(HW.TX(1).ChannelDef) = 20*0.05;    % head Tx coil
HW.TX(1).BlankPostset = 500e-6;
HW.TX(1).BlankOffset = 500e-6;
% warning('ohne Thomann');

HW.TX2RXdeadTime = 55e-3; % Alu-Rohr 3mm
HW.TX2RXdeadTime = 80e-3; % Kurzgeschlossene Prepolspule
HW.TX2RXdeadTime = 70e-3; % bei fLarmor=2 kHz und kurzgeschlossene Prepolspule
HW.TX2RXdeadTime = 20e-3; % bei fLarmor=8 kHz und kurzgeschlossene Prepolspule
if numel(HW.TX) > 1
  HW.TX(2).PaUout2Amplitude(2) = 1;  % DNP amplitudes in Volt; FIXME: Measure coil efficiency!
  HW.TX(2).Max.Amplitude(2) = HW.TX(2).Max.PaUout(2) .* HW.TX(2).PaUout2Amplitude(2);
  HW.TX(2).AmplitudeName = 'TX B2';
  HW.TX(2).AmplitudeUnit = 'Veff';
  HW.TX(2).AmplitudeUnitScale = 1/sqrt(2);
end

%% RX settings
HW.RX(1).VGAGainDef = HW.RX(1).VGAGainDef / HW.RX(1).GainDef * HW.RX(1).Amplitude2Uin / 100e-3;  % 20 mV maximum input voltage
% HW.RX(1).VGAGainDef = HW.RX(1).VGAGainDef / HW.RX(1).GainDef * HW.RX(1).Amplitude2Uin / 200e-3;  % 20 mV maximum input voltage
% HW.RX(1).LNAGain = 500*sqrt(2)/3;
HW.RX(1).LNAGain = 430*sqrt(2)/3; %
% HW.RX(1).Rin = 200e3;
HW.RX(1).Rin = 36e3;
HW.RX(1).Rin = 120e3;

%% DNP amplifier
if numel(HW.TX) > 1
  iDevice = 2; %#ok<NASGU>
  % LoadTxPa_16022017; % PD amplifier
  LoadRF1180_100_01;  % RFPA amplifier 100W (1MHz-180MHz)
end

%% pre-polarization settings
HW.DefSeqValues.Prepol.use = false;
HW.DefSeqValues.Prepol.tPrepol = 6.0;  % pre-polarization time in seconds
HW.DefSeqValues.Prepol.tOffNonAdiabatic = 1e-6;2.6e-3;  % duration of non-adiabatic switch-off ramp in seconds
HW.DefSeqValues.Prepol.tOffAdiabatic = 1e-6;10e-3;60.0e-3;  % duration of adiabatic switch-off ramp in seconds
HW.DefSeqValues.Prepol.tEnd = -80e-3;  % end of prepol before center of excitation pulse (at Seq.tEcho/2) in seconds
HW.DefSeqValues.Prepol.useGradSPI = false;  % use SPI signal to switch prepolarization coil
HW.DefSeqValues.Prepol.useDigitalIO = true;  % use digital output signal to switch prepolarization coil
HW.DefSeqValues.Prepol.digiOutPrepol = 2^(1-1);  % digital output for prepolarization pulse
HW.DefSeqValues.Prepol.digiOutNonAdiabatic = 0;  % digital output for non-adiabatic switch-off ramp
HW.DefSeqValues.Prepol.digiOutAdiabatic = 0;2^(2-1);  % digital output for adiabatic switch-off ramp
HW.DefSeqValues.Prepol.ResetDDS = true;  % reset the DDS at each tRep with prepol pulse
HW.DefSeqValues.Prepol.ClampOnOutput = 0;

HW.TX2RXdeadTime = 13e-3+15e-3; % +15 ms mit großer Aluröhre und Kopfspule
HW.RX(1).ClampCoil.tOffset = 35e-3;
% HW.RX(1).ClampCoil.tOffset = 50e-3;
HW.RX(1).ClampCoil.nSampleOffset = 4;
HW.RX(1).ClampCoil.Enable = false;


HW.RX(1).ClampCoil.DigitalOutputSignal = 2^(2-1);
% HW.RX(1).ClampCoil.tOffset = 30e-3;HW.TX2RXdeadTime/2;  % start clamping signal before start of acquisition in seconds
HW.RX(1).ClampCoil.tPostset = 20e-3;  % stop clamping signal after end of acquisition in seconds
HW.RX(1).ClampCoil.nSamplePostset = 4;

% HW.RX(1).ClampCoil.tOffset = 70e-3;HW.TX2RXdeadTime/2;  % start clamping signal before start of acquisition in seconds
% HW.RX(1).ClampCoil.tPostset = 20e-3;  % stop clamping signal after end of acquisition in seconds

HW.DigitalIO(1).InvertChannelOut = HW.RX(1).ClampCoil.DigitalOutputSignal;

if HW.RX(1).ClampCoil.Enable
  HW.TX2RXdeadTime = max(HW.TX2RXdeadTime, HW.RX(1).ClampCoil.tOffset+1e-3);  % account for ring-down time of TX coil or delay of relay?
  HW.RX2TXdeadTime = max(HW.RX2TXdeadTime, HW.RX(1).ClampCoil.tPostset);
end


%% DNP settings
if numel(HW.TX) > 1
  HW.DefSeqValues.DNP.Device = 2;  % DNP pulse TX Device, e.g. 1
  HW.DefSeqValues.DNP.digitalIO.Device = 2;  % Device that sets the digital output signal
  dBmTX = 40;  % DNP pulse amlitude in dBm
else
  % for local tests only
  warning('second console is not connected or switched off.');
  HW.DefSeqValues.DNP.Device = 1;  % DNP pulse TX Device, e.g. 1
  HW.DefSeqValues.DNP.digitalIO.Device = 1;  % Device that sets the digital output signal
  dBmTX = -40;
end
HW.DefSeqValues.DNP.use = false;
HW.DefSeqValues.DNP.Channel = 2;  % DNP pulse TX Channel, e.g. 2
HW.DefSeqValues.DNP.Frequency = 70e6;  % DNP pulse frequency in Hz
HW.DefSeqValues.DNP.Amplitude = (0.001*10^(dBmTX/10) * ...
  HW.TX(HW.DefSeqValues.DNP.Device).LoadRin(HW.DefSeqValues.DNP.Channel))^0.5 * sqrt(2) * ...
  HW.TX(HW.DefSeqValues.DNP.Device).PaUout2Amplitude(HW.DefSeqValues.DNP.Device)/1;  % amplitude in Volt
HW.DefSeqValues.DNP.tDNP = 2.0;  % duration of DNP pulse in seconds
HW.DefSeqValues.DNP.tEnd = -20e-3;  % end of prepol before center of excitation pulse (at Seq.tEcho/2) in seconds
HW.DefSeqValues.DNP.useExtDDS = true;  % use external DDS
HW.DefSeqValues.DNP.ResetDDS = true;  % reset the DDS at each tRep with prepol pulse
HW.DefSeqValues.DNP.digitalIO.use = false;  % Switch DDS with digital output signal
HW.DefSeqValues.DNP.digitalIO.Channel = 1;  % Digital output channel for the DNP pulse signal
HW.DefSeqValues.DNP.continuous = false;  % keep DNP pulse continuously running during experiment

% For continuous DNP mode, interrupt the DNP pulse for the acquisition windows
HW.RX(1).ClampCoil.Enable = false;  % use clamping signal around acquisitions
HW.RX(1).ClampCoil.DigitalOutputSignal = 2^(HW.DefSeqValues.DNP.digitalIO.Channel-1);  % digital output signal added for clamping
HW.RX(1).ClampCoil.DigitalOutputDevice = HW.DefSeqValues.DNP.digitalIO.Device;  % emit signal on specified MMRT device
HW.RX(1).ClampCoil.SignalOff = true;  % set the signal to off for duration around acquisitions.
HW.RX(1).ClampCoil.nSampleOffset = 3;  % additional samples with signal before the acquisition
HW.RX(1).ClampCoil.nSamplePostset = 3;  % additional samples with signal after the acquisition
HW.RX(1).ClampCoil.tOffset = 9e-3;  % interrupt signal before start of acquisition in seconds
HW.RX(1).ClampCoil.tPostset = 9e-3;  % reset signal after end of acquisition in seconds

%% Gradient system

iDevice = 1; %#ok<NASGU>
% LoadSystem_PrimoGAIA_2ndKernel_GradSystem_small;
LoadSystem_PrimoGAIA_2ndKernel_GradSystem_large;

%% B0 field cage
iDevice = 2;
LoadSystem_PrimoGAIA_2ndKernel_B0Cage;

if numel(HW.Grad) > 1
  %% Second drive-l for B0 field cage

  iDevice = 2;

  %% DNP coil settings
  % FIXME: Adjust those values
  HW.TX(2).CoilPowerDissipation = 1e9;  % power dissipation in Watt
  HW.TX(2).CoilTemperature = 20;  % temperature of coil at start of measurement in degrees C
  HW.TX(2).CoilMaxTemperature = 60;  % maximum temperature (in model) before damage occurs in degrees C
  % thermal capacity in J/K;  FIXME: Adjust values?
  specificHeatCapacity = 385; % J/kg/K (copper wire)
  % volumeCoil = ???;  % volume of copper in coil in m^3
  % mass = volumeCoil*9e3;  % in kg (with density of copper wires)
  mass = 2000e-3;  % mass of copper wires in kg
  HW.TX(2).CoilThermalCapacity = specificHeatCapacity*mass;

end


%% prepolarization coil amplifier settings
if 0
  iDevice = 1;  % device number with the SPI connection
  LoadSystem_SPI_Prepol;


  % properties for coil overheat
  iCoilGroup = 3;

  % properties for "fuse blow" model
  HW.Grad(iDevice).CoilMaxDcCurrent(5) = 250 * 0.9;  % maximum (nominal) DC current of fuse in Ampere
  HW.Grad(iDevice).CoilCurrentSquareTime(5) = 0.9 * 0.9;  % maximum "accumulated heat" in A^2*sec

  HW.DefSeqValues.Prepol.useGradSPI = true;  % use SPI signal to switch prepolarization coil
  HW.DefSeqValues.Prepol.useGrad = false;  % use signal on gradient channel to switch prepolarization coil
  HW.DefSeqValues.Prepol.useDigitalIO = false;  % use digital output signal to switch prepolarization coil
  HW.DefSeqValues.Prepol.SPI.Channel = 5;  % gradient channel for SPI (should be 5)
  HW.DefSeqValues.Prepol.SPI.Device = iDevice;  % device number (1 or 2)
  HW.DefSeqValues.Prepol.SPI.tRampUp = 300e-3;  % ramp up time in seconds
  HW.DefSeqValues.Prepol.SPI.ampPrepol = 20e-3;  % amplitude for prepolarization in T
  HW.DefSeqValues.Prepol.SPI.tRampDownNonAdiabatic = 50e-3;  % non-adiabatic ramp down time in seconds
  HW.DefSeqValues.Prepol.SPI.ampNonAdiabatic = 200e-6;  % amplitude at end of non-adiabatic ramp in T
  HW.DefSeqValues.Prepol.SPI.tRampDownAdiabatic = 5e-3;  % adiabatic ramp down time in seconds
  % FIXME: The postset times of the two IO signals must(!) be different.
  HW.DefSeqValues.Prepol.SPI.TransformerOnOutput = 2^0;  % digital output to turn transformer on
  HW.DefSeqValues.Prepol.SPI.TransformerOnTOffset = 0.5;  % offset time for transformer on signal in seconds
  HW.DefSeqValues.Prepol.SPI.TransformerOnTPostset = 0.5;  % postset time (time between pulse end and transformer off signal) in seconds
  HW.DefSeqValues.Prepol.SPI.ClampOnOutput = 0;  % digital output to clamp the prepolarization coil
  HW.DefSeqValues.Prepol.SPI.ClampOffTOffset = 0.5;  % offset time for clamp off signal in seconds
  HW.DefSeqValues.Prepol.SPI.ClampOffTPostset = 0.320;  % postset time (time between pulse end and clamp on signal) in seconds

  % HW.Function_Before_Measurement = @check_Prepol_Amplifier;
  % HW.DefSeqValues.Prepol.SPI.ErrorInputDevice = iDevice;  % device for amplifier error signal
  % HW.DefSeqValues.Prepol.SPI.ErrorInputChannel = 6;  % digital input channel for amplifier error signal
end


%% external DDS settings
DDSConfigExePath = fullfile(winqueryreg('HKEY_LOCAL_MACHINE', 'SOFTWARE\Pure Devices\Settings', 'LibPath'), ...
  '..', 'PD-DDS', 'DDSConfig.exe');
DDS = PD.DDS.GetInstance(DDSConfigExePath);
% DDS.RDACRef = 10e3;  % reference resistance for DAC
DDS.IDAC2Uout = [DDS.IDACMin; DDS.IDACMax] \ [135e-3; 484e-3];
