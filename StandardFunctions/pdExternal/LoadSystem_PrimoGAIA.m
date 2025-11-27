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
HW.FindPulseDuration.AQSlice.sizeRead = 0.1;
HW.FindPulseDuration.AQSlice.alfa = 0;
HW.FindPulseDuration.AQSlice.phi = 0.5 * pi;
HW.FindPulseDuration.AQSlice.theta = 0;


%% TX settings
HW.TX(1).ChannelDef = 2;
HW.TX(1).Def.Amplitude(HW.TX(1).ChannelDef) = 0.99 * HW.TX(1).Max.Amplitude(HW.TX(1).ChannelDef);  % send with maximum voltage by default
% HW.TX(1).Def.Uout(HW.TX(1).ChannelDef) = HW.TX(1).Max.Uout(HW.TX(1).ChannelDef);  % send with maximum voltage by default
HW.TX(1).Def.Uout(HW.TX(1).ChannelDef) = 3.99;  % send with 2 volts by default

% HW.TX(1).BlankOffset = 200e-6;  % offset of un-blanking signal before rf pulse
% HW.TX(1).BlankPostset = 100e-6;  % additional un-blanking time after rf pulse
% HW.TX(1).BlankOffsetAQ = 5e-6;   % blank of receiver before TX pulse in s
% HW.TX(1).BlankPostsetAQ = 10e-6;  % blank of receiver after TX pulse in s
% HW.TX(1).BlankAQ = 1;               % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
% HW.TX2RXdeadTime = 25e-3;  % antenne RMN EFNMR
% HW.TX2RXdeadTime = 13e-3;
% HW.TX2RXdeadTime = 15e-3;  % antenne RMN gradio rat
HW.TX2RXdeadTime = 18e-3;  % antenne RMN gradio rat
if numel(HW.TX) > 1
  HW.TX(2).PaUout2Amplitude(2) = 1;  % DNP amplitudes in Volt; FIXME: Measure coil efficiency!
  HW.TX(2).Max.Amplitude(2) = HW.TX(2).Max.PaUout(2) .* HW.TX(2).PaUout2Amplitude(2);

  HW.TX(2).AmplitudeName = 'TX B2';
  HW.TX(2).AmplitudeUnit = 'Veff';
  HW.TX(2).AmplitudeUnitScale = 1/sqrt(2);
end


%% RX settings
% HW.RX(1).VGAGainDef = .1*HW.RX(1).VGAGainDef / HW.RX(1).GainDef * HW.RX(1).Amplitude2Uin / 20e-3;  % 20 mV maximum input voltage
% HW.RX(1).VGAGainDef = HW.RX(1).VGAGainDef / HW.RX(1).GainDef * HW.RX(1).Amplitude2Uin / 20e-3;  % 20 mV maximum input voltage
HW.RX(1).VGAGainDef = HW.RX(1).VGAGainDef / HW.RX(1).GainDef * HW.RX(1).Amplitude2Uin / 20e-3;  % 20 mV maximum input voltage
HW.RX(1).LNAGain = 500*sqrt(2)/3;
HW.RX(1).Rin = 200e3;

%% DNP amplifier
if numel(HW.TX) > 1
  iDevice = 2;
  % LoadTxPa_16022017; % PD amplifier
  LoadRF1180_100_01;  % RFPA amplifier 100W (1MHz-180MHz)
end

%% pre-polarization settings
HW.DefSeqValues.Prepol.use = false;
HW.DefSeqValues.Prepol.tPrepol = 3.0;  % pre-polarization time in seconds
HW.DefSeqValues.Prepol.tOffNonAdiabatic = 6.0e-3;  % duration of non-adiabatic switch-off ramp in seconds
HW.DefSeqValues.Prepol.tOffAdiabatic = 60.0e-3;  % duration of adiabatic switch-off ramp in seconds
HW.DefSeqValues.Prepol.tEnd = -40e-3;  % end of prepol before center of excitation pulse (at Seq.tEcho/2) in seconds
HW.DefSeqValues.Prepol.useGradSPI = false;  % use SPI signal to switch prepolarization coil
HW.DefSeqValues.Prepol.useDigitalIO = true;  % use digital output signal to switch prepolarization coil
HW.DefSeqValues.Prepol.digiOutPrepol = 2^(2-1);  % digital output for prepolarization pulse
HW.DefSeqValues.Prepol.digiOutNonAdiabatic = 0;  % digital output for non-adiabatic switch-off ramp
HW.DefSeqValues.Prepol.digiOutAdiabatic = 2^(1-1);  % digital output for adiabatic switch-off ramp
HW.DefSeqValues.Prepol.ResetDDS = true;  % reset the DDS at each tRep with prepol pulse


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
HW.DefSeqValues.DNP.Frequency = 69.5e6;  % DNP pulse frequency in Hz
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

% 4-channel gradient system (S. Lotter, p. 56)
% S1 --> channel 2 (B0 direction)
% S2 --> channel 1 (tube direction)
% S3 --> channel 4 (unused)
% S4 --> channel 3 (perpendicular to tube and B0 direction)

% S2 --> channel 1 (tube direction)
% S1 --> channel 2 (B0 direction)
% S4 --> channel 3 (perpendicular to tube and B0 direction)
% S3 --> channel 4 (unused)


% Channel names and units for sequence plots
HW.Grad(1).Name = {'Grad x', 'Grad y', 'Grad z', 'Grad 4'};
HW.Grad(1).AmpUnit = {'mT/m', 'mT/m', 'mT/m', 'mT/m'};
HW.Grad(1).AmpUnitScale = [1e-3, 1e-3, 1e-3, 1e-3];

maxGrad_Dev1 = 4;
if 0
  %% 3-channel gradient system (not used!)
  % simple correlation of channels to encoding directions
  HW.Grad(1).x = 1;
  HW.Grad(1).y = 3;
  HW.Grad(1).z = 2;
  HW.Grad(1).B = 4;  % 4th gradient channel

  % coil efficiency
  % Measured in (T/m)/A (Tesla per Meter per Ampere);
  % or more generally: [amplitude unit]/[amplifier input unit]
  % HW.Grad(1).LoadIin2Amp(1:3) = [1.1e-3*((37.5+15)/60), 0.89e-3*((20.3+28.1)/60), 1.0e-3*((21.25+36.25)/90)];  % x y z

  HW.Grad(1).LoadIin2Amp(1:3) = [0.9625e-3*((26.6+26.6)/60)*((28.1+32.8)/60)*((33.6+25.8)/60), ...
     0.7179e-3*((38+26)/60)*((49.2+23.4)/60)*((35.6+20.6)/60)*((18.75+28.9)/60), ...
     0.6389e-3*((45+90)/93)*((42.2+60.9)/93)];  % x y z

  HW.Grad(1).AmpCurrentDependent(4) = false;
  HW.Grad(1).LoadUin2Amp(4) = 1;  % gradient amplitude in Volts

  % resistance of coils in Ohm
  % 3-channel gradient system
  % HW.Grad(1).LoadRin = [14.7+2*68, 9.7+2*68, 12.4+2*68, 100];
  HW.Grad(1).LoadRin = [14.7+2*47, 9.7+2*47, 12.4+2*47, 100];

  if 0
    %% settings for gradient channel that switches pre-polarizer
    HW.Grad(1).AmpCurrentDependent(4) = false;
    HW.Grad(1).LoadUin2Amp(4) = 1;  % "gradient" amplitude in Volts
    HW.Grad(1).Name{4} = 'U Prepol';
    HW.Grad(1).AmpUnit{4} = 'V';
    HW.Grad(1).AmpUnitScale(4) = 1;

    maxGrad_Dev1 = 3;

    % properties for coil overheat
    % disable checks for that "coil"
    iCoilGroup = 2;
    HW.Grad(1).CoilThermalGroup(4) = iCoilGroup;
    HW.Grad(1).CoilPowerDissipation(iCoilGroup) = Inf;  % power dissipation at maximum temperature in Watt
    HW.Grad(1).CoilTemperatur(iCoilGroup) = 20;  % resting temperature in degrees C
    HW.Grad(1).CoilMaxTemperature(iCoilGroup) = Inf;  % maximum temperature in degrees C
    HW.Grad(1).CoilThermalCapacity(iCoilGroup) = Inf;  % thermal capacity in J/K
    % properties for "fuse blow" model
    HW.Grad(1).CoilMaxDcCurrent(4) = Inf;  % maximum (nominal) DC current of fuse in Ampere
    HW.Grad(1).CoilCurrentSquareTime(4) = Inf;  % maximum "accumulated heat" in A^2*sec
  end

end


if 1
  %% 4-channel gradient system

  % correlation channels to encoding directions (to match coordinate system of B0 cage)
  HW.Grad(1).x = 1;
  HW.Grad(1).y = 3;
  HW.Grad(1).z = 2;
  HW.Grad(1).B = 4;  % 4th gradient channel

  % coil efficiency
  % Measured in (T/m)/A (Tesla per Meter per Ampere);
  % or more generally: [amplitude unit]/[amplifier input unit]
  % signs without any filters
  % 4th channel not calibrated!
  HW.Grad(1).LoadIin2Amp(1:3) = [0.000896300850236, ... % should be fine
                                 0.000823693937719, ... % ok
                                 0.000487150116648]; % ok % x y z

  % signs for right-handed coordinate system
  HW.Grad(1).xyzBDir = [-1, -1, 1, 1];  % x y z prepol - sign (for inversing direction)

  % resistance of gradient coils (including cabling)
  % HW.Grad(1).LoadRin = [18.4, 13.5, 16.0, 18.4];  % no filters
  % HW.Grad(1).LoadRin = [18.2, 13.1, 16.4, 18.0];
  % HW.Grad(1).LoadRin = [18.2, 13.1, 16.14, 18.0];  % Low Z with 0.23ohm cable   ====> [18.0, 12.9, 16.0, 18.0]  FH cable
  % HW.Grad(1).LoadRin = [15.05, 10.2, 15.9, 14.8];  % Low Z with 0.23ohm cable
  HW.Grad(1).LoadRin = [14.50, 9.70, 14.76, 12.57];  % Low Z with new double cable 2022-03-17

end


% gradients used for auto shimming routine
HW.Grad(1).ShimGradients = [1, 1, 1, 0];

% extent (and location) of the image volume in meters
HW.Grad(1).ImageVol = [-0.06, 0.06, -0.04, 0.04, -0.04, 0.04]*1.5;  % FIXME: Anpassen
HW.Grad(1).ImageVolOffset = [0, 0, 0];  % offset of tube

% FIXME: Alles anpassen
HW.Grad(1).tRamp = 20e-3;  % minimum ramp time (for imaging gradients) in s
HW.Grad(1).SystemTimeDelay(1:3) = [0.00870491, 0.0118726, 0.00852225];  % Time delay of grad amp
HW.Grad(1).tEC = max(HW.Grad(1).SystemTimeDelay(:))*1;  % eddy current time in s (i.e. time until gradient is stable)

% maximum total power (coil heating, model parameters)
% FIXME: Adapt to properties of used coils.
HW.Grad(1).CoilPowerDissipation(1) = 5;  % power dissipation at maximum temperature for each coil group in Watt
HW.Grad(1).CoilTemperatur(1) = 20;  % resting temperature for each coil group in degrees C
HW.Grad(1).CoilMaxTemperature(1) = 60;  % maximum temperature for each coil group in degrees C
HW.Grad(1).CoilThermalCapacity(1) = 0.06 * 0.5 * 0.5 * 3 * 0.8 * 1000;  % thermal capacity for each coil group in J/K

% properties for "fuse blow" model
HW.Grad(1).CoilMaxDcCurrent(1:maxGrad_Dev1) = 0.75 * 0.9;  % maximum (nominal) DC current of fuse in Ampere
HW.Grad(1).CoilCurrentSquareTime(1:maxGrad_Dev1) = 0.9 * 0.9;  % maximum "accumulated heat" in A^2*sec


if numel(HW.Grad) > 1
  %% Second drive-l for B0 field cage

  iDevice = 2;
  LoadGradAmp_DC600_SN_26;

  % X: from left to right (along gradient tube)
  % Y: from front to back
  % Z: from bottom to top

  % For simplicity, use a simple correlation:
  HW.Grad(2).x = 1;  % corresponds to Grad(5) in pulse programs
  HW.Grad(2).y = 2;  % corresponds to Grad(6) in pulse programs
  HW.Grad(2).z = 3;  % corresponds to Grad(7) in pulse programs
  HW.Grad(2).B = 4;  % corresponds to Grad(8) in pulse programs
  % HW.Grad(2).xyzBDir = [-1, -1, -1, 1];  % x y z prepol - sign (for inversing direction)
  % exchange y-z coordinates (9th March 2021)
  HW.Grad(2).xyzBDir = [-1, 1, -1, 1];  % x y z prepol - sign (for inversing direction)

  % Channel names and units for sequence plots
  HW.Grad(2).Name = {'B0 x', 'B0 y', 'B0 z', 'unused'};
  HW.Grad(2).AmpUnit(1:3) = {[char(181) 'T']};
  HW.Grad(2).AmpUnit{4} = 'A';
  HW.Grad(2).AmpUnitScale = [1e-6, 1e-6, 1e-6, 1];

  HW.Grad(2).HoldShim = 1;

  %% coil efficiency
  % Measured in T/A (Tesla per Ampere);
  % or more generally: [amp unit]/[controller input unit]

  % Alternatively, amplitude in Ampere:
  % HW.Grad(2).LoadIin2Amp = [1, 1, 1, 1];  % x y z XXX?

  % HW.Grad(2).LoadIin2Amp = [152e-6, 590e-6*104/100, 1100e-6*98.5/100, 1];  % B0 cage (x y z XXX?) calibrated in T/A with MR on 28/01/2021
  % re-calibrated after modifications (x-coil 32 parallel windings per side)
  % HW.Grad(2).LoadIin2Amp = [7.312639051188890e-05*1.0108, 4.094945768319414e-04*1.004, 3.611585947860855e-04*1.0148, 1];  % x y z XXX?

  % exchange y & z coordinates (9th March 2021)
  % HW.Grad(2).LoadIin2Amp = [7.312639051188890e-05*1.0108, 3.611585947860855e-04*1.0148, 4.094945768319414e-04*1.004, 1];  % x y z XXX? x? y400 z400

  % reduce number of windings in y-coil and z-coil from 400 to 200 respectively
  % HW.Grad(2).LoadIin2Amp = [7.312639051188890e-05*1.0108, 3.611585947860855e-04*1.0148/2, 4.094945768319414e-04*1.004/2, 1];  % x y z XXX? 10.03.2021 x? y200z200
  % re-calibration, 10/03/2021
  % HW.Grad(2).LoadIin2Amp = [7.387719039280947e-05, 1.836016256569566e-04*1.004997501249375, (1.957216097287849e-04*1.050659630606860)/(8730/8718)/(8730/8768)/(8790/8781), 1];  % x y z XXX? 10.03.2021
  % HW.Grad(2).LoadIin2Amp = [7.387719039280947e-05, 1.836016256569566e-04*1.004997501249375, (1.957216097287849e-04*1.050659630606860)/(8790/8810), 1];
  %HW.Grad(2).LoadIin2Amp = [7.387719039280947e-05, 1.836016256569566e-04*1.004997501249375, (1.957216097287849e-04*1.050659630606860)/(8779/8758), 1];
  % HW.Grad(2).LoadIin2Amp = [7.387719039280947e-05, 1.836016256569566e-04*1.004997501249375, (1.957216097287849e-04*1.050659630606860)/(8780/8772), 1];  % gradio rat
  % HW.Grad(2).LoadIin2Amp = [7.387719039280947e-05, 1.836016256569566e-04*1.004997501249375, (1.957216097287849e-04*1.050659630606860)/(8780/8770), 1];  % gradio rat
  HW.Grad(2).LoadIin2Amp = [7.387719039280947e-05, 1.836016256569566e-04*1.004997501249375, (1.957216097287849e-04*1.050659630606860), 1];  % gradio rat
  % HW.Grad(2).LoadIin2Amp = [7.387719039280947e-05, 1.836016256569566e-04*1.004997501249375, (1.957216097287849e-04*1.050659630606860), 1];  % antenne RMN EFNMR
  % HW.Grad(2).LoadIin2Amp = [7.387719039280947e-05, 1.836016256569566e-04*1.004997501249375, (1.957216097287849e-04*1.050659630606860)/(8800/8770), 1];   %salle 4.7T

  %% resistance of coils in Ohm
  % HW.Grad(2).LoadRin = [55.2, 91.2, 79.9, 1];  % gemessen am 28.01.2021

  % x-coil 32 parallel windings per side
  % z-coil windings reduced from 600 to 400
  % HW.Grad(2).LoadRin = [12.8, 91.2, 51.5, 1];  % gemessen am 10.02.2021

  % exchange y-z coordinates (9th March 2021)
  % HW.Grad(2).LoadRin = [13.11, 25.9, 57.8, 1];  % gemessen am 22.02.2021 x? y400 z400
  HW.Grad(2).LoadRin = [13.11, 25.9/2, 57.8/2, 1];  % gemessen am 10.03.2021 x? y200 z200

  % gradients used for auto shimming routine
  HW.Grad(2).ShimGradients(1:4) = 0;

  % extent (and location) of the image volume in meters
  HW.Grad(2).ImageVol = HW.Grad(1).ImageVol;  % FIXME: Would it make sense if we set something different?
  HW.Grad(2).ImageVolOffset = HW.Grad(1).ImageVolOffset;  % FIXME: Would it make sense if we set something different?

  % FIXME: Alles anpassen
  HW.Grad(2).tRamp = 2e-3;  % minimum ramp time (for B0 cage fields) in s
  HW.Grad(2).SystemTimeDelay = [163e-6, 215e-6, 217e-6, 215e-6]-0e-6;  % time delay of grad amp in s (FIXME: Can this be measured?)
  HW.Grad(2).tEC = max(HW.Grad(2).SystemTimeDelay(3))*0.5;  % eddy current time in s (i.e. time until gradient is stable)

  % maximum total power (coil heating, model parameters)
  % FIXME: Adapt to properties of used coils.
  % 4th channel is currently unused. Set some values anyway.
  HW.Grad(2).CoilThermalGroup = 1:4;
  % HW.Grad(2).CoilPowerDissipation(1:4) = 30;  % Watt
  HW.Grad(2).CoilPowerDissipation(1:4) = 45;  % Watt
  HW.Grad(2).CoilTemperatur(1:4) = 20;  % temperature of coil at start of measurement in degrees C
  HW.Grad(2).CoilMaxTemperature(1:4) = 60;  % maximum temperature (in model) before damage occurs in degrees C
  % thermal capacity in J/K;  FIXME: Adjust values?
  specificHeatCapacity = 500; % J/kg/K (estimated for mix of copper wires and aluminum frame)
  volumeCoil = 2*4*1*0.015^2;  % volume of coil in m^3 (2 coils with 4 sides of ~1 m and a cross-section of ~1.5 cm^2)
  mass = volumeCoil*5e3;  % in kg (estimated density for mix of copper wires and aluminum frame)
  HW.Grad(2).CoilThermalCapacity(1:4) = specificHeatCapacity*mass;

  % properties for "fuse blow" model
  HW.Grad(2).CoilMaxDcCurrent = [3.15, 2.0, 2.0, 1.0] * 0.9;  % maximum (nominal) DC current of fuse in Ampere
  HW.Grad(2).CoilCurrentSquareTime = [0.9, 0.9, 0.9, 0.9] * 0.9;  % maximum "accumulated heat" in A^2*sec

  fprintf('B0 field via cage in z-direction %3.1f %cT (in %s).\n', HW.B0*1e6, char(181), mfilename());
  HW.Grad(2).AmpOffsetExtra(3) = HW.B0;
  % HW.Grad(2).AmpOffsetExtra(2) = HW.B0;
  % HW.Grad(2).AmpOffsetExtra(1) = HW.B0;
  % HW.Grad(2).AmpOffsetExtra(1:3) = -cross(HW.Grad(2).AmpOffset(1:3), [1,0,0]) / norm(HW.Grad(2).AmpOffset(1:3),2) * HW.B0;


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
  HW.Grad(iDevice).B = [4, 5];
  HW.Grad(iDevice).Name{5} = 'B0 Prepol';
  HW.Grad(iDevice).AmpUnit{5} = 'mT';
  HW.Grad(iDevice).AmpUnitScale(5) = 1e-3;

  warning('off', 'PD:Grad:PaCurrentControlledChanged');
  HW.Grad(iDevice).PaCurrentControlled(5) = 1;  % current controlled amplifier
  HW.Grad(iDevice).AmpCurrentDependent(5) = 1;  % current dependent amplifier output

  HW.Grad(iDevice).LoadRin(5) = 0.3;  % resistance of coil in Ohm

  HW.Grad(iDevice).DacBits(5) = 18;  % number of DAC bits

  % dummy values (not needed for current controlled amplifiers?)
  HW.Grad(iDevice).Dac2ExtGradUout(5) = 10/2^HW.Grad(1).DacBits(5);  % DAC output efficiency without load
  HW.Grad(iDevice).Dac2ExtGradUout500Ohm(5) = 10/2^HW.Grad(1).DacBits(5);  % DAC output efficiency with 500 Ohm load

  HW.Grad(iDevice).PaRout(5) = 15000;  % output impedance in Ohm
  HW.Grad(iDevice).PaPmaxInt(5) = 100;  % maximum internal power dissipation (ever used???)

  HW.Grad(iDevice).ExtGradOffsetDAC(5) = 0;  % offset in DAC for 0V output
  HW.Grad(iDevice).PaOffsetU(5) = 0; % amplifier offset in V  (set to 0 for current controlled amplifier)
  HW.Grad(iDevice).PaOffsetI(5) = 0; % amplifier offset in A  % needs to be measured!

  HW.Grad(iDevice).PaUin2PaIout(5) = 200/2.5;  % amplification in A/V
  HW.Grad(iDevice).LoadIin2Amp(5) = 20e-3/200;  % coil efficiency in T/A

  HW.Grad(iDevice).PaUoutMax(5) = 100;  % maximum voltage of amplifier driver in V
  HW.Grad(iDevice).Inductance(5) = 1e-3;  % inductance of coil in H

  % properties for coil overheat
  iCoilGroup = 3;
  HW.Grad(iDevice).CoilThermalGroup(5) = iCoilGroup;
  HW.Grad(1).CoilPowerDissipation(iCoilGroup) = 50;  % power dissipation at maximum temperature in Watt;  FIXME: Adjust value
  HW.Grad(1).CoilTemperatur(iCoilGroup) = 20;  % resting temperature in degrees C
  HW.Grad(1).CoilMaxTemperature(iCoilGroup) = 60;  % maximum temperature in degrees C
  heatCapacityAlu = 0.9e3; % J/kg/K
  massAlu = 1.2e3; % kg
  HW.Grad(1).CoilThermalCapacity(iCoilGroup) = heatCapacityAlu*massAlu;  % thermal capacity in J/K;  FIXME: Adjust value
  % properties for "fuse blow" model
  HW.Grad(iDevice).CoilMaxDcCurrent(5) = 250 * 0.9;  % maximum (nominal) DC current of fuse in Ampere
  HW.Grad(iDevice).CoilCurrentSquareTime(5) = 0.9 * 0.9;  % maximum "accumulated heat" in A^2*sec

  HW.DefSeqValues.Prepol.useGradSPI = true;  % use SPI signal to switch prepolarization coil
  HW.DefSeqValues.Prepol.useGrad = false;  % use signal on gradient channel to switch prepolarization coil
  HW.DefSeqValues.Prepol.useDigitalIO = false;  % use digital output signal to switch prepolarization coil
  HW.DefSeqValues.Prepol.SPI.Channel = 5;  % gradient channel for SPI (should be 5)
  HW.DefSeqValues.Prepol.SPI.Device = iDevice;  % device number (1 or 2)
  HW.DefSeqValues.Prepol.SPI.tRampUp = 20e-3;  % ramp up time in seconds
  HW.DefSeqValues.Prepol.SPI.ampPrepol = 20e-3;  % amplitude for prepolarization in T
  HW.DefSeqValues.Prepol.SPI.tRampDownNonAdiabatic = 20e-3;  % non-adiabatic ramp down time in seconds
  HW.DefSeqValues.Prepol.SPI.ampNonAdiabatic = 100e-6;  % amplitude at end of non-adiabatic ramp in T
  HW.DefSeqValues.Prepol.SPI.tRampDownAdiabatic = 4e-3;  % adiabatic ramp down time in seconds
  % FIXME: The postset times of the two IO signals must(!) be different.
  HW.DefSeqValues.Prepol.SPI.TransformerOnOutput = 2^0;  % digital output to turn transformer on
  HW.DefSeqValues.Prepol.SPI.TransformerOnTOffset = 0.15;  % offset time for transformer on signal in seconds
  HW.DefSeqValues.Prepol.SPI.TransformerOnTPostset = 0.15;  % postset time (time between pulse end and transformer off signal) in seconds
  HW.DefSeqValues.Prepol.SPI.ClampOnOutput = 2^1;  % digital output to clamp the prepolarization coil
  HW.DefSeqValues.Prepol.SPI.ClampOffTOffset = 0.12;  % offset time for clamp off signal in seconds
  HW.DefSeqValues.Prepol.SPI.ClampOffTPostset = 10e-3;  % postset time (time between pulse end and clamp on signal) in seconds

  HW.Function_Before_Measurement = @check_Prepol_Amplifier;
  HW.DefSeqValues.Prepol.SPI.ErrorInputDevice = iDevice;  % device for amplifier error signal
  HW.DefSeqValues.Prepol.SPI.ErrorInputChannel = 6;  % digital input channel for amplifier error signal
end

% %% external DDS settings
% DDSConfigExePath = fullfile(winqueryreg('HKEY_LOCAL_MACHINE', 'SOFTWARE\Pure Devices\Settings', 'LibPath'), ...
%   '..', 'PD-DDS', 'DDSConfig.exe');
% DDS = PD.DDS.GetInstance(DDSConfigExePath);
% % DDS.RDACRef = 10e3;  % reference resistance for DAC
% DDS.IDAC2Uout = [DDS.IDACMin; DDS.IDACMax] \ [135e-3; 484e-3];
