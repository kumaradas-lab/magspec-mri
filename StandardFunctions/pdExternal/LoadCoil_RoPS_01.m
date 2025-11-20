%% Load settings for a named coil for the Rocks Profiling System #01
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

HW.GUI.showCoilName = true;  % show name of or selection box for coils

if isempty(HW.TX.CoilName), return; end

if isemptyfield(HW.RecoveryCPMG, 'versionPP'), HW.RecoveryCPMG.versionPP = 1; end

% default values for all coils
HW.TX.Def.PaUout(2) = 70; %98.5;            % default output amplitude in Volts peak
HW.TX.BlankPostset = 80e-9;  % rf amplifier on/off and blank external LNA
HW.TX.DampCoil.AllowShort = true; % set to true to allow a damping pulse to be shortened if necessary

HW.Grad.SystemTimeDelay(1:4) = [5.112e-05, 7.9536e-05, 8.752e-05, 0]; % Time delay of grad amp

% All coils are loop-gap coils.
% The coils have likely been build with 35 um copper foil. The exact height for
% each coil is unknown as of 2024.

% settings specific for each coil
switch HW.TX.CoilName

  case '70mm Slice'
    % damp coil settings
    HW.TX.DampCoil.Enable = true;    % enable coil damping
    HW.TX.DampCoil.DigitalOutputLatency = -0.12e-6;  % latency of damping circuit in seconds
    HW.TX.DampCoil.DigitalOutputDuration = 10e-6; % duration of signal on digital out in seconds
    HW.TX.DampCoil.DampingDuration = 15e-6; % effective duration of damping the coil in seconds
    HW.TX.DampCoil.DigitalOutputChannel = 2; % Digital Output channel
    HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 1e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)

    % Maximum image volume
    HW.Grad.ImageVol = [-0.04, 0.04, -0.010, 0.010, -0.04, 0.04]; % [xmin xmax ymin ymax zmin zmax]
    HW.Grad.ImageVolOffset = [0, 0, 0]; % offset of coil to gradient system

    % coil efficiency
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000003947]; % 2018-11-02T12:07:41, factor voltage amplitude to B1 field strength on the coil output (TX) from CPMG Echo train
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004576]; % 2019-03-04T10:11:01 (p180 = 24.895 탎 @ 103.073 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004408]; % 2019-03-07T12:48:17 (p180 = 21.005 탎 @ 126.817 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004388]; % 2019-03-07T13:54:09 (p180 = 25.000 탎 @ 107.057 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004484]; % 2019-03-07T14:19:51 (p180 = 25.000 탎 @ 104.762 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;

    HW.TX.BlankPostsetAQ = 9.5e-6;  % blank internal LNA
    HW.TX.BlankAQ = 1;
    % reference sample amplitude (100% water)
    % FIXME: Add sensible settings

    % tEchoMin = 80e-6;         % minimum Echo time that is desired for this coil
    HW.RecoveryCPMG.tAQEchoDelay = 6e-6;    % rise time of RX amplitude in s
    % HW.RecoveryCPMG.tFlip180Def = (tEchoMin - 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay))/(1+1.4); % desired pulse length for 180 degrees pulse
    HW.RecoveryCPMG.tFlip180Def = 25e-6; % desired pulse length for 180 degrees pulse
    HW.RecoveryCPMG.tAQEcho = 34e-6;
    HW.RecoveryCPMG.tEchoMin = round((HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + ...
      2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample;

    % coil heating model
    % copper coil
    coilLength = 78e-3;  % estimated length of coil in meter
    foilThickness = 35e-6;  % thickness of conductor in meter
    coilDiameter = 70e-3;  % diameter of coil in meter
    densityCopper = 8.85e3;  % density of copper in kg/m^3
    specificHeatCapacityCopper = 385;  % specific heat capacity of copper in J/kg/K
    massCoil = coilDiameter * pi * coilLength * foilThickness * densityCopper;
    % capacitors
    heatCapacity = 29 * 4000 * (2.8e-3^2*2.5e-3) * 750;  % heat capacity of capacitors in J/K
    % surrounding PTFE
    heatThickness = 1e-3;
    densityPTFE = 2.2e3;  % density of copper in kg/m^3
    massPTFE = coilDiameter * pi * coilLength * heatThickness * densityPTFE;
    specificHeatCapacityPTFE = 1.5e3;  % specific heat capacity of PTFE in J/kg/K
    HW.TX.CoilThermalCapacity = massCoil * specificHeatCapacityCopper + heatCapacity + massPTFE * specificHeatCapacityPTFE;
    HW.RecoveryCPMG.CoilPowerDissipationHigh = 10;  % power dissipation in Watt (with air cooling)
    HW.RecoveryCPMG.CoilPowerDissipationLow = 2;  % power dissipation in Watt (without air cooling)
    HW.TX.CoilPowerDissipation = HW.RecoveryCPMG.CoilPowerDissipationHigh;  % power dissipation in Watt (default)


  case '72mm Slice'
    % damp coil settings
    HW.TX.DampCoil.Enable = true;    % enable coil damping
    HW.TX.DampCoil.DigitalOutputLatency = -1e-6;  % latency of damping circuit in seconds
    HW.TX.DampCoil.DigitalOutputDuration = 8e-6;  % duration of signal on digital out in seconds
    HW.TX.BlankPostsetAQ = 6e-6;  % blank internal LNA
    HW.TX.BlankPostset = 3.8e-6;  % rf amplifier on/off and blank external LNA
    HW.TX.BlankAQ = 1;

    HW.TX.DampCoil.DigitalOutputChannel = 2; % Digital Output channel

    % Maximum image volume
    HW.Grad.ImageVol = [-0.04, 0.04, -0.010, 0.010, -0.04, 0.04]; % [xmin xmax ymin ymax zmin zmax]
    HW.Grad.ImageVolOffset = [0, 0, 0]; % offset of coil to gradient system

    % coil efficiency
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004484]; % 2019-03-07T14:19:51 (p180 = 25.000 탎 @ 104.762 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000003244]; % 2019-04-01T13:50:39 (p180 = 30.000 탎 @ 120.685 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration

    HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;

    % reference sample amplitude (100% water)
    % FIXME: Add sensible settings

    % tEchoMin = 80e-6;         % minimum Echo time that is desired for this coil
    HW.RecoveryCPMG.tAQEchoDelay = 6e-6;    % rise time of RX amplitude in s
    % HW.RecoveryCPMG.tFlip180Def = (tEchoMin - 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay))/(1+1.4); % desired pulse length for 180 degrees pulse

    if HW.RecoveryCPMG.versionPP == 1
      HW.TX.DampCoil.DampingDuration = 9e-6;  % effective duration of damping the coil in seconds
      HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 6e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
      HW.RecoveryCPMG.tFlip180Def = 28e-6;
      HW.RecoveryCPMG.tAQEcho = 34e-6;
      HW.RecoveryCPMG.tEchoMin = round((HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + ...
        2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample;
    elseif HW.RecoveryCPMG.versionPP == 2
      HW.TX.DampCoil.DampingDuration = 10e-6;  % effective duration of damping the coil in seconds
      HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 6e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
      HW.RecoveryCPMG.tFlip180Def = 27e-6;
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/4 *1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 7e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(80e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    elseif HW.RecoveryCPMG.versionPP == 2.5
      HW.RecoveryCPMG.tFlip180Def = 27e-6;
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/2 *1e6)/1e6;
      HW.TX.DampCoil.DigitalOutputDuration = 10e-6;  % duration of signal on digital out in seconds
      HW.TX.DampCoil.DampingDuration = 13e-6;  % effective duration of damping the coil in seconds
      HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 6e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
      HW.RecoveryCPMG.tAQEchoDelay = 8e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(92e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    else
      error('LoadCoil:RecoveryCPMG:UnknownVersionPP', ...
        'Unknown version (%d) for sequence_RecoveryCPMG for selected coil.', HW.RecoveryCPMG.versionPP);
    end

    % coil heating model
    % copper coil
    coilLength = 78e-3;  % estimated length of coil in meter
    foilThickness = 35e-6;  % thickness of conductor in meter
    coilDiameter = 72e-3;  % diameter of coil in meter
    densityCopper = 8.85e3;  % density of copper in kg/m^3
    specificHeatCapacityCopper = 385;  % specific heat capacity of copper in J/kg/K
    massCoil = coilDiameter * pi * coilLength * foilThickness * densityCopper;
    % capacitors
    heatCapacity = 29 * 4000 * (2.8e-3^2*2.5e-3) * 750;  % heat capacity of capacitors in J/K
    % surrounding PTFE
    heatThickness = 1e-3;
    densityPTFE = 2.2e3;  % density of copper in kg/m^3
    massPTFE = coilDiameter * pi * coilLength * heatThickness * densityPTFE;
    specificHeatCapacityPTFE = 1.5e3;  % specific heat capacity of PTFE in J/kg/K
    HW.TX.CoilThermalCapacity = massCoil * specificHeatCapacityCopper + heatCapacity + massPTFE * specificHeatCapacityPTFE;
    HW.RecoveryCPMG.CoilPowerDissipationHigh = 10;  % power dissipation in Watt (with air cooling)
    HW.RecoveryCPMG.CoilPowerDissipationLow = 2;  % power dissipation in Watt (without air cooling)
    HW.TX.CoilPowerDissipation = HW.RecoveryCPMG.CoilPowerDissipationHigh;  % power dissipation in Watt (default)


  case '72mm Slice H-free'
    % damp coil settings
    HW.TX.DampCoil.Enable = true;    % enable coil damping
    HW.TX.DampCoil.DigitalOutputLatency = 0e-6;  % latency of damping circuit in seconds
    HW.TX.DampCoil.DigitalOutputDuration = 10e-6;  % duration of signal on digital out in seconds
    HW.TX.DampCoil.DampingDuration = 13e-6;  % effective duration of damping the coil in seconds
    HW.TX.BlankPostsetAQ = 5e-6;  % blank internal LNA
    HW.TX.BlankPostset = 4e-6;  % rf amplifier on/off and blank external LNA
    HW.TX.BlankAQ = 1;
    HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 9e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)

    HW.TX.DampCoil.DigitalOutputChannel = 2; % Digital Output channel

    % Maximum image volume
    HW.Grad.ImageVol = [-0.04, 0.04, -0.010, 0.010, -0.04, 0.04]; % [xmin xmax ymin ymax zmin zmax]
    HW.Grad.ImageVolOffset = [0, 0, 0]; % offset of coil to gradient system

    % coil efficiency
    % Schirm aus Titanblech
    % Q_ungebremst = 50; Q_gebremst = 6.2;
    % HW.TX.PaUout2Amplitude = [0.000031379, 0.000002705]; % 2021-01-21T10:13:24 (p180 = 33.000 탎 @ 131.560 V) from CPMG Echo train by AutoPulseDuration_CPMG
    % HW.TX.PaUout2Amplitude = [0.000031379, 0.000002672]; % 2021-01-21T11:04:38 (p180 = 31.000 탎 @ 141.797 V) from CPMG Echo train by AutoPulseDuration_CPMG
    % HW.TX.PaUout2Amplitude = [0.000031379, 0.000002652]; % 2021-01-21T12:03:03 (p180 = 31.000 탎 @ 142.822 V) from CPMG Echo train by AutoPulseDuration_CPMG

    % Kupfer-Schirm, Transistor getauscht, Kopplung der Bremse reduziert (Geometrie)
    % Q_ungebremst = 110; Q_gebremst = 16.5;
    % new pulse program
    % HW.TX.PaUout2Amplitude = [0.000031379, 0.000003762];  % 2021-03-30T10:08:56 (p180 = 26.000 탎 @ 120.071 V) from CPMG Echo train by AutoPulseDuration_CPMG
    % old pulse program without slice
    % HW.TX.PaUout2Amplitude = [0.000031379, 0.000003875];  % 2021-03-30T11:24:04 (p180 = 26.000 탎 @ 116.573 V) from CPMG Echo train by AutoPulseDuration_CPMG
    % old pulse program with slice
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004093];  % 2021-03-30T11:24:54 (p180 = 26.000 탎 @ 110.348 V) from CPMG Echo train by AutoPulseDuration_CPMG

    HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;

    % reference sample amplitude (100% water)
    % FIXME: Add sensible settings

    % tEchoMin = 80e-6;         % minimum Echo time that is desired for this coil
    % HW.RecoveryCPMG.tFlip180Def = (tEchoMin - 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay))/(1+1.4); % desired pulse length for 180 degrees pulse
    HW.RecoveryCPMG.tFlip180Def = 26e-6;

    if HW.RecoveryCPMG.versionPP == 1
      HW.RecoveryCPMG.tAQEcho = 26e-6;
      HW.RecoveryCPMG.tAQEchoDelay = 6e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(80e-6, round((HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + ...
        2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample);
    elseif HW.RecoveryCPMG.versionPP == 2
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/4 *1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 6e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(80e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    elseif HW.RecoveryCPMG.versionPP == 2.5
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/2 *1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 6e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(100e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    else
      error('LoadCoil:RecoveryCPMG:UnknownVersionPP', ...
        'Unknown version (%d) for sequence_RecoveryCPMG for selected coil.', HW.RecoveryCPMG.versionPP);
    end

    % coil heating model
    % copper coil
    coilLength = 78e-3;  % estimated length of coil in meter
    foilThickness = 35e-6;  % thickness of conductor in meter
    coilDiameter = 72e-3;  % diameter of coil in meter
    densityCopper = 8.85e3;  % density of copper in kg/m^3
    specificHeatCapacityCopper = 385;  % specific heat capacity of copper in J/kg/K
    massCoil = coilDiameter * pi * coilLength * foilThickness * densityCopper;
    % capacitors
    heatCapacity = 29 * 4000 * (2.8e-3^2*2.5e-3) * 750;  % heat capacity of capacitors in J/K
    % surrounding PTFE
    heatThickness = 1e-3;
    densityPTFE = 2.2e3;  % density of copper in kg/m^3
    massPTFE = coilDiameter * pi * coilLength * heatThickness * densityPTFE;
    specificHeatCapacityPTFE = 1.5e3;  % specific heat capacity of PTFE in J/kg/K
    HW.TX.CoilThermalCapacity = massCoil * specificHeatCapacityCopper + heatCapacity + massPTFE * specificHeatCapacityPTFE;
    HW.RecoveryCPMG.CoilPowerDissipationHigh = 10;  % power dissipation in Watt (with air cooling)
    HW.RecoveryCPMG.CoilPowerDissipationLow = 2;  % power dissipation in Watt (without air cooling)
    HW.TX.CoilPowerDissipation = HW.RecoveryCPMG.CoilPowerDissipationHigh;  % power dissipation in Watt (default)


  case '50mm Slice'
    % damp coil settings
    HW.TX.DampCoil.Enable = true;    % enable coil damping
    HW.TX.DampCoil.DigitalOutputLatency = -0.12e-6;  % latency of damping circuit in seconds
    HW.TX.DampCoil.DigitalOutputDuration = 5e-6; % duration of signal on digital out in seconds
    HW.TX.DampCoil.DampingDuration = 8e-6; % effective duration of damping the coil in seconds
    HW.TX.DampCoil.DigitalOutputChannel = 2; % Digital Output channel
    HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 2e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    HW.TX.BlankPostsetAQ = 4e-6;  % blank internal LNA
    HW.TX.BlankAQ = 1;

    % Maximum image volume
    HW.Grad.ImageVol = [-0.03, 0.03, -0.010, 0.010, -0.03, 0.03]; % [xmin xmax ymin ymax zmin zmax]
    HW.Grad.ImageVolOffset = [0, 0, 0]; % offset of coil to gradient system

    % coil efficiency
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000005284]; % 2018-11-05T17:08:46, factor voltage amplitude to B1 field strength on the coil output (TX)
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000005642]; % 2019-03-07T14:16:24 (p180 = 20.000 탎 @ 104.062 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;

    HW.TX.Def.PaUout(2) = 90;
    % reference sample amplitude (100% water)
    % FIXME: Add sensible settings

    % tEchoMin = 70e-6;         % minimum Echo time that is desired for this coil
    HW.RecoveryCPMG.tAQEchoDelay = 5e-6;    % rise time of RX amplitude in s
    % HW.RecoveryCPMG.tFlip180Def = (tEchoMin - 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay))/(1+pi/2); % desired pulse length for 180 degrees pulse
    HW.RecoveryCPMG.tFlip180Def = 20e-6; % desired pulse length for 180 degrees pulse
    HW.RecoveryCPMG.tAQEcho = 36e-6;
    HW.RecoveryCPMG.tEchoMin = round((HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + ...
      2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample;

    % coil heating model
    % copper coil
    coilLength = 58e-3;  % estimated length of coil in meter
    foilThickness = 35e-6;  % thickness of conductor in meter
    coilDiameter = 50e-3;  % diameter of coil in meter
    densityCopper = 8.85e3;  % density of copper in kg/m^3
    specificHeatCapacityCopper = 385;  % specific heat capacity of copper in J/kg/K
    massCoil = coilDiameter * pi * coilLength * foilThickness * densityCopper;
    % capacitors
    heatCapacity = 21 * 4000 * (2.8e-3^2*2.5e-3) * 750;  % heat capacity of capacitors in J/K
    % surrounding PTFE
    heatThickness = 1e-3;
    densityPTFE = 2.2e3;  % density of copper in kg/m^3
    massPTFE = coilDiameter * pi * coilLength * heatThickness * densityPTFE;
    specificHeatCapacityPTFE = 1.5e3;  % specific heat capacity of PTFE in J/kg/K
    HW.TX.CoilThermalCapacity = massCoil * specificHeatCapacityCopper + heatCapacity + massPTFE * specificHeatCapacityPTFE;
    HW.RecoveryCPMG.CoilPowerDissipationHigh = 10;  % power dissipation in Watt (with air cooling)
    HW.RecoveryCPMG.CoilPowerDissipationLow = 2;  % power dissipation in Watt (without air cooling)
    HW.TX.CoilPowerDissipation = HW.RecoveryCPMG.CoilPowerDissipationHigh;  % power dissipation in Watt (default)


  case '52mm Slice'
    % damp coil settings
    HW.TX.DampCoil.Enable = true;    % enable coil damping
    HW.TX.DampCoil.DigitalOutputLatency = -1e-6;  % latency of damping circuit in seconds
%     HW.TX.DampCoil.DigitalOutputDuration = 5.5e-6; % duration of signal on digital out in seconds
    HW.TX.DampCoil.DigitalOutputDuration = 9e-6; % duration of signal on digital out in seconds
    HW.TX.DampCoil.DampingDuration = 10e-6; % effective duration of damping the coil in seconds
    HW.TX.DampCoil.DigitalOutputChannel = 2; % Digital Output channel
    HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 7e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    HW.TX.BlankPostset = 3e-6;  % rf amplifier on/off and blank external LNA
    HW.TX.BlankPostsetAQ = 4.4e-6;  % blank internal LNA
    HW.TX.BlankAQ = 1;

    % Maximum image volume
    HW.Grad.ImageVol = [-0.03, 0.03, -0.010, 0.010, -0.03, 0.03]; % [xmin xmax ymin ymax zmin zmax]
    HW.Grad.ImageVolOffset = [0, 0, 0]; % offset of coil to gradient system

    % coil efficiency
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000005469]; % 2019-03-19T10:31:19 (p180 = 20.000 탎 @ 107.363 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004762]; % 2019-04-01T09:55:58 (p180 = 20.000 탎 @ 123.308 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004641]; % 2019-04-01T13:57:07 (p180 = 20.000 탎 @ 126.510 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration

    HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;

    HW.TX.Def.PaUout(2) = 90;
    % reference sample amplitude (100% water)
    % FIXME: Add sensible settings

    % tEchoMin = 70e-6;         % minimum Echo time that is desired for this coil
    HW.RecoveryCPMG.tAQEchoDelay = 5e-6;    % rise time of RX amplitude in s
    % HW.RecoveryCPMG.tFlip180Def = (tEchoMin - 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay))/(1+pi/2); % desired pulse length for 180 degrees pulse
    % desired pulse length for 180 degrees pulse
    if HW.RecoveryCPMG.versionPP == 1
      HW.RecoveryCPMG.tFlip180Def = 20e-6;
      HW.RecoveryCPMG.tAQEcho = 26e-6;
      HW.RecoveryCPMG.tEchoMin = round((HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + ...
        2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample;
    elseif HW.RecoveryCPMG.versionPP == 2
      HW.RecoveryCPMG.tFlip180Def = 23.6e-6-80e-9; % FIXME: AQ.fSample?
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/4 *1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 9e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(70e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    elseif HW.RecoveryCPMG.versionPP == 2.5
      HW.RecoveryCPMG.tFlip180Def = 22e-6;
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/2 *1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 9.5e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(80e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    else
      error('LoadCoil:RecoveryCPMG:UnknownVersionPP', ...
        'Unknown version (%d) for sequence_RecoveryCPMG for selected coil.', HW.RecoveryCPMG.versionPP);
    end

    % coil heating model
    % copper coil
    coilLength = 58e-3;  % estimated length of coil in meter
    foilThickness = 35e-6;  % thickness of conductor in meter
    coilDiameter = 52e-3;  % diameter of coil in meter
    densityCopper = 8.85e3;  % density of copper in kg/m^3
    specificHeatCapacityCopper = 385;  % specific heat capacity of copper in J/kg/K
    massCoil = coilDiameter * pi * coilLength * foilThickness * densityCopper;
    % capacitors
    heatCapacity = 21 * 4000 * (2.8e-3^2*2.5e-3) * 750;  % heat capacity of capacitors in J/K
    % surrounding PTFE
    heatThickness = 1e-3;
    densityPTFE = 2.2e3;  % density of copper in kg/m^3
    massPTFE = coilDiameter * pi * coilLength * heatThickness * densityPTFE;
    specificHeatCapacityPTFE = 1.5e3;  % specific heat capacity of PTFE in J/kg/K
    HW.TX.CoilThermalCapacity = massCoil * specificHeatCapacityCopper + heatCapacity + massPTFE * specificHeatCapacityPTFE;
    HW.RecoveryCPMG.CoilPowerDissipationHigh = 10;  % power dissipation in Watt (with air cooling)
    HW.RecoveryCPMG.CoilPowerDissipationLow = 2;  % power dissipation in Watt (without air cooling)
    HW.TX.CoilPowerDissipation = HW.RecoveryCPMG.CoilPowerDissipationHigh;  % power dissipation in Watt (default)


  case '42mm Slice'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% damp coil settings
    HW.TX.DampCoil.Enable = true;    % enable coil damping
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HW.TX.DampCoil.DigitalOutputLatency = -1e-6;  % latency of damping circuit in seconds
    HW.TX.DampCoil.DigitalOutputDuration = 9.2e-6; % duration of signal on digital out in seconds
    HW.TX.DampCoil.DampingDuration = 10e-6; % effective duration of damping the coil in seconds
    HW.TX.DampCoil.DigitalOutputChannel = 2; % Digital Output channel
    HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 9e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    HW.TX.BlankPostsetAQ = 8.4e-6;  % blank internal LNA
    HW.TX.BlankPostset = 5e-6;  % rf amplifier on/off and blank external LNA
    HW.TX.BlankAQ = 1;

    % Maximum image volume
    HW.Grad.ImageVol = [-0.03, 0.03, -0.010, 0.010, -0.03, 0.03]; % [xmin xmax ymin ymax zmin zmax]
    HW.Grad.ImageVolOffset = [0, 0, 0]; % offset of coil to gradient system

    % coil efficiency
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000005568];
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000006180]; % 2019-03-07T15:04:04 (p180 = 18.000 탎 @ 105.561 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000006184]; % 2019-03-07T15:15:13 (p180 = 18.000 탎 @ 105.503 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000005376]; % 2019-04-01T14:04:07 (p180 = 18.000 탎 @ 121.366 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000005317]; % 2019-05-21T10:38:15 (p180 = 18.000 탎 @ 122.713 V) from CPMG Echo train by AutoPulseDuration_CPMG

    HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;
    % reference sample amplitude (100% water)
    % FIXME: Add sensible settings

    % rf amplifier settings
    % HW.TX.Def.PaUout(2) = 70;            % default output amplitude in Volts peak

    % HW.RecoveryCPMG.tEchoMin = 60e-6;  % minimum Echo time that is desired for this coil
    % HW.RecoveryCPMG.tAQEchoDelay = 6e-6;  % rise time of RX amplitude in s
    % HW.RecoveryCPMG.tFlip180Def = (tEchoMin - 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay))/(1+pi/2); % desired pulse length for 180 degrees pulse
    if HW.RecoveryCPMG.versionPP == 1
      HW.RecoveryCPMG.tFlip180Def = 18e-6;  % desired pulse length for 180 degrees pulse
      HW.RecoveryCPMG.tAQEcho = 18e-6;
      HW.RecoveryCPMG.tAQEchoDelay = 7e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = round((HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + ...
        2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample;
    elseif HW.RecoveryCPMG.versionPP == 2
      HW.RecoveryCPMG.tFlip180Def = 18.3e-6;
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/4 *1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 7e-6;5e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(60e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    elseif HW.RecoveryCPMG.versionPP == 2.5
      HW.RecoveryCPMG.tFlip180Def = 18e-6;
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/2 *1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 7.5e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(60e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    else
      error('LoadCoil:RecoveryCPMG:UnknownVersionPP', ...
        'Unknown version (%d) for sequence_RecoveryCPMG for selected coil.', HW.RecoveryCPMG.versionPP);
    end

    % coil heating model
    % copper coil
    coilLength = 38e-3;  % estimated length of coil in meter
    foilThickness = 35e-6;  % thickness of conductor in meter
    coilDiameter = 42e-3;  % diameter of coil in meter
    densityCopper = 8.85e3;  % density of copper in kg/m^3
    specificHeatCapacityCopper = 385;  % specific heat capacity of copper in J/kg/K
    massCoil = coilDiameter * pi * coilLength * foilThickness * densityCopper;
    % capacitors
    heatCapacity = 14 * 4000 * (2.8e-3^2*2.5e-3) * 750;  % heat capacity of capacitors in J/K
    % surrounding PTFE
    heatThickness = 1e-3;
    densityPTFE = 2.2e3;  % density of copper in kg/m^3
    massPTFE = coilDiameter * pi * coilLength * heatThickness * densityPTFE;
    specificHeatCapacityPTFE = 1.5e3;  % specific heat capacity of PTFE in J/kg/K
    HW.TX.CoilThermalCapacity = massCoil * specificHeatCapacityCopper + heatCapacity + massPTFE * specificHeatCapacityPTFE;
    HW.RecoveryCPMG.CoilPowerDissipationHigh = 10;  % power dissipation in Watt (with air cooling)
    HW.RecoveryCPMG.CoilPowerDissipationLow = 2;  % power dissipation in Watt (without air cooling)
    HW.TX.CoilPowerDissipation = HW.RecoveryCPMG.CoilPowerDissipationHigh;  % power dissipation in Watt (default)


  case '42mm Image'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% damp coil settings
    HW.TX.DampCoil.Enable = true;    % enable coil damping
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HW.TX.DampCoil.DigitalOutputLatency = -1e-6;  % latency of damping circuit in seconds
    HW.TX.DampCoil.DigitalOutputDuration = 3.0e-6; % duration of signal on digital out in seconds
    HW.TX.DampCoil.DampingDuration = 9e-6; % effective duration of damping the coil in seconds
    HW.TX.DampCoil.DigitalOutputChannel = 2; % Digital Output channel
    HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 12e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    HW.TX.BlankPostsetAQ = 2.5e-6;  % blank internal LNA
    HW.TX.BlankPostset = 2e-6;  % rf amplifier on/off and blank external LNA

    % Maximum image volume
    HW.Grad.ImageVol = [-0.03, 0.03, -0.035, 0.035, -0.03, 0.03]; % [xmin xmax ymin ymax zmin zmax]
    HW.Grad.ImageVolOffset = [0, 0, 0]; % offset of coil to gradient system

    % coil efficiency
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004629];
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004520]; % 2019-03-07T15:37:04 (p180 = 20.000 탎 @ 129.901 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004520]; % 2019-03-07T15:45:45 (p180 = 20.600 탎 @ 126.107 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000004094]; % 2019-04-02T09:58:14 (p180 = 23.600 탎 @ 121.530 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;
    % reference sample amplitude (100% water)
    % FIXME: Add sensible settings

    % rf amplifier settings
    % HW.TX.Def.PaUout(2) = 98.5;            % default output amplitude in Volts peak

    % HW.RecoveryCPMG.tEchoMin = 60e-6;  % minimum Echo time that is desired for this coil
    if HW.RecoveryCPMG.versionPP == 1
      HW.RecoveryCPMG.tAQEchoDelay = 9e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tFlip180Def = 22e-6;
      HW.RecoveryCPMG.tAQEcho = 14e-6;
      HW.RecoveryCPMG.tEchoMin = round((HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + ...
        2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample;
    elseif HW.RecoveryCPMG.versionPP == 2
      HW.RecoveryCPMG.tAQEchoDelay = 9.5e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tFlip180Def = 20.5e-6;
      HW.RecoveryCPMG.tAQEcho = 14e-6;
      HW.RecoveryCPMG.tAQEchoDelay = 8.4e-6;  % rise time of RX amplitude in s
      HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 11e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
      HW.RecoveryCPMG.tEchoMin = max(60e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample);
    elseif HW.RecoveryCPMG.versionPP == 2.5
      HW.RecoveryCPMG.tAQEchoDelay = 8.5e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tFlip180Def = 20.5e-6;
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/2 * 1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 10e-6;  % rise time of RX amplitude in s
      HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 11e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
      tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/2 * 1e6)/1e6;
      HW.RecoveryCPMG.tEchoMin = max(77e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample);
    else
      error('LoadCoil:RecoveryCPMG:UnknownVersionPP', ...
        'Unknown version (%d) for sequence_RecoveryCPMG for selected coil.', HW.RecoveryCPMG.versionPP);
    end

    % coil heating model
    % copper coil
    coilLength = 38e-3;  % estimated length of coil in meter
    foilThickness = 35e-6;  % thickness of conductor in meter
    coilDiameter = 42e-3;  % diameter of coil in meter
    densityCopper = 8.85e3;  % density of copper in kg/m^3
    specificHeatCapacityCopper = 385;  % specific heat capacity of copper in J/kg/K
    massCoil = coilDiameter * pi * coilLength * foilThickness * densityCopper;
    % capacitors
    heatCapacity = 14 * 4000 * (2.8e-3^2*2.5e-3) * 750;  % heat capacity of capacitors in J/K
    % surrounding PTFE
    heatThickness = 1e-3;
    densityPTFE = 2.2e3;  % density of copper in kg/m^3
    massPTFE = coilDiameter * pi * coilLength * heatThickness * densityPTFE;
    specificHeatCapacityPTFE = 1.5e3;  % specific heat capacity of PTFE in J/kg/K
    HW.TX.CoilThermalCapacity = massCoil * specificHeatCapacityCopper + heatCapacity + massPTFE * specificHeatCapacityPTFE;
    HW.RecoveryCPMG.CoilPowerDissipationHigh = 10;  % power dissipation in Watt (with air cooling)
    HW.RecoveryCPMG.CoilPowerDissipationLow = 2;  % power dissipation in Watt (without air cooling)
    HW.TX.CoilPowerDissipation = HW.RecoveryCPMG.CoilPowerDissipationHigh;  % power dissipation in Watt (default)


  case '32mm Image'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% damp coil settings
    HW.TX.DampCoil.Enable = true;    % enable coil damping
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HW.TX.DampCoil.DigitalOutputLatency = -1e-6;  % latency of damping circuit in seconds
    HW.TX.DampCoil.DigitalOutputDuration = 5.5e-6; % duration of signal on digital out in seconds
    HW.TX.DampCoil.DampingDuration = 7.2e-6; % effective duration of damping the coil in seconds
    HW.TX.DampCoil.DigitalOutputChannel = 2; % Digital Output channel
    % HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 10.3e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 10.3e-6 - 3.6e-6/2; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    HW.TX.BlankPostsetAQ = 4e-6;  % blank internal LNA
    HW.TX.BlankPostset = 2e-6;  % rf amplifier on/off and blank external LNA

    % Maximum image volume
    HW.Grad.ImageVol = [-0.03, 0.03, -0.035, 0.035, -0.03, 0.03]; % [xmin xmax ymin ymax zmin zmax]
    HW.Grad.ImageVolOffset = [0, 0, 0]; % offset of coil to gradient system

    % coil efficiency
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000007617];
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000008017]; % 2019-03-07T15:51:26 (p180 = 18.000 탎 @ 81.376 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000008079]; % 2019-03-07T16:20:22 (p180 = 16.500 탎 @ 88.092 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000007038]; % 2019-04-01T14:23:07 (p180 = 16.600 탎 @ 100.513 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    % HW.TX.PaUout2Amplitude = [0.000031379, 0.000006963]; % 2019-04-17T09:12:42 (p180 = 14.000 탎 @ 120.472 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000006955]; % 2019-04-17T09:26:46 (p180 = 13.000 탎 @ 129.888 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration

    HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;
    % reference sample amplitude (100% water)
    % FIXME: Add sensible settings

    % rf amplifier settings
    % HW.TX.Def.PaUout(2) = 98.5;            % default output amplitude in Volts peak

    % HW.RecoveryCPMG.tEchoMin = 50e-6;  % minimum Echo time that is desired for this coil
    % HW.RecoveryCPMG.tAQEchoDelay = 6e-6;  % rise time of RX amplitude in s
    % HW.RecoveryCPMG.tFlip180Def = (tEchoMin - 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay))/(1+pi/2); % desired pulse length for 180 degrees pulse
    HW.RecoveryCPMG.tAQEchoDelay = 5e-6;    % rise time of RX amplitude in s
    if HW.RecoveryCPMG.versionPP == 1
      HW.RecoveryCPMG.tFlip180Def = 16.6e-6;
      HW.RecoveryCPMG.tAQEcho = 12e-6;
      HW.RecoveryCPMG.tEchoMin = round((HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + ...
        2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample;
      % FIXME: double precision error
      HW.RecoveryCPMG.tEchoMin = 50e-6;
    elseif HW.RecoveryCPMG.versionPP == 2
      HW.RecoveryCPMG.tFlip180Def = 15.3e-6;
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/4 *1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 4e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(50e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    elseif HW.RecoveryCPMG.versionPP == 2.5
      HW.RecoveryCPMG.tFlip180Def = 14e-6;
      HW.TX.DampCoil.DigitalOutputDuration = 6e-6;  % duration of signal on digital out in seconds
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/2 *1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 5.5e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(0*56e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    else
      error('LoadCoil:RecoveryCPMG:UnknownVersionPP', ...
        'Unknown version (%d) for sequence_RecoveryCPMG for selected coil.', HW.RecoveryCPMG.versionPP);
    end

    % coil heating model
    % copper coil
    coilLength = 38e-3;  % estimated length of coil in meter
    foilThickness = 35e-6;  % thickness of conductor in meter
    coilDiameter = 32e-3;  % diameter of coil in meter
    densityCopper = 8.85e3;  % density of copper in kg/m^3
    specificHeatCapacityCopper = 385;  % specific heat capacity of copper in J/kg/K
    massCoil = coilDiameter * pi * coilLength * foilThickness * densityCopper;
    % capacitors
    heatCapacity = 14 * 4000 * (2.8e-3^2*2.5e-3) * 750;  % heat capacity of capacitors in J/K
    % surrounding PTFE
    heatThickness = 1e-3;
    densityPTFE = 2.2e3;  % density of copper in kg/m^3
    massPTFE = coilDiameter * pi * coilLength * heatThickness * densityPTFE;
    specificHeatCapacityPTFE = 1.5e3;  % specific heat capacity of PTFE in J/kg/K
    HW.TX.CoilThermalCapacity = massCoil * specificHeatCapacityCopper + heatCapacity + massPTFE * specificHeatCapacityPTFE;
    HW.RecoveryCPMG.CoilPowerDissipationHigh = 10;  % power dissipation in Watt (with air cooling)
    HW.RecoveryCPMG.CoilPowerDissipationLow = 2;  % power dissipation in Watt (without air cooling)
    HW.TX.CoilPowerDissipation = HW.RecoveryCPMG.CoilPowerDissipationHigh;  % power dissipation in Watt (default)


  case '32mm Combi'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% damp coil settings
    HW.TX.DampCoil.Enable = true;    % enable coil damping
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HW.TX.DampCoil.DigitalOutputLatency = 0e-6;  % latency of damping circuit in seconds
    HW.TX.DampCoil.DigitalOutputDuration = 6e-6;  % duration of signal on digital out in seconds
    HW.TX.DampCoil.DampingDuration = 7.6e-6;  % effective duration of damping the coil in seconds
    HW.TX.DampCoil.DigitalOutputChannel = 2;  % Digital Output channel
    % HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 10.3e-6;  % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    % HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 10.3e-6 - 5e-6/2 - 5.4e-6/2;  % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 4e-6;  % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    HW.TX.BlankPostsetAQ = 4e-6;  % blank internal LNA
    HW.TX.BlankPostset = 2e-6;  % rf amplifier on/off and blank external LNA

    % Maximum image volume
    HW.Grad.ImageVol = [-0.03, 0.03, -0.035, 0.035, -0.03, 0.03];  % [xmin xmax ymin ymax zmin zmax]
    HW.Grad.ImageVolOffset = [0, 0, 0];  % offset of coil to gradient system

    % coil efficiency
    %HW.TX.PaUout2Amplitude = [0.000033000, 0.000007983];  % 2020-11-30T13:05:55 by D_emo_Auto_PulseDuration_CPMG, factor voltage amplitude to B1 field strength on the coil output (TX) from CPMG Echo train
    %HW.TX.PaUout2Amplitude = [0.000033000, 0.000006470];  % 2020-12-17T13:55:03 (p180 = 16.000 탎 @ 113.432 V) from CPMG Echo train by AutoPulseDuration_CPMG
    %HW.TX.PaUout2Amplitude = [0.000033000, 0.000006753];  % 2020-12-17T15:13:10 (p90 = 12.421 탎 @ 70.000 V) from 1d Spin Echo by AutoPulseDuration_imaging

    % nach Umbau (Widerstand in Tune-Spule raus)
    % Q_ungebremst ~ 53
    % Q_gebremst ~ 10
    % HW.TX.PaUout2Amplitude = [38.318000, 10.517611]*1e-6; % 2021-02-05T13:51:33 (tFlip90 = 5.583 탎 @ 100.000 V) from 1d Spin Echo by Find_PulseDuration
    % HW.TX.PaUout2Amplitude = [0.000038318, 0.000007626];  % 2021-02-24T09:38:24 (p180 = 16.000 탎 @ 96.249 V) from CPMG Echo train by AutoPulseDuration_CPMG
    % HW.TX.PaUout2Amplitude = [0.000038318, 0.000007580];  % 2021-02-24T11:44:16 (p180 = 16.000 탎 @ 96.825 V) from CPMG Echo train by AutoPulseDuration_CPMG
    % HW.TX.PaUout2Amplitude = [0.000038318, 0.000007516];  % 2021-02-24T12:18:41 (p180 = 16.600 탎 @ 94.124 V) from CPMG Echo train by AutoPulseDuration_CPMG
    HW.TX.PaUout2Amplitude = [0.000038318, 0.000007730];  % 2021-02-24T12:54:26 (p180 = 18.000 탎 @ 84.394 V) from CPMG Echo train by AutoPulseDuration_CPMG

    HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;
    % reference sample amplitude (100% water)
    % FIXME: Add sensible settings

    % rf amplifier settings
    % HW.TX.Def.PaUout(2) = 98.5;  % default output amplitude in Volts peak

    % HW.RecoveryCPMG.tEchoMin = 50e-6;  % minimum Echo time that is desired for this coil
    % HW.RecoveryCPMG.tAQEchoDelay = 6e-6;  % rise time of RX amplitude in s
    % HW.RecoveryCPMG.tFlip180Def = (tEchoMin - 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay))/(1+pi/2);  % desired pulse length for 180 degrees pulse
    % HW.RecoveryCPMG.tAQEchoDelay = 5e-6;  % rise time of RX amplitude in s
    if HW.RecoveryCPMG.versionPP == 1
      HW.RecoveryCPMG.tFlip180Def = 16.6e-6;
      HW.RecoveryCPMG.tAQEcho = 12e-6;
      HW.RecoveryCPMG.tAQEchoDelay = 3e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = round((HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + ...
        2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample;
      % FIXME: double precision error
      HW.RecoveryCPMG.tEchoMin = 50e-6;
    elseif HW.RecoveryCPMG.versionPP == 2
      HW.RecoveryCPMG.tFlip180Def = 18e-6;
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/4 *1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 4e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(50e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    elseif HW.RecoveryCPMG.versionPP == 2.5
      HW.RecoveryCPMG.tFlip180Def = 18e-6;
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/2 *1e6)/1e6;
      HW.RecoveryCPMG.tAQEchoDelay = 4.2e-6;  % rise time of RX amplitude in s
      HW.RecoveryCPMG.tEchoMin = max(60e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    else
      error('LoadCoil:RecoveryCPMG:UnknownVersionPP', ...
        'Unknown version (%d) for sequence_RecoveryCPMG for selected coil.', HW.RecoveryCPMG.versionPP);
    end

    HW.Grad.SystemTimeDelay(1:3) = [7.7712e-05, 0.000111744, 0.000118216];  % Time delay of grad amp

    % coil heating model
    % copper coil
    coilLength = 38e-3;  % estimated length of coil in meter
    foilThickness = 35e-6;  % thickness of conductor in meter
    coilDiameter = 32e-3;  % diameter of coil in meter
    densityCopper = 8.85e3;  % density of copper in kg/m^3
    specificHeatCapacityCopper = 385;  % specific heat capacity of copper in J/kg/K
    massCoil = coilDiameter * pi * coilLength * foilThickness * densityCopper;
    % capacitors
    heatCapacity = 14 * 4000 * (2.8e-3^2*2.5e-3) * 750;  % heat capacity of capacitors in J/K
    % surrounding PTFE
    heatThickness = 1e-3;
    densityPTFE = 2.2e3;  % density of copper in kg/m^3
    massPTFE = coilDiameter * pi * coilLength * heatThickness * densityPTFE;
    specificHeatCapacityPTFE = 1.5e3;  % specific heat capacity of PTFE in J/kg/K
    HW.TX.CoilThermalCapacity = massCoil * specificHeatCapacityCopper + heatCapacity + massPTFE * specificHeatCapacityPTFE;
    HW.RecoveryCPMG.CoilPowerDissipationHigh = 10;  % power dissipation in Watt (with air cooling)
    HW.RecoveryCPMG.CoilPowerDissipationLow = 2;  % power dissipation in Watt (without air cooling)
    HW.TX.CoilPowerDissipation = HW.RecoveryCPMG.CoilPowerDissipationHigh;  % power dissipation in Watt (default)


  case '22mm Image'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% damp coil settings
    HW.TX.DampCoil.Enable = true;    % enable coil damping
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HW.TX.DampCoil.DigitalOutputLatency = -1e-6;  % latency of damping circuit in seconds
    HW.TX.DampCoil.DigitalOutputDuration = 5e-6; % duration of signal on digital out in seconds
    HW.TX.DampCoil.DampingDuration = 6e-6; % effective duration of damping the coil in seconds
    HW.TX.DampCoil.DigitalOutputChannel = 2; % Digital Output channel
    HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 2.7e-6; % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    HW.TX.BlankPostsetAQ = 4e-6;  % blank internal LNA
    HW.TX.BlankPostset = 2.5e-6;  % rf amplifier on/off and blank external LNA


    % Maximum image volume
    HW.Grad.ImageVol = [-0.03, 0.03, -0.035, 0.035, -0.03, 0.03]; % [xmin xmax ymin ymax zmin zmax]
    HW.Grad.ImageVolOffset = [0, 0, 0]; % offset of coil to gradient system

    % coil efficiency
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000007617];
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000014795]; % 2019-03-07T16:25:43 (p180 = 16.500 탎 @ 48.105 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000014866]; % 2019-03-07T16:38:21 (p180 = 18.000 탎 @ 43.885 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000014638]; % 2019-03-07T16:44:59 (p180 = 16.200 탎 @ 49.521 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration
    HW.TX.PaUout2Amplitude = [0.000031379, 0.000011114]; % 2019-04-01T14:32:32 (p180 = 18.000 탎 @ 58.699 V) from CPMG Echo train by D_emo_Auto_PulseDuration_CPMG_fixedDuration

    HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;
    % reference sample amplitude (100% water)
    % FIXME: Add sensible settings

    % rf amplifier settings
    % HW.TX.Def.PaUout(2) = 98.5;            % default output amplitude in Volts peak

    % tEchoMin = 50e-6;  % minimum Echo time that is desired for this coil
    % HW.RecoveryCPMG.tAQEchoDelay = 4e-6;  % rise time of RX amplitude in s
    % HW.RecoveryCPMG.tFlip180Def = (tEchoMin - 2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay))/(1+pi/2); % desired pulse length for 180 degrees pulse
    HW.RecoveryCPMG.tAQEchoDelay = 3e-6;    % rise time of RX amplitude in s
    if HW.RecoveryCPMG.versionPP == 1
      HW.RecoveryCPMG.tFlip180Def = 16.6e-6;
      HW.RecoveryCPMG.tAQEcho = 22e-6;
      HW.RecoveryCPMG.tEchoMin = round((HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + ...
        2*max(max(0,HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay),HW.RecoveryCPMG.tAQEchoDelay)) * HW.TX.fSample)/HW.TX.fSample;
    elseif HW.RecoveryCPMG.versionPP == 2
      HW.RecoveryCPMG.tFlip180Def = 16e-6; %16.3e-6-40e-9;
      HW.RecoveryCPMG.tAQEchoDelay = 3e-6;    % rise time of RX amplitude in s
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/4 * 1e6)/1e6;
      HW.RecoveryCPMG.tEchoMin = max(50e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    elseif HW.RecoveryCPMG.versionPP == 2.5
      HW.RecoveryCPMG.tFlip180Def = 16e-6;
      HW.RecoveryCPMG.tAQEchoDelay = 5e-6;    % rise time of RX amplitude in s
      HW.RecoveryCPMG.tAQEcho = round(HW.RecoveryCPMG.tFlip180Def*3/2 *1e6)/1e6;
      HW.RecoveryCPMG.tEchoMin = max(0*56e-6, round(max(HW.RecoveryCPMG.tFlip180Def*7/4 + HW.RecoveryCPMG.tFlip180Def, ...
        HW.RecoveryCPMG.tFlip180Def + HW.RecoveryCPMG.tAQEcho + 2*max(max(0, HW.TX.DampCoil.TX2RXdeadTime-HW.RecoveryCPMG.tAQEchoDelay), HW.RecoveryCPMG.tAQEchoDelay+2.5e-6)) * HW.TX.fSample)/HW.TX.fSample);
    else
      error('LoadCoil:RecoveryCPMG:UnknownVersionPP', ...
        'Unknown version (%d) for sequence_RecoveryCPMG for selected coil.', HW.RecoveryCPMG.versionPP);
    end

    % coil heating model
    % copper coil
    coilLength = 28e-3;  % estimated length of coil in meter
    foilThickness = 35e-6;  % thickness of conductor in meter
    coilDiameter = 22e-3;  % diameter of coil in meter
    densityCopper = 8.85e3;  % density of copper in kg/m^3
    specificHeatCapacityCopper = 385;  % specific heat capacity of copper in J/kg/K
    massCoil = coilDiameter * pi * coilLength * foilThickness * densityCopper;
    % capacitors
    heatCapacity = 10 * 4000 * (2.8e-3^2*2.5e-3) * 750;  % heat capacity of capacitors in J/K
    % surrounding PTFE
    heatThickness = 1e-3;
    densityPTFE = 2.2e3;  % density of copper in kg/m^3
    massPTFE = coilDiameter * pi * coilLength * heatThickness * densityPTFE;
    specificHeatCapacityPTFE = 1.5e3;  % specific heat capacity of PTFE in J/kg/K
    HW.TX.CoilThermalCapacity = massCoil * specificHeatCapacityCopper + heatCapacity + massPTFE * specificHeatCapacityPTFE;
    HW.RecoveryCPMG.CoilPowerDissipationHigh = 10;  % power dissipation in Watt (with air cooling)
    HW.RecoveryCPMG.CoilPowerDissipationLow = 2;  % power dissipation in Watt (without air cooling)
    HW.TX.CoilPowerDissipation = HW.RecoveryCPMG.CoilPowerDissipationHigh;  % power dissipation in Watt (default)


  otherwise
    warning('LoadCoil:UnknownCoilName', 'No configuration for coil "%s" found. Using default settings.', HW.TX.CoilName)
    return;

end
