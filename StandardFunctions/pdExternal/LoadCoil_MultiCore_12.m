%% Load settings for a named coil for the NMR-SCA with 32 mm bore
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if ~exist('iDevice', 'var')
  iDevice = 1;
end

if isempty(HW.TX(iDevice).CoilName), return; end

% settings specific for each coil
switch HW.TX(iDevice).CoilName

  case '30mm'
    % 30 mm long bronce rf coil (ProbeID 4589343)

    % damp coil settings
    HW.TX(iDevice).DampCoil.Enable = true;  % enable coil damping
    HW.TX(iDevice).DampCoil.DigitalOutputLatency = 0.6e-6;  % latency of damping circuit in seconds
    HW.TX(iDevice).DampCoil.DigitalOutputChannel = 1;  % digital output channel for coil damping signal
    HW.TX(iDevice).DampCoil.DigitalOutputDuration = HW.TX(iDevice).BlankPostset + 1.6e-6;  % duration of signal on digital out in seconds
    HW.TX(iDevice).DampCoil.DampingDuration = HW.TX(iDevice).DampCoil.DigitalOutputDuration + 0.8e-6;  % effective duration of damping the coil in seconds
    HW.TX(iDevice).DampCoil.TX2RXdeadTime = HW.TX(iDevice).DampCoil.DampingDuration + 3e-6;  % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    HW.TX(iDevice).BlankPostsetAQ = 1.4e-6;  % blank internal LNA
    HW.TX(iDevice).BlankAQ = 1;

    % gradient system
    HW.Grad(iDevice).ImageVol = [-0.02, 0.02, -0.020, 0.020, -0.02, 0.02];  % [xmin xmax ymin ymax zmin zmax]
    HW.Grad(iDevice).ImageVolOffset = [0, 0, 0];  % offset of coil to gradient system
    HW.RX(iDevice).EffectiveCoilLength = 30e-3;

    HW.Grad(iDevice).LoadRin(1:3) = [2.7, 3.4, 2.5];  % 30 mm x y z
    HW.Grad(iDevice).LoadIin2Amp(1:3) = [+0.021669, +0.021521, +0.023990];  % 2025-11-28T17:50:03, GainChangeXYZ = [ +0.064% +0.089% -0.027% ], centerXYZ = [ -0.431 +0.176 +0.111 ] mm, DiameterXYZ = [ +10.006 +10.009  +9.997 ] mm, SphereDiameter = +10.000 mm
    HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [79.290, 96.160, 61.236]*1e-6;  % time delay of gradient amplifier in s - calibrated with Calibrate_GradientSystem_GradDelay on 2025-12-08 16:03:50

    HW.MagnetShim([1,2,3]) = [-0.000112432, -0.000283597, -0.000437707];  % 2025-12-08T16:08:14 by FindShim (T2* = 22.8 ms @ 18.364000 MHz), x y z in T/m

    % rf coil efficiency
    HW.TX(iDevice).PaUout2Amplitude = [8.345843, 8.345843]*1e-6;  % 2025-12-01T09:46:08 (tFlip90 = 6.396 us @ 110.000 V @ 18.364000 MHz) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2Amplitude = [7.495798, 7.495798]*1e-6;  % 2025-12-08T16:08:33 (p180 = 16.000 탎 @ 97.916 V @ 18.364000 MHz) from CPMG Echo train by DailyCheck_solidEcho
    HW.TX(iDevice).PaUout2AmplitudeEstimated = HW.TX(iDevice).PaUout2Amplitude;

    % CPMG and solid echo settings
    HW.RecoveryCPMG.tEchoMin = 250e-6;  % minimum echo time that is desired for this coil
    HW.RecoveryCPMG.tFlip180Def = 16e-6;  % pulse length for 180 degrees pulse in seconds
    HW.TX(iDevice).Def.PaUout(2) = 90;  % 90 Volts correspond to a 180 degrees pulse duration of approximately 16e-6 s

    % rf coil heating model
    HW.TX(iDevice).CoilTemperature = 32;  % coil temperature in degrees at start of measurement
    HW.TX(iDevice).CoilMaxTemperature = 60;  % maximum coil temperature in degrees during measurement
    HW.TX(iDevice).CoilPowerDissipation = 2;  % estimated power dissipation to environment in Watt at 60 degrees
    % thermal capacity of coil in J/K
    % specific heat capacity of bronce at 25 degrees C: 0.377 J/(g*K)
    % density of bronce at room temperature: 8.8 g/cm^3
    % volume of bronce in coil: 30 mm * 100 mm * 50 um
    % specific heat capacity of capacitors: 0.45 J/(g*K) (https://www.murata.com/en-eu/support/faqs/capacitor/ceramiccapacitor/conf/0011)
    % weight of capacitors: 11 * 0.08 g
    % factor for solder: 1.5
    % specific heat capacity of Al2O3 tube: 0.90 J/(g*K)
    % density of Al2O3: 3.9 g/cm^3
    % volume of tube: ((33.8 mm/2)^2 - (32.2 mm/2)^2) * pi * 178 mm
    HW.TX(iDevice).CoilThermalCapacity = ...
      0.377*8.8*(30e-1*100e-1*5e-3) + ...  % bronce coil
      0.45*11*0.08*1.5 + ...  % capacitors
      0.90*3.9*((33.8e-1/2)^2 - (32.2e-1/2)^2) * pi * 178e-1;  % Al2O3 ceramic tube


  case '60mm'
    % 60 mm long bronce rf coil (ProbeID 4588997)

    % damp coil settings
    HW.TX(iDevice).DampCoil.Enable = true;  % enable coil damping
    HW.TX(iDevice).DampCoil.DigitalOutputLatency = 0.6e-6;  % latency of damping circuit in seconds
    HW.TX(iDevice).DampCoil.DigitalOutputChannel = 1;  % digital output channel for coil damping signal
    HW.TX(iDevice).DampCoil.DigitalOutputDuration = HW.TX(iDevice).BlankPostset + 1.6e-6;  % duration of signal on digital out in seconds
    HW.TX(iDevice).DampCoil.DampingDuration = HW.TX(iDevice).DampCoil.DigitalOutputDuration + 2e-6;  % effective duration of damping the coil in seconds
    HW.TX(iDevice).DampCoil.TX2RXdeadTime = HW.TX(iDevice).DampCoil.DampingDuration + 3e-6;  % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)
    HW.TX(iDevice).BlankPostsetAQ = 1.4e-6;  % blank internal LNA
    HW.TX(iDevice).BlankAQ = 1;

    % gradient system
    HW.Grad(iDevice).ImageVol = [-0.02, 0.02, -0.032, 0.032, -0.02, 0.02];  % [xmin xmax ymin ymax zmin zmax]
    HW.Grad(iDevice).ImageVolOffset = [0, 0, 0];  % offset of coil to gradient system
    HW.RX(iDevice).EffectiveCoilLength = 60e-3;

    HW.Grad(iDevice).LoadRin(1:3) = [2.8, 3.5, 2.4];  % 60 mm x y z
    HW.Grad(iDevice).LoadIin2Amp(1:3) = [+0.021680, +0.021557, +0.023995];  % 2025-11-28T15:26:45, GainChangeXYZ = [ +0.034% -0.072% -0.065% ], centerXYZ = [ -0.096 +0.280 +0.316 ] mm, DiameterXYZ = [ +10.003  +9.993  +9.994 ] mm, SphereDiameter = +10.000 mm
    HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [80.018, 94.112, 61.262]*1e-6;  % time delay of gradient amplifier in s - calibrated with Calibrate_GradientSystem_GradDelay on 2025-12-08 15:58:31

    HW.MagnetShim([1,2,3]) = [0.000094201, -0.000284065, -0.000437377];  % 2025-12-08T16:16:34 by FindShim (T2* = 19.1 ms @ 18.364000 MHz), x y z in T/m

    % rf coil efficiency
    % HW.TX(iDevice).PaUout2Amplitude = [7.050108, 7.050108]*1e-6;  % 2025-11-17T11:05:10 (p90 = 8.000 탎 @ 104.106 V @ 18.364000 MHz) from FID amplitude by DailyCheck_solidEcho
    % kurze Probe:
    % HW.TX(iDevice).PaUout2Amplitude = [6.542983, 6.542983]*1e-6;  % 2025-12-01T13:50:25 (p180 = 16.000 탎 @ 112.190 V @ 18.360000 MHz) from CPMG Echo train by DailyCheck_solidEcho
    % lange Probe:
    HW.TX(iDevice).PaUout2Amplitude = [6.305094, 6.305094]*1e-6;  % 2025-12-08T16:16:52 (p180 = 16.000 탎 @ 116.407 V @ 18.364000 MHz) from CPMG Echo train by DailyCheck_solidEcho
    HW.TX(iDevice).PaUout2AmplitudeEstimated = HW.TX(iDevice).PaUout2Amplitude;

    % CPMG and solid echo settings
    HW.RecoveryCPMG.tEchoMin = 250e-6;  % minimum echo time that is desired for this coil
    HW.RecoveryCPMG.tFlip180Def = 16e-6;  % desired pulse length for 180 degrees pulse
    HW.TX(iDevice).Def.PaUout(2) = 110;  % 110 Volts correspond to a 180 degrees pulse duration of approximately 16e-6 s

    % rf coil heating model 30 mm bronce rf coil
    HW.TX(iDevice).CoilTemperature = 32;  % coil temperature in degrees at start of measurement
    HW.TX(iDevice).CoilMaxTemperature = 60;  % maximum coil temperature in degrees during measurement
    HW.TX(iDevice).CoilPowerDissipation = 2;  % estimated power dissipation to environment in Watt at 60 degrees
    % thermal capacity of coil in J/K
    % specific heat capacity of bronce at 25 degrees C: 0.377 J/(g*K)
    % density of bronce at room temperature: 8.8 g/cm^3
    % volume of bronce in coil: 60 mm * 100 mm * 50 um
    % specific heat capacity of capacitors: 0.45 J/(g*K) (https://www.murata.com/en-eu/support/faqs/capacitor/ceramiccapacitor/conf/0011)
    % weight of capacitors: 11 * 0.08 g
    % factor for solder: 1.5
    % specific heat capacity of Al2O3 tube: 0.90 J/(g*K)
    % density of Al2O3: 3.9 g/cm^3
    % volume of tube: ((33.8 mm/2)^2 - (32.2 mm/2)^2) * pi * 178 mm
    HW.TX(iDevice).CoilThermalCapacity = ...
      0.377*8.8*(60e-1*100e-1*5e-3) + ...  % bronce coil
      0.45*11*0.08*1.5 + ...  % capacitors
      0.90*3.9*((33.8e-1/2)^2 - (32.2e-1/2)^2) * pi * 178e-1;  % Al2O3 ceramic tube


  otherwise
    warning('LoadCoil:UnknownCoilName', ...
      'No configuration for coil "%s" found. Using default settings.', ...
      HW.TX(iDevice).CoilName);
    return;

end
