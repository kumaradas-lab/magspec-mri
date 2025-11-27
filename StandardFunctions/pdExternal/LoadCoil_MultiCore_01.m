%% Load settings for a named coil

% HW.GUI.showCoilName = true;  % show name of or selection box for coils

if ~exist('iDevice', 'var'), iDevice = 1; end

if isempty(HW.TX(iDevice).CoilName), return; end

if strcmp(HW.UserName, 'teach')
  % copied from LoadSystem_Specific
  switch HW.TX(iDevice).CoilName
    case {'probe_H1C13_H1', 'probe_H1C13_C13', 'probe_H1C13_dual'}
      % Magnet 130
      HW.B0 = (24.303e6)/HW.Gamma.H1*2*pi;

      % H1C13 Probenkopf
      HW.Grad(iDevice).LoadRin = [6.32, 6.06, 5.45, 5.50];  % x y z B0 gemessen inkl. RJ45-Kabel
      HW.Grad(iDevice).LoadIin2Amp = [0.21254, 0.19832, 0.32516, 0.015207]; % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [7.7808e-05   7.9568e-05   7.5088e-05]; % Time delay of grad amp
      HW.MagnetShim = [-0.001027354, -0.000294694, -0.000563739, 0.000000000]; % 2019-09-19T16:31:24 by FindShim, x y z in T/m and B0 in T

    case {'probe_H1N15_H1', 'probe_H1N15_N15', 'probe_H1N15_dual'}
      % Magnet 130
      HW.B0 = (24.303e6)/HW.Gamma.H1*2*pi;

      % H1N15 Probenkopf
      HW.Grad(iDevice).LoadRin = [6.14, 6.26, 5.55, 5.50];  % x y z B0 gemessen inkl. RJ45-Kabel
      HW.Grad(iDevice).LoadIin2Amp = [0.19325, 0.20272, 0.32883, 0.015207]; % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [6.9592e-05, 6.5984e-05, 7.6272e-05]; % Time delay of grad amp

      HW.MagnetShim = [-0.001027354, -0.000294694, -0.000563739, 0.000000000]; % !!! WERT von probe_H1C13 muss von Kunden bestimmt werden

    case {'probe_H1Xe129_H1', 'probe_H1Xe129_Xe129', 'probe_H1Xe129_dual'}
      % Magnet 130
      HW.B0 = (24.303e6)/HW.Gamma.H1*2*pi;

      % H1Xe129 Probenkopf
      HW.Grad(iDevice).LoadRin  =[6.20, 6.44, 5.97, 5.50]; % x y z B0 gemessen inkl. RJ45-Kabel (23.10.2019)
      HW.Grad(iDevice).LoadIin2Amp = [0.20074, 0.21637, 0.33889, 0.015207]; % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [7.1976e-05, 6.9496e-05, 7.9424e-05]; % Time delay of grad amp

      % Shimwerte von probe_H1 kopiert
      HW.MagnetShim = [-0.000922358, -0.000420482, -0.000421957, 0.000000000]; % 2019-09-19T16:25:50 by FindShim, x y z in T/m and B0 in T

    case {'m224_probe_H1C13_H1',  'm224_probe_H1C13_C13', 'm224_probe_H1C13_dual'}  % Magnet 224
      HW.B0 = (18.964e6)/HW.Gamma.H1*2*pi;

      % H1C13 Probenkopf 20 mm Magnet
      HW.Grad(iDevice).LoadRin = [7.37-1.041163, 7.19-0.557871, 6.08-0.180048, 5.5]; % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.0689807, 0.0679459, 0.0715349, 0.0051]; % x y z B0  (T/(m*A)) Tesla per Meter per Ampere

      HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [127.612, 149.320, 102.380]*1e-6;  % time delay of gradient amplifier in s
      %in DC-600: HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [85.498, 100.262, 64.332]*1e-6;  % time delay of gradient amplifier in s

      HW.MagnetShim([1,2,3]) = [-0.000648419, 0.000530268, -0.000778691];  % 2024-09-27T14:55:00 by FindShim (T2* = 8.4 ms @ 18.963358 MHz), x y z in T/m

    case {'m224_probe_H1N15_H1', 'm224_probe_H1N15_N15', 'm224_probe_H1N15_dual'}  % Magnet 224
      HW.B0 = (18.964e6)/HW.Gamma.H1*2*pi;

      % H1N15 Probenkopf 20 mm Magnet
      HW.Grad(iDevice).LoadRin = [6.61-0.378840, 7.16-0.411752, 5.71-0.109211, 5.5]; % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.0687858, 0.0677294, 0.0715089, 0.0051]; % x y z B0  (T/(m*A)) Tesla per Meter per Ampere

      HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [131.962, 153.270, 108.378]*1e-6;  % time delay of gradient amplifier in s
      % in DC-600: HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [88.580, 106.096, 65.346]*1e-6;  % time delay of gradient amplifier in s

      HW.MagnetShim([1,2,3]) = [-0.000640212, 0.000499943, -0.000772939];  % 2024-09-27T14:11:11 by FindShim (T2* = 8.4 ms @ 18.963109 MHz), x y z in T/m

  end

  if isa(HW, 'PD.HWClass')
    % reload calibration settings
    % This won't work (and is not needed) when this is called while creating the
    % LoadMySystem files for all users. So, skip this block.
    HW.LoadMySystemConfig();
    HW.LoadMagnetShimCal();
  end

else
  if isempty(HW.UserName) || ...
      (all(~strcmp(HW.UserName, HW.UserNameList)) && ~strcmp(HW.TX(iDevice).CoilName, 'probe_H1')) || ...
      (any(strcmp(HW.UserName, HW.UserNameList)) && ~strncmp(HW.UserName, HW.TX(iDevice).CoilName, numel(HW.UserName)))
    error('PD:LoadCoil:NameMismatch', ...
      ['HW.UserName and HW.TX.CoilName cannot be used in this combination.\n', ...
      'Is HW.TX.CoilName in your LoadMySystem set correctly?']);
  end

end


% settings specific for each coil
switch HW.TX(iDevice).CoilName
  case 'probe_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [45.917517, 45.917517] * 1e-6; % 2019-08-30T10:30:52 (tFlip90 = 34.560 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Na23_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [14.686260, 14.686260] * 1e-6; % H1-Pulslänge mit H1Na23 Probenkopf 108.056 µs - 2019-08-29T13:38:56 (tFlip90 = 108.056 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Na23_Na23'
    HW.GammaDef = HW.Gamma.Na23;
    HW.TX(iDevice).PaUout2Amplitude = [0.000013499, 0.000013499] * HW.Gamma.H1/HW.Gamma.Na23; % Na23-Pulslänge mit H1Na23 Probenkopf 117.56 µs - 2019-08-29T13:43:08 by D_emo_Auto_PulseDuration_CPMG_fixedVoltage, factor voltage amplitude to B1 field strength on the coil output (TX) from CPMG Echo train
  case 'probe_H1Xe129_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [21.233065, 21.233065]*1e-6; % H1-Pulslänge mit H1Xe129 Probenkopf 74.739 µs - 2019-10-23T10:44:02 (tFlip90 = 74.739 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Xe129_Xe129'
    HW.GammaDef = HW.Gamma.Xe129;
    HW.TX(iDevice).PaUout2Amplitude = [13.224433, 13.224433]*1e-6 * HW.Gamma.H1/HW.Gamma.Xe129; % Xe129-Pulslänge mit H1Xe129 Probenkopf 120 µs - 2019-10-23T10:46:22 (tFlip90 = 120.000 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Xe129_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.Xe129;
    HW.TX(iDevice).PaUout2Amplitude = [21.233065, 21.233065]*1e-6; % H1-Pulslänge mit H1Xe129 Probenkopf 74.739 µs - 2019-10-23T10:44:02 (tFlip90 = 74.739 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [13.224433, 13.224433]*1e-6 * HW.Gamma.H1/HW.Gamma.Xe129; % Xe129-Pulslänge mit H1Xe129 Probenkopf 120 µs - 2019-10-23T10:46:22 (tFlip90 = 120.000 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeEstimatedX = HW.TX(iDevice).PaUout2AmplitudeX;
  case 'probe_H1N15_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [15.924614, 15.924614]*1e-6; % H1-Pulslänge mit H1N15 Probenkopf 99.653 µs -  2019-10-10T12:48:26 (tFlip90 = 99.653 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1N15_N15'
    HW.GammaDef = HW.Gamma.N15;
    HW.TX(iDevice).PaUout2Amplitude = [12.695456, 12.695456]*1e-6 * HW.Gamma.H1/HW.Gamma.N15; % N15-Pulslänge mit H1N15 Probenkopf 125 µs -  2019-10-11T09:04:03 (tFlip90 = 125.000 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1N15_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.N15;
    HW.TX(iDevice).PaUout2Amplitude = [15.924614, 15.924614]*1e-6;  % H1-Pulslänge mit H1N15 Probenkopf 99.653 µs -  2019-10-10T12:48:26 (tFlip90 = 99.653 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [12.695456, 12.695456]*1e-6 * HW.Gamma.H1/HW.Gamma.N15;  % N15-Pulslänge mit H1N15 Probenkopf 125 µs -  2019-10-11T09:04:03 (tFlip90 = 125.000 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeEstimatedX = HW.TX(iDevice).PaUout2AmplitudeX;
  case 'probe_H1C13_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [16.336914, 16.336914] * 1e-6; % H1-Pulslänge mit H1C13 Probenkopf 97.138 µs - 2019-08-30T08:36:45 (tFlip90 = 97.138 µs @ 3.700
  case 'probe_H1C13_C13'
    HW.GammaDef = HW.Gamma.C13;
    HW.TX(iDevice).PaUout2Amplitude = [13.224433, 13.224433] * 1e-6 * HW.Gamma.H1/HW.Gamma.C13; % C13-Pulslänge mit Ölprobe auf ca. 120 µs bestimmt - 2019-08-30T15:08:32 (tFlip90 = 120.000 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1C13_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.C13;
    HW.TX(iDevice).PaUout2Amplitude = [16.336914, 16.336914] * 1e-6;  % H1-Pulslänge mit H1C13 Probenkopf 97.138 µs - 2019-08-30T08:36:45 (tFlip90 = 97.138 µs @ 3.700
    HW.TX(iDevice).PaUout2AmplitudeX = [13.224433, 13.224433] * 1e-6 * HW.Gamma.H1/HW.Gamma.C13;  % C13-Pulslänge mit Ölprobe auf ca. 120 µs bestimmt - 2019-08-30T15:08:32 (tFlip90 = 120.000 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeEstimatedX = HW.TX(iDevice).PaUout2AmplitudeX;

  case 'm224_probe_H1C13_H1'
    HW.GammaDef = HW.Gamma.H1;
    % HW.TX(iDevice).PaUout2Amplitude = [4.224240, 4.547855]*1e-6;  % 1H an TRx - 2024-09-24T13:26:13 (tFlip90 = 375.673 us @ 3.700 V @ 18.964175 MHz) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2Amplitude = [4.224240, 4.547855]*1e-6;  % 1H mit RF200 #5 - 2024-09-24T13:56:12 (tFlip90 = 43.036 us @ 30.000 V @ 18.964129 MHz) from 1d Spin Echo by Find_PulseDuration
  case 'm224_probe_H1C13_C13'
    HW.GammaDef = HW.Gamma.C13;
    % HW.TX(iDevice).PaUout2Amplitude = [16.033552, 15.270316]*1e-6;  % TRx mit Na23 in Magspec 218 (kalt) 2024-09-25T09:05:15 (tFlip90 = 374.175 us (393.520 us @ C13) @ 3.700 V @ 4.757340 MHz) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2Amplitude = [16.033552, 15.270316]*1e-6;  % RF200#5 mit Na23 in Magspec 218 (kalt) 2024-09-25T09:12:11 (tFlip90 = 48.455 us (50.960 us @ C13) @ 30.000 V @ 4.757231 MHz) from 1d Spin Echo by Find_PulseDuration
  case 'm224_probe_H1C13_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.C13;
    % HW.TX(iDevice).PaUout2Amplitude = [4.224240, 4.547855]*1e-6;  % 1H an TRx - 2024-09-24T13:26:13 (tFlip90 = 375.673 us @ 3.700 V @ 18.964175 MHz) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2Amplitude = [4.224240, 4.547855]*1e-6;  % 1H mit RF200 #5 - 2024-09-24T13:56:12 (tFlip90 = 43.036 us @ 30.000 V @ 18.964129 MHz) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [16.033552, 15.270316]*1e-6;  % mit Na23 in Magspec 218 (kalt) 2024-09-25T09:12:11 (tFlip90 = 48.455 us (50.960 us @ C13) @ 30.000 V @ 4.757231 MHz) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeEstimatedX = HW.TX(iDevice).PaUout2AmplitudeX;

  case 'm224_probe_H1N15_H1'
    HW.GammaDef = HW.Gamma.H1;
    % HW.TX(iDevice).PaUout2Amplitude = [3.200255, 3.257785]*1e-6;  % 1H an TRx - 2024-09-24T14:02:01 (tFlip90 = 495.877 us @ 3.700 V @ 18.964999 MHz) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2Amplitude = [3.200255, 3.257785]*1e-6;  % 1H mit RF200 #5 - 2024-09-24T14:00:07 (tFlip90 = 60.078 us @ 30.000 V @ 18.964671 MHz) from 1d Spin Echo by Find_PulseDuration
  case 'm224_probe_H1N15_N15'
    HW.GammaDef = HW.Gamma.N15;
    % HW.TX(iDevice).PaUout2Amplitude = [17.970001, 22.102659]*1e-6;  % TRx mit mit Na23 in Minispec (kalt) 2024-09-25T17:03:17 (tFlip90 = 333.854 us (884.562 us @ N15) @ 3.700 V @ 1.990046 MHz) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2Amplitude = [17.970001, 22.102659]*1e-6;  % RF200#5 mit mit Na23 in Minispec (kalt) 2024-09-25T17:08:28 (tFlip90 = 33.477 us (87.390 us @ N15) @ 30.000 V @ 1.990046 MHz) from 1d Spin Echo by Find_PulseDuration
  case 'm224_probe_H1N15_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.N15;
    % HW.TX(iDevice).PaUout2Amplitude = [3.200255, 3.257785]*1e-6;  % 1H an TRx - 2024-09-24T14:02:01 (tFlip90 = 495.877 us @ 3.700 V @ 18.964999 MHz) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2Amplitude = [3.200255, 3.257785]*1e-6;  % 1H mit RF200 #5 - 2024-09-24T14:00:07 (tFlip90 = 60.078 us @ 30.000 V @ 18.964671 MHz) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [17.970001, 22.102659]*1e-6;  % RF200#5 mit mit Na23 in Minispec (kalt) 2024-09-25T17:08:28 (tFlip90 = 33.477 us (87.390 us @ N15) @ 30.000 V @ 1.990046 MHz) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeEstimatedX = HW.TX(iDevice).PaUout2AmplitudeX;
  otherwise
    warning('LoadCoil:UnknownCoilName', 'No configuration for coil "%s" found. Using default settings.', HW.TX(iDevice).CoilName);
    return;
end

HW.TX(iDevice).PaUout2AmplitudeEstimated = HW.TX(iDevice).PaUout2Amplitude;
