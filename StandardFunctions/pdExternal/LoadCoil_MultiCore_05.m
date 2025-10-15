%% Load settings for a named coil (or probe if deployed)

% HW.GUI.showCoilName = true;  % show name of or selection box for coils

% FIXME: Support multiple connected devices?
iDevice = 1;

if isempty(HW.TX(iDevice).CoilName), return; end

% Users for dual nucleus configurations have the suffix "_dual". Additionally,
% users might have a numeric suffix for the resistor switch box. Strip those
% suffices before comparing to coil name.
user_no_suffix = regexp(HW.UserName, '(_dual)?(_(?:[0-9]*|DC600))?$', 'split');
if ~isdeployed() && ...
    (isempty(HW.UserName) || ...
     (any(strcmp(HW.UserName, HW.UserNameList)) && ...
      ~strncmp(user_no_suffix{1}, HW.TX(iDevice).CoilName, numel(user_no_suffix{1}))))
  error('PD:LoadCoil:NameMismatch', ...
    ['HW.UserName and HW.TX.CoilName cannot be used in this combination.\n', ...
    'Is HW.TX.CoilName in your LoadMySystem set correctly?']);
end

% settings specific for each coil
switch HW.TX.CoilName
  case 'probe_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [43.544670, 43.544670]*1e-6;  % 2020-11-03T17:12:17 (tFlip90 = 36.444 탎 @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    if isdeployed()
      % resistor switch box is set to 0 Ohms (without DC-600)
      HW.Grad(iDevice).LoadRin = [4.44+0.002625+0.010768, 4.89+0.092672-0.010542, 4.84+0.196801-0.008441, 200];  % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.253149, 0.252638, 0.311733, 0.015207];  % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [20.904, 24.728, 25.816]*1e-6;  % time delay of gradient amplifier in s - 2023-02-02 17:05
      HW.MagnetShim([1,2,3]) = [0.000381048, 0.000129276, 0.001992783];  % 2023-02-02T17:00:04 by FindShim (T2* = 7.6 ms), x y z in T/m and B0 in T
      HW.RX(iDevice).EffectiveCoilLength = 13.5e-3;
    end
  case 'probe_H1_F19'
    HW.GammaDef = HW.Gamma.F19;
    HW.FindFrequencyGamma = HW.Gamma.F19;
    HW.TX(iDevice).PaUout2Amplitude = [4.627016, 4.627016]*1e-6;  % 2020-11-04T10:04:04 (tFlip90 = 364.531 탎 @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    if isdeployed()
      % resistor switch box is set to 0 Ohms (without DC-600)
      HW.Grad(iDevice).LoadRin = [4.44+0.002625+0.010768, 4.89+0.092672-0.010542, 4.84+0.196801-0.008441, 200];  % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.253149, 0.252638, 0.311733, 0.015207];  % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [20.904, 24.728, 25.816]*1e-6;  % time delay of gradient amplifier in s - 2023-02-02 17:05
      HW.MagnetShim([1,2,3]) = [0.000381048, 0.000129276, 0.001992783];  % 2023-02-02T17:00:04 by FindShim (T2* = 7.6 ms), x y z in T/m and B0 in T
      HW.RX(iDevice).EffectiveCoilLength = 13.5e-3;
    end
  case 'probe_F19'
    HW.GammaDef = HW.Gamma.F19;
    HW.FindFrequencyGamma = HW.Gamma.F19;
    % without RF-100
    % HW.TX(iDevice).PaUout2Amplitude = [39.964995, 39.964995]*1e-6;  % 2020-11-04T11:23:03 (tFlip90 = 42.204 탎 @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    % with RF-100 #6
    HW.TX(iDevice).PaUout2Amplitude = [39.964995, 41.804654]*1e-6;  % 2023-02-17T12:12:04 (tFlip90 = 2.986 us @ 50.000 V) from 1d Spin Echo by Find_PulseDuration
    if isdeployed()
      % resistor switch box is set to 0 Ohms (without DC-600)
      HW.Grad(iDevice).LoadRin = [(4.84-0.075054)/(4.125+4.325)*9, (4.71+0.024190)/(4.21875+2.57818)*(4.45313+2.69531), (4.26+0.082375)/(3.65625+4.6875)*9, 200];  % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.254361*(4.5+4.59375)/9, 0.25406, 0.30964*(4.125+5.0625)/9, 0.015207];  % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [3.324e-05, 3.6768e-05, 4.16e-05];  % Time delay of grad amp
      HW.MagnetShim([1,2,3]) = [0.000226409, 0.000112733, 0.001841038];  % 2023-02-17T11:58:00 by FindShim (T2* = 11.8 ms), x y z in T/m and B0 in T
      HW.RX(iDevice).EffectiveCoilLength = 13.5e-3;
    end
  case 'probe_F19_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    % without RF-100
    % HW.TX(iDevice).PaUout2Amplitude = [5.142463, 5.142463]*1e-6;  % 2020-11-04T10:43:44 (tFlip90 = 308.594 탎 @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    % with RF-100 #6
    HW.TX(iDevice).PaUout2Amplitude = [4.961585, 4.484760]*1e-6;  % 2023-02-17T11:26:44 (tFlip90 = 26.185 us @ 50.000 V) from 1d Spin Echo by Find_PulseDuration
    if isdeployed()
      % resistor switch box is set to 0 Ohms (without DC-600)
      HW.Grad(iDevice).LoadRin = [(4.84-0.075054)/(4.125+4.325)*9, (4.71+0.024190)/(4.21875+2.57818)*(4.45313+2.69531), (4.26+0.082375)/(3.65625+4.6875)*9, 200];  % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.254361*(4.5+4.59375)/9, 0.25406, 0.30964*(4.125+5.0625)/9, 0.015207];  % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [3.324e-05, 3.6768e-05, 4.16e-05]; % Time delay of grad amp
      HW.MagnetShim([1,2,3]) = [0.000226409, 0.000112733, 0.001841038];  % 2023-02-17T11:58:00 by FindShim (T2* = 11.8 ms), x y z in T/m and B0 in T
      HW.RX(iDevice).EffectiveCoilLength = 13.5e-3;
    end
  case 'probe_H1C13_H1'
    % Q=13
    HW.GammaDef = HW.Gamma.H1;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    % without RF-100
    HW.TX(iDevice).PaUout2Amplitude = [9.462452, 9.462452]*1e-6;  % 1H - 2023-04-05T11:21:18 (tFlip90 = 167.708 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    % with RF-100 #6 (nicht gemessen!)
    if isdeployed()
      % resistor switch box is set to 0 Ohms (without DC-600)
      HW.Grad(iDevice).LoadRin = [4.82+0.1880, 5.09+0.3540, 4.81+0.3572, 200];  % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.247947, 0.235389, 0.304866, 2];  % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [41.614, 54.732, 49.752]*1e-6;  % time delay of gradient amplifier in s
      HW.MagnetShim([1,2,3]) = [0.001070608, -0.000196937, -0.000192357];  % 2023-04-24T13:22:14 by FindShim (T2* = 19.3 ms), x y z in T/m and B0 in T
      HW.RX(iDevice).EffectiveCoilLength = 13.5e-3;
    end
  case 'probe_H1C13_C13'
    % Q=33
    HW.GammaDef = HW.Gamma.C13;
    HW.FindFrequencyGamma = HW.Gamma.C13;
    % without RF-100
    HW.TX(iDevice).PaUout2Amplitude = [56.050104, 56.050104]*1e-6;  % 13C - 2023-04-05T23:48:59 (tFlip90 = 112.580 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    % with RF-100 #6 (nicht gemessen!)
    if isdeployed()
      % resistor switch box is set to 0 Ohms (without DC-600)
      HW.Grad(iDevice).LoadRin = [4.82+0.1880, 5.09+0.3540, 4.81+0.3572, 200];  % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.247947, 0.235389, 0.304866, 2];  % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [41.614, 54.732, 49.752]*1e-6;  % time delay of gradient amplifier in s
      HW.MagnetShim([1,2,3]) = [0.001070608, -0.000196937, -0.000192357];  % 2023-04-24T13:22:14 by FindShim (T2* = 19.3 ms), x y z in T/m and B0 in T
      HW.RX(iDevice).EffectiveCoilLength = 13.5e-3;
    end
  case 'probe_H1N15_H1'
    % Q=29
    HW.GammaDef = HW.Gamma.H1;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    % without RF-100
    HW.TX(iDevice).PaUout2Amplitude = [9.670681, 9.670681]*1e-6;  % 1H - 2023-04-17T09:31:43 (tFlip90 = 164.097 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    % with RF-100 #6 (nicht gemessen!)
    if isdeployed()
      % resistor switch box is set to 0 Ohms (without DC-600)
      HW.Grad(iDevice).LoadRin = [4.73+0.1799, 4.46+0.3788, 4.94+0.4285, 200];  % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.248666, 0.234669, 0.30511, 2]; % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [54.574, 69.676, 62.390]*1e-6;  % time delay of gradient amplifier in s
      HW.MagnetShim([1,2,3]) = [0.001226034, 0.000058690, 0.000233335];  % 2023-04-25T08:29:49 by FindShim (T2* = 18.5 ms), x y z in T/m and B0 in T
      HW.RX(iDevice).EffectiveCoilLength = 13.5e-3;
    end
  case 'probe_H1N15_N15'
    % Q=22
    HW.GammaDef = HW.Gamma.N15;
    HW.FindFrequencyGamma = HW.Gamma.N15;
    % without RF-100
    HW.TX(iDevice).PaUout2Amplitude = [120.469560, 120.469560]*1e-6;  % 15N mit Demo_FID manuell bestimmt - 2023-04-18T12:10:56 (tFlip90 = 130.000 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    % with RF-100 #6 (nicht gemessen!)
    if isdeployed()
      % resistor switch box is set to 0 Ohms (without DC-600)
      HW.Grad(iDevice).LoadRin = [4.73+0.1799, 4.46+0.3788, 4.94+0.4285, 200];  % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.248666, 0.234669, 0.30511, 2]; % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [54.574, 69.676, 62.390]*1e-6;  % time delay of gradient amplifier in s
      HW.MagnetShim([1,2,3]) = [0.001226034, 0.000058690, 0.000233335];  % 2023-04-25T08:29:49 by FindShim (T2* = 18.5 ms), x y z in T/m and B0 in T
      HW.RX(iDevice).EffectiveCoilLength = 13.5e-3;
    end
  case 'probe_H1P31_H1'
    % Q=30 (?)
    HW.GammaDef = HW.Gamma.H1;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    % without RF-100
    HW.TX(iDevice).PaUout2Amplitude = [10.823194, 10.823194]*1e-6;  % 2023-04-18T16:40:50 (tFlip90 = 146.623 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    % with RF-100 #6 (nicht gemessen!)
    if isdeployed()
      % resistor switch box is set to 0 Ohms (without DC-600)
      HW.Grad(iDevice).LoadRin = [5.94-1.0091, 4.51+0.3803, 5.06+0.3510, 200]; % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.248756, 0.235473, 0.306371, 2]; % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [52.722, 68.790, 58.654]*1e-6;  % time delay of gradient amplifier in s
      HW.MagnetShim([1,2,3]) = [0.001358681, 0.000052989, 0.000186768];  % 2023-04-25T09:53:57 by FindShim (T2* = 19.7 ms), x y z in T/m and B0 in T
      HW.RX(iDevice).EffectiveCoilLength = 13.5e-3;
    end
  case 'probe_H1P31_P31'
    % Q=43
    HW.GammaDef = HW.Gamma.P31;
    HW.FindFrequencyGamma = HW.Gamma.P31;
    % without RF-100
    HW.TX(iDevice).PaUout2Amplitude = [48.686072, 48.686072]*1e-6;  % 31P - 2023-04-18T16:47:37 (tFlip90 = 80.521 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    % with RF-100 #6 (nicht gemessen!)
    if isdeployed()
      % resistor switch box is set to 0 Ohms (without DC-600)
      HW.Grad(iDevice).LoadRin = [5.94-1.0091, 4.51+0.3803, 5.06+0.3510, 200]; % x y z B0 in (Ohm)
      HW.Grad(iDevice).LoadIin2Amp = [0.248756, 0.235473, 0.306371, 2]; % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(1:3) = [52.722, 68.790, 58.654]*1e-6;  % time delay of gradient amplifier in s
      HW.MagnetShim([1,2,3]) = [0.001358681, 0.000052989, 0.000186768];  % 2023-04-25T09:53:57 by FindShim (T2* = 19.7 ms), x y z in T/m and B0 in T
      HW.RX(iDevice).EffectiveCoilLength = 13.5e-3;
    end
  otherwise
    warning('LoadCoil:UnknownCoilName', ...
      'No configuration for coil "%s" found. Using default settings.', HW.TX(iDevice).CoilName);
    return;
end

HW.TX(iDevice).PaUout2AmplitudeEstimated = HW.TX(iDevice).PaUout2Amplitude;
