%% Load settings for a named coil

% HW.GUI.showCoilName = true;  % show name of or selection box for coils

if ~exist('iDevice', 'var'), iDevice = 1; end

if isempty(HW.TX(iDevice).CoilName), return; end

if isempty(HW.UserName) || ...
    (all(~strcmp(HW.UserName, HW.UserNameList)) && ~strcmp(HW.TX(iDevice).CoilName, 'probe_H1')) || ...
    (any(strcmp(HW.UserName, HW.UserNameList)) && ~strncmp(HW.UserName, HW.TX(iDevice).CoilName, numel(HW.UserName)))
  error('PD:LoadCoil:NameMismatch', ...
    ['HW.UserName and HW.TX.CoilName cannot be used in this combination.\n', ...
    'Is HW.TX.CoilName in your LoadMySystem set correctly?']);
end

HW.GammaX = Inf;

% settings specific for each coil
switch HW.TX(iDevice).CoilName
  case 'probe_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [44.198804, 44.255226]*1e-6;  % 2023-09-12T14:46:06 (tFlip90 = 41.424 us @ 3.131 V) (tFlip90 = 2.654 us @ 50.000 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1C13_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [11.704897, 11.872117]*1e-6;  % 2023-09-14T12:35:15 (tFlip90 = 160.938 us @ 3.131 V)(tFlip90 = 9.891 us @ 50.000 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1C13_C13'
    HW.GammaDef = HW.Gamma.C13;
    HW.TX(iDevice).PaUout2Amplitude = [50.505002, 50.505002]*1e-6;  % 13C - 2023-09-18T15:58:15 (tFlip90 = 147.650 us @ 3.131 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1C13_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.C13;
    HW.TX(iDevice).PaUout2Amplitude = [11.999730, 12.209637]*1e-6;  % 2023-09-14T12:35:15 (tFlip90 = 156.285 us @ 3.131 V)(tFlip90 = 9.618 us @ 50.000 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [50.505002, 50.505002]*1e-6;  % 13C - 2023-09-18T15:58:15 (tFlip90 = 147.650 us @ 3.131 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeEstimatedX = HW.TX(iDevice).PaUout2AmplitudeX;
  case 'probe_H1N15_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [12.600941, 12.871815]*1e-6;  % 2023-09-18T08:35:19 (tFlip90 = 148.828 us @ 3.131 V)(tFlip90 = 9.123 us @ 50.000 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1N15_N15'
    HW.GammaDef = HW.Gamma.N15;
    HW.TX(iDevice).PaUout2Amplitude = [82.702745, 82.702745]*1e-6;  % 15N - 2023-09-20T03:02:05 (tFlip90 = 223.785 us @ 3.131 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1N15_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.N15;
    HW.TX(iDevice).PaUout2Amplitude = [12.600941, 12.871815]*1e-6;  % 2023-09-18T08:35:19 (tFlip90 = 148.828 us @ 3.131 V)(tFlip90 = 9.123 us @ 50.000 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [82.702745, 82.702745]*1e-6;  % 15N - 2023-09-20T03:02:05 (tFlip90 = 223.785 us @ 3.131 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeEstimatedX = HW.TX(iDevice).PaUout2AmplitudeX;
  otherwise
    warning('LoadCoil:UnknownCoilName', 'No configuration for coil "%s" found. Using default settings.', HW.TX(iDevice).CoilName);
    return;
end

HW.TX(iDevice).PaUout2AmplitudeEstimated = HW.TX(iDevice).PaUout2Amplitude;
