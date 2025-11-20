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
  case 'm175_probe_H1Na23_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [8.861588, 9.924331]*1e-6; % 2023-09-20T14:08:37 (tFlip90 = 202.691 us @ 2.919 V) from 1d Spin Echo by Find_PulseDuration

  case 'm175_probe_H1Na23_Na23'
    HW.GammaDef = HW.Gamma.Na23;
    HW.TX(iDevice).PaUout2Amplitude = [48.377869, 51.539045]*1e-6; % Na23 - 2023-09-20T14:18:28 (tFlip90 = 124.010 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration

  case 'm175_probe_H1Na23_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.Na23;
    HW.TX(iDevice).PaUout2Amplitude = [8.861588, 9.924331]*1e-6; % 2023-09-20T14:08:37 (tFlip90 = 202.691 us @ 2.919 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [48.377869, 51.539045]*1e-6; % Na23 - 2023-09-20T14:18:28 (tFlip90 = 124.010 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeEstimatedX = HW.TX(iDevice).PaUout2AmplitudeX;

  otherwise
    warning('LoadCoil:UnknownCoilName', 'No configuration for coil "%s" found. Using default settings.', HW.TX(iDevice).CoilName);
    return;
end

HW.TX(iDevice).PaUout2AmplitudeEstimated = HW.TX(iDevice).PaUout2Amplitude;
