%% Load settings for a named coil

% HW.GUI.showCoilName = true;  % show name of or selection box for coils

if isempty(HW.TX.CoilName), return; end

if isempty(HW.UserName) || ...
    (any(strcmp(HW.UserName, HW.UserNameList)) && ~strncmp(HW.UserName, HW.TX.CoilName, numel(HW.UserName)))
  error('PD:LoadCoil:NameMismatch', 'HW.UserName and HW.TX.CoilName cannot be used in this combination.');
end

% settings specific for each coil
switch HW.TX.CoilName
  case 'probe_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    HW.TX.PaUout2Amplitude = [43.544670, 43.544670]*1e-6; % 2020-11-03T17:12:17 (tFlip90 = 36.444 탎 @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1_F19'
    HW.GammaDef = HW.Gamma.F19;
    HW.FindFrequencyGamma = HW.Gamma.F19;
    HW.TX.PaUout2Amplitude = [4.627016, 4.627016]*1e-6; % 2020-11-04T10:04:04 (tFlip90 = 364.531 탎 @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_F19'
    HW.GammaDef = HW.Gamma.F19;
    HW.FindFrequencyGamma = HW.Gamma.F19;
    HW.TX.PaUout2Amplitude = [39.964995, 39.964995]*1e-6; % 2020-11-04T11:23:03 (tFlip90 = 42.204 탎 @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_F19_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    HW.TX.PaUout2Amplitude = [5.142463, 5.142463]*1e-6; % 2020-11-04T10:43:44 (tFlip90 = 308.594 탎 @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  otherwise
    warning('LoadCoil:UnknownCoilName', 'No configuration for coil "%s" found. Using default settings.', HW.TX.CoilName);
    return;
end

HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;