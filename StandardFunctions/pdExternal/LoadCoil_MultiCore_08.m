%% Load settings for a named coil

% for usage with MMRT #207 and magnet 206

% HW.GUI.showCoilName = true;  % show name of or selection box for coils
if ~exist('iDevice', 'var'), iDevice = 1; end

if isempty(HW.TX(iDevice).CoilName), return; end


% settings specific for each coil and/or resonance
switch HW.TX(iDevice).CoilName
  case 'probe_1H13C_1H'
    HW.GammaDef = HW.Gamma.H1;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [9.866643, 9.866643]*1e-6;  % 2023-02-01T10:00:09 (tFlip90 = 164.508 us @ 3.183 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_1H13C_13C'
    HW.GammaDef = HW.Gamma.C13;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [82.064349, 82.064349]*1e-6;  % 2023-01-27T13:45:13 (tFlip90 = 76.892 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  otherwise
    warning('LoadCoil:UnknownCoilName', ...
      'No configuration for coil "%s" found. Using default settings.', ...
      HW.TX(iDevice).CoilName);
    return;
end

HW.TX(iDevice).PaUout2AmplitudeEstimated = HW.TX(iDevice).PaUout2Amplitude;
