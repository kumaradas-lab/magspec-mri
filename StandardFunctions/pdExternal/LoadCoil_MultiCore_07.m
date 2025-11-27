%% Load settings for a named coil

% HW.GUI.showCoilName = true;  % show name of or selection box for coils
if ~exist('iDevice', 'var'), iDevice = 1; end

if isempty(HW.TX.CoilName), return; end


% settings specific for each coil
switch HW.TX.CoilName
  case 'MS187_probe_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [37.889029, 36.961662]*1e-6;  % 2022-11-23T16:15:41 (tFlip90 = 41.884 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'MS187_probe_H1C13_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [21.893960, 21.893960]*1e-6;  % 2022-11-28T14:59:18 (tFlip90 = 72.483 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'MS187_probe_H1C13_C13'
    HW.GammaDef = HW.Gamma.C13;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [13.224433, 13.224433] * 1e-6 * HW.Gamma.H1/HW.Gamma.C13; % C13-Pulslänge mit Ölprobe auf ca. 120 µs bestimmt - 2019-08-30T15:08:32 (tFlip90 = 120.000 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'MS187_probe_H1C13_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.C13;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [21.538002, 21.538002]*1e-6;  % 2022-11-30T11:14:13 (tFlip90 = 73.681 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [13.224433, 13.224433] * 1e-6 * HW.Gamma.H1/HW.Gamma.C13;  % C13-Pulslänge mit Ölprobe auf ca. 120 µs bestimmt - 2019-08-30T15:08:32 (tFlip90 = 120.000 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeEstimatedX = HW.TX(iDevice).PaUout2AmplitudeX;
  case 'MS187_probe_damp_Li7'
    HW.GammaDef = HW.Gamma.Li7;
    HW.FindFrequencyGamma = HW.Gamma.Li7;
    HW.TX(iDevice).PaUout2Amplitude = [40.267104, 40.267104]*1e-6;  % 2023-03-02T15:49:15 (tFlip90 = 101.406 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration

  otherwise
    warning('LoadCoil:UnknownCoilName', 'No configuration for coil "%s" found. Using default settings.', HW.TX.CoilName);
    return;
end

HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;
