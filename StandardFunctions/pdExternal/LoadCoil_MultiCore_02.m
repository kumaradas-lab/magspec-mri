%% Load settings for a named coil

% HW.GUI.showCoilName = true;  % show name of or selection box for coils

if isempty(HW.TX.CoilName), return; end

if isempty(HW.UserName) || ...
    (all(~strcmp(HW.UserName, HW.UserNameList)) && ~strcmp(HW.TX.CoilName, 'probe_H1')) || ...
    (any(strcmp(HW.UserName, HW.UserNameList)) && ~strncmp(HW.UserName, HW.TX.CoilName, numel(HW.UserName)))
  error('PD:LoadCoil:NameMismatch', ...
    ['HW.UserName and HW.TX.CoilName cannot be used in this combination.\n', ...
    'Is HW.TX.CoilName in your LoadMySystem set correctly?']);
end

% settings specific for each coil
switch HW.TX.CoilName
  case 'probe_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX.PaUout2Amplitude = [42.446617, 41.081947]*1e-6; % 2019-12-04T08:12:07 (tFlip90 = 37.387 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
%   case 'probe_H1Na23_H1'
%     HW.GammaDef = HW.Gamma.H1;
%     HW.TX.PaUout2Amplitude = [14.686260, 36.955000] * 1e-6; % H1-Pulslänge mit H1Na23 Probenkopf 108.056 µs - 2019-08-29T13:38:56 (tFlip90 = 108.056 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
%   case 'probe_H1Na23_Na23'
%     HW.GammaDef = HW.Gamma.Na23;
%     HW.TX.PaUout2Amplitude = [0.000013499, 0.000036955] * HW.Gamma.H1/HW.Gamma.Na23; % Na23-Pulslänge mit H1Na23 Probenkopf 117.56 µs - 2019-08-29T13:43:08 by D_emo_Auto_PulseDuration_CPMG_fixedVoltage, factor voltage amplitude to B1 field strength on the coil output (TX) from CPMG Echo train
%   case 'probe_H1Xe129_H1'
%     HW.GammaDef = HW.Gamma.H1;
%     HW.TX.PaUout2Amplitude = [21.233065, 33.000000]*1e-6; % H1-Pulslänge mit H1Xe129 Probenkopf 74.739 µs - 2019-10-23T10:44:02 (tFlip90 = 74.739 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
%   case 'probe_H1Xe129_Xe129'
%     HW.GammaDef = HW.Gamma.Xe129;
%     HW.TX.PaUout2Amplitude = [13.224433, 33.000000]*1e-6 * HW.Gamma.H1/HW.Gamma.Xe129; % Xe129-Pulslänge mit H1Xe129 Probenkopf 120 µs - 2019-10-23T10:46:22 (tFlip90 = 120.000 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  otherwise
    warning('LoadCoil:UnknownCoilName', 'No configuration for coil "%s" found. Using default settings.', HW.TX.CoilName);
    return;
end

HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;
