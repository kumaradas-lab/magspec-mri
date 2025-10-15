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
    HW.FindFrequencyGamma = HW.Gamma.H1;
    HW.FindFrequencySweep.maxTime = 1000;
    HW.TX.PaUout2Amplitude = [40.038079, 40.038079]*1e-6; % 2020-07-09T10:39:09 (tFlip90 = 39.636 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_C13'
    HW.GammaDef = HW.Gamma.C13;
    HW.FindFrequencyGamma = HW.Gamma.C13;
    HW.FindFrequencySweep.maxTime = inf;
    
    % C13-Pulslänge mit Ölprobe auf ca. 45 µs bestimmt 
    HW.TX.PaUout2Amplitude = [35.265155, 35.265155]* 1e-6 * HW.Gamma.H1/HW.Gamma.C13; % 2020-07-09T10:49:24 (tFlip90 = 45.000 µs @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  otherwise
    warning('LoadCoil:UnknownCoilName', 'No configuration for coil "%s" found. Using default settings.', HW.TX.CoilName);
    return;
end

HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;
