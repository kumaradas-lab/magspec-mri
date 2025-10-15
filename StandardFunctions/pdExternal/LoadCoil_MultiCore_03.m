%% Load settings for a named coil

% HW.GUI.showCoilName = true;  % show name of or selection box for coils

if isempty(HW.TX.CoilName), return; end

if isempty(HW.UserName) || ...
    (any(strcmp(HW.UserName, HW.UserNameList)) && ~strncmp(HW.UserName, HW.TX.CoilName, numel(HW.UserName)))
  error('PD:LoadCoil:NameMismatch', ...
    ['HW.UserName and HW.TX.CoilName cannot be used in this combination.\n', ...
    'Is HW.TX.CoilName in your LoadMySystem set correctly?']);
end

% settings specific for each coil
switch HW.TX.CoilName
  case 'probe_H1_10'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX.PaUout2Amplitude = [39.690526, 43.141959]*1e-6; % 2020-03-12 - Channel 1: 39.983 탎 @ 3.700 V - Channel 2: 46.424 탎 @ 2.932 V
  case 'probe_H1_5'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX.PaUout2Amplitude = [61.574457, 67.462923]*1e-6; % 2020-03-12 - Channel 1: 25.773 탎 @ 3.700 V - Channel 2: 29.688 탎 @ 2.932 V
  otherwise
    warning('LoadCoil:UnknownCoilName', 'No configuration for coil "%s" found. Using default settings.', HW.TX.CoilName);
    return;
end

HW.TX.PaUout2AmplitudeEstimated = HW.TX.PaUout2Amplitude;
