%% Load settings for a named coil

% HW.GUI.showCoilName = true;  % show name of or selection box for coils

if ~exist('iDevice', 'var'), iDevice = 1; end

if isempty(HW.TX(iDevice).CoilName), return; end

% if isempty(HW.UserName) || ...
%     (all(~strcmp(HW.UserName, HW.UserNameList)) && ~strcmp(HW.TX(iDevice).CoilName, 'probe_H1')) || ...
%     (any(strcmp(HW.UserName, HW.UserNameList)) && ~strncmp(HW.UserName, HW.TX(iDevice).CoilName, numel(HW.UserName)))
%   error('PD:LoadCoil:NameMismatch', ...
%     ['HW.UserName and HW.TX.CoilName cannot be used in this combination.\n', ...
%     'Is HW.TX.CoilName in your LoadMySystem set correctly?']);
% end
% Users might have a numeric suffix for the resistor switch box. Strip those
% suffices before comparing to coil name.
% user_no_suffix = regexp(HW.UserName, '(_dualX?)?(_(?:[0-9]*|DC600))?$', 'split');
user_no_suffix = regexp(HW.UserName, '(_(?:[0-9]*|DC600))?$', 'split');
if ~isdeployed() && ...
    (isempty(HW.UserName) || ...
     (any(strcmp(HW.UserName, HW.UserNameList)) && ...
      ~strncmp(user_no_suffix{1}, HW.TX(iDevice).CoilName, numel(user_no_suffix{1}))))
  error('PD:LoadCoil:NameMismatch', ...
    ['HW.UserName and HW.TX.CoilName cannot be used in this combination.\n', ...
    'Is HW.TX.CoilName in your LoadMySystem set correctly?']);
end

% settings specific for each coil
switch HW.TX(iDevice).CoilName
  case 'probe_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [41.831845, 41.831845]*1e-6;  % 2022-04-08T13:39:25 (tFlip90 = 37.936 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1_F19'
    HW.GammaDef = HW.Gamma.F19;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    % mit Demo_Auto_ParameterSearch_dual_nuclear und Tri-Fluor-Essigsäure in Wasser bestimmt:
    HW.TX(iDevice).PaUout2Amplitude = [3.979367, 33.000000]*1e-6;  % 2022-11-16T12:43:52 (tFlip90 = 423.894 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1_F19_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.F19;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    % mit Demo_Auto_ParameterSearch_dual_nuclear und Tri-Fluor-Essigsäure in Wasser bestimmt:
    HW.TX(iDevice).PaUout2Amplitude = [39.711046, 33.000000]*1e-6;  % 2022-11-16T12:28:25 (tFlip90 = 39.962 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [3.979367, 33.000000]*1e-6;  % 2022-11-16T12:43:52 (tFlip90 = 423.894 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1_F19_dualX'
    HW.GammaDef = HW.Gamma.F19;
    HW.GammaX = HW.Gamma.H1;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    % mit Demo_Auto_ParameterSearch_dual_nuclear und Tri-Fluor-Essigsäure in Wasser bestimmt:
    HW.TX(iDevice).PaUout2AmplitudeX = [39.711046, 33.000000]*1e-6;  % 2022-11-16T12:28:25 (tFlip90 = 39.962 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2Amplitude = [3.979367, 33.000000]*1e-6;  % 2022-11-16T12:43:52 (tFlip90 = 423.894 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Li7_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [32.390084, 32.390084]*1e-6;  % 2022-04-08T13:46:49 (tFlip90 = 48.994 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Li7_Li7'
    HW.GammaDef = HW.Gamma.Li7;
    HW.TX(iDevice).PaUout2Amplitude = [25.477775, 25.477775]*1e-6;  % 2022-04-07T09:32:15 (tFlip90 = 160.271 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Li7_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.Li7;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    % mit Demo_Auto_ParameterSearch_dual_nuclear und LiCl-Lösung bestimmt:
    HW.TX(iDevice).PaUout2Amplitude = [32.390084, 32.390084]*1e-6;  % 2022-04-08T13:46:49 (tFlip90 = 48.994 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [25.991053, 25.991053]*1e-6;  % 2022-04-07T09:57:33 (tFlip90 = 157.105 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Al27_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [23.284478, 23.284478]*1e-6;  % 2023-01-24T13:05:50 (tFlip90 = 68.154 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Al27_Al27'
    HW.GammaDef = HW.Gamma.Al27;
    % FIXME:
    warning('coil efficiency not yet calibrated');
    HW.TX(iDevice).PaUout2Amplitude = [25.477775, 25.477775]*1e-6;  % 2022-04-07T09:32:15 (tFlip90 = 160.271 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Al27_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.Al27;
    HW.FindFrequencyGamma = HW.Gamma.H1;
    % FIXME:
    warning('coil efficiency not yet calibrated');
    % mit Demo_Auto_ParameterSearch_dual_nuclear und LiCl-Lösung bestimmt:
    HW.TX(iDevice).PaUout2Amplitude = [32.390084, 32.390084]*1e-6;  % 2022-04-08T13:46:49 (tFlip90 = 48.994 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [25.991053, 25.991053]*1e-6;  % 2022-04-07T09:57:33 (tFlip90 = 157.105 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Al27_Na23'
    HW.GammaDef = HW.Gamma.Na23;
    % mit gesättigter Kochsalzlösung bestimmt
    HW.TX(iDevice).PaUout2Amplitude = [28.727829, 28.727829]*1e-6;  % 2023-01-24T13:58:37 (tFlip90 = 208.834 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1Al27_Na23_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.Na23;
    HW.TX(iDevice).PaUout2Amplitude = [23.284478, 23.284478]*1e-6;  % 2023-01-24T13:05:50 (tFlip90 = 68.154 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
    HW.TX(iDevice).PaUout2AmplitudeX = [28.727829, 28.727829]*1e-6;  % 2023-01-24T13:58:37 (tFlip90 = 208.834 us @ 3.700 V) from 1d Spin Echo by Find_PulseDuration
  otherwise
    warning('LoadCoil:UnknownCoilName', 'No configuration for coil "%s" found. Using default settings.', HW.TX.CoilName);
    return;
end

% HW.TX(iDevice).PaUout2AmplitudeEstimated = HW.TX(iDevice).PaUout2Amplitude;
