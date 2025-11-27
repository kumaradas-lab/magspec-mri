%% Load settings for a named coil

% HW.GUI.showCoilName = true;  % show name of or selection box for coils

if ~exist('iDevice', 'var'), iDevice = 1; end

if isempty(HW.TX(iDevice).CoilName), return; end

if strcmp(HW.UserName, 'teach')
  % copied from LoadSystem_Specific
  switch HW.TX(iDevice).CoilName
    case {'probe_H1C13_H1', 'probe_H1C13_C13', 'probe_H1C13_dual'}
      % probe_H1C13
      HW.Grad(iDevice).LoadRin = [2.95, 3.26, 2.35, 5];  % Gradient - 11.11.2024
      HW.Grad(iDevice).LoadIin2Amp = [0.0261231, 0.0258773, 0.0266735, 0.0051];  % x y z B0  (T/(m*A)) Tesla per Meter per Ampere
      HW.Grad(iDevice).SystemTimeDelay(HW.Grad(iDevice).xyzB(1:3)) = [104.144, 104.010, 73.878]*1e-6;  % time delay of gradient amplifier in s

      HW.MagnetShim([1,2,3]) = [0.000295507, 0.000210236, 0.000943888];  % 2024-11-11T14:44:23 by FindShim (T2* = 16.2 ms @ 18.577939 MHz), x y z in T/m

  end

  if isa(HW, 'PD.HWClass')
    % reload calibration settings
    % This won't work (and is not needed) when this is called while creating the
    % LoadMySystem files for all users. So, skip this block.
    HW.LoadMySystemConfig();
    HW.LoadMagnetShimCal();
  end

else
  if isempty(HW.UserName) || ...
      (all(~strcmp(HW.UserName, HW.UserNameList)) && ~strcmp(HW.TX(iDevice).CoilName, 'probe_H1C13_H1')) || ...
      (any(strcmp(HW.UserName, HW.UserNameList)) && ~strncmp(HW.UserName, HW.TX(iDevice).CoilName, numel(HW.UserName)))
    error('PD:LoadCoil:NameMismatch', ...
      ['HW.UserName and HW.TX.CoilName cannot be used in this combination.\n', ...
      'Is HW.TX.CoilName in your LoadMySystem set correctly?']);
  end

end


% settings specific for each coil
switch HW.TX(iDevice).CoilName
  case 'probe_H1C13_H1'
    HW.GammaDef = HW.Gamma.H1;
    HW.TX(iDevice).PaUout2Amplitude = [3.359211, 3.283662]*1e-6;  % 2024-11-11T14:44:49 (tFlip90 = 35.763 us @ 50.000 V @ 18.577939 MHz) from 1d Spin Echo by Find_PulseDuration
  case 'probe_H1C13_C13'
    HW.GammaDef = HW.Gamma.C13;
    % acetic acid (C13)
    HW.TX(iDevice).PaUout2Amplitude = [3.185085, 14.942144]*1e-6;  % 2024-11-19T18:16:33 (tFlip90 = 31.250 us @ 50.000 V) from bulk FID by Find_PulseDurationFID

  case 'probe_H1C13_dual'
    HW.GammaDef = HW.Gamma.H1;
    HW.GammaX = HW.Gamma.C13;
    HW.TX(iDevice).PaUout2Amplitude = [3.359211, 3.283662]*1e-6;  % 2024-11-11T14:44:49 (tFlip90 = 35.763 us @ 50.000 V @ 18.577939 MHz) from 1d Spin Echo by Find_PulseDuration
    % acetic acid (C13)
    HW.TX(iDevice).PaUout2AmplitudeX = [3.185085, 14.942144]*1e-6;  % 2024-11-19T18:16:33 (tFlip90 = 31.250 us @ 50.000 V) from bulk FID by Find_PulseDurationFID
    HW.TX(iDevice).PaUout2AmplitudeEstimatedX = HW.TX(iDevice).PaUout2AmplitudeX;

  otherwise
    warning('LoadCoil:UnknownCoilName', 'No configuration for coil "%s" found. Using default settings.', HW.TX(iDevice).CoilName);
    return;

end

HW.TX(iDevice).PaUout2AmplitudeEstimated = HW.TX(iDevice).PaUout2Amplitude;
