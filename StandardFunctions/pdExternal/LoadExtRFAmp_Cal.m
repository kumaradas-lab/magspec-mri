%% Load the calibration file for the RF-100 amplifier
% This script is called by the LoadRF100_NN scripts.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2019-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

if ~exist('iDevice', 'var'), iDevice = 1; end

if any(HW.TX(iDevice).Uout2PaUout ~= 1)
  switch UseExtRFAmpSwitch
    case 0
      fileName = sprintf('Calibration%s_%02d.mat', ...
        HW.TX(iDevice).ExtRFType, HW.TX(iDevice).ExtRFSN);
    case 1
      fileName = sprintf('Calibration%s_%02d_Switch.mat', ...
        HW.TX(iDevice).ExtRFType, HW.TX(iDevice).ExtRFSN);
    case 2
      fileName = sprintf('Calibration%s_%02d_ExternalSwitch.mat', ...
        HW.TX(iDevice).ExtRFType, HW.TX(iDevice).ExtRFSN);
  end
  % search for calibration files only in pdExternal folder
  fileName = fullfile(HW.RootPath, 'StandardFunctions', 'pdExternal', fileName);
  if exist(fileName, 'file')
    load(fileName);
    if isfield(CalibrationRfAmp, 'TX')
      HW.TX(iDevice).CalibrationRfAmp = CalibrationRfAmp.TX;
    else
      HW.TX(iDevice).CalibrationRfAmp = CalibrationRfAmp;
    end
    clear CalibrationRfAmp
    if numel(HW.TX(iDevice).CalibrationRfAmp) >= HW.TX(iDevice).ChannelDef && ...
        ~isemptyfield(HW.TX(iDevice).CalibrationRfAmp(HW.TX(iDevice).ChannelDef), 'withSwitch') && ...
        HW.TX(iDevice).CalibrationRfAmp(HW.TX(iDevice).ChannelDef).withSwitch ~= UseExtRFAmpSwitch
      warning('PD:LoadExtRFAmp_Cal:InconsistentSwitch', ...
        ['UseExtRFAmpSwitch was set to %d. ', ...
        'But HW.TX.CalibrationRfAmp.withSwitch in the calibration file is %d.'], ...
        UseExtRFAmpSwitch, HW.TX(iDevice).CalibrationRfAmp(HW.TX(iDevice).ChannelDef).withSwitch);
    end
  else
    warning('PD:LoadExtRFAmp_Cal:NoCalibrationFound', ...
      ['No calibration file for the rf amplifier was found. ', ...
      'Please, contact Pure Devices.']);
  end
end
if isemptyfield(HW.TX(iDevice), {'CalibrationRfAmp', 'withSwitch'})
  HW.TX(iDevice).CalibrationRfAmp(HW.TX(iDevice).ChannelDef).withSwitch = UseExtRFAmpSwitch;
end
clear UseExtRFAmpSwitch

if isa(HW, 'PD.HWClass')
  % This check works here only if HW is an object.
  % For a HW structure, the used field are calculated later on and cannot yet be
  % used here.
  if any(isnan(HW.TX(iDevice).Uout2PaUoutCalibrationGainMax))
    warning('PD:LoadExtRFAmp_Cal:InvalidUout2PaUoutCalibrationGainMax', ...
      ['Unable to calculate "HW.TX.Uout2PaUoutCalibrationGainMax".\n', ...
      'Are you sure the amplifier is rated for the set voltage?']);
  end
  if any(isnan(HW.TX(iDevice).Uout2PaUoutCalibrationGainDef))
    warning('PD:LoadExtRFAmp_Cal:InvalidUout2PaUoutCalibrationGainDef', ...
      ['Unable to calculate "HW.TX.Uout2PaUoutCalibrationGainDef".\n', ...
      'Are you sure the amplifier is rated for the set voltage?']);
  end
end
