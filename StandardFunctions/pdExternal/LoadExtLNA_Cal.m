%% Load the calibration file for an external low-noise amplifier
% This script is called by the LoadTxRxLNA_NN scripts.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2022-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

if ~exist('iDevice', 'var'), iDevice = 1; end

if any(HW.TX(iDevice).Uout2PaUout ~= 1)
  fileName = sprintf('CalibrationLNA_%02d.mat', HW.RX(iDevice).LnaSN);
  % search for calibration files only in pdExternal folder
  fileName = fullfile(HW.RootPath, 'StandardFunctions', 'pdExternal', fileName);
  if exist(fileName, 'file')
    load(fileName);
    % Use the same mechanism that is used for rf amplifiers for the frequency
    % dependent attenuation of the TX channel with the switch in the LNA.
    HW.TX(iDevice).CalibrationRfAmp = CalibrationLNA.TX;
    clear CalibrationLNA;
  else
    warning('PD:LoadExtLNA_Cal:NoCalibrationFound', ...
      ['No calibration file for the low-noise amplifier #%02d was found. ', ...
      'Please, contact Pure Devices.'], HW.RX(iDevice).LnaSN);
  end
  clear fileName;
end

if isa(HW, 'PD.HWClass')
  % This check works here only if HW is an object.
  % For a HW structure, the used fields are calculated later on and cannot yet
  % be used here.
  if any(isnan(HW.TX(iDevice).Uout2PaUoutCalibrationGainMax))
    warning('PD:LoadExtLNA_Cal:InvalidUout2PaUoutCalibrationGainMax', ...
      ['Unable to calculate "HW.TX.Uout2PaUoutCalibrationGainMax".\n', ...
      'Are you sure the low-noise amplifier is rated for the set voltage?']);
  end
  if any(isnan(HW.TX(iDevice).Uout2PaUoutCalibrationGainDef))
    warning('PD:LoadExtLNA_Cal:InvalidUout2PaUoutCalibrationGainDef', ...
      ['Unable to calculate "HW.TX.Uout2PaUoutCalibrationGainDef".\n', ...
      'Are you sure the low-noise amplifier is rated for the set voltage?']);
  end
end
