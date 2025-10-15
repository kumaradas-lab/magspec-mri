%% Load the calibration file for the RF-100 amplifier
% This script is called by the LoadRF100_NN scripts.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2019-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

if UseRF100Switch
  if any(HW.TX.Uout2PaUout~=1)
    fileName = ['CalibrationRfAmp_' num2str(HW.TX.ExtRFSN, '%02d') 'Switch.mat'];
    if exist(fileName, 'file')
      load(fileName);
      HW.TX.CalibrationRfAmp = CalibrationRfAmp;
      clear CalibrationRfAmp
    end
  end
else
  if any(HW.TX.Uout2PaUout~=1)
    fileName = ['CalibrationRfAmp_' num2str(HW.TX.ExtRFSN, '%02d') '.mat'];
    if exist(fileName, 'file')
      load(fileName);
      HW.TX.CalibrationRfAmp = CalibrationRfAmp;
      clear CalibrationRfAmp
    end
  end
end
clear UseRF100Switch

if isa(HW, 'PD.HW')
  % This check works here only if HW is an object.
  % For a HW structure, the used field are calculated later on and cannot yet be
  % used here.
  if any(isnan(HW.TX.Uout2PaUoutCalibrationGainMax))
    warning('PD:LoadRF100_Cal:InvalidUout2PaUoutCalibrationGainMax', ...
      ['Unable to calculate "HW.TX.Uout2PaUoutCalibrationGainMax".\n', ...
      'Are you sure the amplifier is rated for the set voltage?']);
  end
  if any(isnan(HW.TX.Uout2PaUoutCalibrationGainDef))
    warning('PD:LoadRF100_Cal:InvalidUout2PaUoutCalibrationGainDef', ...
      ['Unable to calculate "HW.TX.Uout2PaUoutCalibrationGainDef".\n', ...
      'Are you sure the amplifier is rated for the set voltage?']);
  end
end
