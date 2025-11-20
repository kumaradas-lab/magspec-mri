function totalMagnPol = get_CalibrationSampleTotalMagneticPolarization(HW, T_C, iDevice)
%% Calculate expected total magnetic polarization for selected calibration sample
%
%   totalMagnPol = get_CalibrationSampleTotalMagneticPolarization(HW, T_C, iDevice)
%
%
% INPUT:
%
%   HW
%       HW object or structure
%
%   T_C
%       Thermodynamic equilibrium temperature in degrees Celsius.
%       (Default: Temperature "TempCal" during calibration of the sample if
%       available or 30 otherwise)
%
%   iDevice
%       Index in case multiple devices are connected.
%       (Default: 1)
%
%
% OUTPUT:
%
%   totalMagnPol
%       Total magnetic polarization in T*m^3 of the currently selected
%       calibration sample (i.e., HW.RX.CalibrationSampleName) at the current
%       temperature (T_C) and magnetic field amplitude (HW.B0) in thermodynamic
%       equilibrium.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input parameters
if nargin < 2
  T_C = [];
end
if nargin < 3 || isempty(iDevice)
  iDevice = 1;
end


%% read calibration values from file
CalibrationSample = LoadCalibrationSample(HW, iDevice);

% get line with data for currently selected calibration sample
selectedCalibrationSample = ...
  CalibrationSample(strcmp(CalibrationSample.UID, HW.RX(iDevice).CalibrationSampleName), :);

if isempty(selectedCalibrationSample)
  error('PD:get_CalibrationSampleMagneticPolarization:UnknownSample', ...
    'No calibration data for sample with UID "%s" could be found.', ...
    HW.RX(iDevice).CalibrationSampleName);
end

% FIXME: Do we need to check whether all required fields are available?


%% get expected total magnetic polarization amplitude
if isempty(T_C)
  if any(strcmp(selectedCalibrationSample.Properties.VariableNames, 'TempCal'))
    % assume that the measurement is at the same temperature as the calibration
    T_C = selectedCalibrationSample.TempCal;
  else
    T_C = 30;
  end
end

totalMagnPol = get_TotalMagneticPolarization(HW, selectedCalibrationSample.mol, ...
  T_C, HW.B0, selectedCalibrationSample.gamma, selectedCalibrationSample.I);

end
