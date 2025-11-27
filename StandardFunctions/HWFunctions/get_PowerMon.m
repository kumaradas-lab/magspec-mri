function powerMon = get_PowerMon(doPlot, average, tightTolerances, iDevice)
%% Get power status of MMRT
%
%     powerMon = get_PowerMon(doPlot, average, tightTolerances, iDevice)
%   alternatively:
%     powerMon = get_PowerMon(devStatus)
%
% This function queries and returns voltage and current readings from various
% sensors in the console. It can be used to isolate potential hardware issues.
%
%
% INPUT:
%
% For the first calling form, the sensor values are read from the device:
%
%   doPlot
%       Boolean value to indicate whether a plot with the measured voltage and
%       current values should be shown.
%       (Default: false)
%
%   average
%       Number of repeated readings from the voltage and current sensors.
%       (Default: average)
%
%   tightTolerances
%       Boolean value to indicate whether to use tight tolerances (true) or
%       relaxed tolerances (false) for the voltage and current readings.
%       (Default: true)
%
%   iDevice
%       Index for the devices (in case multiple devices are connected
%       simultaneously).
%       (Default: 1)
%
% For the alternative calling form, a structure with voltage and current values
% can be provided (without reading from the device in this function):
%
%   devStatus
%       Structure with values for each sensor. For this calling form, loose
%       tolerances are used.
%
%
% OUTPUT:
%
%   powerMon
%       Structure containing value with the readings for all sensors and derived
%       properties. The returned fields include:
%
%     average
%         Value of the input parameter "average".
%
%     RUNTIME_FPGA
%         Time in seconds since the intialization of the device.
%
%     U_V_Batt, U_1V8, U_1212, U_3V8, I_Vin2, U_55, I_Vin, U_Vin, U_3V3, U_ADC,
%     U_1V0, U_2V5, T_FPGA
%         [average x 1] column vectors with the readings from the sensors.
%
%     I_System, P_System
%         [average x 1] column vectors with derived values.
%
%     mean
%         Structure with the average value of all readings and derived values.
%
%     P_System_OK, U_Vin_OK, U_internal_OK, I_OK
%         Boolean values indicating whether all measured values are within
%         tolerance.
%
%     U_internal_1Pro_OK
%         Boolean value indicating whether all measured voltages deviate less
%         than one percent from their nominal values.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2025 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function PD.Talker

