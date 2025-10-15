classdef MagnetBase < handle
%% Class with the base properties for communication with the magnet (USB)
%
% A user should not create an instance of this class manually. Instead, it
% will be created automatically as HW.Magnet if HW.Magnet.useBase is set to
% `true` in the LoadMySystem.m script file.
% Among other things, it can be used to check the temperature of magnets with
% temperature stabilization unit.
%
%
% PROPERTIES:
%
%   version
%       Read-only. String with version information.
%
%   macAddress
%       Read-only. String with MAC address of the connected magnet
%       (identifier).
%
%   ipAddress
%       Read-only. String with the IP address of the connected magnet.
%
%   minRuntime
%       Minimum runtime in seconds before the magnet's temperature can be
%       considered as stable. (Default: 400)
%
%   maxCurrentTemperatureDev
%       Tolerance for the temperature of the magnet in degrees C for
%       considering it stable. (Default: 0.05)
%
%   maxTemperatureStdev
%       Maximum allowed standard deviation of the temperature of the magnet in
%       degrees C for considering is stable. (Default: 0.05)
%
%   useBase
%       Boolean value to indicate whether the communication link with the
%       temperature stabilization unit in the magnet should be established
%       during "LoadSystem".
%       (Default: false)
%
%   useSampleHeater
%       Boolean value to indicate wether the communication link with the
%       sample heater unit in the magnet should be established during
%       "LoadSystem". This requires that `useBase` is also set to `true`. See
%       also: PD.MagnetSampleHeater (Default: true)
%
%   enableMagnetHeating
%       Boolean value to indicate whether the PID controller for the magnet
%       heating is enabled (true) or paused (false). This requires that
%       `useBase` is set to `true` and the magnet firmware is new enough.
%       (Default: true)
%
%
% METHODS:
%
%   runtime = HW.Magnet.getRuntime()
%       Query the current runtime of the magnet in seconds.
%
%   targetTemp = HW.Magnet.getTemperatureTarget()
%       Query the current target temperature of the magnet in degrees C.
%
%   dispTemp = HW.Magnet.getTemperatureDisplay()
%       Query the temperature in degrees C that is currently displayed on the
%       display of the magnet.
%
%   tempStd = HW.Magnet.getTemperatureStdev()
%       Query the current standard deviation of the magnet's temperature
%       (derived from recently measured values of the temperature sensors).
%
%   tempState = checkTemperature()
%       Check if the temperature of the magnet can be considered stable.
%       Returns an enumeration object of type PD.MagnetReadyState.
%
%   HW.Magnet.setVoltageDAC(value, chIdx, wait)
%       Set the voltage for the AUX output channels.
%
%   voltageDAC = HW.Magnet.getVoltageDAC()
%       Get the voltage that is currently set for the AUX output channels.
%
%   voltageADC = HW.Magnet.getVoltageADC()
%       Get voltages at input channels of LTC2493.
%
%
% ----------------------------------------------------------------------------
% (C) Copyright 2021-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ----------------------------------------------------------------------------

end


%#function PD.MagnetReadyState
%#function isemptyfield

