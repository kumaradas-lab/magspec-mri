classdef SampleHeaterBase < handle
%% Class for communication with directly connected sample heater (USB)
%
% A user should not create an instance of this class manually. Instead, it
% will be created automatically as HW.SampleHeater if HW.SampleHeater.use
% is set to `true` in the LoadMySystem.m script file.
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
%   use
%       Boolean value to indicate wether the communication link with a
%       sample heater unit that is directly connected to the PC should be
%       established during "LoadSystem".
%       (Default: false)
%
%   targetTemperature
%       Set or get target temperature of sample heater in degrees C. Setting
%       the target temperature to NaN turns off the controller. Setting it to
%       a finite value in range, starts the controller to reach the set
%       temperature.
%
%   targetTemperatureMin, targetTemperatureMax
%       Read-only limits for allowed target temperature of sample heater
%       in degrees C. If a value outside this range is set for
%       sampleHeaterTargetTemperature, an error is thrown.
%
%   targetTemperatureDevMax
%       Maximum deviation of the current temperature to the target temperature
%       in degrees C for the test whether the target temperature was reached.
%
%   targetTemperatureStdevMax
%       Maximum standard deviation of the recent temperature measurements in
%       degrees C for the test whether the target temperature was reached.
%
%   targetTemperatureWaitTime
%       Additional time after the above criterion are reached until
%       waitforSampleHeaterTemperature returns.
%
%   targetTemperatureCallback
%       Callback function that is executed just before updating the target
%       temperature of the sample heater. If it is non-empty, it must be a
%       function handle with the following prototype:
%           temperature_out = @fcn(sampleHeater, temperature_in)
%       where sampleHeater is an object of type PD.MagnetSampleHeater,
%       temperature_in is the new temperature as set by the user, and
%       temperature_out is the target temperature that will actually be set
%       for the sample heater. (Default: [])
%
%   PIDParameters
%       3-element vector with the coefficients for the PID controller.
%
%   enableCooling
%       The Peltier cooler is only enabled if this property is set to true
%       explicitly. Make sure to use dry air when enabling cooling below room
%       temperature.
%
%
% METHODS:
%
%   runtime = HW.SampleHeater.getRuntime()
%       Query the current runtime of the magnet in seconds.
%
%   HW.SampleHeater.checkTemperature()
%       Returns true if the target temperature is reached (see
%       sampleHeaterTargetTemperatureDevMax and
%       sampleHeaterTargetTemperatureStdevMax).
%
%   HW.SampleHeater.waitForTemperature(verbose)
%       Repeatedly check if the target temperature was reached (see
%       checkSampleHeaterTemperature) and return when it is stable after an
%       additional time (see sampleHeaterTargetTemperatureWaitTime).
%
%   isOn = HW.SampleHeater.isOnAirflow()
%       Returns true if the air flow is switched on. Returns false otherwise.
%
%   airPressureValue = HW.SampleHeater.getAirPressureValue()
%       Get the value for the air pressure in mbar.
%
%   tempAvg = HW.SampleHeater.getTemperatureAverage()
%       Get mean temperature (of the last N samplings, ~25 s) in degrees C.
%
%   tempStd = HW.SampleHeater.getTemperatureStdev()
%       Get standard deviation of recently measured temperatures (~25 s) in
%       degrees C.
%
%   tempCurrent = HW.SampleHeater.getTemperaturePort1()
%       Returns the temperature that is measured at the thermo couple
%       connected to port 1. Returns NaN if no thermo couple is connected.
%
%   tempCurrent = HW.SampleHeater.getTemperaturePort2()
%       Returns the temperature that is measured at the thermo couple
%       connected to port 2. Returns NaN if no thermo couple is connected.
%
%
% ----------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ----------------------------------------------------------------------------

end


%#function isemptyfield

