classdef MagnetSampleHeater < PD.MagnetBase
%% Communication with magnet including sample heater
%
% This class is used to control the sample heater and query data from it.
%
% PROPERTIES:
%
%   sampleHeaterTargetTemperature
%       Set or get target temperature of sample heater in degrees C. Setting
%       the target temperature to NaN turns off the air flow. Setting it to a
%       finite value in range, switches the airflow on and start the
%       controller to reach the set temperature.
%
%   sampleHeaterAirflowTemperatureMax
%       Set or get the maximum temperature of the airflow in degrees C. This
%       temperature must be higher than sampleHeaterTargetTemperature for the
%       controller to reach the target temperature. It puts an upper limit to
%       the airflow temperature the controller uses.
%
%   sampleHeaterAirflow
%       Set or get the airflow of the sample heater in percent of maximum.
%
%   sampleHeaterTargetTemperatureMin, sampleHeaterTargetTemperatureMax
%       User-definable limits for allowed target temperature of sample heater
%       in degrees C. If a value outside this range is set for
%       sampleHeaterTargetTemperature, an error is thrown
%
%   sampleHeaterTargetTemperatureDevMax
%       Maximum deviation of the current temperature to the target temperature
%       in degrees C for the test whether the target temperature was reached.
%
%   sampleHeaterTargetTemperatureStdevMax
%       Maximum standard deviation of the recent temperature measurements in
%       degrees C for the test whether the target temperature was reached.
%
%   sampleHeaterTargetTemperatureWaitTime
%       Additional time after the above criterion are reached until
%       waitforSampleHeaterTemperature returns.
%
%   sampleHeaterTargetTemperatureCallback
%       Callback function that is executed just before updating the target
%       temperature of the sample heater. If it is non-empty, it must be a
%       function handle with the following prototype:
%           temperature_out = @fcn(sampleHeater, temperature_in)
%       where sampleHeater is an object of type PD.MagnetSampleHeater,
%       temperature_in is the new temperature as set by the user, and
%       temperature_out is the target temperature that will actually be set
%       for the sample heater. (Default: [])
%
%   sampleHeaterAirflowTemperatureMax
%       Scalar value with the maximum temperature of the air flow in degrees
%       Celcius that the sample heater can use.
%
%   sampleHeaterAirflowTemperatureMaxLimits
%       Two-element vector with increasing values that are used as the limits
%       for sampleHeaterAirflowTemperatureMax. (Default: [0, 120])
%
%   sampleHeaterPIDParameters
%       3-element vector with the coefficients for the PID controller.
%
%
% METHODS:
%
%   checkSampleHeaterTemperature()
%       Returns true if the target temperature is reached (see
%       sampleHeaterTargetTemperatureDevMax and
%       sampleHeaterTargetTemperatureStdevMax).
%
%   waitforSampleHeaterTemperature(verbose)
%       Repeatedly check if the target temperature was reached (see
%       checkSampleHeaterTemperature) and return when it is stable after an
%       additional time (see sampleHeaterTargetTemperatureWaitTime).
%
%   checkMagnetTemperature()
%       Check if any of the temperature sensors exceed their respective
%       limits. In case they do, the sample heating is automatically switched
%       off and this function throws an error.
%
%   getSampleHeaterTemperatureOpticalFiber(channel)
%       Returns the temperature that is measured at the optical fiber port
%       with number "channel". If omitted or empty, channel defaults to 1.
%
%
% ----------------------------------------------------------------------------
% (C) Copyright 2021-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ----------------------------------------------------------------------------

end


%#function PD.MagnetBase
%#function isemptyfield

