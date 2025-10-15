function temperature = SampleHeater_TargetTemperature_Callback(sampleHeater, temperature)
%% Callback function that is executed when setting the target temperature of the sample heater
%
%   temperature = SampleHeater_TargetTemperature_Callback(sampleHeater, temperature)
%
%
% INPUT:
%
%   sampleHeater
%       Object of type PD.MagnetSampleHeater
%
%   temperature
%       Temperature in degrees Celcius as set by the user for the sample heater
%       target.
%
%
% OUTPUT:
%
%   temperature
%       Possibly adjusted target temperature for the sample heater in degrees
%       Celcius.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


% PID Limits
control_max_p = round(temperature/100*50, 1);
sampleHeater.setSampleHeaterPIDLimitsP([-80, control_max_p]);

% control_max_i = round(max(3, (5/16)*temperature-5), 1);  % @ Ki=0.15
control_max_i = round(max(3, (7/20)*temperature-5), 1);  % @ Ki=0.2
control_min_i = round(min(-10, control_max_i-10), 1);
sampleHeater.setSampleHeaterPIDLimitsI([control_min_i, control_max_i]);


% limit maximum temperature of air flow
user_max_air_temperature = min(395, round(temperature + 1.5*control_max_i, 1));
sampleHeater.setSampleHeaterAirflowTemperatureMax(user_max_air_temperature);

end
