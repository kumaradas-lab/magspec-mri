function coilTemp = check_CoilTemperature(coilTemp, Seq)
%% Model temperature of coil (gradient system or rf)
%
%   coilTemp = check_CoilTemperature(coilTemp)
%
% Used model:
% * Coil is heated according to coilTemp.Work and coilTemp.Time.
% * Uniform warming of the probe head with coilTemp.CoilThermalCapacity.
% * Coils are damaged if the temperature exceeds coilTemp.CoilMaxTemperature.
% * Heat dissipates with coilTemp.Dissipation at coilTemp.CoilMaxTemperature.
%   No heat dissipation at coilTemp.CoilTemperature.  Linear model in between.
% * The modeled temperature is returned in coilTemp.Temperature.
%
%
% Note: This is an internal function without any input validation.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2021 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

end
