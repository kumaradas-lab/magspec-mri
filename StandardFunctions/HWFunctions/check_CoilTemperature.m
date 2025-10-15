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
% The following differential equation is used to model the accumulated heat in
% the coil:
%
% ODE:   d/dt T(t) = heat - D/CT * (T(t) - T0) / (Tmax-T0)
% solution:
%         ==> T(t) = T0 + heat*tau * (1-exp(-(t-t1)/tau)) + (T1-T0) * exp(-(t-t1)/tau)
% with
%   tau = CT/D * (Tmax-T0)
% where
%   T(t): temperature of the coil in dec C
%   heat: heating during the step in K/s
%   T0:   temperature without any heating in deg C
%   Tmax: temperature with full dissipation in deg C
%   D:    dissipation at Tmax in J/s
%   CT:   thermal capacity of the coil in J/K
%   T1:   temperature at the beginning of the step in deg C
%   t1:   time at the beginning of the step in s
%
%
% Note: This is an internal function without any input validation.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2022 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

end
