function use_ExternalClock(tf, iDevices, force)
%% Switch to extern clock
%
%   use_ExternalClock(tf, iDevices, force)
%
%
% INPUT:
%
%   tf
%       Boolean. If true, switch to external 10MHz clock generator on Clk port.
%       If false, use internal clock. (Default: true)
%
%
%   iDevices
%       Indices of device(s) to be switched. (Default: 1:numDevices)
%
%   force
%       Set the external clock usage even if it looks like it didn't change from
%       the last function call. (Default: false)
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2022-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function PD.Talker

