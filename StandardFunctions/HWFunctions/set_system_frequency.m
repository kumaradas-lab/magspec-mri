function fSystem = set_system_frequency(HW, fSystem, force)
%% Change the system frequency the MRT device is using
%
%     fSystem = set_system_frequency(HW, fSystem, force)
%
% Change the system frequency used by all connected devices. By default, a
% system frequency of 125 MHz is used (corresponding to a 8 ns clock).
%
%
% INPUT:
%
%   HW
%       PD.HWClass object
%
%   fSystem
%       New system frequency in Hz. The value is rounded to the closest possible
%       value. (Default: 125e6)
%
%   force
%       Set the system frequency even if it looks like the divisors do not
%       change from the values stored in HW. (Default: false)
%
%
% OUTPUT:
%
%   fSystem
%       Actually used system frequency in Hz.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2022-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

end


%#function PD.Talker

