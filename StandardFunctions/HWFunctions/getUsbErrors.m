function [usbErrors, timeSinceLastCall] = getUsbErrors(talker)
%% Check USB packet collisions
%
%     [usbErrors, timeSinceLastCall] = getUsbErrors(talker)
%
% INPUT:
%
%   talker
%           object of type MTalker.
%
% OUTPUT:
%
%   usbErrors
%           Number of USB errors since last call (possible uint16 overflow) if
%           supported. Otherwise, the function returns -1.
%
%   timeSinceLastCall
%           Time since last call of the function in seconds.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end
