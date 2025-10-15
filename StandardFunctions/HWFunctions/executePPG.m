function [data, errorLevel] = executePPG(commandArray, config, HW)
%% Start pulse program command array and wait for received data
%
%   [data, errorLevel] = executePPG(commandArray, config, HW)
%
% INPUT:
%
%   commandArray
%       uint64 vector with commands for pulse program
%
%   config
%       optional structure with the following fields
%     checkRunning
%         Handle to a function that returns true if the execution should
%         continue or false if should be interrupted.
%     getDigitalInputTime
%         Boolean to indicate if timestamps for the digital inputs should be
%         querried.
%
%   HW
%       Object of type PD.HWClass
%
% OUTPUT:
%
%   data
%       Structure with the following fields.
%       NOTE: The layout of the data in all of these fields might differ from
%             the layout returned by get_data.
%     data
%         Complex data received at the end of the pulse program.
%     overflow
%         Boolean matrix in which samples are marked for which the receiver was
%         saturated for at least one of the non-downsampled input. It has the
%         same size as data.
%     timestamp
%         Structure with fields of the digital input timestamps.
%
%   errorLevel
%       Integer value indicating if an error occured. 0 means no error or
%       warning.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function PD.Talker

