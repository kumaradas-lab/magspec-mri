function isAutoStarted = StartTcpServer(mySave, iDevice)
%% Start TCP Server for MMRT connection
%
%     isAutoStarted = StartTcpServer(mySave, iDevice)
%
%
% INPUT:
%
%   mySave
%       Structure with information about the TCP server connection.
%       `mySave.reg.BinPath` must have been set to the root path where the
%       "tcpServer" folder with the files for the TCP server executable are
%       installed.
%
%
%   iDevice
%       Index of the device for which the server should be started (for multiple
%       connected MMRT devices).
%       (Default: 1)
%
%
% OUTPUT:
%
%   isAutoStarted
%       Boolean value that indicates whether the TCP server was automatically
%       started.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

end


%#function getRunningTasks
%#function isemptyfield

