function tasks = getRunningTasks(procName)
%% Get structure with a list of all running tasks
%
%       tasks = getRunningTasks(procName)
%
%
% INPUT:
%
%   procName
%       Name of the process to search in the list of running tasks. If empty,
%       a list of all running tasks is returned. Default: ''.
%
%
% OUTPUT:
%
%   tasks
%       Structure with all running tasks as returned by the operating system
%       with the following fields:
%
%     name
%         Name of the task as seen in the Task Manager on Windows, or the full
%         command that was used to start the process on macOS or Linux.
%
%     PID
%         Process ID of the task.
%
%     sessionName
%         Name of the session in which the task is running on Windows, or empty
%         on macOS or Linux.
%
%     sessionID
%         Number of the session in which the task is running.
%
%     mem
%         Memory currently used by the process in KB on Windows or Linux, or
%         empty on macOS.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end
