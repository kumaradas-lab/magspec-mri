function res = setProcessPriority(priority)  %#ok<STOUT,INUSD>
%% Set the scheduling priority of the currently running process
%
%     res = setProcessPriority(priority)
%
% This function has any effect only on Windows. On other platforms, it returns
% -1 without actually changing any system settings.
%
% INPUT:
%   priority
%       (Optional) integer between 0 and 5 selecting the scheduling priority:
%       0: realtime (use with care!)
%       1: high
%       2: above normal
%       3: normal (default)
%       4: below normal
%       5: idle
%
% OUTPUT:
%   res
%       If the function succeeds, the return value is nonzero.
%       If the function fails, the return value is zero.
%
% For more details see:
% https://docs.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-setpriorityclass
%
% ------------------------------------------------------------------------------
% (C) Copyright 2020-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com

end
