function checkEvalBuffer(buf)
%% Check type of buf and evaluate as a function
%
%     checkEvalBuffer(buf)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if isempty(buf)
  return;
end

if iscell(buf)
  if ~strcmp(func2str(buf{1}), 'localModeWindowButtonUpFcn')
    % do not execute standard buttonupfcn in zoom,pan,rotate modes
    feval(buf{1}, [], [], buf{2:end});
  end
elseif ischar(buf)
  eval(buf);
elseif isa(buf, 'function_handle')
  buf();
end

end
