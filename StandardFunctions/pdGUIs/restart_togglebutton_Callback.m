function restart_togglebutton_Callback(h, scriptname)
%% Callback function for "Run" togglebutton that re-starts a script
%
%     restart_togglebutton_Callback(h, scriptname)
%
% INPUT:
%
%   h
%       Handle to the pushbutton.
%
%   scriptname
%       String with the name of the script/function that is executed in the base
%       workspace when the toggle button is (re-)activated.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

persistent isRunning

if get(h, 'Value')
  % button is activated

  if isRunning
    % avoid nesting callbacks if button was pressed too fast repeatedly
    return;
  end

  % execute script
  isRunning = true;  %#ok<NASGU>
  evalin('base', scriptname);
  isRunning = false;
else
  % button is de-activated

  % Disable button to avoid re-activating while the script is still running.
  % That can still happen occasionally if the button is pressed again before
  % this callback function has been executed after de-activating.
  % It must be (re-)enabled in the script or function.

  set(h, 'Enable', 'off');
end

end
