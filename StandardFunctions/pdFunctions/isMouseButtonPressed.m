function buttonPressed = isMouseButtonPressed()
%% Check if mouse button is pressed at the moment
%
%   buttonPressed = isMouseButtonPressed()
%
% Returns the sum of the following int32 values:
% 0: no button is pressed
% 1: left mouse button is pressed
% 2: right mouse button is pressed
% 4: middle mouse button is pressed (not scrolling)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2019 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

%% return early if not on Windows
buttonPressed = 0;
if ~ispc()
  warning('PD:isMouseButtonPressed:NoWindows', ...
    'isMouseButtonPressed only works on Windows.');
  return;
end

%% initialize
if ~libisloaded('user32_GetAsyncKeyState')
  loadlibrary('user32', @isMouseButtonPressed_user32_proto64, ...
    'alias', 'user32_GetAsyncKeyState');
end

%% get mouse button state
if libisloaded('user32_GetAsyncKeyState')
  calllib('user32_GetAsyncKeyState', 'GetAsyncKeyState', int32(1));
  calllib('user32_GetAsyncKeyState', 'GetAsyncKeyState', int32(2));
  calllib('user32_GetAsyncKeyState', 'GetAsyncKeyState', int32(4));
  % Call the functions twice
  left = calllib('user32_GetAsyncKeyState', 'GetAsyncKeyState', int32(1)) < 0;
  right = calllib('user32_GetAsyncKeyState', 'GetAsyncKeyState', int32(2)) < 0;
  middle = calllib('user32_GetAsyncKeyState', 'GetAsyncKeyState', int32(4)) < 0;
  buttonPressed = sum(int32([1 2 4]) .* int32([left right middle]), 'native');
end

end
