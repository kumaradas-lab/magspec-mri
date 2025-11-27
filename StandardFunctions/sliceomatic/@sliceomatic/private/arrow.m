function hArrow = arrow(hSliderAxes, arrowDir, pos, hArrow)
%% Draw an arrow in hSliderAxes at position pos
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%  -1 0 1   -1 0 1     12345     12345
% 5 *---*  5   *     1 *-*     1   *-*
% 4 |   |  4  / \      |  \       /  |
% 3 *   *  3 *   *   0 |   *   0 *   |
% 2  \ /   2 |   |     |  /       \  |
% 1   *    1 *---*  -1 *-*    -1   *-*

switch arrowDir
  case 'down'
    pts = [ 0 1; -1 3; -1 5; 1 5; 1 3 ];
    mp = 'SOM leftright';
  case 'up'
    pts = [ 0 5; 1 3; 1 1; -1 1; -1 3 ];
    mp = 'SOM leftright';
  case 'right'
    pts = [ 5 0; 3 -1; 1 -1; 1 1; 3 1 ];
    mp = 'SOM topbottom';
  case 'left'
    pts = [ 1 0; 3 1; 5 1; 5 -1; 3 -1 ];
    mp = 'SOM topbottom';
end

f = [1:5, 1];

% Modify the arrows to look reasonable no matter the data aspect ratio.
% FIXME: The resize function should take care of the aspect ratio of the arrows.
if any(strcmp(arrowDir, {'up', 'down'}))
  lim = get(hSliderAxes, 'XLim');
  fivep = abs(lim(1)-lim(2))/15/5;
  pts(:,1) = pts(:,1)*fivep + pos;
else
  lim = get(hSliderAxes, 'YLim');
  fivep = abs(lim(1)-lim(2))/15/5;
  pts(:,2) = pts(:,2)*fivep + pos;
end

if nargin > 3
  % use existing arrow
  set(hArrow, 'Vertices', pts);
else
  % create patch and add app data
  hArrow = patch('Vertices', pts, 'Faces', f, ...
    'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'k', ...
    'LineWidth', 1, 'Parent', hSliderAxes, 'Tag', 'sliceomaticarrow');
  setappdata(hArrow, 'motionpointer', mp);
end
setappdata(hArrow, 'arrowcenter', pos);

end
