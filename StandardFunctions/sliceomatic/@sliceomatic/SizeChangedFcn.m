function SizeChangedFcn(this)
%% Adapt positions of axes in parent object
%
% ------------------------------------------------------------------------------
% (C) Copyright 2017-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% Get size of parent in pixels
oldUnits = get(this.hParent, 'Units');
set(this.hParent, 'Units', 'pixels');
parentPosition = get(this.hParent, 'Position');
set(this.hParent, 'Units', oldUnits);

if isempty(parentPosition)
  % figure is being constructed and not yet ready
  return;
end

if this.isSaved
  % additional tasks to re-install necessary callbacks when loading from file
  d = getappdata(this.hParent, 'sliceomatic');
  this.SetupRestoreAxes(d, this.xmesh, this.ymesh, this.zmesh, true);
  this.isSaved = false;
end

%% sliders
sliderWidth = 20;
horzBorderTB = 100;
widthTB = parentPosition(3)-2*horzBorderTB;
if widthTB < 40
  widthTB = 40;
end
horzBorderLR = 40;
buttomBorderLR = 50;
heightLR = parentPosition(4)-130;
if heightLR < 40
  heightLR = 40;
end
vertBorderTop = 120;
vertBorderButtom = 50;
set(this.hSliderX, ...
  'Position', [horzBorderTB, parentPosition(4)-vertBorderTop+sliderWidth, widthTB, sliderWidth]);
set(this.hSliderY, ...
  'Position', [horzBorderLR, buttomBorderLR, sliderWidth, heightLR]);
set(this.hSliderZ, ...
  'Position', [parentPosition(3)-horzBorderLR-sliderWidth, buttomBorderLR, sliderWidth, heightLR]);
set(this.hSliderIso, ...
  'Position', [horzBorderTB, vertBorderButtom, widthTB, sliderWidth]);

%% main axes
horzBorder = horzBorderLR+sliderWidth;
if strcmp(get(this.hSliderX, 'Visible'), 'off')
  vertBorderTop = 80;
  vertBorderButtom = 0;
  horzBorder = 0;
end
widthMain = parentPosition(3)-2*horzBorder;
if widthMain < 160
  widthMain = 160;
end
heightMain = parentPosition(4)-vertBorderTop-vertBorderButtom-sliderWidth;
if heightMain < 120
  heightMain = 120;
end
set(this.hAxes, 'OuterPosition', [horzBorder, vertBorderButtom+sliderWidth+10, widthMain, heightMain])

%% uipanels
horzBorderPanel = 10;
heightPanel = 50;
widthPanel = 80;
set(this.colormapPanel, ...
  'Position', [horzBorderPanel, parentPosition(4)-heightPanel, widthPanel, heightPanel]);
% set(this.alphamapPanel, ...
%   'Position', [parentPosition(3)-horzBorderPanel-widthPanel, parentPosition(4)-heightPanel, widthPanel heightPanel]);
set(this.alphamapPanel, ...
  'Position', [2*horzBorderPanel+widthPanel, parentPosition(4)-heightPanel, widthPanel heightPanel]);
set(this.orientationPanel, ...
  'Position', [parentPosition(3)-horzBorderPanel-widthPanel, parentPosition(4)-heightPanel, widthPanel heightPanel]);

end
