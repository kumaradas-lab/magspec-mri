function d = figmenus(this, d)
%% Set up sliceomatic's GUI menues
%
%   d = this.figmenus(d)
%
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% Main Figure Menu
if this.hFigure == this.hParent
  if ~isappdata(this.hFigure, 'sliceomatic_oldmenubar')
    setappdata(this.hFigure, 'sliceomatic_oldmenubar', get(this.hFigure, 'Menubar'));
  end
  set(this.hFigure, 'Menubar', 'none');
end

% File menu
currentUiMenues = findobj(this.hFigure, 'Type', 'uimenu', '-depth', 1);
if isempty(currentUiMenues) || ...
    ~any(~cellfun(@isempty, strfind(lower(get(currentUiMenues, 'Label')), 'file', 'ForceCellOutput', true)))
  % Only create file menu if there is not yet one with a similar name
  this.fileMenu = uimenu(this.hFigure, 'Label', 'File');
  % Copy the main axes to a new static figure
  d.fcopy = uimenu(this.fileMenu, ...
    'Label',    'Copy Main Axes', ...
    'Callback', @(hObject, eventData) this.callbacks('copy'));
  if ispc() % Windows only
    % Copy the main axes and sliders to the clipboard
    d.fclipboard = uimenu(this.fileMenu, ...
      'Label',    'Copy to Clipboard', ...
      'Callback', @(h,e) this.CopyToClipboard());
  end
  d.fprint  = uimenu(this.fileMenu, ...
    'Label',    'Print...', ...
    'Callback', @(hObject, eventData) this.callbacks('print'));
  d.faviexport  = uimenu(this.fileMenu, ...
    'Label',    'Export as Video...', ...
    'Callback', @(hObject, eventData) this.ExportAvi());
  if exist('Image3D2VolView', 'file') == 2
    d.fvolviewexport  = uimenu(this.fileMenu, ...
      'Label',    'Export to VolView...', ...
      'Callback', @(hObject, eventData) this.ExportVolView());
  end
  d.fsaveprefs = uimenu(this.fileMenu, ...
    'Label',    'Save Preferences', ...
    'Callback', {@SavePrefs, this});

  % How to get these props onto the print figure?
  %d.fprints = uimenu(this.fileMenu,'label','Print Setup...','callback','printdlg -setup');
  % ---
  d.fexit = uimenu(this.fileMenu, ...
    'Label',    'Close', ...
    'Callback', 'closereq', ...
    'Separator', 'on');
end

% Controls Menu
this.controlsMenu = uimenu(this.hFigure, ...
  'Label',    'Controls', ...
  'Callback', {@controlmenu, this});
if exist('uitoolfactory', 'file') == 2
  d.anntoolbar = uimenu(this.controlsMenu, ...
    'Label',    'Annotations toolbar', ...
    'Callback', @(hObject, eventData) this.callbacks('annotationtoolbar'));
end
d.camtoolbar = uimenu(this.controlsMenu, ...
  'Label',    'Camera toolbar', ...
  'Callback', @(hObject, eventData) this.callbacks('cameratoolbar'));
d.dcalpha = uimenu(this.controlsMenu, 'Label', 'Controls Transparency');
d.dcalpha1 = uimenu(d.dcalpha, ...
  'Label',    '1.0', ...
  'Callback', @(hObject, eventData) this.callbacks('controlalpha', 1));
d.dcalpha8 = uimenu(d.dcalpha, ...
  'Label',    '0.8', ...
  'Callback', @(hObject, eventData) this.callbacks('controlalpha', .8));
d.dcalpha6 = uimenu(d.dcalpha, ...
  'Label',    '0.6', ...
  'Callback', @(hObject, eventData) this.callbacks('controlalpha', .6));
d.dcalpha5 = uimenu(d.dcalpha, ...
  'Label',    '0.5', ...
  'Callback', @(hObject, eventData) this.callbacks('controlalpha', .5));
d.dcalpha4 = uimenu(d.dcalpha, ...
  'Label',    '0.4', ...
  'Callback', @(hObject, eventData) this.callbacks('controlalpha', .4));
d.dcalpha2 = uimenu(d.dcalpha, ...
  'Label',    '0.2', ...
  'Callback', @(hObject, eventData) this.callbacks('controlalpha', .2));
d.dcalpha0 = uimenu(d.dcalpha, ...
  'Label',    '0.0', ...
  'Callback', @(hObject, eventData) this.callbacks('controlalpha', 0));
d.dcanimstep = uimenu(this.controlsMenu, ...
  'Label',    'Animation', ...
  'Callback', @(hObject, eventData) this.callbacks('toggleanimation'));
d.dclabels = uimenu(this.controlsMenu, ...
  'Label',    'Tick Labels', ...
  'Callback', @(hObject, eventData) this.callbacks('controllabels'));
d.dcvis = uimenu(this.controlsMenu, ...
  'Label',    'Visible', ...
  'Callback', @(hObject, eventData) this.callbacks('controlvisible'));
% d.dsetrange= uimenu(this.controlsMenu','label','Set Range','callback','@setvolumerange');
% d.dcslice = uimenu(this.controlsMenu,'label','Slice Controls','callback','sliceomatic useslicecontrols');
% d.dciso   = uimenu(this.controlsMenu,'label','Iso Surface Control','callback','sliceomatic useisocontrols','separator','on');

% Remove this once we have more controls to enable and disable.
%  set(this.controlsMenu,'vis','off');

% Default for new slices menu
this.objectDefaultsMenu = uimenu(this.hFigure, ...
  'Label',    'Object Defaults', ...
  'Callback', {@defaultmenu, this.hParent});
d.dfacet  = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Slice Color Faceted', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultfaceted'));
d.dflat   = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Slice Color Flat', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultflat'));
d.dinterp = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Slice Color Interp', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultinterp'));
d.dtex    = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Slice Color Texture', ...
  'Callback', @(hObject, eventData) this.callbacks('defaulttexture'));
d.dcnone  = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Slice Color None', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultcolornone'));
d.dtnone  = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Slice Transparency None', ...
  'Callback', @(hObject, eventData) this.callbacks('defaulttransnone'), ...
  'Separator', 'on');
d.dtflat  = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Slice Transparency Flat', ...
  'Callback', @(hObject, eventData) this.callbacks('defaulttransflat'));
d.dtinterp= uimenu(this.objectDefaultsMenu, ...
  'Label',    'Slice Transparency Interp', ...
  'Callback', @(hObject, eventData) this.callbacks('defaulttransinterp'));
d.dttex   = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Slice Transparency Texture', ...
  'Callback', @(hObject, eventData) this.callbacks('defaulttranstexture'));
d.dlflat  = uimenu(this.objectDefaultsMenu, ...
  'Label',    'IsoSurface Lighting Flat', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultlightflat'), ...
  'Separator', 'on');
d.dlsmooth= uimenu(this.objectDefaultsMenu, ...
  'Label',    'IsoSurface Lighting Smooth', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultlightsmooth'));
d.disot02 = uimenu(this.objectDefaultsMenu, ...
  'Label',    'IsoSurface Transparency 0.2', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultIsoAlpha', .2), ...
  'Separator', 'on');
d.disot05 = uimenu(this.objectDefaultsMenu, ...
  'Label',    'IsoSurface Transparency 0.5', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultIsoAlpha', .5));
d.disot08 = uimenu(this.objectDefaultsMenu, ...
  'Label',    'IsoSurface Transparency 0.8', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultIsoAlpha', .8));
d.disot10 = uimenu(this.objectDefaultsMenu, ...
  'Label',    'IsoSurface Transparency None', ...
  'Callback', @(hObject, eventData) this.callbacks( 'defaultIsoAlpha', 1));
%d.dcsmooth= uimenu(this.objectDefaultsMenu,'label','Contour Line Smoothing','callback','sliceomatic defaultcontoursmooth');
d.dcflat  = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Contour Color Flat', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultcontourflat'), ...
  'Separator', 'on');
d.dcinterp= uimenu(this.objectDefaultsMenu, ...
  'Label',    'Contour Color Interp', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultcontourinterp'));
d.dcblack = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Contour Color Black', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultcontourblack'));
d.dcwhite = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Contour Color White', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultcontourwhite'));
d.dclinew = uimenu(this.objectDefaultsMenu, ...
  'Label',    'Contour Line Width');
d.dcl1    = uimenu(d.dclinew, ...
  'Label',    '1', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultcontourlinewidth', 1));
d.dcl2    = uimenu(d.dclinew, ...
  'Label',    '2', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultcontourlinewidth', 2));
d.dcl3    = uimenu(d.dclinew, ...
  'Label',    '3', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultcontourlinewidth', 3));
d.dcl4    = uimenu(d.dclinew, ...
  'Label',    '4', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultcontourlinewidth', 4));
d.dcl5    = uimenu(d.dclinew, ...
  'Label',    '5', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultcontourlinewidth', 5));
d.dcl6    = uimenu(d.dclinew, ...
  'Label',    '6', ...
  'Callback', @(hObject, eventData) this.callbacks('defaultcontourlinewidth', 6));

d.defcolor = 'texture';
%   d.defalpha = 'texture';
d.defalpha = 'none';
d.deflight = 'smooth';
d.defisoalpha = 0.5;
d.defcontourcolor = 'black';
d.defcontourlinewidth = 1;
% This exposes an unpleasant R14 bug
d.defcontoursmooth = 'off';

% investigate hardware opengl.
inc = 0;
try  %#ok<TRYNC>
  if verLessThan('Matlab', '9.6')
    % The function "opengl" is deprecated in Matlab R2019a and later.
    od = opengl('data');
    if isfield(od, 'Software')
      % R14 version of MATLAB
      if ~od.Software
        inc = 10;
      end
    else
      % Older version of MATLAB
      if ~(strcmp(od.Renderer, 'Mesa X11') || ...
          strcmp(od.Renderer, 'GDI Generic'))
        inc = 10;
      end
    end
  else
    % Use "rendererinfo" in Matlab R2019a and later.
    od = rendererinfo(this.hAxes);
    if ~isfield(od, 'GraphicsRenderer') || ~strcmp(od.GraphicsRenderer, 'OpenGL Hardware')
      inc = 10;
    end
  end
end

d.animincrement = inc;

d = OverrideStickyUserPreferences(this, d);

% Set props for all slices menu
this.allSlicesMenu = uimenu(this.hFigure, ...
  'Label', 'All Slices');
uimenu(this.allSlicesMenu, ...
  'Label',    'Color Faceted', ...
  'Callback', @(hObject, eventData) this.callbacks('allfacet'));
uimenu(this.allSlicesMenu, ...
  'Label',    'Color Flat', ...
  'Callback', @(hObject, eventData) this.callbacks('allflat'));
uimenu(this.allSlicesMenu, ...
  'Label',    'Color Interp', ...
  'Callback', @(hObject, eventData) this.callbacks('allinterp'));
uimenu(this.allSlicesMenu, ...
  'Label',    'Color Texture', ...
  'Callback', @(hObject, eventData) this.callbacks('alltex'));
uimenu(this.allSlicesMenu, ...
  'Label',    'Color None', ...
  'Callback', @(hObject, eventData) this.callbacks('allnone'));
uimenu(this.allSlicesMenu, ...
  'Label',    'Transparency None', ...
  'Callback', @(hObject, eventData) this.callbacks('alltnone'), ...
  'Separator', 'on');
uimenu(this.allSlicesMenu, ...
  'Label',    'Transparency 0.5', ...
  'Callback', @(hObject, eventData) this.callbacks('alltp5'));
uimenu(this.allSlicesMenu, ...
  'Label',    'Transparency Flat', ...
  'Callback', @(hObject, eventData) this.callbacks('alltflat'));
uimenu(this.allSlicesMenu, ...
  'Label',    'Transparency Interp', ...
  'Callback', @(hObject, eventData) this.callbacks('alltinterp'));
uimenu(this.allSlicesMenu, ...
  'Label',    'Transparency Texture', ...
  'Callback', @(hObject, eventData) this.callbacks('allttex'));

if 0
% Setup Help style options
this.helpMenu = uimenu(this.hFigure, 'Label', 'Help');
uimenu(this.helpMenu, ...
  'Label',    'Help', ...
  'Callback', 'doc sliceomatic/sliceomatic');
uimenu(this.helpMenu, ...
  'Label',    'Check for Updates', ...
  'Callback', 'web http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=764&objectType=FILE');
uimenu(this.helpMenu, ...
  'Label',    'About Author', ...
  'Callback', 'web http://www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do?objectId=803709&objectType=author');
end

%% Context Menus
% Slice Context Menu
this.sliceContextMenu = uicontextmenu('Parent', this.hFigure, 'Callback', {@slicecontextmenu, this.hParent}, ...
  'UserData', this, 'Tag', 'sliceomaticSliceContextMenu');
d.vistog = uimenu(this.sliceContextMenu, ...
  'Label',    'Visibility', ...
  'Callback', @(hObject, eventData) this.callbacks('togglevisible'));
d.uicdelete = uimenu(this.sliceContextMenu, ...
  'Label',    'Delete', ...
  'Callback', @(hObject, eventData) this.callbacks('deleteslice'));
d.smcolorm  = uimenu(this.sliceContextMenu, ...
  'Label',    'Color', ...
  'Separator', 'on');
d.smfacet   = uimenu(d.smcolorm, ...
  'Label',    'Color Faceted', ...
  'Callback', @(hObject, eventData) this.callbacks('setfaceted'));
d.smflat    = uimenu(d.smcolorm, ...
  'Label',    'Color Flat', ...
  'Callback', @(hObject, eventData) this.callbacks('setflat'));
d.sminterp  = uimenu(d.smcolorm, ...
  'Label',    'Color Interp', ...
  'Callback', @(hObject, eventData) this.callbacks('setinterp'));
d.smtex     = uimenu(d.smcolorm, ...
  'Label',    'Color Texture', ...
  'Callback', @(hObject, eventData) this.callbacks('settexture'));
d.smnone    = uimenu(d.smcolorm, ...
  'Label',    'Color None', ...
  'Callback', @(hObject, eventData) this.callbacks('setnone'));
d.smtransm  = uimenu(this.sliceContextMenu, 'Label', 'Transparency');
d.smtnone   = uimenu(d.smtransm, ...
  'Label',    'Transparency None', ...
  'Callback', @(hObject, eventData) this.callbacks('setalphanone'));
d.smtp5     = uimenu(d.smtransm, ...
  'Label',    'Transparency 0.5', ...
  'Callback', @(hObject, eventData) this.callbacks('setalphapoint5'));
d.smtflat   = uimenu(d.smtransm, ...
  'Label',    'Transparency Flat', ...
  'Callback', @(hObject, eventData) this.callbacks('setalphaflat'));
d.smtinterp = uimenu(d.smtransm, ...
  'Label',    'Transparency Interp', ...
  'Callback', @(hObject, eventData) this.callbacks('setalphainterp'));
d.smttex    = uimenu(d.smtransm, ...
  'Label',    'Transparency Texture', ...
  'Callback', @(hObject, eventData) this.callbacks('setalphatexture'));
d.smcontour = uimenu(this.sliceContextMenu, ...
  'Label',    'Add Contour', ...
  'Separator', 'on');
d.smcont0   = uimenu(d.smcontour, ...
  'Label',    'Auto (Slice)', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontour'));
d.smcont0v  = uimenu(d.smcontour, ...
  'Label',    'Auto (Volume)', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourfullauto'));
d.smcont1   = uimenu(d.smcontour, ...
  'Label',    'Select Levels', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontour_select'), ...
  'Separator', 'on');
d.smcsetauto= uimenu(this.sliceContextMenu, ...
  'Label',    'Set Auto Levels (Slice)', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontour_setauto'));
d.smcsetav  = uimenu(this.sliceContextMenu, ...
  'Label',    'Set Auto Levels (Volume)', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontour_setfullauto'));
d.smclevels = uimenu(this.sliceContextMenu, ...
  'Label',    'Set Levels', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontour_setlevels'));
d.smrcontour= uimenu(this.sliceContextMenu, ...
  'Label',    'Remove Contour', ...
  'Callback', @(hObject, eventData) this.callbacks('deleteslicecontour'));
d.smccm     = uimenu(this.sliceContextMenu, 'Label', 'Contour Colors');
d.smcflat   = uimenu(d.smccm, ...
  'Label',    'Contour Flat', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourflat'));
d.smcinterp = uimenu(d.smccm, ...
  'Label',    'Contour Interp', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourinterp'));
d.smcblack  = uimenu(d.smccm, ...
  'Label',    'Contour Black', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourblack'));
d.smcwhite  = uimenu(d.smccm, ...
  'Label',    'Contour White', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourwhite'));
d.smccolor  = uimenu(d.smccm, ...
  'Label',    'Contour Color', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourcolor'));
d.smcsmooth = uimenu(this.sliceContextMenu, ...
  'Visible',  'off', ...
  'Label',    'Smooth Contour Lines', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontoursmooth'));
d.smclinew  = uimenu(this.sliceContextMenu, 'Label', 'Contour Line Width');
d.smcl1     = uimenu(d.smclinew, ...
  'Label',    '1', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourlinewidth', 1));
d.smcl2     = uimenu(d.smclinew, ...
  'Label',    '2', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourlinewidth', 2));
d.smcl3     = uimenu(d.smclinew, ...
  'Label',    '3', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourlinewidth', 3));
d.smcl4     = uimenu(d.smclinew, ...
  'Label',    '4', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourlinewidth', 4));
d.smcl5     = uimenu(d.smclinew, ...
  'Label',    '5', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourlinewidth', 5));
d.smcl6     = uimenu(d.smclinew, ...
  'Label',    '6', ...
  'Callback', @(hObject, eventData) this.callbacks('slicecontourlinewidth', 6));

% Isosurface Context Menu
this.isoContextMenu = uicontextmenu('Callback', {@isocontextmenu, this.hParent});
d.vistogiso = uimenu(this.isoContextMenu, ...
  'Label',    'Visible', ...
  'Callback', @(hObject, eventData) this.callbacks('isotogglevisible'));
d.isodelete = uimenu(this.isoContextMenu, ...
  'Label',    'Delete', ...
  'Callback', @(hObject, eventData) this.callbacks('isodelete'));
d.isoflatlight = uimenu(this.isoContextMenu, ...
  'Label',    'Lighting Flat', ...
  'Callback', @(hObject, eventData) this.callbacks('isoflatlight'), ...
  'Separator', 'on');
d.isosmoothlight = uimenu(this.isoContextMenu, ...
  'Label',    'Lighting Smooth', ...
  'Callback', @(hObject, eventData) this.callbacks('isosmoothlight'));
d.isocolor = uimenu(this.isoContextMenu, ...
  'Label',    'Change Color', ...
  'Callback', @(hObject, eventData) this.callbacks('isocolor'), ...
  'Separator', 'on');
d.isoalpha = uimenu(this.isoContextMenu, 'Label', 'Change Transparency');
d.isoalpha02 = uimenu(d.isoalpha, ...
  'Label',    '0.2', ...
  'Callback', @(hObject, eventData) this.callbacks('isoalpha', .2));
d.isoalpha05 = uimenu(d.isoalpha, ...
  'Label',    '0.5', ...
  'Callback', @(hObject, eventData) this.callbacks('isoalpha', .5));
d.isoalpha08 = uimenu(d.isoalpha, ...
  'Label',    '.8', ...
  'Callback', @(hObject, eventData) this.callbacks('isoalpha', .8));
d.isoalpha10 = uimenu(d.isoalpha, ...
  'Label',    'None', ...
  'Callback', @(hObject, eventData) this.callbacks('isoalpha', 1));
d.isocap = uimenu(this.isoContextMenu, ...
  'Label',    'Add IsoCaps', ...
  'Callback', @(hObject, eventData) this.callbacks('isocaps'), ...
  'Separator', 'on');
d.isoExportStl = uimenu(this.isoContextMenu, ...
  'Label',    'Export to .stl...', ...
  'Callback', @(hObject, eventData) ExportStlLocal(this), ...
  'Separator', 'on');

end


function controlmenu(hFig, eventData, this)
%% Handle doing things to the CONTROLS menu

d = getappdata(this.hParent, 'sliceomatic');

if cameratoolbar('getvisible')
  set(d.camtoolbar, 'Checked', 'on');
else
  set(d.camtoolbar, 'Checked', 'off');
end

if exist('uitoolfactory', 'file') == 2
  if propcheck(d.toolbar, 'Visible', 'on')
    set(d.anntoolbar, 'Checked', 'on');
  else
    set(d.anntoolbar, 'Checked', 'off');
  end
end

set([d.dcalpha1 d.dcalpha8 d.dcalpha6 d.dcalpha5 d.dcalpha6 d.dcalpha2 d.dcalpha0...
  d.dclabels d.dcvis ],...
  'Checked', 'off');

switch get(d.pxx,'facealpha')
  case 1,  set(d.dcalpha1, 'Checked', 'on');
  case .8, set(d.dcalpha8, 'Checked', 'on');
  case .6, set(d.dcalpha6, 'Checked', 'on');
  case .5, set(d.dcalpha5, 'Checked', 'on');
  case .4, set(d.dcalpha4, 'Checked', 'on');
  case .2, set(d.dcalpha2, 'Checked', 'on');
  case 0,  set(d.dcalpha0, 'Checked', 'on');
end

if d.animincrement == 0
  set(d.dcanimstep, 'Checked', 'off');
else
  set(d.dcanimstep, 'Checked', 'on');
end

if ~isempty(get(this.hSliderX, 'XTickLabel'))
  set(d.dclabels, 'Checked', 'on');
end

if strcmp(get(this.hSliderX, 'Visible'), 'on')
  set(d.dcvis, 'Checked', 'on');
end

if 0
  xt = get(get(this.hSliderX, 'Title'), 'String');
  switch xt
    case 'X Slice Controller'
      set(d.dcslice, 'Checked', 'on');
  end

  xt = get(get(this.hSliderIso, 'Title'), 'String');
  switch xt
    case 'Iso Surface Controller'
      set(d.dciso, 'Checked', 'on');
  end
end

end


function defaultmenu(hFig, eventData, hParent)
%% Handle toggling bits on the slice defaults menu

d = getappdata(hParent, 'sliceomatic');

set([d.dfacet d.dflat d.dinterp d.dtex d.dcnone ...
  d.dtnone d.dtflat d.dtinterp d.dttex ...
  d.dlflat d.dlsmooth ...
  d.disot02 d.disot05 d.disot08 d.disot10 ...
  d.dcflat d.dcinterp d.dcblack d.dcwhite ...
  d.smcl1 d.smcl2 d.smcl3 d.smcl4 d.smcl5 d.smcl6 ], 'Checked', 'off');
switch d.defcolor
  case 'faceted'
    set(d.dfacet,'checked','on');
  case 'flat'
    set(d.dflat,'checked','on');
  case 'interp'
    set(d.dinterp,'checked','on');
  case 'texture'
    set(d.dtex,'checked','on');
  case 'none'
    set(d.dcnone,'checked','on');
end
switch d.defalpha
  case 'none'
    set(d.dtnone,'checked','on');
  case 'flat'
    set(d.dtflat,'checked','on');
  case 'interp'
    set(d.dtinterp,'checked','on');
  case 'texture'
    set(d.dttex,'checked','on');
end
switch d.deflight
  case 'flat'
    set(d.dlflat,'checked','on');
  case 'smooth'
    set(d.dlsmooth,'checked','on');
end
switch d.defisoalpha
  case .2, set(d.disot02, 'Checked', 'on');
  case .5, set(d.disot05, 'Checked', 'on');
  case .8, set(d.disot08, 'Checked', 'on');
  case 1., set(d.disot10, 'Checked', 'on');
end
switch d.defcontourcolor
  case 'flat'
    set(d.dcflat,'checked','on');
  case 'interp'
    set(d.dcinterp,'checked','on');
  case 'black'
    set(d.dcblack,'checked','on');
  case 'white'
    set(d.dcwhite,'checked','on');
end
%set(d.dcsmooth,'checked',d.defcontoursmooth);
switch d.defcontourlinewidth
  case 1, set(d.dcl1,'checked','on');
  case 2, set(d.dcl2,'checked','on');
  case 3, set(d.dcl3,'checked','on');
  case 4, set(d.dcl4,'checked','on');
  case 5, set(d.dcl5,'checked','on');
  case 6, set(d.dcl6,'checked','on');
end

end


function slicecontextmenu(hFig, eventData, hParent)
%% Context menu state for slices

d = getappdata(hParent, 'sliceomatic');

[a, s] = getarrowslice();
set([d.smfacet d.smflat d.sminterp d.smtex d.smtnone d.smtp5 ...
  d.smtflat d.smtinterp d.smttex d.smnone d.smcsmooth
  ],'checked','off');
set(d.vistog, 'Checked', get(s,'visible'));

if propcheck(s, 'EdgeColor', [0 0 0])
  set(d.smfacet, 'Checked', 'on');
elseif propcheck(s, 'FaceColor', 'flat')
  set(d.smflat, 'Checked', 'on');
end
if propcheck(s, 'FaceColor', 'interp')
  set(d.sminterp, 'Checked', 'on');
end
if propcheck(s, 'FaceColor', 'texturemap')
  set(d.smtex, 'Checked', 'on');
end
if propcheck(s, 'FaceColor', 'none')
  set(d.smnone, 'Checked', 'on');
end
if propcheck(s, 'FaceAlpha', 1)
  set(d.smtnone, 'Checked', 'on');
end
if propcheck(s, 'FaceAlpha', .5)
  set(d.smtp5, 'Checked', 'on');
end
if propcheck(s, 'FaceAlpha', 'flat')
  set(d.smtflat, 'Checked', 'on');
end
if propcheck(s, 'FaceAlpha', 'interp')
  set(d.smtinterp, 'Checked', 'on');
end
if propcheck(s, 'FaceAlpha', 'texturemap')
  set(d.smttex, 'checked', 'on');
end
cm = [d.smcflat d.smcinterp d.smcblack d.smcwhite d.smccolor ...
  d.smcl1 d.smcl2 d.smcl3 d.smcl4 d.smcl5 d.smcl6 ];
set(cm, 'Checked', 'off');
if isempty(getappdata(s, 'contour'))
  set(d.smcontour,'enable','on');
  set(d.smcsetauto,'enable','off');
  set(d.smcsetav,'enable','off');
  set(d.smclevels,'enable','off');
  set(d.smrcontour,'enable','off');
  set(d.smcsmooth,'enable','off');
  set(cm,'enable','off');
else
  set(d.smcontour,'enable','off')
  set(d.smcsetauto,'enable','on');
  set(d.smcsetav,'enable','on');
  set(d.smclevels,'enable','on');
  set(d.smrcontour,'enable','on')
  set(d.smcsmooth,'enable','on');
  set(cm,'enable','on')
  c = getappdata(s,'contour');
  if propcheck(c,'linesmoothing','on')
    set(d.smcsmooth,'checked','on');
  end
  ec = get(c,'edgecolor');
  if isa(ec, 'char')
    switch ec
      case 'flat'
        set(d.smcflat,'checked','on');
      case 'interp'
        set(d.smcinterp,'checked','on');
    end
  else
    if all(ec == [ 1 1 1 ])
      set(d.smcwhite,'checked','on');
    elseif all(ec == [ 0 0 0 ])
      set(d.smcblack,'checked','on');
    else
      set(d.smccolor,'checked','on');
    end
  end
  clw = get(c, 'Linewidth');
  switch clw
    case 1, set(d.smcl1,'checked','on');
    case 2, set(d.smcl2,'checked','on');
    case 3, set(d.smcl3,'checked','on');
    case 4, set(d.smcl4,'checked','on');
    case 5, set(d.smcl5,'checked','on');
    case 6, set(d.smcl6,'checked','on');
  end
end

end


function isocontextmenu(hAxes, eventData, hParent)
%% Context menu state for isosurfaces

d = getappdata(hParent, 'sliceomatic');

[a, s] = getarrowslice();
if propcheck(s,'facelighting','flat')
  set(d.isoflatlight,'checked','on');
  set(d.isosmoothlight,'checked','off');
else
  set(d.isoflatlight,'checked','off');
  set(d.isosmoothlight,'checked','on');
end
set(d.vistogiso,'checked',get(s,'visible'));
if ~isempty(getappdata(s,'isosurfacecap'))
  set(d.isocap,'checked','on');
else
  set(d.isocap,'checked','off');
end
set([d.isoalpha10 d.isoalpha08 d.isoalpha05 d.isoalpha02], 'Checked', 'off');
isoalpha = get(s, 'FaceAlpha');
switch isoalpha
  case .2, set(d.isoalpha02, 'Checked', 'on');
  case .5, set(d.isoalpha05, 'Checked', 'on');
  case .8, set(d.isoalpha08, 'Checked', 'on');
  case 1., set(d.isoalpha10, 'Checked', 'on');
end

end


function ExportStlLocal(this)
%% Callback in the iso surface context menu

[a, ~] = getarrowslice();
pos = getappdata(a, 'arrowcenter');  % the line the arrow points at is the iso value

this.ExportStl(pos);

end


function SavePrefs(hObject, eventData, this)

%appdata structure knows everything about the implementation
d = getappdata(this.hParent, 'sliceomatic');

%extract only preferences that need to be sticky
prefs.anntoolbar_Checked = get(d.toolbar, 'Visible');
prefs.defcolor = d.defcolor;
prefs.defalpha = d.defalpha;
prefs.deflight = d.deflight;
prefs.defcontourcolor = d.defcontourcolor;
prefs.defcontourlinewidth = d.defcontourlinewidth;
prefs.defcontoursmooth = d.defcontoursmooth;

prefs.camtoolbar_checked = cameratoolbar('getvisible');
prefs.ticklabels = get(this.hSliderX, 'XTickLabelMode');
prefs.animincrement = d.animincrement;
prefs.controlalpha = get(d.pxx, 'FaceAlpha');

%store mini structure (locally where Slice-O-Matic installed)
fileName = UserStickyPrefsFileName;
save(fileName, 'prefs')
disp([ 'Saved: ' fileName])

end


function d = OverrideStickyUserPreferences(this, d)

%characteristic prefs file (stored locally where Slice-O-Matic installed)
fileName = UserStickyPrefsFileName;

%override particular field values (if file exists)
if exist(fileName,'file')
  load(fileName)
  set(d.toolbar,'visible', prefs.anntoolbar_Checked)
  if prefs.camtoolbar_checked
    cameratoolbar('show');
  else
    cameratoolbar('hide');
  end
  d.defcolor = prefs.defcolor;
  d.defalpha = prefs.defalpha;
  d.deflight = prefs.deflight;
  d.defcontourcolor = prefs.defcontourcolor;
  d.defcontourlinewidth = prefs.defcontourlinewidth;
  d.defcontoursmooth = prefs.defcontoursmooth;

  if strcmp('auto', prefs.ticklabels)
    set([this.hSliderX,  this.hSliderIso], 'XTickLabelMode', 'auto');
    set([this.hSliderY, this.hSliderZ], 'YTickLabelMode', 'auto');
  else
    set([this.hSliderX, this.hSliderIso], 'XTickLabel', []);
    set([this.hSliderY, this.hSliderZ], 'YTickLabel', []);
  end
  d.animincrement = prefs.animincrement;
  set([d.pxx d.pxy d.pxz] , 'FaceAlpha', prefs.controlalpha);
  iso = findobj(this.hSliderIso, 'type', 'image');
  set(iso, 'AlphaData', prefs.controlalpha);

  disp('Sticky preferences loaded.')
end

end


function fileName = UserStickyPrefsFileName
localPath = fileparts(which(mfilename));
fileName = fullfile(localPath, 'Sliceomatic.Prefs.mat');

end
