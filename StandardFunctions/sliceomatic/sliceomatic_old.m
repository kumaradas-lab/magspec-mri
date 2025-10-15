function varargout = sliceomatic(varargin)
%% SLICEOMATIC - Slice and isosurface volume exploration GUI
%
% SLICEOMATIC(DATA) - Use 3D double matrix DATA as a volume data.
% SLICEOMATIC(DATA, X, Y, Z) - X,Y,Z are vectors with mesh.
% SLICEOMATIC(HPARENT, ...) - HPARENT is the handle to the parent (can be a
%                             handle to a figure or to a uipanel).
% SLICEOMATIC(HAXES, ...) - HAXES is the handle to the main sliceomatic axes.
%                           The original axes handle will be invalid after the
%                           call to SLICEOMATIC. Use the returned axes handle
%                           instead.
% [HAXES, HPARENT, HFIGURE] = SLICEOMATIC(...)
%
% Example:
%
%       [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
%       v = x .* exp(-x.^2 - y.^2 - z.^2);
%       sliceomatic(v)
%
% Using SLICEOMATIC with no arguments is equivalent to the above
% example.
%
% SLICEOMATIC(DATA, X, Y, Z) - Run sliceomatic using the specified data
% coordinates for the volume DATA. X, Y, and Z are the vectors over which
% DATA is defined.
%
% ex:
%       x = -2:.2:2; y = -2:.25:2; z = -2:.16:2;
%       [X,Y,Z] = meshgrid(x,y,z);
%       v = X .* exp(-X.^2 - Y.^2 - Z.^2);
%       sliceomatic(v,x,y,z)
%
%
% Using the GUI:
% -------------
%
% Create/Delete slices:
%
% The white bars on the top, left, and right allow insertion of new slices on
% the X, Y, and Z planes.  Click in an empty area to add a new slice or surface.
% Right click on a control arrow to reconfigure or delete that slice.
%
% Create/Delete isosurfaces:
%
% The colored bar at the bottom is used to place and position an isosurface.
% The color in the bar indicates a position (as seen in the slice) where the
% isosurface will go.  Right click on a control arrow to reconfigure or delete
% the isosurface.
%
% Orientation of the view:
%
% When the rotate camera button is on, the popup menu will control the camera.
% Turn off camera rotation in order to get individual control over properties of
% the slices and isosurfaces.
%
% Changing Defaults:
%
% The defaults menu provides default features of newly created slices and
% surfaces.  The AllSlices menu controls properties of all the slices and
% surfaces in the scene.  Use popup menus on the objects themselves, or the
% control arrows to change individual properties.
%
% Color & Alpha Maps:
%
% The Colormap popdown controls the currently active colormap.
% This map is used to color the slices.  The Alphamap popdown controls the
% alphamap used on the slices.
%
% Use the color or alpha maps to change how areas of your data are highlighted.
%
% Controls Control:
%
% The Controls menu allows you to adjust how the controls look.  An important
% feature here is the "Animate" item.  You can enable or disable an animation
% when some changes are made.  Since this is decorative, it may be important to
% turn this off for large data sets.
%
% Doing Cool Stuff:
% ----------------
%
% Exploration:
% You can get a quick feel of the current data set by adding a slice using the
% ColorTexture option.  Such a slice can be dragged through the data very
% quickly.
%
% Highlight an Area:
% If certain values in your data are interesting (very large, very small, or
% very median values) you can use transparency to make parts of your slices
% disappear.  Choose AlphaTexture options from the defaults, and sweep some
% slices across your data.  Use the AlphaMap to pick out certain data sets.  The
% example given here looks best with the `vdown' alphamap.
%
% Contours on slices:
% You can add a contour onto a slice to further extract shapes from the data you
% are exploring.  Auto-selecting contour limits will choose contours on a per
% slice basis.  Auto-selecting contour lines from a volume arbitrarily specifies
% 10 levels based on the limits of the volume.
%
% Hidden Shapes:
% Use the isosurface control bar to create an isosurface.  Be patient with this
% control.  It often takes a while for the surface to be created.  Click and
% hold the mouse button down until the first surface appears.  Drag the surface
% through the values until you get something you like, then let go.  If your
% data set is very large, you will need to wait while the new and more accurate
% isosurface is created.
%
% Volumes:
% You can simulate a volume object by creating lots of stacked slices.  Simply
% use the proper Alphamap and transparent textures to highlight the correct
% data, and a large stack of slices will let you simulate a volume object.
%
% Customized Graphics:
% -------------------
%
% To add your own graphics into the sliceomatic display, whatever that may be,
% you can use the following technique:
%
% 1) click on a control arrow
% 2) use gco to get the data for that object
%    slice = getappdata(gco, 'arrowslice')
% 3) use GET to get the cdata and position data which you can use to add your
%    own graphics. 
%
% Setting Default Values:
% ----------------------
%
% If you want to change some default setup feature of sliceomatic, use the "Save
% Preferences" menu item.  Future sliceomatic sessions will then retrieve those
% settings.
%
%
% BUGS:
% ----
%
% 1) Inaccurate Slice
%    Sliceomatic does not use the `slice' command.  All slices are
%    created by explicitly extracting data from the volume.  As such,
%    only slices at integer values are allowed.
%
% 2) Crashes MATLAB
%    Sliceomatic uses the default OpenGL setup.  If you encounter
%    frequent crashes you can start by enabling software OpenGL
%    rendering.  This should always fix the problem, and would
%    likely slow things down too.  You should also update your
%    graphics driver for your video card.  On Windows, in
%    particular, drivers are updated frequently.  For detail on how
%    to overcome these problems, visit this web page:
%  http://www.mathworks.com/support/tech-notes/1200/1201.html
%
% See Also: SLICE, ISOSURFACE, ISOCAPS, CONTOURC, COLORMAP, SURFACE

% This is version 2.3 of sliceomatic.
%
% Sliceomatic is a tool I wrote for fun.  There are no warrenties expressed or
% implied.

% Written by Eric Ludlam <eludlam@mathworks.com>
% Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008
%           The MathWorks Inc
%
% Modified by Emiliano Spezi <emiliano.spezi@physics.org> on 22 May 2003
% Added capability: axes limits control, slice leveling controls,
% and contour level specification.
%
% Patch from David Schwartz on Sep 4, 2008
% Warning fix, and isosurface color fixes.
%
% Patch from Gerhard Stoeckel on Nov 17, 2008
% Fix colormap setting to rand.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2018 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% Check input
if nargin>0 && isscalar(varargin{1}) && ...
    (ishghandle(varargin{1}, 'figure') || ... % check whether first argument is a figure handle
    (varargin{1}>0 && abs(round(varargin{1})-varargin{1})<eps) || ... % or a positive integer
    ishghandle(varargin{1}, 'uipanel') || ... % or a uipanel handle
    ishghandle(varargin{1}, 'axes') || isa(varargin{1}, 'sliceomatic_axes')) % or an axes handle
  if isa(varargin{1}, 'sliceomatic_axes')
    d.axmain_obj = varargin{1};
    d.axmain = d.axmain_obj.hAxes;
    varargin{1} = d.axmain;
  elseif ishghandle(varargin{1}, 'axes')
    d.axmain = varargin{1};
    d.handleParent = get(d.axmain, 'Parent');
  else
    d.handleParent = varargin{1};
  end
  d.handleFigure = ancestor(varargin{1}, 'figure');
  if isempty(d.handleFigure)
    d.handleFigure = varargin{1};
  end
  varargin(1) = [];
end

if numel(varargin)==0
  x = -2:.2:2;
  y = -2:.25:2;
  z = -2:.16:2;
  [XX,YY,ZZ] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
  v = XX .* exp(-XX.^2 - YY.^2 - ZZ.^2);
  if exist('d', 'var')
    [hAxes, hParent, hFigure] = sliceomatic(d.handleParent, v, x, y+3, z);
  else
    [hAxes, hParent, hFigure] = sliceomatic(v, x, y+3, z);
  end
  if nargout > 0
    varargout{1} = hAxes;
  end
  if nargout > 1
    varargout{2} = hParent;
  end
  if nargout > 2
    varargout{3} = hFigure;
  end
  return
end

if ~isnumeric(varargin{1})
  disp('sliceomatic data must be DOUBLE');
  hFigure = [];
  if nargout > 0
    varargout{1} = hFigure;
  end
  if nargout > 1
    varargout{2} = hFigure;
  end
  if nargout > 2
    varargout{3} = hFigure;
  end
  return
end

%% Set-up new sliceomatic
try
  d.data = varargin{1};
  varargin(1) = []; % remove argument from varargin
  haveMesh = false;
  if ~isempty(varargin)
    if all(cellfun(@isnumeric, varargin(1:3)))
      xmesh = varargin{1};
      ymesh = varargin{2};
      zmesh = varargin{3};
      varargin(1:3) = []; % remove arguments from varargin
      haveMesh = true;
      % check if size of mesh fits data
      if ~isvector(xmesh) || length(xmesh) ~= size(d.data, 2)
        error('sliceomatic:SizeMesh', '"X" must be a vector of size(data, 2)');
      elseif iscolumn(xmesh)
        xmesh = xmesh.';
      end
      if ~isvector(ymesh) || length(ymesh) ~= size(d.data, 1)
        error('sliceomatic:SizeMesh', '"Y" must be a vector of size(data, 1)');
      elseif iscolumn(ymesh)
        ymesh = ymesh.';
      end
      if ~isvector(zmesh) || length(zmesh) ~= size(d.data, 3)
        error('sliceomatic:SizeMesh', '"Z" must be a vector of size(data, 3)');
      elseif iscolumn(zmesh)
        zmesh = zmesh.';
      end
    end
    if ~isempty(varargin) && (isnumeric(varargin{1}) || ishghandle(varargin{1}))
      d.handleFigure = varargin{1};
      varargin(1) = []; % remove argument from varargin
      warning('sliceomatic:FigureHandleLast', ...
        'Passing the figure handle as the last argument is deprecated.');
    end
    if ~isempty(varargin)
      % property value pairs
      if ~(numel(varargin)==1 && isstruct(varargin{1})) || (~mod(numel(varargin), 2) && ~all(cellfun(@ischar, varargin(1:2:end))))
        error('sliceomatic:InvalidTrailingInput', ...
          'Trailing arguments to sliceomatic must be a structure or property-value-pairs');
      end
      if (numel(varargin)==1 && isstruct(varargin{1}))
        props = fieldnames(varargin{1});
        vals = struct2cell(varargin{1});
      else
        props = varargin(1:2:end);
        vals = varargin(2:2:end);
      end
      for iPropVal = 1:numel(props)
        d.(lower(props{iPropVal})) = vals{iPropVal};
      end
    end
  end
  % default settings
  if ~isfield(d, 'showtoolbars'), d.showtoolbars = 'on'; end

  if haveMesh
    d = sliceomaticfigure(d, xmesh, ymesh, zmesh);
    d = sliceomaticsetdata(d, xmesh, ymesh, zmesh);
  else
    d = sliceomaticfigure(d);
    d = sliceomaticsetdata(d);
  end
  setappdata(d.handleParent, 'sliceomatic', d);
  hFigure = d.handleFigure;

  set(hFigure, 'CurrentAxes', d.axmain);
catch ME
  warning('sliceomatic:CaughtError', ...
    'The following error was caught while executing sliceomatic. Trying to continue:\n%s', ...
    getReport(ME));
  if ~isfield(d, 'axmain_obj'), d.axmain_obj = []; end
  if ~isfield(d, 'handleParent'), d.handleParent = []; end
  if ~exist('hFigure', 'var'), hFigure = []; end
end

if nargout > 0
  varargout{1} = d.axmain_obj;
end
if nargout > 1
  varargout{2} = d.handleParent;
end
if nargout > 2
  varargout{3} = hFigure;
end

end
