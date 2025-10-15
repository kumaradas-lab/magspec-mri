classdef sliceomatic < handle
  %% Slice and isosurface volume exploration GUI
  % SLICEOMATIC(DATA) - Use 3D double matrix DATA as a volume data.
  % SLICEOMATIC(DATA, X, Y, Z) - X,Y,Z are vectors with mesh.
  % SLICEOMATIC(HPARENT, ...) - HPARENT is the handle to the parent (can be a
  %                             positive integer double valior or ahandle to a
  %                             figure or to a uipanel).
  % SLICEOMATIC(HAXES, ...) - HAXES is the handle to the main sliceomatic axes
  %                           or a sliceomatic object.
  %                           The original axes handle will be invalid after the
  %                           call to SLICEOMATIC. Use the returned axes handle
  %                           instead.
  % [HAXES, HPARENT, HFIGURE] = SLICEOMATIC(...)
  %
  % Example:
  %
  %       [x, y, z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
  %       v = x .* exp(-x.^2 - y.^2 - z.^2);
  %       sliceomatic(v)
  %
  % SLICEOMATIC(DATA, X, Y, Z) - Run sliceomatic using the specified data
  % coordinates for the volume DATA. X, Y, and Z are the vectors over which
  % DATA is defined.
  %
  % Example:
  %       x = -2:.2:2; y = -2:.25:2; z = -2:.16:2;
  %       [X, Y, Z] = meshgrid(x, y, z);
  %       v = X .* exp(-X.^2 - Y.^2 - Z.^2);
  %       sliceomatic(v, x, y+3, z)
  %
  % Using SLICEOMATIC with no arguments is equivalent to the above
  % example.
  %
  %
  % Using the GUI:
  % -------------
  %
  % Create/Delete slices:
  %
  % The white bars on the top, left, and right allow insertion of new slices on
  % the X, Y, and Z planes.  Click in an empty area to add a new slice or
  % surface.  Right click on a control arrow to reconfigure or delete that
  % slice.
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
  % Turn off camera rotation in order to get individual control over properties
  % of the slices and isosurfaces.
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
  % Use the color or alpha maps to change how areas of your data are
  % highlighted.
  %
  % Controls Control:
  %
  % The Controls menu allows you to adjust how the controls look.  An important
  % feature here is the "Animate" item.  You can enable or disable an animation
  % when some changes are made.  Since this is decorative, it may be important
  % to turn this off for large data sets.
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
  % slices across your data.  Use the AlphaMap to pick out certain data sets.
  % The example given here looks best with the `vdown' alphamap.
  %
  % Contours on slices:
  % You can add a contour onto a slice to further extract shapes from the data
  % you are exploring.  Auto-selecting contour limits will choose contours on a
  % per slice basis.  Auto-selecting contour lines from a volume arbitrarily
  % specifies 10 levels based on the limits of the volume.
  %
  % Hidden Shapes:
  % Use the isosurface control bar to create an isosurface.  Be patient with
  % this control.  It often takes a while for the surface to be created.  Click
  % and hold the mouse button down until the first surface appears.  Drag the
  % surface through the values until you get something you like, then let go.
  % If your data set is very large, you will need to wait while the new and more
  % accurate isosurface is created.
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
  % If you want to change some default setup feature of sliceomatic, use the
  % "Save Preferences" menu item.  Future sliceomatic sessions will then
  % retrieve those settings.
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
  %
  % This is version 3.0 of sliceomatic.
  %
  % The original sliceomatic from Matlab File Exchange that was used as base for
  % this tool contained the following disclaimer:
  % Sliceomatic is a tool I wrote for fun.  There are no warrenties expressed or
  % implied.
  %
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
  % ----------------------------------------------------------------------------
  % (C) Copyright 2016-2024 Pure Devices GmbH, Wuerzburg, Germany
  % www.pure-devices.com
  % ----------------------------------------------------------------------------


  properties
    FadeInIso = true  % slowly fade in the transparency of iso surfaces
  end


  properties (SetAccess = private)
    hAxes       % handle to the main axes
    hParent     % handle to the parent of sliceomatic
    hFigure     % handle to the ancestor figure

    hSliderX    % handle to the axes with the X slider
    hSliderY    % handle to the axes with the Y slider
    hSliderZ    % handle to the axes with the Z slider
    hSliderIso  % handle to the axes with the iso-surface slider

    hSliderIsoImage  % handle to the image in the iso-surface slider
  end


  properties (SetAccess = private, GetAccess = private)
    xmesh
    ymesh
    zmesh

    colormapPanel
    alphamapPanel
    orientationPanel

    fileMenu
    controlsMenu
    objectDefaultsMenu
    allSlicesMenu
    helpMenu
    sliceContextMenu
    isoContextMenu

    draggedArrow    % handle to the arrow that is currently dragged

    isSaved = false
  end


  methods

    function [this, hParent, hFigure] = sliceomatic(varargin)
      %% Contructor

      %% Check for handle input
      maybeRetainSlices = false;
      if nargin>0 && isscalar(varargin{1}) && ...
          (ishghandle(varargin{1}, 'figure') || ... % check whether first argument is a figure handle
          (varargin{1}>0 && abs(round(varargin{1})-varargin{1})<eps) || ... % or a positive integer
          ishghandle(varargin{1}, 'uipanel') || ... % or a uipanel handle
          ishghandle(varargin{1}, 'axes') || isa(varargin{1}, 'sliceomatic')) % or an axes handle
        if isa(varargin{1}, 'sliceomatic')
          this = varargin{1};
          varargin{1} = this.hAxes;
          maybeRetainSlices = true;
          % FIXME: Store current slice positions and iso values if it makes
          % sense to add them later to the new data.
        elseif ishghandle(varargin{1}, 'axes')
          this.hAxes = varargin{1};
          this.hParent = get(this.hAxes, 'Parent');
        else
          this.hParent = varargin{1};
        end
        if this.hParent>0 && abs(round(this.hParent)-this.hParent)<eps && ...
          ~ishghandle(this.hParent, 'figure')
          this.hParent = figure(this.hParent);
        end
        d.handleFigure = ancestor(varargin{1}, 'figure');
        if isempty(d.handleFigure)
          d.handleFigure = varargin{1};
        end
        varargin(1) = [];
      end

      %% Display demo
      if numel(varargin)==0
        % show example data
        x = -2:.2:2;
        y = -2:.25:2;
        z = -2:.16:2;
        [XX, YY, ZZ] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
        v = XX .* exp(-XX.^2 - YY.^2 - ZZ.^2);
        if ishghandle(this.hParent)
          this = sliceomatic(this, v, x, y+3, z);
        else
          this = sliceomatic(v, x, y+3, z);
        end
        hParent = this.hParent;
        hFigure = this.hFigure;
        return;
      end

      %% Check input
      if ~isnumeric(varargin{1})
        disp('sliceomatic data must be numeric');  % FIXME: Should this be an error?
        hParent = [];
        hFigure = [];
        return;
      end

      d.data = varargin{1};
      varargin(1) = []; % remove argument from varargin
      xmesh = NaN;
      ymesh = NaN;
      zmesh = NaN;
      if ~isempty(varargin)
        if numel(varargin) > 2 && all(cellfun(@isnumeric, varargin(1:3)))
          xmesh = varargin{1};
          ymesh = varargin{2};
          zmesh = varargin{3};
          varargin(1:3) = []; % remove arguments from varargin

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
        if ~isempty(varargin) ...
            && (isnumeric(varargin{1}) ...
                || (numel(varargin{1})==1 && ishghandle(varargin{1})))
          d.handleFigure = varargin{1};
          varargin(1) = []; % remove argument from varargin
          warning('sliceomatic:FigureHandleLast', ...
            'Passing the figure handle as the last argument is deprecated.');
        end
        if ~isempty(varargin)
          % property value pairs
          if ~(numel(varargin)==1 && isstruct(varargin{1})) && ...
              ~(~mod(numel(varargin), 2) && all(cellfun(@ischar, varargin(1:2:end))))
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

      %% Create sliceomatic GUI and (potentially) restore previous slices
      this.SetupRestoreAxes(d, xmesh, ymesh, zmesh, maybeRetainSlices);
      hFigure = this.hFigure;

      % Call SizeChangedFcn to trigger a layout update
      % (Necessary in newer Matlab versions starting somewhen between
      % Matlab R2012a and Matlab R2020b.)
      this.SizeChangedFcn();

      if nargout > 1
        hParent = this.hParent;
      end
    end


    function delete(this)
      %% Destructor

      % let the DeleteFcn of MainAxes handle the cleanup
      if ishghandle(this.hAxes)
        delete(this.hAxes);
      end
    end


    function sl_new = copyobj(this, newParent)
      %% Copying and re-parenting

      % get sliceomatic data from parent
      hOrigParent = get(this.hAxes, 'Parent');
      d = getappdata(hOrigParent, 'sliceomatic');

      % create new sliceomatic object
      sl_new = sliceomatic(newParent);

      % set references to new objects
      d.axmain = sl_new.hAxes;
      d.handleParent = get(d.axmain, 'Parent');
      d.handleFigure = ancestor(d.axmain, 'figure');
      d.axmain_obj = sl_new;

      % update copy with data from source
      if all(isfinite(this.xmesh)) || all(isfinite(this.ymesh)) || all(isfinite(this.zmesh))
        d = sl_new.SetupAxes(d, this.xmesh, this.ymesh, this.zmesh);
        d = sl_new.sliceomaticsetdata(d, this.xmesh, this.ymesh, this.zmesh);
      else
        d = sl_new.SetupAxes(d);
        d = sl_new.sliceomaticsetdata(d);
      end

      setappdata(d.handleParent, 'sliceomatic', d);

      % copy axes properties from original sliceomatic to new one
      hIso = this.GetSliderIso();
      hIso_new = sl_new.GetSliderIso();
      hX = this.GetSliderX();
      hX_new = sl_new.GetSliderX();
      hY = this.GetSliderY();
      hY_new = sl_new.GetSliderY();
      hZ = this.GetSliderZ();
      hZ_new = sl_new.GetSliderZ();
      % use creator functions for these axes properties
      % (to avoid re-parenting of the original objects)
      axesPropertyObjects = {...
        @xlabel, 'XLabel'; ...
        @ylabel, 'YLabel'; ...
        @zlabel, 'ZLabel'; ...
        @title, 'Title'};  % FIXME: what else?
      for iProp = 1:size(axesPropertyObjects, 1)
        % FIXME: Copy other properties apart from 'String'?
        axesPropertyObjects{iProp, 1}(sl_new.hAxes, ...
          get(get(this.hAxes, axesPropertyObjects{iProp, 2}), 'String'));
        axesPropertyObjects{iProp, 1}(hIso_new, ...
          get(get(hIso, axesPropertyObjects{iProp, 2}), 'String'));
        axesPropertyObjects{iProp, 1}(hX_new, ...
          get(get(hX, axesPropertyObjects{iProp, 2}), 'String'));
        axesPropertyObjects{iProp, 1}(hY_new, ...
          get(get(hY, axesPropertyObjects{iProp, 2}), 'String'));
        axesPropertyObjects{iProp, 1}(hZ_new, ...
          get(get(hZ, axesPropertyObjects{iProp, 2}), 'String'));
      end
      % copy the remainder
      axesProperties = {...
        'XLim', 'YLim', 'ZLim', 'CLim', ...
        'XLimMode', 'YLimMode', 'ZLimMode', 'CLimMode', ...
        'Units', 'Position', 'OuterPosition', 'ActivePositionProperty', ...
        'DataAspectRatio', 'DataAspectRatioMode', ...
        'PlotBoxAspectRatio', 'PlotBoxAspectRatioMode', ...
        'CameraPosition', 'CameraPositionMode', ...
        'CameraTarget', 'CameraTargetMode', ...
        'CameraUpVector', 'CameraUpVectorMode', ...
        'CameraViewAngle', 'CameraViewAngleMode'};  % FIXME: what else?
      for iProp = 1:numel(axesProperties)
        set(sl_new.hAxes, axesProperties{iProp}, get(this.hAxes, axesProperties{iProp}));
        set(hIso_new, axesProperties{iProp}, get(hIso, axesProperties{iProp}));
        set(hX_new, axesProperties{iProp}, get(hX, axesProperties{iProp}));
        set(hY_new, axesProperties{iProp}, get(hY, axesProperties{iProp}));
        set(hZ_new, axesProperties{iProp}, get(hZ, axesProperties{iProp}));
      end

      % copy slices and iso surfaces
      sl_new.AddSliceX(this.GetAllSlicesPosX());
      sl_new.AddSliceY(this.GetAllSlicesPosY());
      sl_new.AddSliceZ(this.GetAllSlicesPosZ());
      sl_new.AddIsoSurface(this.GetAllIsoValues());

      % FIXME: Copy over more settings?
    end


    function this = saveobj(this)
      %% Save object to file

      warn_state = warning('off', 'MATLAB:structOnObject');
      reset_warn_state = onCleanup(@() warning(warn_state));
      stru = struct(this);
      clear reset_warn_state;

      % Set property to select different code path when loading from file.
      stru.isSaved = true;
    end


    function this = loadobj(this)
      %% Load object from file

      warning('Sliceomatic object should be saved as part of the figure');
    end


    function disp(this)
      %% Hide that this is a wrapper object

      if isscalar(this)
        disp(this.hAxes)
      else
        this.hAxes % FIXME: This displays "ans =" twice
      end
    end


    function set(this, varargin)
      %% Forward "set" to hAxes

      if isa(this, 'sliceomatic')
        if mod(numel(varargin), 2) ~= 0
          error('Arguments must be property-value-pairs')
        end
        % get current parent for accessing sliceomatic data
        hOrigParent = get(this.hAxes, 'Parent');
        for iProp = 1:2:numel(varargin)
          if strcmpi(varargin{iProp}, 'CLim')
            clim = varargin{iProp+1};
            % adapt range of ISO slider
            set(this.hSliderIso, 'XLim', clim);
            % adapt size of background image in ISO slider
            isoimg = findobj(this.hSliderIso, 'Type', 'image', ...
              'Tag', 'sliceomaticisocontrolimage');
            set(isoimg, 'XData', clim);
            % change colors of existing ISO surfaces
            d = getappdata(hOrigParent, 'sliceomatic');
            cm_str = get(d.colormapdropdown, 'String');
            cm_fcn_str = cm_str{get(d.colormapdropdown, 'Value')};
            if ~any(strcmp(cm_fcn_str, {'rand', 'custom'})) % FIXME: Add support for those
              cm = feval(cm_fcn_str);
              isosurfs = findobj(this.hAxes, 'Type', 'patch', ...
                'Tag', 'sliceomaticisosurface');
              clen = clim(2) - clim(1);
              for iIso = 1:length(isosurfs)
                value = getappdata(isosurfs(iIso), 'isosurfacevalue');
                idx = fix((value - clim(1))*length(cm)/clen) + 1;
                if idx < 1 || idx > size(cm, 1)
                  % remove iso surface and arrow that are out of range
                  arrow = getappdata(isosurfs(iIso), 'controlarrow');
                  delete(isosurfs(iIso));
                  delete(arrow);
                else
                  % map color
                  set(isosurfs(iIso), 'FaceColor', cm(idx,:));
                end
              end
            end

            % Finally set 'CLim' of the main axes
            set(this.hAxes, varargin{iProp:iProp+1});

          elseif any(strcmpi(varargin{iProp}, {'XLim', 'YLim', 'ZLim', 'XDir', 'YDir', 'ZDir', 'XScale', 'YScale', 'ZScale'}))
            newVal = varargin{iProp+1};

            % adapt orientation of slice slider
            if strcmpi(varargin{iProp}(1), 'X')
              set(this.hSliderX, ['X' varargin{iProp}(2:end)], newVal);
            elseif strcmpi(varargin{iProp}(1), 'Y')
              set(this.hSliderY, ['Y' varargin{iProp}(2:end)], newVal);
            elseif strcmpi(varargin{iProp}(1), 'Z')
              set(this.hSliderZ, ['Y' varargin{iProp}(2:end)], newVal);
            end

            if any(strcmpi(varargin{iProp}, {'XLim', 'YLim', 'ZLim', 'XDir', 'YDir', 'ZDir'}))
              % update value for motion meta-slice
              d = getappdata(hOrigParent, 'sliceomatic');
              d.(lower(varargin{iProp})) = newVal;
              setappdata(hOrigParent, 'sliceomatic', d);
            end

            % Finally set property of the main axes
            set(this.hAxes, varargin{iProp}, newVal);

          elseif any(strcmpi(varargin{iProp}, {'XLimMode', 'YLimMode', 'ZLimMode', 'CLimMode'}))
            if strcmpi(varargin{iProp+1}, 'auto')
              oldLim = get(this.hAxes, varargin{iProp}(1:4));
              set(this.hAxes, varargin{iProp:iProp+1});
              newLim = get(this.hAxes, varargin{iProp}(1:4));
              if ~isequal(oldLim, newLim)
                % Let our [XYZC]Lim update function do its job
                set(this, varargin{iProp}(1:4), newLim);
                % re-set the mode
                set(this.hAxes, varargin{iProp:iProp+1});
              end
            end

          else
            d = getappdata(hOrigParent, 'sliceomatic');
            if isprop(this.hAxes, varargin{iProp})
              set(this.hAxes, varargin{iProp}, varargin{iProp+1});
            elseif isfield(d, lower(varargin{iProp}))
              d.(lower(varargin{iProp})) = varargin{iProp+1};
              setappdata(hOrigParent, 'sliceomatic', d);
            else
              error('PD:sliceomatic:NoSuchProperty', ...
                'No such property "%s" for sliceomatic axes', varargin{iProp});
            end
          end
        end
      else
        if mod(numel(varargin), 2) ~= 0
          error('Arguments must be property-value-pairs')
        end
        for iProp = 1:2:numel(varargin)
          if isa(varargin{iProp+1}, 'sliceomatic')
            set(this, varargin{iProp}, varargin{iProp+1}.hAxes);
          else
            set(this, varargin{iProp}, varargin{iProp+1});
          end
        end
      end
    end


    function val = get(this, prop)
      %% Forward "get" to hAxes

      if nargin < 2
        % FIXME: This doesn't display the "special" sliceomatic properties
        val = get(this.hAxes);
      elseif isprop(this.hAxes, prop)
        val = get(this.hAxes, prop);
      else
        % FIXME: Eventually switch over to the object properties
        d = getappdata(this.hParent, 'sliceomatic');
        if isfield(d, lower(prop))
          val = d.(lower(prop));
        else
          error('PD:sliceomatic:NoSuchProperty', ...
            'No such property "%s" for sliceomatic axes', prop);
        end
      end
    end


    function val = ancestor(this, varargin)
      %% Forward "ancestor" to hAxes

      val = ancestor(this.hAxes, varargin{:});
    end

    function AddSliceX(this, X)
      %% Add slice at position X

      % input check
      if ~isnumeric(X) || ~isreal(X)
        error('PD:sliceomatic:AddSliceX:InvalidInput', ...
          'Input X must be real numbers.');
      end

      d = getappdata(this.hParent, 'sliceomatic');

      for iSlice = 1:numel(X)
        if X(iSlice) < d.xlim(1) || X(iSlice) > d.xlim(2)
          warning('PD:sliceomatic:AddSliceX:OutOfRange', ...
            'New slice at X = %.4g is out of range.', X(iSlice));
          continue;
        end

        % add new arrow and slice
        this.hFigure = ancestor(this.hAxes, 'figure');  % FIXME: Can be removed if nested re-parenting works.
        newa = arrow(this.hSliderX, 'up', X(iSlice));
        set(this.hFigure, 'CurrentAxes', this.hAxes);
        new = this.DrawLocalSlice(d.data, X(iSlice), [], []);
        setappdata(new, 'controlarrow', newa);
        setappdata(newa, 'arrowslice', new);
        set(new, 'AlphaData', get(new, 'CData'), 'AlphaDataMapping', 'scaled');
        set(newa, 'ButtonDownFcn', @(hObject, eventData) this.callbacks('Xmove'));
        set([new, newa], 'UIContextMenu', this.sliceContextMenu);

        this.draggedArrow = newa;
      end

    end


    function AddSliceY(this, Y)
      %% Add slice at position Y

      % input check
      if ~isnumeric(Y) || ~isreal(Y)
        error('PD:sliceomatic:AddSliceY:InvalidInput', ...
          'Input Y must be real numbers.');
      end

      d = getappdata(this.hParent, 'sliceomatic');

      for iSlice = 1:numel(Y)
        if Y(iSlice) < d.ylim(1) || Y(iSlice) > d.ylim(2)
          warning('PD:sliceomatic:AddSliceY:OutOfRange', ...
            'New slice at Y = %.4g is out of range.', Y(iSlice));
          continue;
        end

        % add new arrow and slice
        this.hFigure = ancestor(this.hAxes, 'figure');  % FIXME: Can be removed if nested re-parenting works.
        newa = arrow(this.hSliderY, 'left', Y(iSlice));
        set(this.hFigure, 'CurrentAxes', this.hAxes);
        new = this.DrawLocalSlice(d.data, [], Y(iSlice), []);
        setappdata(new, 'controlarrow', newa);
        setappdata(newa, 'arrowslice', new);
        set(new, 'AlphaData', get(new, 'CData'), 'AlphaDataMapping', 'scaled');
        set(newa, 'ButtonDownFcn', @(hObject, eventData) this.callbacks('Ymove'));
        set([new, newa], 'UIContextMenu', this.sliceContextMenu);

        this.draggedArrow = newa;
      end

    end


    function AddSliceZ(this, Z)
      %% Add slice at position Z

      % input check
      if ~isnumeric(Z) || ~isreal(Z)
        error('PD:sliceomatic:AddSliceZ:InvalidInput', ...
          'Input Z must be real numbers.');
      end

      d = getappdata(this.hParent, 'sliceomatic');

      for iSlice = 1:numel(Z)
        if Z(iSlice) < d.zlim(1) || Z(iSlice) > d.zlim(2)
          warning('PD:sliceomatic:AddSliceZ:OutOfRange', ...
            'New slice at Z = %.4g is out of range.', Z(iSlice));
          continue;
        end

        % add new arrow and slice
        this.hFigure = ancestor(this.hAxes, 'figure');  % FIXME: Can be removed if nested re-parenting works.
        newa = arrow(this.hSliderZ, 'right', Z(iSlice));
        set(this.hFigure, 'CurrentAxes', this.hAxes);
        new = this.DrawLocalSlice(d.data, [], [], Z(iSlice));
        setappdata(new, 'controlarrow', newa);
        setappdata(newa, 'arrowslice', new);
        set(new, 'AlphaData', get(new, 'CData'), 'AlphaDataMapping', 'scaled');
        set(newa, 'ButtonDownFcn', @(hObject, eventData) this.callbacks('Zmove'));
        set([new, newa], 'UIContextMenu', this.sliceContextMenu);

        this.draggedArrow = newa;
      end

    end


    function AddIsoSurface(this, V)
      %% Add iso surface at value V

      % input check
      if ~isnumeric(V) || ~isreal(V)
        error('PD:sliceomatic:AddIsoSurface:InvalidInput', ...
          'Input V must be real numbers.');
      end

      clim = get(this.hSliderIso, 'XLim');
      d = getappdata(this.hParent, 'sliceomatic');

      for iIso = 1:numel(V)
        if V(iIso) < clim(1) || V(iIso) > clim(2)
          warning('PD:sliceomatic:AddIsoSurface:OutOfRange', ...
            'New iso surface at V = %.4g is out of range.', V(iIso));
          continue;
        end

        % add new arrow and iso surface
        this.hFigure = ancestor(this.hAxes, 'figure');  % FIXME: Can be removed if nested re-parenting works.
        newa = arrow(this.hSliderIso, 'down', V(iIso));
        set(this.hFigure, 'CurrentAxes', this.hAxes);
        new = this.DrawLocalIsoSurface(d.reducelims, d.reduce, d.reducesmooth, V(iIso));
        if this.FadeInIso
          slowset(new, 'FaceAlpha', d.defisoalpha);
        else
          set(new, 'FaceAlpha', d.defisoalpha);
        end
        set([newa, new], 'UIContextMenu', this.isoContextMenu);
        setappdata(new, 'controlarrow', newa);
        setappdata(new, 'reduced', 1);
        setappdata(newa, 'arrowiso', new);
        set(newa, 'ButtonDownFcn', @(hObject, eventData) this.callbacks('ISOmove'));

        this.draggedArrow = newa;
      end

    end


    function DeleteSlice(~, hSl)
      %% Delete slices including their controls

      for iSlice = 1:numel(hSl)
        % input check
        if ~ishghandle(hSl(iSlice)) ...
            || (~isappdata(hSl(iSlice), 'controlarrow') && ~isappdata(hSl(iSlice), 'arrowslice'))
          error('PD:sliceomatic:DeleteSlice:InvalidInput', ...
            'Input must be a handle to a slice or the corresponding arrow.');
        end

        [a, s] = getarrowslice(hSl(iSlice));
        if ~strcmp(get(s, 'Type'), 'surface')
          error('PD:sliceomatic:DeleteSlice:InvalidType', ...
            'Input must be a handle to a slice or the corresponding arrow.');
        end

        if ~isempty(getappdata(s, 'contour'))
          delete(getappdata(s, 'contour'));
        end
        delete(s);
        delete(a);
      end
    end


    function DeleteIso(~, hIso)
      %% Delete iso surfaces including their controls

      for iIso = 1:numel(hIso)
        % input check
        if ~ishghandle(hIso(iIso)) ...
            || (~isappdata(hIso(iIso), 'controlarrow') && ~isappdata(hIso(iIso), 'arrowiso'))
          error('PD:sliceomatic:DeleteSlice:InvalidInput', ...
            'Input must be a handle to an iso surface or the corresponding arrow.');
        end

        [a, s] = getarrowslice(hIso(iIso));
        if ~strcmp(get(s, 'Type'), 'patch')
          error('PD:sliceomatic:DeleteSlice:InvalidType', ...
            'Input must be a handle to an iso surface or the corresponding arrow.');
        end

        hArrows = getappdata(hIso(iIso), 'Arrows');
        hArrows(hArrows == a) = [];
        setappdata(hIso(iIso), 'Arrows', hArrows);
        cap = getappdata(s, 'sliceomaticisocap');
        if ~isempty(cap)
          delete(cap);
        end
        delete(s);
        delete(a);
      end
    end


    function UpdateData(this, data)
      %% Update data in existing sliceomatic
      %
      %   this.UpdateData(data)
      %
      % This method replaces the underlying data in an existing sliceomatic
      % without re-creating the axes. Only the axes content is re-drawn. That
      % can avoid flickering.

      % input check
      if ~isnumeric(data) || ~isreal(data)
        error('PD:sliceomatic:UpdateData:InvalidInput', ...
          'Input data must be real numbers.');
      end

      d = getappdata(this.hParent, 'sliceomatic');
      if ~isequal(size(d.data), size(data))
        % FIXME: Can we work around that?
        error('PD:sliceomatic:UpdateData:InvalidSize', ...
          'Size of updated data must match existing data.');
      end

      % FIXME: Should we also allow to (optionally) update the x-y-z meshes?

      % retain positions of existing slices and iso surface
      oldX = this.GetAllSlicesPosX();
      oldY = this.GetAllSlicesPosY();
      oldZ = this.GetAllSlicesPosZ();
      oldIso = this.GetAllIsoValues();

      % remove existing slices and iso surfaces
      this.DeleteSlice(this.GetAllSlices());
      this.DeleteIso(this.GetAllIsos());

      % update data in structure
      d.data = data;

      % adapt limits
      lim = [min(d.data(isfinite(d.data))), max(d.data(isfinite(d.data)))];
      if isempty(lim)
        lim = [-1, 1];
      end
      if lim(1) == lim(2)
        if lim(1) > 0
          lim(1) = 0.9*lim(1);
          lim(2) = 1.1*lim(2);
        else
          lim(1) = 1.1*lim(1);
          lim(2) = 0.9*lim(2);
        end
      end
      set(this.hAxes, 'CLim', lim, 'ALim', lim);
      if ishghandle(this.hSliderIso)
        set(this.hSliderIso, 'XLim', lim);
        set(this.hSliderIsoImage, 'XData', lim);
      end

      % update data in sliceomatic
      d = this.sliceomaticsetdata(d, this.xmesh, this.ymesh, this.zmesh);

      % store updated data in figure
      setappdata(this.hParent, 'sliceomatic', d);

      % add slices and iso surfaces at previous positions
      this.AddSliceX(oldX);
      this.AddSliceY(oldY);
      this.AddSliceZ(oldZ);
      this.AddIsoSurface(oldIso);
    end


    function ExportStl(this, V, stlFileStr, zoomFactor)
      %% Export iso surface with value V to .stl file (3d printer)
      %
      %   hsl.ExportStl(V, stlFileStr, zoomFactor)
      %
      % INPUT:
      %
      %   V
      %       scalar value for the iso surface.
      %
      %   stlFileStr
      %       string with the file name for the .stl file. If omitted or empty,
      %       a dialog is shown to select a file name.
      %
      %   zoomFactor
      %       scalar value with the scaling factor from the units on the axes to
      %       the exported .stl. If omitted or empty, a scaling factor of 1 is
      %       used. Depending on the used program for 3d printing, .stl files
      %       are imported in mm, cm or inches.

      if ~exist('stlwrite', 'file')
        error('PD:sliceomatic:ExportStl:NoStlWrite', ...
          'Function "stlwrite" does not exist.');
      end

      d = getappdata(this.hParent, 'sliceomatic');
      hp = this.DrawLocalIsoSurface(d.reducelims, d.reduce, d.reducesmooth, V);
      set(hp, 'Visible', 'off');

      faces = get(hp, 'Faces');
      vertices = get(hp, 'Vertices');
      delete(hp);

      if nargin < 3 || isempty(stlFileStr)
        docDir = winqueryreg('HKEY_CURRENT_USER', ...
          'Software\Microsoft\Windows\CurrentVersion\Explorer\User Shell Folders', ...
          'Personal');
        [stlFileStr, stlPathStr] = uiputfile('*.stl', 'Export to .stl', ...
          fullfile(docDir, 'iso.stl'));
        if isnumeric(stlFileStr) && stlFileStr == 0
          return;
        end
      end

      if nargin < 4 || isempty(zoomFactor)
        zoomFactor = inputdlg('Scaling factor (from values on axes) for export', ...
          'Zoom factor', 1, {'1'});
        if isempty(zoomFactor)
          return;
        end
        zoomFactor = str2double(zoomFactor{1});
        if ~isnumeric(zoomFactor)
          zoomFactor = 1;
        end
      end

      stlwrite(fullfile(stlPathStr, stlFileStr), faces, vertices*zoomFactor, ...
        'Title', sprintf('%s: Exported from sliceomatic (V=%s)', datestr(now), num2str(V)));

    end


    function posX = GetAllSlicesPosX(this)
      %% Get x value of all x slices currently in the axes

      posX = this.GetAllSlicesPosType('X');
    end


    function posY = GetAllSlicesPosY(this)
      %% Get y value of all y slices currently in the axes

      posY = this.GetAllSlicesPosType('Y');
    end


    function posZ = GetAllSlicesPosZ(this)
      %% Get z value of all z slices currently in the axes

      posZ = this.GetAllSlicesPosType('Z');
    end


    function isoValues = GetAllIsoValues(this)
      %% Get values of all iso-surfaces currently in the axes

      hIsos = this.GetAllIsos();
      isoValues = arrayfun(@(x) getappdata(x, 'isosurfacevalue'), hIsos, 'UniformOutput', true);
    end


    function hSliderX = GetSliderX(this)
      %% Get handle to the slider for the x-axis

      hSliderX = this.hSliderX;
    end


    function hSliderY = GetSliderY(this)
      %% Get handle to the slider for the y-axis

      hSliderY = this.hSliderY;
    end


    function hSliderZ = GetSliderZ(this)
      %% Get handle to the slider for the z-axis

      hSliderZ = this.hSliderZ;
    end


    function hSliderIso = GetSliderIso(this)
      %% Get handle to the slider for the iso-surfaces

      hSliderIso = this.hSliderIso;
    end


    function CopyToClipboard(this)
      %% Copy the content of the sliceomatic parent to the clipboard

      % show waitbar
      hwb = waitbar(0, 'Copying graphics to clipboard...');
      hwb_protect = onCleanup(@() close(hwb));

      % make progress indeterminate
      hwbKids = allchild(hwb);
      jProgressBar = hwbKids(ishghandle(hwbKids, 'hgjavacomponent')).JavaPeer;
      jProgressBar.setIndeterminate(1);

      % Keep same aspect ratio as currently displayed on screen
      hf = figure('Visible', 'off', ...
        'Units', get(this.hParent, 'Units'), 'Position', get(this.hParent, 'Position'));
      hf_protect = onCleanup(@() close(hf));

      % copy main axes and sliders
      copyobj([this.hAxes, this.hSliderX, this.hSliderY, this.hSliderZ, this.hSliderIso], hf);

      % print to clipboard
      print(hf, '-dmeta');
    end

  end


  methods (Access = private)

    % Create the figure window to be used by the sliceomatic GUI
    d = SetupAxes(this, d, xmesh, ymesh, zmesh);

    % Create the data used for sliceomatic in the appdata D
    d = sliceomaticsetdata(this, d, xmesh, ymesh, zmesh);


    function d = SetupRestoreAxes(this, d, xmesh, ymesh, zmesh, maybeRetainSlices)
      %% Create sliceomatic GUI and (potentially) restore previous slices

      %% get old slice position and iso values
      if maybeRetainSlices
        %  If the meshes match, store the current slice positions and iso values
        %  and try to restore them later with the updated data.
        d_old = getappdata(this.hParent, 'sliceomatic');
        if ~isequaln(xmesh, this.xmesh) || ~isequaln(ymesh, this.ymesh) || ~isequaln(zmesh, this.zmesh) ...
            || isempty(d_old) || ~isstruct(d_old) || ~isfield(d_old, 'axmain') ...
            || ~isequal(size(d.data), size(d_old.data))
          maybeRetainSlices = false;
        else
          oldPosX = this.GetAllSlicesPosX();
          oldPosY = this.GetAllSlicesPosY();
          oldPosZ = this.GetAllSlicesPosZ();
          oldIsoVals = this.GetAllIsoValues();
          if isvalid(d_old.axmain)
            oldAxes = d_old.axmain;
            oldCameraPosition = get(oldAxes, 'CameraPosition');
            oldCameraTarget = get(oldAxes, 'CameraTarget');
            oldCameraUpVector = get(oldAxes, 'CameraUpVector');
            oldCameraViewAngle = get(oldAxes, 'CameraViewAngle');
          else
            oldAxes = [];
          end
        end
      end

      %% Set-up new sliceomatic
      try
        d = this.SetupAxes(d, xmesh, ymesh, zmesh);
        d = this.sliceomaticsetdata(d, xmesh, ymesh, zmesh);
        d.axmain_obj = this;
        setappdata(this.hParent, 'sliceomatic', d);
        this.hFigure = d.handleFigure;

        if maybeRetainSlices
          if ~isempty(oldAxes)
            set(this.hAxes, 'CameraPosition', oldCameraPosition);
            set(this.hAxes, 'CameraTarget', oldCameraTarget);
            set(this.hAxes, 'CameraUpVector', oldCameraUpVector);
            set(this.hAxes, 'CameraViewAngle', oldCameraViewAngle);
          end
          this.AddSliceX(oldPosX);
          this.AddSliceY(oldPosY);
          this.AddSliceZ(oldPosZ);
          this.AddIsoSurface(oldIsoVals);
        end

        set(this.hFigure, 'CurrentAxes', this.hAxes);
      catch ME
        warning('sliceomatic:CaughtError', ...
          'The following error was caught while executing sliceomatic. Trying to continue:\n%s', ...
          getReport(ME));
        if ~exist('hFigure', 'var'), this.hFigure = []; end
      end
    end


    % Adapt positions of axes in parent object
    SizeChangedFcn(this);

    % Executes on deletion of main axes.
    AxesDeleteFcn(this, hAxes, eventData);

    % Executes on deletion of parent figure.
    FigureDeleteFcn(this, hFigure, eventData);

    % Set up sliceomatic's GUI menues
    d = figmenus(this, d);

    % set up slider for iso controls
    SetupIsoControls(this, onoff);

    % Callback functions for sliceomatic
    callbacks(this, varargin);

    % Slice management
    s = DrawLocalSlice(this, data, X, Y, Z, oldslice);

    % Create a contour on slice
    DrawLocalContour(this, slice, oldcontour, levels);


    function slicePos = GetSlicePos(this, coor, type)
      %% Get position of slice from cursor in slider

      meshField = [type, 'mesh'];

      if isnan(this.(meshField))
        slicePos = round(coor);
      else
        diffCoor = abs(this.(meshField) - coor);
        slice_number = find(diffCoor == min(diffCoor) & ...
          diffCoor <= abs(this.(meshField)(1) - this.(meshField)(2))/2, 1);
        if isequal(get(this, 'zdir'), 'reverse')
          slice_number = length(this.(meshField)) - slice_number + 1;
        end
        slicePos = this.(meshField)(slice_number);
      end

    end


    % Draw or refresh iso surface
    hIso = DrawLocalIsoSurface(this, volume, data, datanormals, value, oldiso);

    % Draw or refresh caps for iso surfaces
    hIsoCap = DrawLocalIsoCaps(this, hIsoSurface, oldIsoCap);

    % Prepare graphics objects for dragging an arrow
    ArrowDragPrepare(this);

    % Finish dragging an arrow. Finalize graphics objects
    ArrowDragFinish(this);


    function hSlices = GetAllSlices(this)
      %% Get an array with the handles to all slices

      hSlices = findobj(this.hParent, 'Type', 'surface', 'Tag', 'sliceomaticslice');
    end


    function hIsos = GetAllIsos(this)
      %% Get an array with the handles to all iso surfaces

      hIsos = findobj(this.hParent, 'Type', 'patch', 'Tag', 'sliceomaticisosurface');
    end


    function hCaps = GetAllCaps(this)
      %% Get an array with the handles to all iso surface caps

      % FIXME: This function seems to be unused.
      hCaps = findobj(this.hParent, 'Type', 'patch', 'Tag', 'sliceomaticisocap');
    end


    function pos = GetAllSlicesPosType(this, type)
      %% Get the position of all slices of given type

      hSlices = this.GetAllSlices();
      slicetype = arrayfun(@(x) getappdata(x, 'slicetype'), hSlices, 'UniformOutput', false);
      hSlices = hSlices(strcmp(slicetype, type));
      if isempty(hSlices)
        pos = [];
        return;
      end
      data = get(hSlices, [type, 'Data']);
      if iscell(data)
        pos = cellfun(@(x) x(1,1), data, 'UniformOutput', true);
      else
        pos = data(1,1);
      end
    end

  end

end
