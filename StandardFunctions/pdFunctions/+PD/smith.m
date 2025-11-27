classdef smith < handle
%% Draw a blank Z- or Y-Smith chart
%
%     hsm = PD.smith(ca, 'prop', value, ...);
%
% INPUT:
%   ca      handle to an axes object in which to draw the Smith chart. If empty,
%           the Smith chart is drawn in the current axes.
%   property-value-pairs
%           Pairs of properties and values for the Smith chart graphics
%           properties (see below).
%
% OUTPUT:
%   hsm     handle to a Smith chart object mimicking an hghandle, i.e. you can
%           use the set (and get) function with property-value-pairs to change
%           (or query) the Smith chart graphics properties or set (and get) them
%           directly. Some functions that work on hghandles might not work or
%           work incorrectly on a Smith chart handle.
%
% Smith chart graphics properties:
%   General properties
%     Parent            Handle to the parent of the smitch chart, e.g. an axes
%                       or hggroup.
%     Type              Type of smith chart ('z', 'y', 'yz'), default: 'z'
%                       If 'z' or 'y', the properties for the major Smith chart
%                       are used. If 'yz', the properties for the major Smith
%                       chart are used for the Z-Smith chart, and the properties
%                       for the minor Smith chart are used for the Y-Smith
%                       chart.
%   Axis label properties
%     LabelColor        Color of the axis labels (any valid Matlab color value),
%                       default: [0 0 0]
%     LabelSize         Fontsize of the axis labels, default: 10
%     LabelVisible      Visibility of the axis labels ('on' or 'off'),
%                       default: 'on'
%   Properties for major Z- or Y-Smith chart
%     Color             Color of the major grid lines, default: [0.4 0.4 0.4]
%     LineStyle         LineStyle of the major grid lines, default: '-'
%     LineWidth         LineWidth of the major grid lines, default: 0.5
%     MinorColor        Color of the minor grid lines, default: [0.7 0.7 0.7]
%     MinorLineStyle    LineStyle of the minor grid lines, default: '--'
%     MinorLineWidth    LineWidth of the minor grid lines, default: 0.5
%     MinorVisible      Visibility of the minor grid lines, default: 'off'
%     Visible           Visibility of the major grid lines, default: 'on'
%     Value             Values of the major grid lines. The first row is the
%                       actual positions. The corresponding value in the second
%                       row indicates the start of that grid line. Additionally,
%                       a closed line at max(Value(2,:)) is drawn. Default:
%                       [ 0.2 0.5 1.0 2.0  5.0; ...
%                         1.0 2.0 5.0 5.0 30.0]
%     MinorValue        Values of the minor grid lines. The first row is the
%                       actual positions. The corresponding value in the second
%                       row indicates the start of that grid line. Additionally,
%                       a closed line at max(MinorValue(2,:)) is drawn.
%                       Default:
%                       [ 0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.75  1.0 1.5  2.0  3.5  5.0 10.0; ...
%                         1.0  2.0 1.0  1.0 5.0 1.0 2.0 5.0  10.0 5.0 10.0 30.0 30.0 30.0]
%   Properties for minor Y-Smith chart (only applicable if 'Type' is 'yz')
%     SubColor            Color of the major grid lines, default: [0.8 0.8 0.8]
%     SubLineStyle        LineStyle of the major grid lines, default: ':'
%     SubLineWidth        LineWidth of the major grid lines, default: 0.5
%     SubMinorColor       Color of the minor grid lines, default: [0.9 0.9 0.9]
%     SubMinorLineStyle   LineStyle of the minor grid lines, default: '-.'
%     SubMinorLineWidth   LineWidth of the minor grid lines, default: 0.5
%     SubMinorVisible     Visibility of the minor grid lines, default: 'off'
%     SubVisible          Visibility of the major grid lines, default: 'on'
%     SubValue            Values of the major grid lines. The first row is the
%                         actual positions. The corresponding value in the
%                         second row indicates the start of that grid line.
%                         Additionally, a closed line at max(SubValue(2,:)) is
%                         drawn. Only used if 'SubValueMode' is 'manual'.
%                         Default:
%                         [ 0.2 0.5 1.0 2.0  5.0; ...
%                           1.0 2.0 5.0 5.0 30.0]
%     SubValueMode        Mode ('auto' or 'manual'). If 'auto', the same values
%                         are used in the minor Smith chart like in the major
%                         Smith chart. If 'manual', the values in 'SubValue' are
%                         used. Default: 'auto'
%     SubMinorValue       Values of the minor grid lines. The first row is the
%                         actual positions. The corresponding value in the
%                         second row indicates the start of that grid line.
%                         Additionally, a closed line at max(SubMinorValue(2,:))
%                         is drawn. Only used if 'SubMinorValueMode' is
%                         'manual'. Default:
%                         [ 0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.75  1.0 1.5  2.0  3.5  5.0 10.0; ...
%                           1.0  2.0 1.0  1.0 5.0 1.0 2.0 5.0  10.0 5.0 10.0 30.0 30.0 30.0]
%     SubMinorValueMode   Mode ('auto' or 'manual'). If 'auto', the same values
%                         are used in the minor Smith chart like in the major
%                         Smith chart. If 'manual', the values in
%                         'SubMinorValue' are used. Default: 'auto'
%
% Methods:
%     draw()              Re-draw the grid with the current properties.
%     clear()             Clear all children of smith chart and reset color
%                         order and line style order.
%     clearSelf()         Clear all grid lines and labels.
%     ish = ishghandle(hsm, type)
%                         Overloaded method. See ishghandle.
%     smithValue = createGrid(gridTicks)    % (Static)
%                         Create 'Value' array from groups of tick values.
%                         INPUT:
%                           gridTicks     cell array with groups of grid ticks.
%                         OUTPUT:
%                           smithValue    2×N array which can be used to set the
%                                         '*Value' properties. The lines of each
%                                         group end at the maximum value in each
%                                         group.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2017-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

  properties
    % general properties
    Parent
    Type              = 'z'
    % label properties
    LabelColor        = [0 0 0]
    LabelSize         = 10
    LabelVisible      = 'on'
    % Z-smith properties
    Color             = [0.4 0.4 0.4]
    LineStyle         = '-'
    LineWidth         = 0.5
    MinorColor        = [0.7 0.7 0.7]
    MinorLineStyle    = '--'
    MinorLineWidth    = 0.5
    MinorVisible      = 'off';
    Visible           = 'on';
    Value             = [ 0.2 0.5 1.0 2.0  5.0; ...
                          1.0 2.0 5.0 5.0 30.0]
    MinorValue        = [ 0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.75  1.0 1.5  2.0  3.5  5.0 10.0; ...
                          1.0  2.0 1.0  1.0 5.0 1.0 2.0 5.0  10.0 5.0 10.0 30.0 30.0 30.0]
    % Y-smith properties
    SubColor          = [0.8 0.8 0.8]
    SubLineStyle      = ':'
    SubLineWidth      = 0.5
    SubMinorColor     = [0.9 0.9 0.9]
    SubMinorLineStyle = '-.'
    SubMinorLineWidth = 0.5
    SubMinorVisible   = 'off';
    SubVisible        = 'on';
    SubValue          = []
    SubValueMode      = 'auto'
    SubMinorValue     = []
    SubMinorValueMode = 'auto'
  end


  properties (SetAccess = protected)
    Name              = 'Smith chart'
  end


  properties (Access = protected)
    holdState
    axesListener

    hGroup

    hBackground
    hBaseLines
    hRLabelsZ
    hILabelsZ
    hRLinesZ
    hILinesZ
    hRMinorLinesZ
    hIMinorLinesZ
    hRLabelsY
    hILabelsY
    hRLinesY
    hILinesY
    hRMinorLinesY
    hIMinorLinesY

    doDraw = false
    drawPending = false
    recCount = 0
  end


  methods

    function this = smith(varargin)
      if nargin >= 1 && ishghandle(varargin{1})
        ca = varargin{1};
        varargin(1) = [];
      else
        ca = gca;
      end
      this.Parent = ca;

      %% parse input
      p = inputParser;
      myProps = properties(this);
      if verLessThan('matlab', '8.5.0')
        for iProp = 1:length(myProps)
          p.addParamValue(myProps{iProp}, this.(myProps{iProp}));
        end
      else
        for iProp = 1:length(myProps)
          p.addParameter(myProps{iProp}, this.(myProps{iProp}));
        end
      end
      p.parse(varargin{:});

      fNames = fieldnames(p.Results);
      for iFName = 1:length(fNames)
        this.(fNames{iFName}) = p.Results.(fNames{iFName});
      end

      if strcmp(this.SubValueMode, 'auto')
        this.SubValue = this.Value;
      end
      if strcmp(this.SubMinorValueMode, 'auto')
        this.SubMinorValue = this.MinorValue;
      end

      %% draw smith chart
      this.doDraw = true;
      this.draw();
    end


    function clearSelf(this)
      % clear all lines and labels
      if ~isempty(this.hGroup) && ishghandle(this.hGroup)
        set(this.hGroup, 'DeleteFcn', []);
        delete(this.hGroup(ishghandle(this.hGroup)))
        this.hGroup=[];
      end
      delete(this.hBackground(ishghandle(this.hBackground)))
      this.hBackground=[];
      delete(this.hBaseLines(ishghandle(this.hBaseLines)))
      this.hBaseLines=[];
      delete(this.hRLabelsZ(ishghandle(this.hRLabelsZ)))
      this.hRLabelsZ=[];
      delete(this.hILabelsZ(ishghandle(this.hILabelsZ)))
      this.hILabelsZ=[];
      delete(this.hRLinesZ(ishghandle(this.hRLinesZ)))
      this.hRLinesZ=[];
      delete(this.hILinesZ(ishghandle(this.hILinesZ)))
      this.hILinesZ=[];
      delete(this.hRMinorLinesZ(ishghandle(this.hRMinorLinesZ)))
      this.hRMinorLinesZ=[];
      delete(this.hIMinorLinesZ(ishghandle(this.hIMinorLinesZ)))
      this.hIMinorLinesZ=[];
      delete(this.hRLabelsY(ishghandle(this.hRLabelsY)))
      this.hRLabelsY=[];
      delete(this.hILabelsY(ishghandle(this.hILabelsY)))
      this.hILabelsY=[];
      delete(this.hRLinesY(ishghandle(this.hRLinesY)))
      this.hRLinesY=[];
      delete(this.hILinesY(ishghandle(this.hILinesY)))
      this.hILinesY=[];
      delete(this.hRMinorLinesY(ishghandle(this.hRMinorLinesY)))
      this.hRMinorLinesY=[];
      delete(this.hIMinorLinesY(ishghandle(this.hIMinorLinesY)))
      this.hIMinorLinesY=[];
    end


    function clear(this)
      kids = get(this.Parent, 'Children');
      hNonSmith = findall(kids, 'flat', '-not', 'Tag', 'smith');
      delete(hNonSmith);
      % reset ColorOrder and LineStyleOrder
      ha = ancestor(this, 'axes');
      set(ha, 'ColorOrderIndex', 1);
      set(ha, 'LineStyleOrderIndex', 1);
    end

    function ish = ishghandle(this, type)
      ish = (nargin < 2 || strcmpi(type, 'smith')) && ishghandle(this.hGroup);
    end


    function draw(this)
      if ~this.doDraw
        this.drawPending = true;
        return
      end
      this.clearSelf();

      ha = ancestor(this.Parent, 'axes');
      axis(ha, 'equal', 'off')
      % hggroup
      this.hGroup = hggroup('Parent', this.Parent, 'Tag', 'smith');
      setappdata(this.hGroup, 'smith_object', this);

      % find and clear other smith objects in parent
      kids = get(this.Parent, 'Children');
      hSmith = findall(kids, 'flat', 'Tag', 'smith');
      for iSmith = 1:length(hSmith)
        if hSmith(iSmith) ~= this.hGroup
          otherSmith = getappdata(hSmith(iSmith), 'smith_object');
          if isa(otherSmith, class(this))
            otherSmith.clearSelf();
          end
          if ishghandle(hSmith(iSmith))
            delete(hSmith(iSmith));
          end
        end
      end

      %% Grid
      % background
      z = exp(1i*2*pi*(0:.01:1));
      this.hBackground = patch(real(z), imag(z), 'w', ...
        'FaceColor', ones(1,3)*0.999, 'EdgeColor', 'none', 'Parent', this.hGroup);

      % draw potential sub-grid below main grid
      if any(lower(this.Type) == 'y')
        %% draw minor grid (below major grid)
        if strcmpi(this.Type, 'yz')
          v = this.SubMinorValue;
        else
          v = this.MinorValue;
        end
        % constant real part
        for iLine = 1:size(v, 2)
          z = v(1,iLine) + ...
            1i*[-linspace(-1, 0, 100).^2, linspace(0, 1, 100).^2]*v(2,iLine);
          Z = (z-1)./(z+1);
          this.hRMinorLinesY(iLine) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData', imag(Z));
        end
        z = max(v(2,:)) + ...
          1i*[-linspace(-1, 0, 100).^2, linspace(0, 1, 100).^2]*500;
        Z = (z-1)./(z+1);
        this.hRMinorLinesY(end+1) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData', imag(Z));

        % constant imaginary part
        for iLine = 1:size(v, 2)
          z = linspace(0, 1, 100).^2*v(2,iLine) - 1i*v(1,iLine);
          Z = (z-1)./(z+1);
          this.hIMinorLinesY(2*iLine-1) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData', imag(Z));
          this.hIMinorLinesY(2*iLine) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData', -imag(Z));
        end
        z = linspace(0, 1, 100).^2*500 - 1i*max(v(2,:));
        Z = (z-1)./(z+1);
        this.hIMinorLinesY(end+1) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData', imag(Z));
        this.hIMinorLinesY(end+1) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData', -imag(Z));
        if strcmpi(this.Type, 'yz')
          set([this.hRMinorLinesY, this.hIMinorLinesY], ...
              'Color', this.SubMinorColor, ...
              'LineStyle', this.SubMinorLineStyle, ...
              'LineWidth', this.SubMinorLineWidth, ...
              'Visible', this.SubMinorVisible);
        else
          set([this.hRMinorLinesY, this.hIMinorLinesY], ...
              'Color', this.MinorColor, ...
              'LineStyle', this.MinorLineStyle, ...
              'LineWidth', this.MinorLineWidth, ...
              'Visible', this.MinorVisible);
        end

        %% draw major grid
        if strcmpi(this.Type, 'yz')
          v = this.SubValue;
        else
          v = this.Value;
        end
        % constant real part
        this.hRLinesY = [];
        for iLine = 1:size(v, 2)
          z = v(1,iLine) + ...
            1i*[-linspace(-1, 0, 100).^2, linspace(0, 1, 100).^2]*v(2,iLine);
          Z = (z-1)./(z+1);
          this.hRLinesY(iLine) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData', imag(Z));
        end
        z = max(v(2,:)) + ...
          1i*[-linspace(-1, 0, 100).^2, linspace(0, 1, 100).^2]*500;
        Z = (z-1)./(z+1);
        this.hRLinesY(end+1) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData', imag(Z));

        % constant imaginary part
        for iLine = 1:size(v, 2)
          z = linspace(0, 1, 100).^2*v(2,iLine) - 1i*v(1,iLine);
          Z = (z-1)./(z+1);
          this.hILinesY(2*iLine-1) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData', imag(Z));
          this.hILinesY(2*iLine) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData',- imag(Z));
        end
        z = linspace(0, 1, 100).^2*500 - 1i*max(v(2,:));
        Z = (z-1)./(z+1);
        this.hILinesY(end+1) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData', imag(Z));
        this.hILinesY(end+1) = line('Parent', this.hGroup, 'XData', -real(Z), 'YData',- imag(Z));
        if strcmpi(this.Type, 'yz')
          set([this.hRLinesY, this.hILinesY], ...
              'Color', this.SubColor, ...
              'LineStyle', this.SubLineStyle, ...
              'LineWidth', this.SubLineWidth, ...
              'Visible', this.Visible);
        else
          set([this.hRLinesY, this.hILinesY], ...
              'Color', this.Color, ...
              'LineStyle', this.LineStyle, ...
              'LineWidth', this.LineWidth, ...
              'Visible', this.Visible);
        end
      end

      % draw main grid above potential sub-grid
      if any(lower(this.Type) == 'z')
        %% draw minor grid below major grid
        % constant real part
        for iLine = 1:size(this.MinorValue, 2)
          z = this.MinorValue(1,iLine) + ...
            1i*[-linspace(-1, 0, 100).^2, linspace(0, 1, 100).^2]*this.MinorValue(2,iLine);
          Z = (z-1)./(z+1);
          this.hRMinorLinesZ(iLine) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', imag(Z));
        end
        z = max(this.MinorValue(2,:)) + ...
          1i*[-linspace(-1, 0, 100).^2, linspace(0, 1, 100).^2]*500;
        Z = (z-1)./(z+1);
        this.hRMinorLinesZ(end+1) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', imag(Z));

        % constant imaginary part
        for iLine = 1:size(this.MinorValue, 2)
          z = linspace(0, 1, 100).^2*this.MinorValue(2,iLine) - 1i*this.MinorValue(1,iLine);
          Z = (z-1)./(z+1);
          this.hIMinorLinesZ(2*iLine-1) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', imag(Z));
          this.hIMinorLinesZ(2*iLine) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', -imag(Z));
        end
        z = linspace(0, 1, 100).^2*500 - 1i*max(this.MinorValue(2,:));
        Z = (z-1)./(z+1);
        this.hIMinorLinesZ(end+1) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', imag(Z));
        this.hIMinorLinesZ(end+1) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', -imag(Z));
        set([this.hRMinorLinesZ, this.hIMinorLinesZ], ...
            'Color', this.MinorColor, ...
            'LineStyle', this.MinorLineStyle, ...
            'LineWidth', this.MinorLineWidth, ...
            'Visible', this.MinorVisible);

        %% draw major grid
        % constant real part
        for iLine = 1:size(this.Value, 2)
          z = this.Value(1,iLine) + ...
            1i*[-linspace(-1, 0, 100).^2, linspace(0, 1, 100).^2]*this.Value(2,iLine);
          Z = (z-1)./(z+1);
          this.hRLinesZ(iLine) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', imag(Z));
        end
        z = max(this.Value(2,:)) + ...
          1i*[-linspace(-1, 0, 100).^2, linspace(0, 1, 100).^2]*500;
        Z = (z-1)./(z+1);
        this.hRLinesZ(end+1) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', imag(Z));

        % constant imaginary part
        for iLine = 1:size(this.Value, 2)
          z = linspace(0, 1, 100).^2*this.Value(2,iLine) - 1i*this.Value(1,iLine);
          Z = (z-1)./(z+1);
          this.hILinesZ(2*iLine-1) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', imag(Z));
          this.hILinesZ(2*iLine) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', -imag(Z));
        end
        z = linspace(0, 1, 100).^2*500 - 1i*max(this.Value(2,:));
        Z = (z-1)./(z+1);
        this.hILinesZ(end+1) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', imag(Z));
        this.hILinesZ(end+1) = line('Parent', this.hGroup, 'XData', real(Z), 'YData', -imag(Z));
        set([this.hRLinesZ, this.hILinesZ], ...
            'Color', this.Color, ...
            'LineStyle', this.LineStyle, ...
            'LineWidth', this.LineWidth, ...
            'Visible', this.Visible);
      end

      %% baseLines
      z = exp(1i*2*pi*(0:.01:1));
      this.hBaseLines(1) = line('Parent', this.hGroup, 'XData', real(z), 'YData', imag(z));
      this.hBaseLines(2) = line('Parent', this.hGroup, 'XData', [-1 1], 'YData', [0 0]);
      set(this.hBaseLines, 'Color', this.Color, ...
          'LineStyle', this.LineStyle, ...
          'LineWidth', this.LineWidth);
      % draw dot at (0,1)
      this.hBaseLines(3) = line('Parent', this.hGroup, 'XData', 0, 'YData', 0, ...
        'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'MarkerSize', 5);

      %% draw labels (tick marks) above grid
      if any(lower(this.Type) == 'z')
        % real part labels
        this.hRLabelsZ(1) = text(-1, 0, '0.0 ', ...
          'HorizontalAlignment', 'right', ...
            'Parent', this.hGroup);
        xR = real((this.Value(1,:)-1)./(this.Value(1,:)+1));
        RLabelsZ = cellfun(@(v) sprintf(' %.*f', -floor(log10(abs(v)))+(v>=1), v), num2cell(this.Value(1,:)), 'UniformOutput', false);
        this.hRLabelsZ = [this.hRLabelsZ, ...
          text(xR, zeros(size(xR)), RLabelsZ, 'Parent', this.hGroup, ...
               'Rotation', 90, 'VerticalAlignment', 'bottom')'];
        this.hRLabelsZ(end+1) = text(1, 0, ' \infty', ...
          'HorizontalAlignment', 'left', 'Interpreter', 'tex', ...
            'Parent', this.hGroup);
        if strcmpi(this.Type, 'yz')
          set(this.hRLabelsZ([1 end]), 'VerticalAlignment', 'bottom');
        end

        % imaginary part labels
        zI = 1i*this.Value(1,:);
        xI = 1.13 * real((zI-1)./(zI+1));
        yI = 1.08 * imag((zI-1)./(zI+1));
        ILabelsZ = cellfun(@(v) sprintf('+j%.*f', -floor(log10(abs(v)))+(v>=1), v), num2cell(this.Value(1,:)), 'UniformOutput', false);
        this.hILabelsZ = text(xI, yI, ILabelsZ, 'Parent', this.hGroup)';
        ILabelsZ = cellfun(@(v) sprintf('-j%.*f', -floor(log10(abs(v)))+(v>=1), v), num2cell(this.Value(1,:)), 'UniformOutput', false);
        this.hILabelsZ = [this.hILabelsZ, ...
          text(xI, -yI, ILabelsZ, 'Parent', this.hGroup)'];
        set(this.hILabelsZ, 'HorizontalAlignment', 'center');

        set([this.hRLabelsZ, this.hILabelsZ], ...
          'FontSize', this.LabelSize, ...
          'Visible', this.LabelVisible);
      end

      if any(lower(this.Type) == 'y')
        % real part labels
        this.hRLabelsY(1) = text(-1, 0, '\infty ', ...
          'HorizontalAlignment', 'right', 'Interpreter', 'tex', ...
            'Parent', this.hGroup);
        xR = -real((this.Value(1,:)-1)./(this.Value(1,:)+1));
        RLabelsY = cellfun(@(v) sprintf(' %.*f', -floor(log10(abs(v)))+(v>=1), v), num2cell(this.Value(1,:)), 'UniformOutput', false);
        this.hRLabelsY = [this.hRLabelsY, ...
          text(xR, zeros(size(xR)), RLabelsY, 'Parent', this.hGroup, ...
               'Rotation', -90, 'VerticalAlignment', 'bottom')'];
        this.hRLabelsY(end+1) = text(1, 0, ' 0.0', ...
          'HorizontalAlignment', 'left', ...
            'Parent', this.hGroup);
        set(this.hRLabelsY([1 end]), 'Rotation', 0)
        if strcmpi(this.Type, 'yz')
          set(this.hRLabelsY([1 end]), 'VerticalAlignment', 'top');
        end

        % imaginary part labels
        zI = 1i*this.Value(1,:);
        if strcmpi(this.Type, 'yz')
          offsetMultX = -1.3;
          offsetMultY = 1.2;
        else
          offsetMultX = -1.13;
          offsetMultY = 1.08;
        end
        xI = offsetMultX * real((zI+1)./(zI-1));
        yI = offsetMultY * imag((zI+1)./(zI-1));
        ILabelsY = cellfun(@(v) sprintf('+j%.*f', -floor(log10(abs(v)))+(v>=1), v), num2cell(this.Value(1,:)), 'UniformOutput', false);
        this.hILabelsY = text(xI, yI, ILabelsY, 'Parent', this.hGroup)';
        ILabelsY = cellfun(@(v) sprintf('-j%.*f', -floor(log10(abs(v)))+(v>=1), v), num2cell(this.Value(1,:)), 'UniformOutput', false);
        this.hILabelsY = [this.hILabelsY, ...
          text(xI, -yI, ILabelsY, 'Parent', this.hGroup)'];
        set(this.hILabelsY, 'HorizontalAlignment', 'center');

        set([this.hRLabelsY, this.hILabelsY], ...
          'FontSize', this.LabelSize, ...
          'Visible', this.LabelVisible);
      end
      set(this.allHandles, 'Tag', 'smith');

      if ~isempty(this.axesListener)
        delete(this.axesListener{1})
        this.axesListener = [];
      end
      ca = ancestor(this.hGroup, 'axes');
      if strcmpi(this.Type, 'yz')
        xlim(ca, [-1.15 1.15])
        ylim(ca, [-1.1 1.25])
      else
        xlim(ca, [-1.1 1.1])
        ylim(ca, [-1 1.15])
      end

      % sort existing lines to the front
      kids = get(ca, 'Children');
      hNonSmith = findall(kids, 'flat', '-not', 'Tag', 'smith');
      set(ca, 'Children', [hNonSmith; this.hGroup])

      % add listener for redrawing if overwritten
      set(this.hGroup, 'DeleteFcn', @(src,evt)this.onDelete);
      hold(ca, 'on');

      this.drawPending = false;
    end

  end


  methods (Static = true)

    function smithValue = createGrid(gridTicks)
      %% create Value array from groups of tick values
      % let lines end at max value of each group
      gridValues = cellfun(@(s) [s;max(s)*ones(size(s))], gridTicks, 'UniformOutput', false);
      Values = horzcat(gridValues{:});
      % eliminate lines that are drawn on top of each other
      [uniqueValues, ia] = unique(Values(1,end:-1:1), 'stable');
      smithValue = [uniqueValues; Values(2,size(Values,2)+1-ia)];
    end

  end

  methods (Access = private)

    function ah = allHandles(this)
      ah = [...
        this.hBackground, ...
        this.hBaseLines, ...
        this.hRLinesZ, ...
        this.hILinesZ, ...
        this.hRLabelsZ, ...
        this.hILabelsZ, ...
        this.hRLinesY, ...
        this.hILinesY, ...
        this.hRLabelsY, ...
        this.hILabelsY];
      ah(~ishghandle(ah)) = [];
    end


    function onDelete(this)
      if isempty(this.hGroup) || ~ishghandle(this.hGroup)
        return
      end
      ca = ancestor(this.hGroup, 'axes');
      if isempty(ca) || ~ishghandle(ca) || ...
          (strcmp(get(ca, 'Visible'), 'off') && ishghandle(this.hGroup))
        return
      end
      this.draw();
    end
  end


  %% set and get methods for "graphics properties"
  methods

    function set(this, varargin)
      this.recCount = this.recCount + 1;
      this.doDraw = false;
      iArgin = 1;
      try
        while iArgin <= numel(varargin)
          if isstruct(varargin{iArgin})
            % convert to property-value-pairs
            fields = fieldnames(varargin{iArgin});
            values = cellfun(@(x) varargin{iArgin}.(x), fields, ...
                             'UniformOutput', false);
            fv = [fields, values]';
            this.set(fv{:});
            iArgin = iArgin+1;
          elseif ischar(varargin{iArgin})
            myProperties = properties(this);
            isProp = strcmpi(myProperties, varargin{iArgin});
            if sum(isProp) == 1
              this.(myProperties{isProp}) = varargin{iArgin+1};
              iArgin = iArgin+2;
            else
              error('smith: Invalid property name "%s".', varargin{iArgin});
            end
          else
            error('smith: invalid set')
          end
        end
      catch ME
        this.recCount = 0;
        this.doDraw = true;
        rethrow(ME)
      end
      this.recCount = this.recCount - 1;
      if this.recCount < 1 && this.drawPending
        this.doDraw = true;
        this.draw();
      end
    end

    function value = get(this, prop)
      myProperties = properties(this);
      isProp = strcmpi(myProperties, prop);
      if sum(isProp) == 1
        value = this.(myProperties{isProp});
      else
        error('smith: Invalid property name "%s".', prop);
      end
    end


    % label set functions
    function set.LabelColor(this, newLabelColor)
      if isequal(newLabelColor, this.LabelColor)
        return
      end
      allLabels = [...
        this.hRLabelsZ, ...
        this.hILabelsZ, ...
        this.hRLabelsY, ...
        this.hILabelsY];
      set(allLabels, 'Color', newLabelColor)
      this.LabelColor = newLabelColor;
    end

    function set.LabelSize(this, newLabelSize)
      if isequal(newLabelSize, this.LabelSize)
        return
      end
      allLabels = [...
        this.hRLabelsZ, ...
        this.hILabelsZ, ...
        this.hRLabelsY, ...
        this.hILabelsY];
      set(allLabels, 'FontSize', newLabelSize)
      this.LabelSize = newLabelSize;
    end

    function set.LabelVisible(this, newLabelVisible)
      if isequal(newLabelVisible, this.LabelVisible)
        return
      end
      allLabels = [...
        this.hRLabelsZ, ...
        this.hILabelsZ, ...
        this.hRLabelsY, ...
        this.hILabelsY];
      set(allLabels, 'Visible', newLabelVisible)
      this.LabelVisible = newLabelVisible;
    end

    % set functions for major Z- or Y-Smith chart
    function set.Color(this, newColor)
      if isequal(newColor, this.Color)
        return
      end
      set([this.hRLinesZ, this.hILinesZ], 'Color', newColor)
      if any(lower(this.Type) == 'z')
        set(this.hBaseLines, 'Color', newColor)
      end
      this.Color = newColor;
    end

    function set.LineStyle(this, newLineStyle)
      if isequal(newLineStyle, this.LineStyle)
        return
      end
      set([this.hRLinesZ, this.hILinesZ], 'LineStyle', newLineStyle)
      if any(lower(this.Type) == 'z')
        set(this.hBaseLines, 'LineStyle', newLineStyle)
      end
      this.LineStyle = newLineStyle;
    end

    function set.LineWidth(this, newLineWidth)
      if isequal(newLineWidth, this.LineWidth)
        return
      end
      set([this.hRLinesZ, this.hILinesZ], 'LineWidth', newLineWidth)
      if any(lower(this.Type) == 'z')
        set(this.hBaseLines, 'LineWidth', newLineWidth)
      end
      this.LineWidth = newLineWidth;
    end

    function set.MinorColor(this, newColor)
      if isequal(newColor, this.MinorColor)
        return
      end
      set([this.hRMinorLinesZ, this.hIMinorLinesZ], 'Color', newColor)
      this.MinorColor = newColor;
    end

    function set.MinorLineStyle(this, newLineStyle)
      if isequal(newLineStyle, this.MinorLineStyle)
        return
      end
      set([this.hRMinorLinesZ, this.hIMinorLinesZ], 'LineStyle', newLineStyle)
      this.MinorLineStyle = newLineStyle;
    end

    function set.MinorLineWidth(this, newLineWidth)
      if isequal(newLineWidth, this.MinorLineWidth)
        return
      end
      set([this.hRMinorLinesZ, this.hIMinorLinesZ], 'LineWidth', newLineWidth)
      this.MinorLineWidth = newLineWidth;
    end

    function set.MinorVisible(this, newMinorVisible)
      if isequal(newMinorVisible, this.MinorVisible)
        return
      end
      set([this.hRMinorLinesZ, this.hIMinorLinesZ], 'Visible', newMinorVisible)
      this.MinorVisible = newMinorVisible;
    end

    % set functions for minor Z- or Y-Smith chart
    function set.SubColor(this, newSubColor)
      if isequal(newSubColor, this.SubColor)
        return
      end
      set([this.hRLinesY, this.hILinesY], 'Color', newSubColor)
      if strcmpi(this.Type, 'y')
        set(this.hBaseLines, 'Color', newSubColor)
      end
      this.SubColor = newSubColor;
    end

    function set.SubLineStyle(this, newSubLineStyle)
      if isequal(newSubLineStyle, this.SubLineStyle)
        return
      end
      set([this.hRLinesY, this.hILinesY], 'LineStyle', newSubLineStyle)
      if strcmpi(this.Type, 'y')
        set(this.hBaseLines, 'LineStyle', newSubLineStyle)
      end
      this.SubLineStyle = newSubLineStyle;
    end

    function set.SubLineWidth(this, newSubLineWidth)
      if isequal(newSubLineWidth, this.SubLineWidth)
        return
      end
      set([this.hRLinesY, this.hILinesY], 'LineWidth', newSubLineWidth)
      if strcmpi(this.Type, 'y')
        set(this.hBaseLines, 'LineWidth', newSubLineWidth)
      end
      this.SubLineWidth = newSubLineWidth;
    end

    function set.SubMinorColor(this, newSubColor)
      if isequal(newSubColor, this.SubMinorColor)
        return
      end
      set([this.hRMinorLinesY, this.hIMinorLinesY], 'Color', newSubColor)
      this.SubMinorColor = newSubColor;
    end

    function set.SubMinorLineStyle(this, newSubLineStyle)
      if isequal(newSubLineStyle, this.SubMinorLineStyle)
        return
      end
      set([this.hRMinorLinesY, this.hIMinorLinesY], 'LineStyle', newSubLineStyle)
      this.SubMinorLineStyle = newSubLineStyle;
    end

    function set.SubMinorLineWidth(this, newSubLineWidth)
      if isequal(newSubLineWidth, this.SubMinorLineWidth)
        return
      end
      set([this.hRMinorLinesY, this.hIMinorLinesY], 'LineWidth', newSubLineWidth)
      this.SubMinorLineWidth = newSubLineWidth;
    end

    function set.SubMinorVisible(this, newMinorVisible)
      if isequal(newMinorVisible, this.SubMinorVisible)
        return
      end
      set([this.hRMinorLinesY, this.hIMinorLinesY], 'Visible', newMinorVisible)
      this.SubMinorVisible = newMinorVisible;
    end

    % general set functions
    function set.Type(this, newType)
      if strcmpi(newType, this.Type)
        return
      end
      if ~any(strcmpi(newType, {'z', 'y', 'yz'}))
        error('smith: ''Type'' must be either ''z'', ''y'' or ''yz''.')
      end
      % clear all lines and labels
      this.clearSelf();

      % redraw with new type
      this.Type = lower(newType);
      this.draw();
    end

    function set.Value(this, newValue)
      if isequal(newValue, this.Value)
        return
      end
      sizeNewValue = size(newValue);
      if length(sizeNewValue) ~=2 || sizeNewValue(1) ~= 2 || ~isnumeric(newValue)
        error('smith: ''Value'' must be a numeric 2×N matrix.')
      end
      % clear all lines and labels
      this.clearSelf();

      % redraw with new values
      this.Value = newValue;
      this.draw();
    end

    function set.MinorValue(this, newValue)
      if isequal(newValue, this.MinorValue)
        return
      end
      sizeNewValue = size(newValue);
      if length(sizeNewValue) ~=2 || sizeNewValue(1) ~= 2 || ~isnumeric(newValue)
        error('smith: ''MinorValue'' must be a numeric 2×N matrix.')
      end
      % clear all lines and labels
      this.clearSelf();

      % redraw with new values
      this.MinorValue = newValue;
      this.draw();
    end

    function set.SubValue(this, newValue)
      if isequal(newValue, this.SubValue)
        this.SubValueMode = 'manual';
        return
      end
      sizeNewValue = size(newValue);
      if length(sizeNewValue) ~=2 || sizeNewValue(1) ~= 2 || ~isnumeric(newValue)
        error('smith: ''SubValue'' must be a numeric 2×N matrix.')
      end
      % clear all lines and labels
      this.clearSelf();

      % redraw with new values
      this.SubValueMode = 'manual';
      this.SubValue = newValue;
      this.draw();
    end

    function value = get.SubValue(this)
      if strcmpi(this.SubValueMode, 'auto')
        value = this.Value;
      else
        value = this.SubValue;
      end
    end

    function set.SubValueMode(this, newValue)
      if isequal(newValue, this.SubValueMode)
        return
      end
      if strcmpi(newValue, 'auto')
        this.SubValueMode = lower(newValue);
        this.draw();
      elseif strcmpi(newValue, 'manual')
        this.SubValueMode = lower(newValue);
        if isempty(this.SubValue)
          this.SubValue = this.Value;
        end
      else
        error('smith: ''SubValueMode'' must either be ''auto'' or ''manual''.')
      end
    end

    function set.SubMinorValue(this, newValue)
      if isequal(newValue, this.SubMinorValue)
        this.SubMinorValueMode = 'manual';
        return
      end
      sizeNewValue = size(newValue);
      if length(sizeNewValue) ~=2 || sizeNewValue(1) ~= 2 || ~isnumeric(newValue)
        error('smith: ''SubMinorValue'' must be a numeric 2×N matrix.')
      end
      % clear all lines and labels
      this.clearSelf();

      % redraw with new values
      this.SubMinorValueMode = 'manual';
      this.SubMinorValue = newValue;
      this.draw();
    end

    function value = get.SubMinorValue(this)
      if strcmpi(this.SubMinorValueMode, 'auto')
        value = this.MinorValue;
      else
        value = this.SubMinorValue;
      end
    end

    function set.SubMinorValueMode(this, newValue)
      if isequal(newValue, this.SubMinorValueMode)
        return
      end
      if strcmpi(newValue, 'auto')
        this.SubMinorValueMode = lower(newValue);
        this.draw();
      elseif strcmpi(newValue, 'manual')
        this.SubMinorValueMode = lower(newValue);
        if isempty(this.SubMinorValue)
          this.SubMinorValue = this.MinorValue;
        end
      else
        error('smith: ''SubMinorValueMode'' must either be ''auto'' or ''manual''.')
      end
    end

    function set.Parent(this, value)
      if isequal(value, this.Parent)
        return
      end
      this.clearSelf();
      this.Parent = value;
      this.draw();
    end

  end

end
