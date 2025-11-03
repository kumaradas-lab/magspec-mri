function callbacks(this, varargin)
%% Callback functions for sliceomatic
%
%     this.callbacks(varargin)
%
% Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008
%           The MathWorks Inc
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% Check input
if ~ischar(varargin{1})
  warning('sliceomatic:WrongCallback', 'Sliceomatic callback must be a string');
  return;
end

%% Interpret commands
d = getappdata(this.hParent, 'sliceomatic');
hFigure = ancestor(this.hAxes, 'figure');
if ~isequal(hFigure, this.hFigure)
  % re-ancestoring
  % FIXME: This should be done in a listener to the ancestor figure.
  %        How can this be done?
  this.FigureDeleteFcn();
  this.hFigure = hFigure;
end

try
  switch varargin{1}
    case 'Xnew'
      if strcmp(get(hFigure, 'SelectionType'), 'normal')
        pt = get(gcbo(), 'CurrentPoint');
        X = pt(1,1);
        this.AddSliceX(X);

        % Make sure whatever buttonupfcn on the figure is run now to "turn
        % off" whatever was going on before we got our callback on the
        % arrow.
        buf = get(hFigure, 'WindowButtonUpFcn');
        checkEvalBuffer(buf);
        this.ArrowDragPrepare();
        setpointer(hFigure, 'SOM leftright');
        set(d.motionmetaslice, 'Visible', 'off');

        % store data for tip at arrow
        % FIXME: This is pretty similar to what is done in "motion". Take care
        % to use the same position offsets.
        tipdata.parentaxes = this.hSliderX;
        tipdata.value = this.GetSlicePos(X, 'x');
        tipdata.position = [X, -1];
        tipdata.verticalalign = 'top';
        tipdata.horizontalalign = 'center';
        setappdata(this.draggedArrow, 'tipdata', tipdata);
      end

    case 'Ynew'
      if strcmp(get(hFigure, 'SelectionType'), 'normal')
        pt = get(gcbo(), 'CurrentPoint');
        Y = pt(1,2);
        this.AddSliceY(Y);

        % Make sure whatever buttonupfcn on the figure is run now to "turn
        % off" whatever was going on before we got our callback on the
        % arrow.
        buf = get(hFigure, 'WindowButtonUpFcn');
        checkEvalBuffer(buf);
        this.ArrowDragPrepare();
        setpointer(hFigure, 'SOM topbottom');
        set(d.motionmetaslice, 'Visible', 'off');

        % store data for tip at arrow
        % FIXME: This is pretty similar to what is done in "motion". Take care
        % to use the same position offsets.
        tipdata.parentaxes = this.hSliderY;
        tipdata.value = this.GetSlicePos(Y, 'y');
        tipdata.position = [7, Y];
        tipdata.verticalalign = 'middle';
        tipdata.horizontalalign = 'left';
        setappdata(this.draggedArrow, 'tipdata', tipdata);
      end

    case 'Znew'
      if strcmp(get(hFigure, 'SelectionType'), 'normal')
        pt = get(gcbo(), 'CurrentPoint');
        Y = pt(1,2);
        this.AddSliceZ(Y);

        % Make sure whatever buttonupfcn on the figure is run now to "turn
        % off" whatever was going on before we got our callback on the
        % arrow.
        buf = get(hFigure, 'WindowButtonUpFcn');
        checkEvalBuffer(buf);
        this.ArrowDragPrepare();
        setpointer(hFigure, 'SOM topbottom');
        set(d.motionmetaslice, 'Visible', 'off');

        % store data for tip at arrow
        % FIXME: This is pretty similar to what is done in "motion". Take care
        % to use the same position offsets.
        tipdata.parentaxes = this.hSliderZ;
        tipdata.value = this.GetSlicePos(Y, 'z');
        tipdata.position = [-1, Y];
        tipdata.verticalalign = 'middle';
        tipdata.horizontalalign = 'right';
        setappdata(this.draggedArrow, 'tipdata', tipdata);
      end

    case 'ISO'
      if strcmpi(get(hFigure, 'SelectionType'), 'normal')
        pt = get(gcbo(), 'CurrentPoint');
        V = pt(1,1);
        this.AddIsoSurface(V);

        % Make sure whatever buttonupfcn on the figure is run now to "turn
        % off" whatever was going on before we got our callback on the
        % arrow.
        buf = get(hFigure, 'WindowButtonUpFcn');
        checkEvalBuffer(buf);
        this.ArrowDragPrepare();
        setpointer(hFigure, 'SOM leftright');

        % store data for tip at arrow
        % FIXME: This is pretty similar to what is done in "motion". Take care
        % to use the same position offsets.
        tipdata.parentaxes = this.hSliderIso;
        tipdata.value = V;
        tipdata.position = [V, 7];
        tipdata.verticalalign = 'bottom';
        tipdata.horizontalalign = 'center';
        setappdata(this.draggedArrow, 'tipdata', tipdata);
      end

    case 'Xmove'
      if strcmp(get(hFigure, 'SelectionType'), 'normal')
        setpointer(hFigure, 'SOM leftright');
        this.draggedArrow = getarrowslice();
        this.ArrowDragPrepare();
      end
    case 'Ymove'
      if strcmp(get(hFigure, 'SelectionType'), 'normal')
        this.draggedArrow = getarrowslice();
        this.ArrowDragPrepare();
      end
    case 'Zmove'
      if strcmp(get(hFigure, 'SelectionType'), 'normal')
        this.draggedArrow = getarrowslice();
        this.ArrowDragPrepare();
      end
    case 'ISOmove'
      if strcmp(get(hFigure, 'SelectionType'), 'normal')
        setpointer(hFigure, 'SOM leftright');
        this.draggedArrow = getarrowslice();
        this.ArrowDragPrepare();
      end
    case 'motion'
      % Make sure our cursor is ok
      a = this.draggedArrow;      % The arrow being dragged
      s = getappdata(a, 'arrowslice');  % The slice to 'move'
      if isempty(s)
        s = getappdata(a, 'arrowiso');  % or the isosurface
      end
      aa = get(a, 'Parent');    % arrow's parent axes
      pos = getappdata(a, 'arrowcenter');  % the line the arrow points at.
      apos = get(aa, 'CurrentPoint');

      % Bind the axes position to the limits of that axes.
      xlimits = get(aa, 'XLim');
      ylimits = get(aa, 'YLim');

      if apos(1,1) < xlimits(1)
        apos(1,1) = xlimits(1);
      elseif apos(1,1) > xlimits(2)
        apos(1,1) = xlimits(2);
      end

      if apos(1,2) < ylimits(1)
        apos(1,2) = ylimits(1);
      elseif apos(1,2) > ylimits(2)
        apos(1,2) = ylimits(2);
      end

      if aa==this.hSliderX || aa==this.hSliderIso
        % We are moving an X slice or iso surface
        xdiff = apos(1,1) - pos;
        v = get(a, 'Vertices');
        v(:,1) = v(:,1) + xdiff;
        set(a, 'Vertices', v);
        np = apos(1,1);
        % This might be a slice, or an isosurface!
        if aa == this.hSliderIso
          new = this.DrawLocalIsoSurface(d.reducelims, d.reduce, d.reducesmooth, ...
            apos(1,1), s);
          setappdata(new, 'reduced',1);
          movetipforarrow(this.hParent, a, aa, apos(1,1), [apos(1,1), 7], ...
            'bottom', 'center');
        else
          % disp([ 'apos = ' num2str(apos(1,1))])
          % disp([ 'pos  = ' num2str(pos(1,1))])
          % disp([ 'change=' num2str(this.GetSlicePos(apos(1,1), 'x') ~= this.GetSlicePos(pos(1,1), 'x'))]);
          slicePos = this.GetSlicePos(apos(1,1), 'x');
          if slicePos ~= this.GetSlicePos(pos, 'x')
            this.DrawLocalSlice(d.data, apos(1,1), [], [], s);
          end
          movetipforarrow(this.hParent, a, aa, slicePos, [apos(1,1), -1], ...
            'top', 'center')
        end
      else
        % We are moving a Y or Z slice
        ydiff = apos(1,2) - pos;
        v = get(a, 'Vertices');
        v(:,2) = v(:,2) + ydiff;
        set(a, 'Vertices', v);
        np = apos(1,2);
        if aa == this.hSliderY
          slicePos = this.GetSlicePos(apos(1,2), 'y');
          if slicePos ~= this.GetSlicePos(pos, 'y')
            this.DrawLocalSlice(d.data, [], apos(1,2), [], s);
          end
          movetipforarrow(this.hParent, a, aa, slicePos, [7, apos(1,2)], ...
            'middle', 'left');
        else
          slicePos = this.GetSlicePos(apos(1,2), 'z');
          if slicePos ~= this.GetSlicePos(pos, 'z')
            this.DrawLocalSlice(d.data, [], [], apos(1,2), s);
          end
          movetipforarrow(this.hParent, a, aa, slicePos, [-1, apos(1,2)], ...
            'middle', 'right');
        end
      end
      setappdata(a, 'arrowcenter', np);
      % Question: Is anyone dependent on versions of MATLAB
      % that do not support Java Figures?
      %
      %% This improves animation speed in Java Figures/R14
      %% The Rule: Java Figures don't want this drawnow.
      %try
      %  if isempty(get(gcf,'javaframe'))
      %    drawnow;
      %  end
      %catch
      %  drawnow;
      %end
    %
    % IsoSurface context menu items
    %
    case 'isotogglevisible'
      [a, s] = getarrowslice();
      if propcheck(s,'visible','on')
        set(s,'visible','off');
      else
        set(s,'visible','on');
      end
    case 'isodelete'
      [a, s] = getarrowslice();
      hArrows = getappdata(gco, 'Arrows');
      hArrows(hArrows == a) = [];
      setappdata(gco, 'Arrows', hArrows);
      cap = getappdata(s, 'sliceomaticisocap');
      if ~isempty(cap)
        delete(cap);
      end
      delete(s);
      delete(a);
    case 'isoflatlight'
      [a, s] = getarrowslice();
      set(s,'facelighting','flat');
    case 'isosmoothlight'
      [a, s] = getarrowslice();
      set(s,'facelighting','phong');
    case 'isocolor'
      [a, s] = getarrowslice();
      c=uisetcolor(get(s,'facecolor'));
      slowset(s,'facecolor',c,d.animincrement);
    case 'isoalpha'
      [a, s] = getarrowslice();
      if length(varargin) ~= 2
        error('Not enough arguments to sliceomatic.');
      end
      slowset(s, 'facealpha', varargin{2}, d.animincrement);
    case 'isocaps'
      [a, s] = getarrowslice();
      cap = getappdata(s,'isosurfacecap');
      if isempty(cap)
        new = this.DrawLocalIsoCaps(s);
        set(new, 'UIContextMenu', this.isoContextMenu);
      else
        delete(cap);
        setappdata(s, 'isosurfacecap', []);
      end
    %
    % Slice context menu items
    %
    case 'togglevisible'
      [a, s] = getarrowslice();
      switch get(s, 'Visible')
        case 'on'
          set(s, 'Visible', 'off');
          pushset(a, 'FaceAlpha', .2);
        case 'off'
          set(s, 'Visible', 'on');
          popset(a, 'FaceAlpha');
      end
    case 'setfaceted'
      [a, s] = getarrowslice();
      set(s, 'EdgeColor', 'k', 'FaceColor', 'flat');
      if ischar(get(s, 'FaceAlpha')) && strcmp(get(s, 'FaceAlpha'), 'texturemap')
        set(s, 'FaceAlpha', 'flat');
      end
      textureizeslice(s, 'off');
    case 'setflat'
      [a, s] = getarrowslice();
      set(s, 'EdgeColor', 'none', 'FaceColor', 'flat');
      if ischar(get(s, 'FaceAlpha')) && strcmp(get(s, 'FaceAlpha'), 'texturemap')
        set(s, 'FaceAlpha', 'flat');
      end
      textureizeslice(s, 'off');
    case 'setinterp'
      [a, s] = getarrowslice();
      set(s, 'EdgeColor', 'none', 'FaceColor', 'interp');
      if ischar(get(s, 'FaceAlpha')) && strcmp(get(s, 'FaceAlpha'), 'texturemap')
        set(s, 'FaceAlpha', 'interp');
      end
      textureizeslice(s, 'off');
    case 'settexture'
      [a, s] = getarrowslice();
      set(s, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
      if ischar(get(s, 'FaceAlpha'))
        set(s, 'FaceAlpha', 'texturemap');
      end
      textureizeslice(s, 'on');
    case 'setnone'
      [a, s] = getarrowslice();
      set(s, 'FaceColor', 'none', 'EdgeColor', 'none');
      textureizeslice(s, 'off');
    case 'setalphanone'
      [a, s] = getarrowslice();
      slowset(s,'facealpha',1,d.animincrement);
    case 'setalphapoint5'
      [a, s] = getarrowslice();
      slowset(s,'facealpha',.5,d.animincrement);
    case 'setalphaflat'
      [a, s] = getarrowslice();
      set(s, 'FaceAlpha', 'flat');
      if ischar(get(s,'FaceColor')) && strcmp(get(s, 'FaceColor'), 'texturemap')
        set(s, 'FaceColor', 'flat');
        textureizeslice(s, 'off');
      end
    case 'setalphainterp'
      [a, s] = getarrowslice();
      set(s, 'FaceAlpha', 'interp');
      if ischar(get(s, 'FaceColor')) && strcmp(get(s, 'FaceColor'), 'texturemap')
        set(s, 'FaceColor', 'interp');
        textureizeslice(s, 'off');
      end
    case 'setalphatexture'
      [a, s] = getarrowslice();
      set(s, 'FaceAlpha', 'texturemap');
      if ischar(get(s, 'FaceColor'))
        set(s, 'FaceColor', 'texturemap');
        textureizeslice(s, 'on');
      end
    case 'slicecontour'
      [a, s] = getarrowslice();
      this.DrawLocalContour(s, getappdata(s, 'contour'));
    case 'slicecontourfullauto'
      [a, s] = getarrowslice();
      minmax = get(this.hSliderIso, 'XLim');
      levels = minmax(1):(minmax(2)-minmax(1))/10:minmax(2);
      setappdata(s, 'contourlevels', levels);
      this.DrawLocalContour(s, getappdata(s, 'contour'), levels);
    case 'slicecontour_setauto'
      [a, s] = getarrowslice();
      setappdata(s, 'contourlevels', []);
      this.DrawLocalContour(s, getappdata(s, 'contour'));
    case 'slicecontour_setfullauto'
      [a, s] = getarrowslice();
      minmax = get(this.hSliderIso, 'XLim');
      levels = minmax(1):(minmax(2)-minmax(1))/10:minmax(2);
      setappdata(s, 'contourlevels', levels);
      this.DrawLocalContour(s, getappdata(s, 'contour'), levels);
    case 'slicecontour_select'
      [a, s] = getarrowslice();
      xl = get(this.hSliderIso, 'XLim');
      levels = selectcontourlevels(get(s,'cdata'), xl(1), xl(2));
      setappdata(s, 'contourlevels', levels);
      this.DrawLocalContour(s, getappdata(s, 'contour'), levels);
    case 'slicecontour_setlevels'
      [a, s] = getarrowslice();
      xl = get(this.hSliderIso, 'XLim');
      levels = selectcontourlevels(get(s, 'CData'), xl(1), xl(2));
      setappdata(s, 'contourlevels', levels);
      this.DrawLocalContour(s, getappdata(s, 'contour'), levels);
    case 'deleteslice'
      [a, s] = getarrowslice();
      if ~isempty(getappdata(s,'contour'))
        delete(getappdata(s,'contour'));
      end
      delete(s);
      delete(a);
    case 'deleteslicecontour'
      [a, s] = getarrowslice();
      if ~isempty(getappdata(s,'contour'))
        delete(getappdata(s,'contour'));
      end
      temp=getappdata(s);
      try
        temp.contourlevels;
        setappdata(s,'contourlevels',[]);
      end
      setappdata(s,'contour',[]);
    case 'slicecontourflat'
      [a, s] = getarrowslice();
      c = getappdata(s,'contour');
      if ~isempty(c)
        set(c,'edgecolor','flat');
      end
    case 'slicecontourinterp'
      [a, s] = getarrowslice();
      c = getappdata(s,'contour');
      if ~isempty(c)
        set(c,'edgecolor','interp');
      end
    case 'slicecontourblack'
      [a, s] = getarrowslice();
      c = getappdata(s,'contour');
      if ~isempty(c)
        set(c, 'EdgeColor', 'black');
      end
    case 'slicecontourwhite'
      [a, s] = getarrowslice();
      c = getappdata(s,'contour');
      if ~isempty(c)
        set(c, 'EdgeColor', 'white');
      end
    case 'slicecontoursmooth'
      [a, s] = getarrowslice();
      c = getappdata(s,'contour');
      onoff = get(gcbo, 'Checked');
      switch onoff
        case 'off'
          set(c, 'LineSmoothing', 'on');
        case 'on'
          set(c, 'LineSmoothing', 'off');
      end
    case 'slicecontourcolor'
      [a, s] = getarrowslice;
      c = getappdata(s,'contour');
      if ~isempty(c)
        inputcolor = get(c, 'EdgeColor');
        if ischar(inputcolor)
          inputcolor=[ 1 1 1 ];
        end
        slowset(c,'edgecolor',uisetcolor(inputcolor),d.animincrement);
      end
    case 'slicecontourlinewidth'
      if ischar(varargin{2})
        val = str2double(varargin{2});
      else
        val = varargin{2};
      end
      [a, s] = getarrowslice();
      c = getappdata(s,'contour');
      if ~isempty(c)
        slowset(c, 'linewidth', val, d.animincrement);
      end
    %
    % Menu All Slices
    %
    case 'allfacet'
      s = this.GetAllSlices();
      set(s, 'FaceColor', 'flat', 'EdgeColor', 'k');
      textureizeslice(s, 'off');
    case 'allflat'
      s = this.GetAllSlices();
      set(s, 'FaceColor', 'flat', 'EdgeColor', 'none');
      textureizeslice(s, 'off');
    case 'allinterp'
      s = this.GetAllSlices();
      set(s, 'FaceColor', 'interp', 'EdgeColor', 'none');
      textureizeslice(s, 'off');
    case 'alltex'
      s = this.GetAllSlices();
      set(s, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
      textureizeslice(s, 'on');
    case 'allnone'
      s = this.GetAllSlices();
      set(s, 'FaceColor', 'none', 'EdgeColor', 'none');
      textureizeslice(s, 'off');
    case 'alltnone'
      s = this.GetAllSlices();
      set(s, 'FaceAlpha', 1);
      textureizeslice(s, 'off');
    case 'alltp5'
      s = this.GetAllSlices();
      set(s, 'FaceAlpha', .5);
      textureizeslice(s, 'off');
    case 'alltflat'
      s = this.GetAllSlices();
      set(s, 'FaceAlpha', 'flat');
      textureizeslice(s, 'off');
    case 'alltinterp'
      s = this.GetAllSlices();
      set(s, 'FaceAlpha', 'interp');
      textureizeslice(s, 'off');
    case 'allttex'
      s = this.GetAllSlices();
      set(s, 'FaceAlpha', 'texturemap');
      textureizeslice(s, 'on');
    %
    % Menu Defaults callbacks
    %
    case 'defaultfaceted'
      d.defcolor='faceted';
    case 'defaultflat'
      d.defcolor='flat';
    case 'defaultinterp'
      d.defcolor='interp';
    case 'defaulttexture'
      d.defcolor='texture';
      if strcmp(d.defalpha,'flat') || strcmp(d.defalpha,'interp')
        d.defalpha='texture';
      end
    case 'defaultcolornone'
      d.defcolor='none';
    case 'defaulttransnone'
      d.defalpha='none';
    case 'defaulttransflat'
      d.defalpha='flat';
    case 'defaulttransinterp'
      d.defalpha='interp';
    case 'defaulttranstexture'
      d.defalpha='texture';
      d.defcolor='texture';
    case 'defaultlightflat'
      d.deflight='flat';
    case 'defaultlightsmooth'
      d.deflight='smooth';
    case 'defaultIsoAlpha'
      if ischar(varargin{2})
        d.defisoalpha = str2double(varargin{2});
      else
        d.defisoalpha = varargin{2};
      end
    case 'defaultcontoursmooth'
      d.defaultcontoursmooth='on';
    case 'defaultcontourflat'
      d.defcontourcolor='flat';
    case 'defaultcontourinterp'
      d.defcontourcolor='interp';
    case 'defaultcontourblack'
      d.defcontourcolor='black';
    case 'defaultcontourwhite'
      d.defcontourcolor='white';
    case 'defaultcontourlinewidth'
      if ischar(varargin{2})
        d.defcontourlinewidth = str2double(varargin{2});
      else
        d.defcontourlinewidth = varargin{2};
      end
    %
    % Camera toolbar Toggling
    %
    case 'cameratoolbar'
      cameratoolbar('Toggle');
    case 'annotationtoolbar'
      if propcheck(d.toolbar, 'Visible', 'on')
        set(d.toolbar, 'Visible', 'off');
      else
        set(d.toolbar, 'Visible', 'on');
      end
    %
    % Controller Preferences
    %
    case 'controlalpha'
      if ischar(varargin{2})
        val = str2double(varargin{2});
      else
        val = varargin{2};
      end
      iso = findobj(this.hSliderIso, 'Type', 'image');
      if val == 0
        set([d.pxx, d.pxy, d.pxz, iso], 'Visible', 'off');
      else
        set([d.pxx, d.pxy, d.pxz, iso], 'Visible', 'on');
        slowset([d.pxx, d.pxy, d.pxz], 'FaceAlpha', val, d.animincrement);
        slowset(iso, 'AlphaData', val, d.animincrement);
      end
    case 'toggleanimation'
      if d.animincrement == 0
        d.animincrement = 10;
      else
        d.animincrement = 0;
      end
    case 'controllabels'
      l = get(this.hSliderX, 'XTickLabel');
      if isempty(l)
        set([this.hSliderX, this.hSliderIso], 'XTickLabelMode', 'auto');
        set([this.hSliderY, this.hSliderZ], 'YTickLabelMode', 'auto');
      else
        set([this.hSliderX, this.hSliderIso], 'XTickLabel', []);
        set([this.hSliderY, this.hSliderZ], 'YTickLabel', []);
      end
    case 'controlvisible'
      objs = findobj([this.hSliderIso, this.hSliderX, this.hSliderY, this.hSliderZ]);
      if strcmp(get(this.hSliderX, 'Visible'), 'on')
        set(objs, 'Visible', 'off');
      else
        set(objs, 'Visible', 'on');
      end
      this.SizeChangedFcn();
    %
    % UICONTROL callbacks
    %
    case 'colormap'
      str = get(gcbo(), 'String');
      val = str{get(gcbo(), 'Value')};
      if strcmp(val, 'custom')
        colormapeditor;
        % FIXME: Execution should wait until colormapeditor was closed.
        val = colormap(this.hAxes);
      else
        if strcmp(val, 'rand')
          cm = get(this.hFigure, 'Colormap');
          val = feval(@rand, size(cm));
        else
          val = feval(val);
        end
        slowset(this.hFigure, 'colormap', val, d.animincrement);
      end
      % update colors of ISO surfaces
      isosurfs = this.GetAllIsos();
      clim = get(this.hAxes, 'CLim');
      clen = clim(2) - clim(1);
      for iIso = 1:length(isosurfs)
        value = getappdata(isosurfs(iIso), 'isosurfacevalue');
        idx = fix((value - clim(1))*length(val)/clen) + 1;
        set(isosurfs(iIso), 'FaceColor', val(idx,:));
      end
    case 'alphamap'
      str = get(gcbo(), 'String');
      str = str{get(gcbo(), 'Value')};
      if strcmp(str, 'rand')
        val = rand(size(alphamap));
      else
        val = alphamap(str);
      end
      slowset(this.hFigure, 'alphamap', val, d.animincrement);
    case 'orientation'
      str = get(gcbo(), 'String');
      str = str{get(gcbo(), 'Value')};
      rotateView(this.hAxes, str);
    %
    % Commands
    %
    case 'copy'
      hax_copy = copyobj(this.hAxes, figure);
      set(hax_copy, 'Units', 'normalized', 'Position', [.05, .05, .9, .9]);
    case 'print'
      newf = figure('Visible', 'off', 'Renderer', get(this.hFigure, 'Renderer'));
      hax_copy = copyobj(this.hAxes, newf);
      set(hax_copy, 'Units', 'normalized', 'Position', [.05, .05, .9, .9]);
      printdlg(newf);
      close(newf);
    otherwise
      error('Bad slice-o-matic command.');
  end
catch ME
  disp(ME.message)
  disp(ME.stack)
  disp(getReport(ME));
end
setappdata(this.hParent, 'sliceomatic', d);

end


function movetipforarrow(hParent, arrow, ax, value, position, va, ha)
% Setup the current data tip for a slice arrow, and show it's
% control value

tipdata.parentaxes = ax;
tipdata.value = value;
tipdata.position = position;
tipdata.verticalalign = va;
tipdata.horizontalalign = ha;

setappdata(arrow, 'tipdata', tipdata);

showarrowtip(hParent, arrow);
% Put it onto this.hSliderIso so that
% it always appears on top.
%set(t, 'Parent', this.hSliderIso);
end


%{
function working(onoff)

  ax=getappdata(gcf,'workingaxis');
  d = getappdata(gcf, 'sliceomatic');

  if isempty(ax)
    ax = axes('Units', 'normalized', 'Position', [.3 .4 .4 .2], ...
              'Box', 'on', 'YTick', [], 'XTick', [], ...
              'XLim', [-1 1], 'YLim', [-1 1], ...
              'Color', 'none', 'HandleVisibility', 'off', ...
              'Parent', d.handleParent);
    text('Parent', ax, 'String', 'Working...', 'FontSize', 64, ...
         'Position', [0 0], ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'EraseMode', 'xor');
    setappdata(gcf,'workingaxis',ax);
  end

  disp(['Working...' onoff]);
  set([ax get(ax, 'Children')], 'Visibility', onoff);
end
%}


function rotateView(hAxes, orientation)
%% Rotate the view of the main axes according to keyword

switch orientation
  case 'default'
    set(hAxes, 'View', [-37.5, 30]);
  case 'X-Y'
    set(hAxes, 'View', [0, 90]);
  case 'X-Z'
    set(hAxes, 'View', [0, 0]);
  case 'Y-Z'
    set(hAxes, 'View', [90, 0]);
  otherwise
    error('sliceomatic:rotateView:Unknown', 'Unknown orientation identifier');
end

end
