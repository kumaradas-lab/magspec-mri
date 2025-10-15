function DrawLocalContour(this, slice, oldcontour, levels)
%% Create a contour on SLICE
% When OLDCONTROUR, recycle that contour patch.
% This does not use the CONTOURSLICE command, but instead uses a
% specialized slice created for sliceomatic.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

cdata = get(slice, 'CData');
st = getappdata(slice, 'slicetype');

% Calculate the new contour for CDATA's values.
if nargin < 4
  if isnan(this.zmesh)==1
    c = contourc(double(cdata));
  else
    switch st
      case 'X'
        c = contours(this.zmesh, this.ymesh, cdata);
      case 'Y'
        c = contours(this.zmesh, this.xmesh, cdata);
      case 'Z'
        c = contours(this.xmesh, this.ymesh, cdata);
    end
  end
else
  if isnan(this.zmesh)==1
    c = contourc(double(cdata), levels);
  else
    switch st
      case 'X'
        c = contours(this.zmesh, this.ymesh, cdata, levels);
      case 'Y'
        c = contours(this.zmesh, this.xmesh, cdata, levels);
      case 'Z'
        c = contours(this.xmesh, this.ymesh, cdata, levels);
    end
  end
end

newvertices = [];
newfaces = {};
longest = 1;
cdata = [];

limit = size(c,2);
i_contour = 1;
while (i_contour < limit)
  z_level = c(1,i_contour);
  npoints = c(2,i_contour);
  nexti = i_contour+npoints+1;

  xdata = c(1,i_contour+1:i_contour+npoints);
  ydata = c(2,i_contour+1:i_contour+npoints);

  switch st
    case 'X'
      xv = get(slice, 'XData');
      lzdata = xv(1,1) + 0*xdata;
      vertices = [lzdata.', ydata.', xdata.'];
    case 'Y'
      yv = get(slice, 'YData');
      lzdata = yv(1,1) + 0*xdata;
      vertices = [ydata.', lzdata.', xdata.'];
    case 'Z'
      zv = get(slice, 'ZData');
      lzdata = zv(1,1) + 0*xdata;
      vertices = [xdata.', ydata.', lzdata.'];
  end

  faces = 1:size(vertices, 1);
  faces = faces + size(newvertices, 1);

  longest = max(longest, size(faces, 2));

  newvertices = [newvertices; vertices];
  newfaces{end+1} = faces;

  tcdata = (z_level + 0*xdata).';

  cdata = [cdata; tcdata]; % need to be same size as faces

  i_contour = nexti;
end

% Glom a NaN on the end for loop-breaking
newvertices = [newvertices; NaN, NaN, NaN];
cdata = [cdata; NaN];

vertmax = size(newvertices, 1);

% Fix up FACES, which is a cell array.
faces = [];
for i_face = 1:size(newfaces, 2)
  faces = [ faces; ...
    newfaces{i_face}, ones(1, longest-size(newfaces{i_face},2))*vertmax, vertmax ];
end

if isempty(oldcontour)
  d = getappdata(this.hParent, 'sliceomatic');
  oldcontour = patch('FaceColor', 'none', 'EdgeColor', d.defcontourcolor, ...
    'LineWidth', d.defcontourlinewidth, 'Parent', this.hAxes);
  try
    set(oldcontour, 'linesmoothing', d.defcontoursmooth);
  catch
  end
  setappdata(slice, 'contour', oldcontour);
end

set(oldcontour, 'Vertices', newvertices, ...
  'Faces', faces, ...
  'FaceVertexCData', cdata);

end
