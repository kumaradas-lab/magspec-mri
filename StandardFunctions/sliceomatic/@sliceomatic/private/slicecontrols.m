function slicecontrols(hParent, onoff, xmesh, ymesh, zmesh, xdir, ydir, zdir)
%% Convert figure to contain controls for manipulating slices.

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

%% Check variables
narginchk(2, 8)

d = getappdata(hParent, 'sliceomatic');

if onoff
  if nargin < 3
    % If there is no user supplied mesh, make one up.
    xmesh(1) = 1;
    xmesh(2) = size(d.data,2);
    ymesh(1) = 1;
    ymesh(2) = size(d.data,1);
    zmesh(1) = 1;
    zmesh(2) = size(d.data,3);
    
    xdir = 'normal';
    ydir = 'normal';
    zdir = 'normal';
  end

  set([d.axx d.axy d.axz], 'HandleVisibility', 'on');

  set(d.axx, 'XLim', [xmesh(1) xmesh(end)], ...
    'YLim', [1 5]);
  set(d.pxx, 'Vertices', [ xmesh(1) xmesh(1) -1; xmesh(end) xmesh(1) -1; xmesh(end) 5 -1; xmesh(1) 5 -1], ...
    'Faces', [1 2 3 ; 1 3 4]);

  activelabel(d.axx, 'title', 'X');

  set(d.axy, 'XLim', [1 5], ...
    'YLim', [ymesh(1) ymesh(end)]);
  set(d.pxy, 'Vertices', [ymesh(1) ymesh(1) -1; ymesh(1) ymesh(end) -1; 5 ymesh(end) -1; 5 ymesh(1) -1], ...
    'Faces', [1 2 3 ; 1 3 4]);
  activelabel(d.axy, 'title', 'Y');

  set(d.axz, 'XLim', [1 5], ...
    'YLim', [zmesh(1) zmesh(end)]);
  set(d.pxz, 'Vertices', [zmesh(1) zmesh(1) -1; zmesh(1) zmesh(end) -1; 5 zmesh(end) -1; 5 zmesh(1) -1], ...
    'Faces', [1 2 3 ; 1 3 4]);
  activelabel(d.axz, 'title', 'Z');

  set([d.axx d.axy d.axz], 'HandleVisibility', 'off');

  set(d.axx, 'XDir', xdir);
  set(d.axy, 'YDir', ydir);
  set(d.axz, 'ZDir', zdir);
else
  % Disable these controls.  Perhaps hide all slices?
end

end
