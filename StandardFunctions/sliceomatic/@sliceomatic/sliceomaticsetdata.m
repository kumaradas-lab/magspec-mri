function d = sliceomaticsetdata(this, d, xmesh, ymesh, zmesh)
%% Create the data used for sliceomatic in the appdata D
%
%   d = this.sliceomaticsetdata(d, xmesh, ymesh, zmesh)

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% Check variables
narginchk(5, 5);

%% Simplify the isonormals
% as smooth3 is taking the most time for large volumes
if size(d.data, 1) * size(d.data, 2) * size(d.data, 3) > 10^6
    disp('Smoothing for IsoNormals... SKIPPED!');
    disp('Do not use ISO surface controller for this large volume...');
    d.smooth = d.data;
else
    fprintf('Smoothing for IsoNormals...');
    d.smooth = smooth3(d.data);  % ,'box',5);
    fprintf(' done\n');
end
d.reducenumbers = [floor(size(d.data,2)/20), ...
                   floor(size(d.data,1)/20), ...
                   floor(size(d.data,3)/20)];
d.reducenumbers(d.reducenumbers==0) = 1;


if ~isnan(xmesh)
  % Set axis orientation
  xdir = 'normal';
  ydir = 'normal';
  zdir = 'normal';
  if ~issorted(xmesh)
    xmesh = flip(xmesh, 2);
    xdir = 'reverse';
    d.xlim = [xmesh(1), xmesh(end)];
    xmesh = flip(xmesh, 2);
  else
    d.xlim = [xmesh(1), xmesh(end)];
  end
  if ~issorted(ymesh)
    ymesh = flip(ymesh,2);
    ydir = 'reverse';
    d.ylim = [ymesh(1), ymesh(end)];
    ymesh = flip(ymesh,2);
  else
    d.ylim = [ymesh(1), ymesh(end)];
  end
  % This should not be the case for medical images
  if ~issorted(zmesh)
    zmesh = flip(zmesh, 2);
    zdir = 'reverse';
    d.zlim = [zmesh(1), zmesh(end)];
    zmesh = flip(zmesh, 2);
  else
    d.zlim = [zmesh(1), zmesh(end)];
  end

  % Vol vis suite takes numbers in X/Y form.
  ly = 1:d.reducenumbers(1):size(d.data,2);
  lx = 1:d.reducenumbers(2):size(d.data,1);
  lz = 1:d.reducenumbers(3):size(d.data,3);

  for i = 1:length(ly)
    ly(i) = xmesh(ly(i));
  end
  for i = 1:length(lx)
    lx(i) = ymesh(lx(i));
  end
  for i = 1:length(lz)
    lz(i) = zmesh(lz(i));
  end

  d.reducelims = {ly, lx, lz};
  fprintf('Generating reduction volume...');
  d.reduce = reducevolume(d.data, d.reducenumbers);
  d.reducesmooth = smooth3(d.reduce, 'box', 5);
  fprintf(' done\n');

  % Set axis
  %d.xlim = [xmesh(1) xmesh(end)];
  %d.ylim = [ymesh(1) ymesh(end)];
  %d.zlim = [zmesh(1) zmesh(end)];
  this.xmesh = xmesh;
  this.ymesh = ymesh;
  this.zmesh = zmesh;
  d.xdir = xdir;
  d.ydir = ydir;
  d.zdir = zdir;

else
  % Vol vis suite takes numbers in X/Y form.
  ly = 1:d.reducenumbers(1):size(d.data,2);
  lx = 1:d.reducenumbers(2):size(d.data,1);
  lz = 1:d.reducenumbers(3):size(d.data,3);

  d.reducelims = {ly, lx, lz};
  fprintf('Generating reduction volume...');
  d.reduce = reducevolume(d.data, d.reducenumbers);
  d.reducesmooth = smooth3(d.reduce, 'box', 5);
  fprintf(' done\n');

  d.xlim = [1, size(d.data,2)];
  d.ylim = [1, size(d.data,1)];
  d.zlim = [1, size(d.data,3)];
  this.xmesh = NaN;
  this.ymesh = NaN;
  this.zmesh = NaN;
  d.xdir = 'normal';
  d.ydir = 'normal';
  d.zdir = 'normal';
end

end
