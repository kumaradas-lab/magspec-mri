function [a, s] = getarrowslice(hobj)
%% Return the Arrow and Slice based on the input handle or GCO

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc
% Copyright 2019-2023 Pure Devices GmbH

a = [];
s = [];

if nargin < 1
  hobj = gco();
end

if isempty(hobj)
  return;
end

if isempty(getappdata(hobj, 'controlarrow')) && ...
    isempty(getappdata(hobj, 'isosurface'))
  a = hobj;
  s = getappdata(a, 'arrowslice');
  if isempty(s)
    s = getappdata(a, 'arrowiso');
  end
else
  s = hobj;
  if ~isempty(getappdata(s, 'isosurface'))
    s = getappdata(s, 'isosurface');
  end
  a = getappdata(s, 'controlarrow');
end

end
