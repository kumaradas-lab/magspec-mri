function setpointer(hf, ptr)
%% Set the pointer on the current figure to PTR
%
% This function also provides several specialized SOM (SliceOMatic) pointers.
%
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc
%
% ------------------------------------------------------------------------------
% (C) Copyright 2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

pd_r = [ ...
  nan nan nan nan nan nan nan nan 1   nan nan nan nan nan nan nan
  nan nan nan nan nan nan nan nan 1   1   nan nan nan nan nan nan
  nan nan nan nan nan nan nan nan 1   1   1   nan nan nan nan nan
  nan nan nan nan nan nan nan nan 1   2   1   1   nan nan nan nan
  1   1   1   1   1   1   1   1   1   2   2   1   1   nan nan nan
  1   2   2   2   2   2   2   2   2   2   2   2   1   1   nan nan
  1   2   2   2   2   2   2   2   2   2   2   2   2   1   1   nan
  1   2   2   2   2   2   2   2   2   2   2   2   2   2   1   1
  1   2   2   2   2   2   2   2   2   2   2   2   2   2   1   1
  1   2   2   2   2   2   2   2   2   2   2   2   2   1   1   nan
  1   2   2   2   2   2   2   2   2   2   2   2   1   1   nan nan
  1   1   1   1   1   1   1   1   1   2   2   1   1   nan nan nan
  nan nan nan nan nan nan nan nan 1   2   1   1   nan nan nan nan
  nan nan nan nan nan nan nan nan 1   1   1   nan nan nan nan nan
  nan nan nan nan nan nan nan nan 1   1   nan nan nan nan nan nan
  nan nan nan nan nan nan nan nan 1   nan nan nan nan nan nan nan ];

pd_lr = [ ...
  nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan
  nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan
  nan nan nan nan nan 1   nan nan nan nan 1   nan nan nan nan nan
  nan nan nan nan 1   1   nan nan nan nan 1   1   nan nan nan nan
  nan nan nan 1   1   1   nan nan nan nan 1   1   1   nan nan nan
  nan nan 1   1   2   1   1   1   1   1   1   2   1   1   nan nan
  nan 1   1   2   2   2   2   2   2   2   2   2   2   1   1   nan
  1   1   2   2   2   2   2   2   2   2   2   2   2   2   1   1
  1   1   2   2   2   2   2   2   2   2   2   2   2   2   1   1
  nan 1   1   2   2   2   2   2   2   2   2   2   2   1   1   nan
  nan nan 1   1   2   1   1   1   1   1   1   2   1   1   nan nan
  nan nan nan 1   1   1   nan nan nan nan 1   1   1   nan nan nan
  nan nan nan nan 1   1   nan nan nan nan 1   1   nan nan nan nan
  nan nan nan nan nan 1   nan nan nan nan 1   nan nan nan nan nan
  nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan
  nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan ];

switch ptr
  case 'SOM left'
    pd = fliplr(pd_r);
    set(hf, 'PointerShapeCData', pd, ...
            'PointerShapeHotSpot', [9, 1], ...
            'Pointer', 'custom');

  case 'SOM right'
    pd = pd_r;
    set(hf, 'PointerShapeCData', pd, ...
            'PointerShapeHotSpot', [9, 16], ...
            'Pointer', 'custom');

  case 'SOM bottom'
    pd = pd_r.';
    set(hf, 'PointerShapeCData', pd, ...
            'PointerShapeHotSpot', [16, 9], ...
            'Pointer', 'custom');

  case 'SOM top'
    pd = rot90(pd_r);
    set(hf, 'PointerShapeCData', pd, ...
            'PointerShapeHotSpot', [1, 9], ...
            'Pointer', 'custom');

  case 'SOM leftright'
    pd = pd_lr;
    set(hf, 'PointerShapeCData', pd, ...
            'PointerShapeHotSpot', [9, 9], ...
            'Pointer', 'custom');

  case 'SOM topbottom'
    pd = pd_lr.';
    set(hf, 'PointerShapeCData', pd, ...
            'PointerShapeHotSpot', [9, 9], ...
            'Pointer', 'custom');

  otherwise
    % Set it to the string passed in
    set(hf, 'Pointer', ptr);

end

end
