function slowset(handle, prop, value, increment)
% SLOWSET(H, 'PROPERTY', VALUE, INCREMENT)
%
% Like SET, except that the property is set against the original  value over
% several steps such that the value morphs from the starting value to the end
% value.
%
% H can be a vector of handles, but they must all accept PROPERTY.
%
% Only one property is supported at this time.
%
% Optional fourth argument INCREMENT specifies how many steps of animation to
% use.

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

if nargin == 4
  if increment == 0
    increment = 1;
  end
else
  increment = 10;
end

proplist = fieldnames(get(handle(1)));
tprop = proplist(strncmpi(prop, proplist, length(prop)));
prop = tprop{1};

hp = [];

for i = 1:length(handle)
  hp(i).handle = handle(i);
  hp(i).start = get(hp(i).handle, prop);
  hp(i).end = value;
  if isnumeric(hp(i).end) && isnumeric(hp(i).start)
    hp(i).values = VectorCalc(hp(i), increment);
  else
    set(hp(i).handle, prop, value);
    hp(i).values = [];
  end
end

for inc = 1:increment
  for i = 1:length(handle)
    if ~isempty(hp(i).values)
      newval = reshape(hp(i).values(inc,:,:,:), ...
                       size(hp(i).start,1), ...
                       size(hp(i).start,2));

      set(hp(i).handle, prop, newval);
    end
  end
  pause(.05)
end

end

function values =  VectorCalc(hp, increment)
% Do nothing but go to end value.

values = ones(increment, size(hp.end, 1), size(hp.end, 2), size(hp.end, 3));

for c = 1:numel(hp.end)
  newval = linspace(hp.start(c), hp.end(c), increment);
  values(:,c) = newval';
end

values = reshape(values, increment, size(hp.end, 1), size(hp.end, 2), ...
                 size(hp.end, 3));
end
