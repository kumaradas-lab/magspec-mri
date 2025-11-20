function unwrapped = unwrap3Dmiddle(tounwrap)
%% Unwrap the first 3 non-singleton dimensions of phase matrix (up to 6d)
%
%   unwrapped = unwrap3Dmiddle(tounwrap)
%
% The unwrapping is done along the first three non-singleton dimensions. It aims
% to preserve the phase of the center voxel.
% It is done in three steps: In the first step, all columns of the (squeezed)
% matrix are unwrapped independent from each other. In the second step, all
% center rows are unwrapped. All other rows are unwrapped identically to their
% respective center rows. In the third step, all slices are unwrapped smoothly
% for a "ray" along the center of the slices.
%
% See also: unwrap2Dmiddle
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2021 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

[s(1), s(2), s(3), s(4), s(5), s(6)] = size(tounwrap);
% remove singleton dimensions
tounwrap = squeeze(tounwrap);
unwrapped12 = zeros(size(tounwrap));
unwrapped = zeros(size(tounwrap));
for n = 1:size(tounwrap, 3)
  unwrapped12(:,:,n) = unwrap2Dmiddle(tounwrap(:,:,n));
end
unwrapped3m = unwrap(squeeze(unwrapped12(floor(end/2)+1,floor(end/2)+1,:))) - ...
  squeeze(unwrapped12(floor(end/2)+1,floor(end/2)+1,:));
for n = 1:size(tounwrap, 3)
  unwrapped(:,:,n) = unwrapped12(:,:,n) + unwrapped3m(n);
end
% Set phase of center voxel
unwrapped = unwrapped + ...
  tounwrap(floor(end/2)+1,floor(end/2)+1,floor(end/2)+1) - ...
  unwrapped(floor(end/2)+1,floor(end/2)+1,floor(end/2)+1);
% undo squeeze, restore original singleton dimensions
unwrapped = reshape(unwrapped, s(1), s(2), s(3), s(4), s(5), s(6));

end
