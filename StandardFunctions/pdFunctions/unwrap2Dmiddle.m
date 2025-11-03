function unwrapped = unwrap2Dmiddle(tounwrap)
%% Unwrap a 2d phase matrix
%
%   unwrapped = unwrap2Dmiddle(tounwrap)
%
% The unwrapping is done in two steps: In the first step, all columns of the
% matrix are unwrapped independent from each other. In the second step, the
% center row is unwrapped. All other rows are unwrapped identically to the
% center row.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2020 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

u1 = unwrap(tounwrap, [], 1);
u1m = unwrap(u1(round(end/2),:)) - u1(round(end/2),:);
u2 = u1 + ones(size(u1,1),1)*u1m;
unwrapped = u2 - (u2(round(end/2),round(end/2)) - tounwrap(round(end/2),round(end/2)));

end