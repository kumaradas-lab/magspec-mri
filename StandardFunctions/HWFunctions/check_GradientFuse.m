function [Seq, Grad] = check_GradientFuse(HW, Seq, Grad, iGrad, iDevice)
%% Check the sequence whether it might blow a fuse of the gradient amplifier
%
%   [Seq, Grad] = check_GradientFuse(HW, Seq, Grad, iGrad, iDevice)
%
% Used model:
% * Phase transition at current greater than HW.Grad.CoilMaxDcCurrent.
% * Fuse blows if I^2*t with I > HW.Grad.CoilMaxDcCurrent in sequence exceeds
%   HW.Grad.CoilCurrentSquareTime.
% * I^2*t dissipates with HW.Grad.CoilMaxDcCurrent^2*t.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2025 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

end
