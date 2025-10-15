function [fLarmorOut, tEcho, df] = get_fLarmorFitTotEcho(HW, fLarmorIn, tEcho)
%% Calculate a Larmor frequency whose period has an integer multiple to tEcho
%
%     [fLarmorOut, tEcho, df] = get_fLarmorFitTotEcho(HW, fLarmorIn, tEcho)
%
% INPUT:
%   HW          HW structure or object
%   fLarmorIn   Larmor frequency in Hz. This frequency is changed to the closest
%               that has an integer multiple to the (fixed) Echo time (see
%               below).
%   tEcho       Echo time in seconds. This time is rounded to the closest grid
%               defined by the system frequency of the controller.
%
% OUTPUT:
%   fLarmorOut  Larmor frequency in Hz that fulfills the requirement.
%   tEcho       Echo time in seconds rounded to the system frequency grid.
%   df          Difference in Hz of the returned Larmor frequency to the input
%               Larmor frequency.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
%
% SEE ALSO: get_tEchoFitTofLarmor

%%
tGrid = 1/HW.TX.fSample;
fGrid = HW.TX.fSample/(2^HW.TX.DdsPicBits-1);

% fix tEcho to the system frequency
nGridEcho = round(tEcho/tGrid);
tEcho = nGridEcho*tGrid;

% fix fLarmor to tEcho
nPeriod = round(tEcho.*fLarmorIn);
fLarmorOut = nPeriod./tEcho;

% round to a frequency that can be synthesized
fLarmorOut = round(fLarmorOut/fGrid)*fGrid;

df = fLarmorOut - fLarmorIn;

end
