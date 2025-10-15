function [fLarmorOut, tEcho, df] = get_tEchoFitTofLarmor(HW, fLarmorIn, tEcho)
%% Calculate an Echo time that is an integer multiple of the Larmor period
%
%      [fLarmorOut, tEcho, df] = get_tEchoFitTofLarmor(HW, fLarmorIn, tEcho)
%
% INPUT:
%   HW          HW structure or object
%   fLarmorIn   Larmor frequency in Hz. This frequency is slightly changed to
%               find a tEcho that is an integer multiple of the (adapted) Larmor
%               period.
%   tEcho       Echo time in seconds. This time is an integer multiple of the
%               adapted Larmor frequency and lies on the grid defined by the
%               system frequency of the controller.
%               If tEcho is a scalar, it can be varied by +/- 1000/HW.TX.fSample
%               or 5% of its value (whichever is smaller).
%               Otherwise, the minimum and maximum value in tEcho are used to
%               define the span in which the Echo time might vary.
%
% OUTPUT:
%   fLarmorOut  Larmor frequency in Hz that fulfills the requirements.
%   tEcho       Echo time in seconds that fulfills the requirements.
%   df          Difference in Hz of the returned Larmor frequency to the input
%               Larmor frequency.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
%
% SEE ALSO: get_fLarmorFitTotEcho

%%

tGrid = 1/HW.TX.fSample;
fGrid = HW.TX.fSample/(2^HW.TX.DdsPicBits-1);

% define span of Echo times
if isscalar(tEcho)
  % either 5% of tEcho or 1000*tGrid
  steps = min([round(tEcho/tGrid/20), 1000]);
  dev_tGrid = (-steps:steps)*tGrid;
elseif numel(tEcho)
  % min and max define span
  mintEcho = min(tEcho(:));
  maxtEcho = max(tEcho(:));
  tEcho = mean([mintEcho, maxtEcho]);
  dev_tGrid = (round((mintEcho-tEcho)/tGrid) : round((maxtEcho-tEcho)/tGrid))*tGrid;
end

% round to a frequency that can be synthesized
fLarmorIn = round(fLarmorIn/fGrid)*fGrid;

% fix tEcho to fLarmor
tEcho_fLarmor = round(fLarmorIn*tEcho)./fLarmorIn;

% create vector of possible Echo times (system frequency)
tEcho_Grid = round(tEcho_fLarmor/tGrid)*tGrid + dev_tGrid;

% calculate Larmor frequencies that match each Echo time
fLarmor_tEcho = round(fLarmorIn.*tEcho_Grid)./tEcho_Grid;

% use combination for which the deviation to the desired Larmor frequency is smallest
dev_fLarmor = abs(fLarmorIn-round(fLarmor_tEcho/fGrid)*fGrid);
[~, ifLarmor] = min(dev_fLarmor);
tEcho = tEcho_Grid(ifLarmor);

% round to a frequency that can be synthesized
fLarmorOut = round(fLarmor_tEcho(ifLarmor)/fGrid)*fGrid;

df = fLarmorOut - fLarmorIn;

end
