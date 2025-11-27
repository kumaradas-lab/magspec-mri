function AQ = get_PowerSyncFrequencyMPS(HW, AQ, Frequency, iDevice)
%% Get power synchronization frequency for DC-DC converter on MMRT
%
%     AQ = get_PowerSyncFrequencyMPS(HW, AQ, Frequency, iDevice)
%
% Harmonics of the power synchronization frequency of the DC-DC converter of the
% console are visible in the acquired signal. This function selects a power
% synchronization frequency such that these harmonics are least close (on
% average) to a selection of frequencies for each tRep in the pulse program.
%
%
% INPUT:
%
%   HW
%       Object of class PD.HWClass or HW structure.
%
%   AQ
%       AQ structure at least with the field AutoPowerSyncFrequency which
%       contains a logical vector (with the size of tRep) that indicates for
%       which tRep the frequency should be calculated.
%
%   Frequency
%       Matrix with the same number of columns as the number of tReps in the
%       pulse program. Each column contains the frequencies which should not be
%       disturbed by harmonics of the power synchronization frequency. If the
%       number of frequencies which should not be disturbed is different between
%       tReps, remaining elements must be filled with NaN.
%
%   iDevice
%       Index of the MMRT device which will run the acquisition (if multiple
%       MMRT are connected).
%       (Default: 1)
%
%
% OUTPUT:
%
%   AQ
%       Same as input AQ but with the selected power synchronization frequency
%       in the field "AQ.PowerSyncFrequency" for the tReps for which
%       "AQ.AutoPowerSyncFrequency" is set to true.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function isemptyfield

