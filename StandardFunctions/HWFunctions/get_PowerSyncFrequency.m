function AQ = get_PowerSyncFrequency(HW, AQ, iDevice)
%% Get power synchronization frequency for DC-DC converter on MMRT
%
%     AQ = get_PowerSyncFrequency(HW, AQ, iDevice)
%
% Harmonics of the power synchronization frequency of the DC-DC converter of the
% console are visible in the MR spectrum. This function selects a power
% synchronization frequency such that these harmonics are moved away from the
% center of the acquisition spectrum. The (approximate) relative distance of the
% closest harmonic to the center of the acquisition spectrum is set with
% HW.MMRT.fPowerSyncDist.
%
%
% INPUT:
%
%   HW
%       Object of class PD.HWClass or HW structure.
%
%   AQ
%       AQ structure.
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
% (C) Copyright 2024-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function isemptyfield

