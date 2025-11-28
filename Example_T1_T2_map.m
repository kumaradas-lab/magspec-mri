%% Acquire a T1-T2-map of a sample
%
% Saturation or inversion pulse with increasing delay to a following
% CPMG echo train.

%% preparation
LoadSystem;  % load system parameters (reset to defaults: HW Seq AQ TX Grad)

[HW, mySave] = Find_Frequency_Sweep(HW, mySave,1);  % find magnet frequency


%% measurement settings
Seq.T2Estimated = 0.15;  % estimated T2 time of the sample in seconds
Seq.T1Estimated = 0.15;  % estimated T1 time of the sample in seconds

Seq.Recovery = 'Inversion';  % recovery type: 'Saturation' or 'Inversion'
Seq.SeqAverage.average = 2;  % set to 2 for phase cycling

Seq.plotSeq = 1:3;  % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)

HW.Grad(1).Slice.channel = HW.Grad(1).xyzB(2);  % gradient to select slice in y-direction
Seq.thicknessSlice = 5e-3;  % thickness of slice selected by continuous gradient in m
Seq.Slice.offAfterPreparation = 1;  % switch gradient off after preparation pulse

Seq.fitExpCorr.RingFilter = 1;  % use sliding average over 2 echoes in CPMG echo train

% actual measurement
[data, SeqOut, mySave] = sequence_RecoveryCPMG(HW, Seq, mySave);


%% evaluation
if (SeqOut.nTau1 > 2) && (SeqOut.nEcho > 1)
  [data_iL2D, SeqOutiL] = get_iLaplace2D(data(1), SeqOut);
end

%% -----------------------------------------------------------------------------
% (C) Copyright 2020-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
