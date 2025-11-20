function [HW, mySave, SeqOut] = Find_Frequency_FID(HW, mySave, minTime, iterations, tAQ, tAQStart)
%% Search the Larmor frequency and store it in HW
%
%   [HW, mySave, SeqOut] = Find_Frequency_FID(HW, mySave, minTime, iterations, tAQ, tAQStart)
%
% This function searches the MR frequency and stores it in HW and mySave.
% It is fast and can be adjusted to loose only a part of the magnetization M0
% per iteration  (using the default pulse angle divisor
% HW.FindFrequencyFID.tPulseDivider).
% It uses the FID after a (short) excitation pulse at the Larmor frequency
% corresponding to HW.FindFrequencyGamma to determine HW.B0 (in optionally
% several iterations):
% It calculates the offset frequency to the mixer frequency by the average
% gradient of the phase of the acquired signal weighted by its amplitude. For
% the optional next iterations, this frequency is used as the new frequency for
% the excitation pulse and mixer. The final frequency is used to calculate the
% magnitude of the B0 field (HW.B0) using HW.FindFrequencyGamma. The frequency
% for the default nucleus (HW.GammaDef) is saved in HW.fLarmor.
% Optionally, the B0 field can be adjusted instead of the Larmor frequency by
% using (optional) B0 coils. For this, HW.FindFrequencyFID.shiftB0 can be set to
% true (see below).
%
% After the frequency is determined, the execution is paused for
% HW.FindFrequencyPauseFID seconds.
%
%
% INPUT:
%
% All input parameters but HW are optional. If they are omitted or empty,
% default values are used.
%
%   HW
%       HW structure or object. For settings that affect this measurement, see
%       in particular the structure in HW.FindFrequencyFID.
%
%   mySave
%       mySave structure (necessary for minTime)
%
%   minTime
%       Minimal time in seconds since the last time Find_Frequency_Sweep was
%       executed before new value for the Larmor frequency is searched again.
%       (Default: HW.FindFrequencyFID.maxTime)
%
%   iterations
%       Number of iterations used for frequency determination.
%       (Default: HW.FindFrequencyFID.iterations)
%
%   tAQ
%       Acquisition time in seconds. (Default: HW.FindFrequencyFID.tAQ)
%
%   tAQStart
%       Acquisition start time in seconds. If set to 0, the minimum time is used
%       taking the rf pulse length and deadtime into account.
%       (Default: HW.FindFrequencyFID.tAQStart)
%
% Additionally, the behavior of this function can be changed with the following
% fields in the structure HW.FindFrequencyFID:
%
%   fSample
%       Sampling frequency of the acquisition window(s) in Hertz.
%
%   tPulseDivider
%       Divisor for the flip angle of the rf pulses. A value of 1 corresponds to
%       90 degrees flip angles for each rf pulse.
%
%   average
%       Number of averages.
%
%   shiftB0
%       Boolean value to select whether the default frequency (HW.fLarmor) is
%       adapted to match the current Larmor frequency of the magnet (false), or
%       (optional) B0 coils are used to shift the magnetic field to HW.B0Target
%       (true).
%
%
% OUTPUT:
%
%   HW
%       HW structure or object
%
%   mySave
%       mySave structure
%
%   SeqOut
%       Structure with sequence parameters and results.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if (nargin < 2),                        mySave     = []; end
if (nargin < 3) || isempty(minTime),    minTime    = HW.FindFrequencyFID.maxTime; end
if (nargin < 4) || isempty(iterations), iterations = HW.FindFrequencyFID.iterations; end
if (nargin < 5) || isempty(tAQ),        tAQ        = HW.FindFrequencyFID.tAQ; end
if (nargin < 6) || isempty(tAQStart),   tAQStart   = HW.FindFrequencyFID.tAQStart; end


%% initialization
% first time clear mySave
if isempty(mySave), mySave = struct(); end
if isemptyfield(mySave, 'lastTime'), mySave.lastTime = 0; end
if isemptyfield(mySave, 'HW'), mySave.HW = struct(); end
if isemptyfield(mySave.HW, 'B0'), mySave.HW.B0 = HW.B0; end

% FIXME: Handle measurements on multiple devices
iDevice = 1;

Dummy = false;
if ~isemptyfield(mySave, 'DummySerial')
  if mySave.DummySerial(min(iDevice, numel(mySave.DummySerial))) > 0
    Dummy = true;
    mySave.HW.B0 = HW.B0;
  end
end

SeqOut = struct();


%% actual measurement
% If time since last frequency search has exceeded "minTime", search again.
if ~Dummy && (now*24*3600-mySave.lastTime > minTime)
  fprintf('Searching frequency...\n');

  % define parameters
  fOffsetFID         = 0;                                  % init frequency offset var
  start_Frequency    = HW.B0 * HW.FindFrequencyGamma/2/pi; % start frequency search here
  mySave.HW.tRepInit = HW.tRepInit;                        % store state of HW.tRepInit

  for n = 1:iterations

    if isempty(HW.FindFrequencyFID.tPulse90)
      Seq.tPulse = HW.TX(iDevice).Amp2FlipPiIn1Sec/2 / ...
        HW.TX(iDevice).AmpDef / HW.FindFrequencyFID.tPulseDivider;
    else
      Seq.tPulse = HW.FindFrequencyFID.tPulse90 / HW.FindFrequencyFID.tPulseDivider;
    end

    % Seq
    if isemptyfield(HW.FindFrequencyFID, 'tRep')
      Seq.tRep = 10e-3+tAQ+tAQStart+Seq.tPulse;
    else
      Seq.tRep = HW.FindFrequencyFID.tRep;
    end
    Seq.CLTime = 50e-6;
    Seq.average = HW.FindFrequencyFID.average;

    % AQ
    AQ.fSample      = HW.RX(iDevice).fSample / round(HW.RX(iDevice).fSample/HW.FindFrequencyFID.fSample);
    AQ.Start        = max([Seq.tPulse + get_DeadTimeTX2RX(HW, AQ.fSample(1), iDevice), tAQStart]);
    AQ.nSamples     = max([round(tAQ * AQ.fSample), 2]);
    AQ.Frequency    = start_Frequency - fOffsetFID;
    AQ.Phase        = 0;
    AQ.ResetPhases  = 1;
    AQ.Device       = iDevice;

    % TX
    TX(1).Duration  = Seq.tPulse;
    TX(1).Start     = -Seq.tPulse/2;
    TX(1).Amplitude = HW.TX(iDevice).AmpDef;
    TX(1).Frequency = AQ.Frequency(1);
    TX(1).Phase     = 0;
    TX(1).Device    = iDevice;

    % no gradients
    [Grad(1:HW.Grad(iDevice).n).Time] = deal(NaN);
    [Grad(1:HW.Grad(iDevice).n).Amp]  = deal(0);

    % Plot pulse program
    Seq.plotSeq = [];
    Seq.plotSeqTR = [];

    % run measurement
    [~, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);
    % don't reinitialize the hardware after each loop
    if (HW.FindFrequencyPauseFID <= 2), HW.tRepInit = 0.05; end

    % plot measured data
    if HW.FindFrequencyFID.doPlot
      plot_data_1D(HW, data_1D);
    end

    [MeanPhaseDiffWeighted, phaseDiffWeightedStd] = ...
      get_MeanPhaseDiffWeighted(data(1).data(1:SeqOut.AQ(1).nSamples(1),1,1));
    % calculate frequency offset
    fOffsetFID = double(MeanPhaseDiffWeighted)*SeqOut.AQ(1).fSample(1)/2/pi + fOffsetFID;
    fOffsetFIDStd = phaseDiffWeightedStd * SeqOut.AQ(1).fSample(1)/2/pi;
    SeqOut.fLarmorStdev = fOffsetFIDStd;  % standard uncertainty of measurement

    if HW.FindFrequencyPauseFID >= 1
      fprintf('Waiting %2d seconds after finding Frequency... ', round(HW.FindFrequencyPauseFID));
      sleep(HW.FindFrequencyPauseFID);
      fprintf('done\n');
    elseif HW.FindFrequencyPauseFID > 0
      sleep(HW.FindFrequencyPauseFID);
    end

    % adjust frequency or B0 field amplitude
    Seq.fLarmor = start_Frequency - fOffsetFID;

    if ~isemptyfield(HW.FindFrequencyFID, 'shiftB0') ...
        && HW.FindFrequencyFID.shiftB0
      if isempty(HW.B0Target) || ~isfinite(HW.B0Target)
        warning('PD:Find_Frequency_FID:NoB0Target', ...
          ['HW.FindFrequencySweep.shiftB0 is set to true, but HW.B0Target is not set.\n', ...
          'B0 amplitude cannot be shifted without a target. Adjusting HW.fLarmor instead.']);
      else
        newB0shift = HW.Grad(iDevice).AmpOffset(4) + (HW.B0Target - HW.B0);
        newB0 = HW.B0Target;
        if abs(newB0shift) > HW.Grad(iDevice).HoldShimNormMax(4)*HW.Grad.MaxAmp(4)
          desiredB0shift = newB0shift;
          newB0shift = sign(newB0shift) * HW.Grad(iDevice).HoldShimNormMax(4) * HW.Grad(iDevice).MaxAmp(4);
          newB0 = newB0 - (desiredB0shift - newB0shift);
          warning('PD:Find_Frequency_FID:B0ShiftExceeded', ...
            ['The required B0 shift (%.1f %cT) is larger than the maximum shift (%.1f %cT).\n', ...
            'Consider setting a different HW.B0Target. Adjusting HW.fLarmor to stay in range.'], ...
            desiredB0shift*1e6, 181, HW.Grad(iDevice).HoldShimNormMax(4)*HW.Grad(iDevice).MaxAmp(4)*1e6, 181);
        end
        HW.Grad(iDevice).AmpOffset(4) = newB0shift;
        Seq.fLarmor = newB0 * HW.FindFrequencyGamma/(2*pi);
        fOffsetFID = Seq.fLarmor - start_Frequency;
      end
    end
  end

  % restore state of HW.tRepInit
  HW.tRepInit = mySave.HW.tRepInit;

  % store frequency, B0 field amplitude and shift to mySave structure
  mySave.HW.B0 = Seq.fLarmor * 2*pi/HW.FindFrequencyGamma;
  mySave.HW.fLarmor = mySave.HW.B0 * HW.GammaDef/2/pi;
  mySave.HW.B0shift = HW.Grad(iDevice).AmpOffset(4);
  mySave.lastTime = now*24*3600;
  fprintf('New frequency: %12.3f Hz (SD %6.3f Hz). Calibrated from FID\n', mySave.HW.fLarmor, fOffsetFIDStd);
end

% save the new frequency, B0 field amplitude and shift to HW parameters
HW.B0 = mySave.HW.B0;
HW.fLarmor = mySave.HW.B0 * HW.GammaDef/2/pi;
HW.Grad(iDevice).AmpOffset(4) = mySave.HW.B0shift;

end
