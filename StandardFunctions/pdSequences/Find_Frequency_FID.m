function [HW, mySave] = Find_Frequency_FID(HW, mySave, minTime, iterations, tAQ, tAQStart)
%% Search the Larmor frequency and store it in HW
%
%   [HW, mySave] = Find_Frequency_FID(HW, mySave, minTime, iterations, tAQ)
%
% This function searches the MR frequency and stores it in HW and mySave.
% It is fast and looses only about 10% of the magnetization M0 per iteration
% (using the default pulse angle divisor HW.FindFrequencyFID.tPulseDivider).
% It uses the FID after a short excitation pulse at the Larmor frequency
% corresponding to HW.FindFrequencyGamma to determine HW.B0 (in optionally
% several iterations):
% It calculates the offset frequency to the mixer frequency by the average
% gradient of the phase of the acquired signal weighed by its amplitude. For the
% optional next iterations, this frequency is used as the new frequency for the
% excitation pulse and mixer. The final frequency is used to calculate the
% magnitude of the B0 field (HW.B0) using HW.FindFrequencyGamma. The frequency
% for the default nucleus (HW.GammaDef) is saved in HW.fLarmor.
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
%   HW        HW structure or object. See in particular the structure in
%             HW.FindFrequencyFID.
%
%   mySave    mySave structure (necessary for minTime)
%
%   minTime   Minimal time in seconds since the last time Find_Frequency_Sweep
%             was executed before new value for the Larmor frequency is searched
%             again. (default: 1000)
%
%   iterations
%             number of iterations used for freqeuency determination (default:
%             1)
%
%   tAQ       acquisition time in seconds (default: 0.5e-3)
%
%   tAQStart  acquisition start time in seconds (default: 0 (minimum))
%
%
% OUTPUT:
%
%   HW        HW structure or object
%
%   mySave    mySave structure
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2020 Pure Devices GmbH, Wuerzburg, Germany
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

Dummy = false;
if ~isemptyfield(mySave, 'DummySerial')
  if mySave.DummySerial > 0
    Dummy = true;
    mySave.HW.B0 = HW.B0;
  end
end

% FIXME: Handle measurements on multiple devices
iDevice = 1;

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
    if HW.FindFrequencyPauseFID >= 1
      fprintf('Waiting %2d seconds after finding Frequency... ', round(HW.FindFrequencyPauseFID));
      sleep(HW.FindFrequencyPauseFID);
      fprintf('done\n');
    elseif HW.FindFrequencyPauseFID > 0
      sleep(HW.FindFrequencyPauseFID);
    end
  end

  % restore state of HW.tRepInit
  HW.tRepInit = mySave.HW.tRepInit;

  % store frequency and B0 field
  Seq.fLarmor = start_Frequency-fOffsetFID;
  HW.B0 = Seq.fLarmor * 2*pi/HW.FindFrequencyGamma;

  % store frequency and B0 field to mysave
  mySave.HW.B0 = Seq.fLarmor * 2*pi/HW.FindFrequencyGamma;
  mySave.HW.fLarmor = mySave.HW.B0 * HW.GammaDef/2/pi;
  mySave.lastTime = now*24*3600;
  fprintf('New frequency: %12.3f Hz (SD %6.3f Hz). Calibrated from FID\n', mySave.HW.fLarmor, fOffsetFIDStd);
end

% save the new frequency to HW parameters
HW.B0 = mySave.HW.B0;
HW.fLarmor = mySave.HW.B0 * HW.GammaDef/2/pi;

end
