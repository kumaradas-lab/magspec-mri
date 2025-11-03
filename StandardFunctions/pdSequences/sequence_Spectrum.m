function [HW, mySave, SeqLoop, SliceSelect] = sequence_Spectrum(HW, mySave, Seq, SliceSelect)
%% Acquire an FID for spectroscopy
%
% The default parameters are set for water.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2021 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

if nargin < 4
  if exist(HW.Spectro.ReferenceFidPath, 'file')
    if HW.Spectro.useSliceSelect  % use slice gradient
      Seq.useSliceSelect = 1;
      SliceSelect = HW.Spectro.SliceSelect;
      Seq.SlicePulse = HW.Spectro.SlicePulse;
    else
      Seq.useSliceSelect = 0;
      SliceSelect = [];
      Seq.FlipPulse = HW.Spectro.FlipPulse;
    end
  end
end

if isemptyfield(Seq, 'loops'), Seq.loops = 1; end  % 1,2,3....
if isemptyfield(Seq, 'PreLoops'), Seq.PreLoops = 1; end  % 1,2,3....
if isemptyfield(Seq, 'plot'), Seq.plot = 1; end  % 1,2,3....
if isemptyfield(Seq, 'nEchos'), Seq.nEchos = 0; end  % 0,1,2,3....
if isemptyfield(Seq, 'tFID'), Seq.tFID = 1; end  % duration of FID acquisition in spectroscopy experiment
if isemptyfield(Seq, 'fSample'), Seq.fSample = HW.FindShim.fSample; end  % HW.RX.fSample = 125e6 Hz
if isemptyfield(Seq, 'AQFIDGrid'), Seq.AQFIDGrid = 1; end  % AQ Grid starts with center of excitation
if isemptyfield(Seq, 'RepetitionTime'), Seq.RepetitionTime = 10; end  % HW.RX.fSample = 125e6 Hz
if isemptyfield(Seq, 'Find_Frequency_interval'), Seq.Find_Frequency_interval = 0; end  % Find Frequency interval
Seq.nEchos = 0;          % 0,1,2,3....
Seq.tEcho = Seq.tFID*2;  % echo time (if Seq.AQFID = 2 -> FID AQ Time)
Seq.AQFID = 1;           % relativer Anteil der Echozeit/2 in dem Daten aufgenommen werden]0...1[
Seq.AQEcho = 0.5;        % relativer Anteil der Echozeit in dem Daten aufgenommen werden [0...1[, bei 0 wird 1 Sample aufgenommen

if isemptyfield(Seq, 'Spectro'), Seq.Spectro = []; end
if isemptyfield(Seq.Spectro, 'fOffsetTimeStart'), Seq.Spectro.fOffsetTimeStart = 0.03; end
if isemptyfield(Seq.Spectro, 'fOffsetTimeStop'), Seq.Spectro.fOffsetTimeStop = 0.09; end
if isemptyfield(Seq.Spectro, 'PhaseOffsetTimeStart'), Seq.Spectro.PhaseOffsetTimeStart = 0.03; end
if isemptyfield(Seq.Spectro, 'PhaseOffsetTimeStop'), Seq.Spectro.PhaseOffsetTimeStop = 0.04; end


if Seq.nEchos
  Seq.Spectro.PhaseOffsetTimeStart = Seq.Spectro.PhaseOffsetTimeStart + Seq.tEcho*Seq.nEchos;
  Seq.Spectro.PhaseOffsetTimeStop = Seq.Spectro.PhaseOffsetTimeStop + Seq.tEcho*Seq.nEchos;
  Seq.Spectro.fOffsetTimeStart = Seq.Spectro.fOffsetTimeStart + Seq.tEcho*Seq.nEchos;
  Seq.Spectro.fOffsetTimeStop = Seq.Spectro.fOffsetTimeStop + Seq.tEcho*Seq.nEchos;
end

iDevice = 1;  % FIXME: Support multiple MMRT devices
Channel = 1;  % FIXME: Add support for multiple AQ channels?


prevHoldShim = HW.Grad(iDevice).HoldShim;
prevFindFrequencyPause = HW.FindFrequencyPause;
hwGuard = onCleanup(@() protectHW(HW, iDevice, prevHoldShim, prevFindFrequencyPause));
HW.Grad(iDevice).HoldShim = 1;
HW.FindFrequencyPause = Seq.RepetitionTime;
% update approximate magnet frequency if 1000s have passed since last frequency update
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 1000);
% repeat frequency sweep with improved accuracy
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, Seq.Find_Frequency_interval, ...
  0, 1, HW.tFlip90Def*HW.FindFrequencySweep.tPulseDivider, 1, 1024);

if Seq.RepetitionTime<=30 && (Seq.loops+Seq.PreLoops>1)
  Seq.TimeToNextSequence = Seq.RepetitionTime;
end

for loop = (-Seq.PreLoops+1):Seq.loops
  if Seq.loops+Seq.PreLoops>1
    disp(['loops to run ' num2str(Seq.loops-loop+1) ' (' num2str((Seq.loops-loop+1)*Seq.RepetitionTime) ' s)' ]);
  end

  [data, SeqOut] = sequence_EchoStandard(HW, Seq, SliceSelect);
  tr = SeqOut.nEchos+1;
  Seq.plotSeq = [];
  Seq.plotSeqTR = [];

  if Seq.RepetitionTime <= 30
    Seq.Reinitialize = 0;
    Seq.IgnoreTimingError = 1;
    Seq.TimeFromLastSequence = SeqOut.RepetitionTime - sum(SeqOut.SequenceTime);
    Seq.TimeToNextSequence = Seq.TimeFromLastSequence;
  else
    Seq.StartSequenceTime = SeqOut.StartSequenceTime+SeqOut.RepetitionTime;
    HW.tRepInit = 0.5;
  end

  iAQ = find([SeqOut.AQ(:).Channel] == Channel & [SeqOut.AQ(:).Device] == iDevice, 1, 'first');

  if loop == 1
    SeqLoop = SeqOut;
    SeqLoop.loopdata.fft1_data = zeros(size(data(iAQ).fft1_data ,1), size(data(iAQ).fft1_data,2), size(data(iAQ).fft1_data,3), Seq.loops);
    SeqLoop.loopdata.data = zeros(size(data(iAQ).data,1), size(data(iAQ).data,2), size(data(iAQ).data,3), Seq.loops);
    SeqLoop.loopdata.time_all = zeros(size(data(iAQ).data,1), size(data(iAQ).data,2), size(data(iAQ).data,3), Seq.loops);
    SeqLoop.loopdata.phaseOffset = zeros(Seq.loops, 1);
    SeqLoop.loopdata.fOffset = zeros(Seq.loops, 1);
    SeqLoop.loopdata.StartTime = zeros(Seq.loops, 1);
    Seq.StartSequenceTimeOffset = SeqOut.StartSequenceTime;
  end
  if Seq.plot
    fh6 = figure(6);
    ax6(1) = subplot(2,1,1, 'Parent', fh6);
    plot(ax6(1),...
      repmat(data(iAQ).time_all(:,1,tr), 1, 3), ...
      [abs(data(iAQ).data(:,1,tr)), real(data(iAQ).data(:,1,tr)), imag(data(iAQ).data(:,1,tr))]*data(iAQ).Amplitude2Uin(1)*1e6);
    title(ax6(1), 'Acquired signal');
    ylabel(ax6(1), ['Amplitude in ' char(181) 'V']);
    xlabel(ax6(1), 'Time in s');
    ax6(2) = subplot(2,1,2, 'Parent', fh6);
    plot(ax6(2), ...
      ((data(iAQ).f_fft1_data(:,1,tr))/SeqOut.AQ(iAQ).Frequency(tr)-1)*1e6, ...
      abs(data(iAQ).fft1_data(:,1,tr)).*data(iAQ).Amplitude2Uin(1)*1e6);
    xlabel(ax6(2), 'Frequency in ppm');  % offset ppm
    ylabel(ax6(2), ['Amplitude in ' char(181) 'V']);
    title(ax6(2), 'FFT of acquired signal');
    xlim(ax6(2), [-10, 10]);
    set(ax6(2), 'XDir', 'reverse');
  end
  if loop >= 1
    SeqLoop.loopdata.StartTime(loop) = sum([SeqOut.StartSequenceTime, -Seq.StartSequenceTimeOffset]);
    SeqLoop.loopdata.fft1_data(:,:,:,loop) = data(iAQ).fft1_data;
    SeqLoop.loopdata.time_all(:,:,:,loop) = data(iAQ).time_all;
    SeqLoop.loopdata.data(:,:,:,loop) = data(iAQ).data;
    temp = squeeze(data(iAQ).data(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr));
    SeqLoop.loopdata.phaseOffset(loop) = ...
      get_MeanPhaseWeighted(temp(find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.PhaseOffsetTimeStart,1,'first')   :   find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.PhaseOffsetTimeStop,1,'first')));
    SeqLoop.loopdata.fOffset(loop) = ...
      get_MeanPhaseDiffWeighted(temp(find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStart,1,'first')   :   find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStop,1,'first')))*SeqOut.AQ(iAQ).fSample(1)/2/pi;
    HW.fLarmor = HW.fLarmor - double(SeqLoop.loopdata.fOffset(loop));
  else
    temp = squeeze(data(iAQ).data(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr));
    fOffset = ...
      get_MeanPhaseDiffWeighted(temp(find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStart,1,'first')   :   find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStop,1,'first')))*SeqOut.AQ(iAQ).fSample(1)/2/pi;
    HW.fLarmor = HW.fLarmor - double(fOffset);
  end
end

end

function protectHW(HW, iDevice, prevHoldShim, prevFindFrequencyPause)
%% restore previous HW settings when function finishes or is aborted

HW.Grad(iDevice).HoldShim = prevHoldShim;
HW.FindFrequencyPause = prevFindFrequencyPause;

end
