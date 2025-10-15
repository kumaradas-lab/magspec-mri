function [HW, mySave, SeqLoop, SliceSelect] = sequence_Spectrum(HW, mySave, Seq, SliceSelect)
%% Acquire an FID for spectroscopy
%
% The default parameters are set for water.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2023 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

if nargin < 4
  if exist(HW.Spectro.ReferenceFidPath, 'file') ...
      && HW.Spectro.useSliceSelect  % use slice gradient
    Seq.useSliceSelect = 1;
    SliceSelect = HW.Spectro.SliceSelect;
    Seq.SlicePulse = HW.Spectro.SlicePulse;
  else
    Seq.useSliceSelect = 0;
    SliceSelect = [];
    Seq.FlipPulse = HW.Spectro.FlipPulse;
  end
end

if isemptyfield(Seq, 'Loops')
  % number of loops for averages
  if isemptyfield(Seq, 'loops')
    Seq.Loops = 1;
  else
    % copy from legacy spelling "loops"
    Seq.Loops = Seq.loops;
  end
end
if isemptyfield(Seq, 'PreLoops')
  % number of pre-shots that are discarded in averages
  Seq.PreLoops = 1;
end
if isemptyfield(Seq, 'plot')
  % show plots with data of each single average step
  Seq.plot = 1;
end
if isemptyfield(Seq, 'plotAverage')
  % show plots with averaged data at each step
  Seq.plotAverage = Seq.plot && (Seq.Loops>1);
end
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

if Seq.RepetitionTime<=30 && (Seq.Loops+Seq.PreLoops>1)
  Seq.TimeToNextSequence = Seq.RepetitionTime;
end

for loop = (-Seq.PreLoops+1):Seq.Loops
  if Seq.Loops + Seq.PreLoops > 1
    fprintf('loops to run %d (%.3f s)\n', Seq.Loops-loop+1, (Seq.Loops-loop+1)*Seq.RepetitionTime);
  end

  [data, SeqOut] = sequence_EchoStandard(HW, Seq, SliceSelect);
  tr = SeqOut.nEchos+1;
  Seq.plotSeq = [];
  Seq.plotSeqTR = [];

  if Seq.RepetitionTime <= 30
    Seq.Reinitialize = 0;
    Seq.IgnoreTimingError = 1;
    Seq.TimeFromLastSequence = SeqOut.RepetitionTime - sum(SeqOut.SequenceTime) + SeqOut.CLTime(1) + SeqOut.tOffset(1);
    Seq.TimeToNextSequence = Seq.TimeFromLastSequence;
  else
    Seq.StartSequenceTime = SeqOut.StartSequenceTime+SeqOut.RepetitionTime;
    HW.tRepInit = 0.5;
  end

  iAQ = find([SeqOut.AQ(:).Channel] == Channel & [SeqOut.AQ(:).Device] == iDevice, 1, 'first');
  iData = find([data(:).channel] == Channel & [data(:).device] == iDevice);

  if loop == 1
    SeqLoop = SeqOut;
    SeqLoop.loopdata.fft1_data = zeros(...
      size(data(iData(1)).fft1_data, 1), ...
      size(data(iData(1)).fft1_data, 2), ...
      size(data(iData(1)).fft1_data, 3), ...
      Seq.Loops, ...
      numel(iData));
    SeqLoop.loopdata.data = SeqLoop.loopdata.fft1_data;
    SeqLoop.loopdata.time_all = SeqLoop.loopdata.fft1_data;
    for iAQX = 1:numel(iData)
      SeqLoop.loopdata.f_fft1_data(:,1,tr,1,iAQX) = data(iData(iAQX)).f_fft1_data(:,1,tr);
      SeqLoop.loopdata.Amplitude2Uin(iAQX) = data(iData(iAQX)).Amplitude2Uin(1);
    end
    SeqLoop.loopdata.phaseOffset = zeros(Seq.Loops, 1);
    SeqLoop.loopdata.fOffset = zeros(Seq.Loops, 1);
    SeqLoop.loopdata.StartTime = zeros(Seq.Loops, 1);

    Seq.StartSequenceTimeOffset = SeqOut.StartSequenceTime;
  end

  if Seq.plot
    for iAQX = 1:numel(iData)
      h_fig_data = figure(5+iAQX);
      h_ax_data(1) = subplot(2,1,1, 'Parent', h_fig_data);
      plot(h_ax_data(1),...
        repmat(data(iData(iAQX)).time_all(:,1,tr), 1, 3), ...
        [abs(data(iData(iAQX)).data(:,1,tr)), real(data(iData(iAQX)).data(:,1,tr)), imag(data(iData(iAQX)).data(:,1,tr))]*data(iAQ).Amplitude2Uin(1)*1e6);
      title(h_ax_data(1), 'Acquired signal');
      ylabel(h_ax_data(1), ['Amplitude in ' char(181) 'V']);
      xlabel(h_ax_data(1), 'Time in s');
      grid(h_ax_data(1), 'on');
      h_ax_data(2) = subplot(2,1,2, 'Parent', h_fig_data);
      if iAQX == 1
        freq = SeqOut.AQ(iAQ).Frequency(tr);
      else
        freq = SeqOut.AQ(iAQ).FrequencyX(tr);
      end
      plot(h_ax_data(2), ...
        ((data(iData(iAQX)).f_fft1_data(:,1,tr))/freq-1)*1e6, ...
        abs(data(iData(iAQX)).fft1_data(:,1,tr)).*data(iData(iAQX)).Amplitude2Uin(1)*1e6);
      xlabel(h_ax_data(2), 'Frequency in ppm');  % offset ppm
      ylabel(h_ax_data(2), ['Amplitude in ' char(181) 'V']);
      title(h_ax_data(2), 'FFT of acquired signal');
      xlim(h_ax_data(2), [-10, 10]);
      set(h_ax_data(2), 'XDir', 'reverse');
      grid(h_ax_data(2), 'on');
      drawnow();
    end
  end

  iAQ_ref = 1;  % use GammaDef signal for frequency drift
  if loop >= 1
    SeqLoop.loopdata.StartTime(loop) = sum([SeqOut.StartSequenceTime, -Seq.StartSequenceTimeOffset]);
    for iAQX = 1:numel(iData)
      SeqLoop.loopdata.fft1_data(:,:,:,loop,iAQX) = data(iData(iAQX)).fft1_data;
      SeqLoop.loopdata.time_all(:,:,:,loop,iAQX) = data(iData(iAQX)).time_all;
      SeqLoop.loopdata.data(:,:,:,loop,iAQX) = data(iData(iAQX)).data;
    end
    temp = squeeze(data(iData(iAQ_ref)).data(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr));
    SeqLoop.loopdata.phaseOffset(loop) = ...
      get_MeanPhaseWeighted(temp(find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.PhaseOffsetTimeStart,1,'first')   :   find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.PhaseOffsetTimeStop,1,'first')));
    SeqLoop.loopdata.fOffset(loop) = ...
      get_MeanPhaseDiffWeighted(temp(find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStart,1,'first')   :   find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStop,1,'first')))*SeqOut.AQ(iAQ).fSample(1)/2/pi;
    HW.fLarmor = HW.fLarmor - double(SeqLoop.loopdata.fOffset(loop));

    if Seq.plotAverage
      for iAQX = 1:numel(iData)
        h_fig_avg_data = figure(15+iAQX); clf(h_fig_avg_data);
        h_ax_avg_data(1) = subplot(2,1,1, 'Parent', h_fig_avg_data);
        % naive avaraging (without independent frequency shifts)
        avgData = mean(SeqLoop.loopdata.data(:,:,:,1:loop,iAQX), 4);
        plot(h_ax_avg_data(1),...
          repmat(data(iData(iAQX)).time_all(:,1,tr), 1, 3), ...
          [abs(avgData(:,1,tr)), real(avgData(:,1,tr)), imag(avgData(:,1,tr))]*data(iAQ).Amplitude2Uin(1)*1e6);
        title(h_ax_avg_data(1), 'Averaged signal');
        ylabel(h_ax_avg_data(1), ['Amplitude in ' char(181) 'V']);
        xlabel(h_ax_avg_data(1), 'Time in s');
        grid(h_ax_avg_data(1), 'on');
        h_ax_avg_data(2) = subplot(2,1,2, 'Parent', h_fig_avg_data);
        avgData_fft = fftshift(ifft(avgData(:,1,tr), [], 1));
        if iAQX == 1
          freq = SeqOut.AQ(iAQ).Frequency(tr);
        else
          freq = SeqOut.AQ(iAQ).FrequencyX(tr);
        end
        plot(h_ax_avg_data(2), ...
          ((data(iData(iAQX)).f_fft1_data(:,1,tr))/freq-1)*1e6, ...
          abs(avgData_fft).*data(iData(iAQX)).Amplitude2Uin(1)*1e6);
        xlabel(h_ax_avg_data(2), 'Frequency in ppm');  % offset ppm
        ylabel(h_ax_avg_data(2), ['Amplitude in ' char(181) 'V']);
        title(h_ax_avg_data(2), 'FFT of averaged signal');
        xlim(h_ax_avg_data(2), [-10, 10]);
        set(h_ax_avg_data(2), 'XDir', 'reverse');
        grid(h_ax_avg_data(2), 'on');
        drawnow();
      end
    end

  else
    temp = squeeze(data(iData(iAQ_ref)).data(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr));
    fOffset = ...
      get_MeanPhaseDiffWeighted(temp(find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStart,1,'first')   :   find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStop,1,'first')))*SeqOut.AQ(iAQ).fSample(1)/2/pi;
    HW.fLarmor = HW.fLarmor - double(fOffset);
  end
end

end

function protectHW(HW, iDevice, prevHoldShim, prevFindFrequencyPause)
%% restore previous HW settings when function finishes or is aborted

HW.Grad(iDevice).HoldShim = prevHoldShim;
HW.FindFrequencyPause = prevFindFrequencyPause;

end
