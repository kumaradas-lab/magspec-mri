function [HW, mySave, SeqOut] = Find_Frequency_Sweep(HW, mySave, minTime, span, doPlot, tPulse90, nMeasurements, nSamples)
%% Search the Larmor frequency and store it in HW
%
%   [HW, mySave, SeqOut] = Find_Frequency_Sweep(HW, mySave, minTime, span, doPlot, tPulse90, nMeasurements, nSamples)
%
% This function searches the Larmor frequency around the currently set B0 field
% amplitude (HW.B0) and saves it in HW and mySave.
% It uses a series of interleaved TX pulses and AQ windows with increasing
% (rf or mixer) frequency:
% The offset frequency to the mixer frequency is calculated in each of the
% acquisition windows from the gradient of the phase of the acquired signal. The
% updated Larmor frequency is selected from the AQ window with the lowest
% standard error for the determined offset frequency.
% This frequency is used to calculate the current magnitude of the B0 field
% (HW.B0). The frequency for the default nucleus (HW.GammaDef) is saved in
% HW.fLarmor.
% Optionally, the B0 field can be adjusted instead of the Larmor frequency by
% using (optional) B0 coils. For this, HW.FindFrequencySweep.shiftB0 can be set
% to true (see below).
%
% After the frequency is determined, the execution is paused for
% HW.FindFrequencyPause seconds.
%
%
% INPUT:
%
% All input parameters but HW are optional. If they are omitted or empty, the
% default values as set in the structure HW.FindFrequencySweep are used.
%
%   HW
%       HW structure or object
%
%   mySave
%       mySave structure (necessary for minTime)
%
%   minTime
%       Minimum time in seconds since the last time Find_Frequency_Sweep was
%       executed before new value for the Larmor frequency is searched again.
%       For this to work correctly, make sure to use "mySave" in the input and
%       output arguments of all respective functions.
%       (Default: HW.FindFrequencySweep.maxTime)
%
%   span
%       Frequency span in Hz around fLarmor (HW.B0*HW.FindFrequencyGamma/2/pi)
%       in which the the current Larmor frequency is scanned.
%       (Default: HW.FindFrequencySweep.span)
%
%   doPlot
%       Plot data (bool) or handle to the parent figure or uipanel in which the
%       data is displayed.
%       (Default: HW.FindFrequencySweep.doPlot)
%
%   tPulse90
%       Duration of a 90 degrees pulse at HW.TX.AmpDef in seconds. This value
%       essentially selects the flip angle for the rf pulses in the experiment.
%       That flip angle is divided by HW.FindFrequencySweep.tPulseDivider for
%       the pulses in the frequency sweep. The actually used pulse durations
%       might be increased to match the bandwidth of the acquisition windows.
%       (Default: HW.FindFrequencySweep.tPulse90
%                 or if that is empty: HW.tFlip90Def)
%
%   nMeasurements
%       Number of measurements (i.e., the number of TX pulses followed by an
%       acquisition).
%       (Default: HW.FindFrequencySweep.nMeasurements)
%
%   nSamples
%       Number of samples in each AQ window. The sampling frequency can be set
%       with HW.FindFrequencySweep.fSample.
%       (Default: HW.FindFrequencySweep.nSamples)
%
% Additionally, the behavior of this function can be changed with the following
% fields in the structure HW.FindFrequencySweep:
%
%   fSample
%       Sampling frequency of the acquisition window(s) in Hertz.
%
%   tPulseDivider
%       Divisor for the flip angle of the rf pulses. A value of 1 corresponds to
%       90 degrees flip angles for each rf pulse.
%
%   fOffsetFIDsStdMaxValue
%       If the standard deviation of the frequency offset of the acquisition
%       with the lowest standard devition is higher than this value in Hertz,
%       the measurement is discarded and the Larmor frequency is not updated.
%
%   shiftB0
%       Boolean value to select whether the default frequency (HW.fLarmor) is
%       adapted to match the current Larmor frequency of the magnet (false), or
%       (optional) B0 coils are used to shift the magnetic field to HW.B0Target
%       (true).
%
%   Seq
%       Structure with settings for the sequence. The function might override
%       some fields of this structure. If Seq.TimeToNextSequence is non-empty,
%       the execution is *not* paused by this function for HW.FindFrequencyPause
%       seconds (as described above). Instead, the caller is responsible to
%       setup the correct timing to potentially following experiments.
%       (Default: structure with no fields)
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
% ------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------


%% default input
if (nargin < 2)
  mySave = [];
end
if (nargin < 3) || isempty(minTime)
  minTime = HW.FindFrequencySweep.maxTime;
end
if (nargin < 4) || isempty(span)
  % search span of frequency around fLarmor (HW.B0*HW.FindFrequencyGamma/2/pi)
  span = HW.FindFrequencySweep.span;
end
if (nargin < 5) || isempty(doPlot)
  doPlot = HW.FindFrequencySweep.doPlot;
end
if (nargin < 6) || isempty(tPulse90)
  if isempty(HW.FindFrequencySweep.tPulse90)
    tPulse90 = HW.tFlip90Def;
  else
    tPulse90 = HW.FindFrequencySweep.tPulse90;
  end
end
if (nargin < 7) || isempty(nMeasurements)
  nMeasurements = HW.FindFrequencySweep.nMeasurements;
end
if (nargin < 8) || isempty(nSamples)
  nSamples = HW.FindFrequencySweep.nSamples;
end


%% initialization
% FIXME: Handle measurements on multiple devices
iDevice = 1;

% first time clear mySave
if isempty(mySave)
  mySave = struct();
end
if isemptyfield(mySave, 'lastTime')
  mySave.lastTime = 0;
end
if isemptyfield(mySave, 'HW')
  mySave.HW = struct();
end
if isemptyfield(mySave.HW, 'B0')
  mySave.HW.B0 = HW.B0;
end
if isemptyfield(mySave.HW, 'B0shift')
  mySave.HW.B0shift = HW.Grad(iDevice).AmpOffset(4);
end

Dummy = false;
if ~isemptyfield(mySave, 'DummySerial') ...
    && mySave.DummySerial(min(iDevice, numel(mySave.DummySerial))) > 0
  Dummy = true;
  mySave.HW.B0 = HW.B0;
end

SeqOut = struct();

if isemptyfield(HW.FindFrequencySweep, 'Seq')
  Seq = struct();
else
  Seq = HW.FindFrequencySweep.Seq;
end

Seq.plotSeq = 0;


%% actual measurement
% If time since last frequency search has exceeded "minTime", search again.
if ~Dummy && (now*24*3600-mySave.lastTime > minTime - eps(mySave.lastTime))

  fprintf('Searching frequency...\n');

  % define parameters
  Seq.fLarmor = HW.B0 * HW.FindFrequencyGamma/2/pi;  % center frequency for fLarmor sweep
  Seq.nMeasurements = nMeasurements;  % number of different frequency pulses in the sweep
  Seq.nSamples = nSamples;  % number of samples per AQ window
  Seq.tPulse = tPulse90 / HW.FindFrequencySweep.tPulseDivider;  % use small flip angle depending on "tPulse90"
  % selected flip angle in radians
  selectedFlipAngle = HW.FindFrequencyGamma * Seq.tPulse * HW.TX(iDevice).AmpDef;

  % Generate a sequence in the following order:
  % TX AQ TX AQ TX AQ ...
  % the transmitting and the mixing frequency is increased step by step

  % AQ bandwidth
  AQ.fSample = HW.FindFrequencySweep.fSample;  % sampling rate (should be twice as high as the frequency stepsize)
  if Seq.nMeasurements == 1, span = 0; end
  % mixing Frequency is increased step by step (to match rf pulse center frequency)
  AQ.Frequency = linspace(-span/2, span/2, Seq.nMeasurements)' + Seq.fLarmor;

  % Check for gaps in the frequency coverage.
  if Seq.nMeasurements > 1 && diff(AQ.Frequency(1:2))/2 > AQ.fSample/3
    warning('PD:Find_Frequency_Sweep:IncompleteCoverage', ...
      ['The acquisition bandwidth (%.3f kHz) is not high enough to cover ', ...
       'the selected frequency span (%.3f kHz) without gaps with the ', ...
       'selected number of measurements (%d).\n', ...
       'Consider reducing "span" or increasing "nMeasurements".'], ...
      AQ.fSample*1e-3, span*1e-3, Seq.nMeasurements);
  end

  % TX amplitude for the selected frequencies
  TX(1).Frequency = AQ.Frequency;  % rf pulse frequency matches mixing frequency
  TX(1).Channel = HW.TX(iDevice).ChannelDef;
  if isa(HW, 'PD.HWClass')
    % apply limits from HW structure at this frequency
    TX(1).Amplitude = get_TX_Amplitude(HW, 'Frequency', TX(1).Frequency, ...
      'Norm', HW.TX(iDevice).Def.Norm(TX(1).Channel), ...
      'Mmrt', HW.TX(iDevice).Def.Mmrt(TX(1).Channel), ...
      'Uout', HW.TX(iDevice).Def.Uout(TX(1).Channel), ...
      'PaUout', HW.TX(iDevice).Def.PaUout(TX(1).Channel), ...
      'Amplitude', HW.TX(iDevice).Def.Amplitude(TX(1).Channel), ...
      'Channel', TX(1).Channel);
  else
    % fall back to amplitude at "primary" frequency
    TX(1).Amplitude = HW.TX(iDevice).AmpDef * ones(Seq.nMeasurements, 1);
  end

  % Adjust pulse duration to keep the previously selected flip angle.
  Seq.tPulse = selectedFlipAngle ./ TX(1).Amplitude / HW.FindFrequencyGamma;

  % If necessary, reduce bandwidth of TX pulse to match bandwidth of AQ window.
  Seq.tPulse = max(Seq.tPulse, 1/AQ.fSample);

  % That might require a lower TX pulse amplitude to keep the previously
  % selected flip angle.
  TX(1).Amplitude = selectedFlipAngle ./ Seq.tPulse / HW.FindFrequencyGamma;


  % remaining properties for AQ windows
  AQ.Start = cumsum(Seq.tPulse + get_DeadTimeTX2RX(HW, AQ.fSample(1)) + ...
                    Seq.nSamples/AQ.fSample(1) + ...
                    get_DeadTimeRX2TX(HW, AQ.fSample(1))) - ...
             (Seq.nSamples/AQ.fSample(1) + get_DeadTimeTX2RX(HW, AQ.fSample(1), iDevice));
  AQ.nSamples = max(5, Seq.nSamples);  % each window has the same number of samples
  AQ.Device = iDevice;

  % remaining properties for TX pulses
  TX(1).Duration = Seq.tPulse;
  TX(1).Start = AQ.Start - get_DeadTimeTX2RX(HW, AQ.fSample(1), iDevice) - TX(1).Duration;
  TX(1).Device = iDevice;


  % Seq
  Seq.tRep = AQ.Start(end) + AQ.nSamples/AQ.fSample + ...
    max(10 / AQ.fSample, ...
        HW.RX(iDevice).ClampCoil.Enable * HW.RX(iDevice).ClampCoil.tPostset) + 0.001;

  % don't plot pulse program
  Seq.plotSeq = [];
  Seq.plotSeqTR = [];

  % no gradients
  [Grad(1:sum([HW.Grad(:).n])).Time] = deal(NaN);
  [Grad(1:sum([HW.Grad(:).n])).Amp]  = deal(0);

  % run measurement
  [~, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);

  % plot acquired measurement data
  if ~(isnumeric(doPlot) && doPlot == 1) && ...  % tread "1" as Boolean true (instead of figure(1))
      (ishghandle(doPlot, 'figure') || ishghandle(doPlot, 'uipanel'))
    conf.hParent = doPlot;
    plot_data_1D(HW, data_1D, [], true, conf);
  elseif doPlot
    plot_data_1D(HW, data_1D);
    drawnow();
  end


  %% evaluation
  % frequency offset to mixing frequency (derived from phase of NMR signal)

  % fOffsetFIDsSamples = diff(unwrap(angle(double(data(1).data(1:SeqOut.AQ(1).nSamples(1),:,1))))) * ...
  %   SeqOut.AQ(1).fSample(1)/2/pi;

  [maxA, maxI] = max(abs(data(1).data), [], 1);

  % FIXME: It would be nice if we could distinguish the "systematic"
  %        contribution to the standard error of the phase slope from the
  %        "stochastic" contribution. The former could be caused by a NMR signal
  %        that "shifts" in frequency during acquisition (e.g., due to chemical
  %        composition or B0 inhomogeneity). The latter should be dominated by
  %        the noise level.
  %        Maybe, try to estimate the noise level, compare it to the signal
  %        amplitude, and estimate the error due to that "calculated" SNR?

  % rmsCIC = sqrt(mean(data.cic_corr(:,1).^(-2)));  % rms power loss of CIC filter
  % PnoiseRinDensity = HW.Constant.Boltzmann * HW.TempR50 * 1 * rmsCIC;  % rms power / Hz
  % UnoiseRinDensity = sqrt(PnoiseRinDensity * HW.RX(iDevice).Rin);  % rms voltage / sqrt(Hz)
  %
  % ampNoiseR50 = UnoiseRinDensity * 2^0.5 * sqrt(SeqOut.AQ(1).fSample(1)) ...
  %   / (data(1).Amplitude2Uin(1)/HW.RX(iDevice).LNAGain);
  %
  % highAmp = double(bsxfun(@le, max(maxA/4, ampNoiseR50*10), abs(data(1).data)));

  highAmp = double(bsxfun(@le, maxA/4, abs(data(1).data)));
  tt = 1:numel(maxA);
  firstI = NaN(size(maxA));
  secondI = NaN(size(maxA));
  nhighAmp = NaN(size(maxA));
  fOffsetFIDs = NaN(size(maxA));
  fOffsetFIDsStd = NaN(size(maxA));
  for t = tt(~isnan(maxA))
    % Try to avoid "zero-crossing" of the (complex) signal. That could happen
    % for samples with two (or more) main frequency contributors. At a
    % "zero-crossing", the phase jumps by pi. The mean slope across such a jump
    % doesn' give the correct frequency offset.
    % Take only the samples after the first maximum until they reach 1/4 of the
    % maximum amplitude (to stay away from the "zero-crossing").
    if maxI(t) > 1
      firstI(t) = find(highAmp(1:maxI(t),t), 1, 'first');
    else
      firstI(t) = 1;
    end
    secondIt = find(~highAmp(maxI(t):end,t), 1, 'first') + maxI(t) - 1;
    if isempty(secondIt)
      secondI(t) = size(highAmp, 1);
    else
      secondI(t) = secondIt;
    end
    nhighAmp(t) = secondI(t) - firstI(t) + 1;

    % calculate weighted phase slope for selected samples
    [phaseDiffWeightedMean, phaseDiffWeightedStd] = ...
      get_MeanPhaseDiffWeighted(double(data(1).data(firstI(t):(secondI(t)-1),t,1)), 1, 'omitnan');

    fOffsetFIDs(t) = phaseDiffWeightedMean * SeqOut.AQ(1).fSample(1) / (2*pi);
    fOffsetFIDsStd(t) = phaseDiffWeightedStd * SeqOut.AQ(1).fSample(1) / (2*pi);

    % alternative: estimate standard error by standard deviation of frequency
    % offset calculated from each pair of phase values (diff).
    % if nhighAmp(t) > 3
    %   % standard deviation of frequency offset of each AQ window
    %   fOffsetFIDsStd(t) = std(fOffsetFIDsSamples(firstI(t):(secondI(t)-1),t), 0, 1)/nhighAmp(t)^0.5;
    % end

  end

  % Do not use signals that are close to the Nyquist frequency.
  % Depending on T1, these signals could have pretty high amplitude.
  fOffsetFIDs(abs(fOffsetFIDs) > AQ.fSample/3) = NaN;
  % discard very low amplitudes (compared to the highest remaining signal)
  fOffsetFIDs(maxA < max(maxA(~isnan(fOffsetFIDs)))/10) = NaN;

  fOffsetFIDsStd(isnan(fOffsetFIDs)) = NaN;

  % find the frequency
  % find index of the window with the lowest standard deviation
  [~, fOffsetFIDsStdIndx] = min(fOffsetFIDsStd, [], 2);
  fOffsetFIDsStdMinValue = fOffsetFIDsStd(fOffsetFIDsStdIndx);
  SeqOut.fLarmorStdev = fOffsetFIDsStdMinValue;  % standard uncertainty of measurement
  % mean of offset frequency of each AQ window
  fOffsetFID = fOffsetFIDs - SeqOut.AQ(1).Frequency' + Seq.fLarmor;
  newfLarmor = Seq.fLarmor - double(fOffsetFID(fOffsetFIDsStdIndx)); % calculate new Larmor frequency
  if fOffsetFIDsStdMinValue < HW.FindFrequencySweep.fOffsetFIDsStdMaxValue ...
      && newfLarmor > 0
    % save the new frequency to mySave
    Seq.fLarmor = newfLarmor;
    HW.B0 = Seq.fLarmor * 2*pi/HW.FindFrequencyGamma;  % calculate the new magnetic field strength
    if ~isemptyfield(HW.FindFrequencySweep, 'shiftB0') ...
        && HW.FindFrequencySweep.shiftB0
      if isempty(HW.B0Target) || ~isfinite(HW.B0Target)
        warning('PD:Find_Frequency_Sweep:NoB0Target', ...
          ['HW.FindFrequencySweep.shiftB0 is set to true, but HW.B0Target is not set.\n', ...
          'B0 amplitude cannot be shifted without a target. Adjusting HW.fLarmor instead.']);
      else
        newB0shift = HW.Grad(iDevice).AmpOffset(4) + (HW.B0Target - HW.B0);
        newB0 = HW.B0Target;
        if abs(newB0shift) > HW.Grad(iDevice).HoldShimNormMax(4)*HW.Grad(iDevice).MaxAmp(4)
          desiredB0shift = newB0shift;
          newB0shift = sign(newB0shift) * HW.Grad(iDevice).HoldShimNormMax(4) * HW.Grad(iDevice).MaxAmp(4);
          newB0 = newB0 - (desiredB0shift - newB0shift);
          warning('PD:Find_Frequency_Sweep:B0ShiftExceeded', ...
            ['The required B0 shift (%.1f %cT) is larger than the maximum shift (%.1f %cT).\n', ...
            'Consider setting a different HW.B0Target. Adjusting HW.fLarmor to stay in range.'], ...
            desiredB0shift*1e6, 181, HW.Grad(iDevice).HoldShimNormMax(4)*HW.Grad(iDevice).MaxAmp(4)*1e6, 181);
        end
        HW.Grad(iDevice).AmpOffset(4) = newB0shift;
        HW.B0 = newB0;
        Seq.fLarmor = newB0 * HW.FindFrequencyGamma/(2*pi);
      end
    end
    mySave.HW.B0 = Seq.fLarmor * 2*pi/HW.FindFrequencyGamma;  % save B0 to mySave
    mySave.HW.fLarmor = mySave.HW.B0 * HW.GammaDef/2/pi;  % save fLarmor of GammaDef (actual nuclei) to mySave
    mySave.HW.B0shift = HW.Grad(iDevice).AmpOffset(4);  % save B0 shift to mySave
    mySave.lastTime = now*24*3600;  % set lastTime to now in seconds

    fprintf('New frequency: %12.3f Hz (SD %6.3f Hz)\n', ...
      mySave.HW.fLarmor, SeqOut.fLarmorStdev);

    % add entry to log file in temp folder
    if ispc()
       % Windows
      logFileDir = fullfile(getenv('LOCALAPPDATA'), 'Pure Devices', 'OpenMatlab', 'logs');
    elseif ismac()
      % macOS
      logFileDir = fullfile('~', 'Libary', 'Logs', 'Pure Devices', 'OpenMatlab');
    else
      % Linux
      logFileDir = fullfile('~', 'Pure Devices', 'OpenMatlab', 'logs');
    end
    if ~isdir(logFileDir)
      mkdir(logFileDir); pause(0.1);
    else
      logFiles = dir(fullfile(logFileDir, '*_fLarmor.log'));
      if numel(logFiles) > 9
        % delete oldest files
        [~, iDate] = sort([logFiles.datenum], 'descend');
        delLogFiles = fullfile(logFileDir, {logFiles(iDate(10:end)).name});
        delete(delLogFiles{:});
      end
    end
    logFilePath = fullfile(logFileDir, sprintf('%s_fLarmor.log', datestr(now, 'yyyymmdd')));
    if isemptyfield(HW.FindFrequencySweep, 'fidLog') || ...
        ~strcmp(fopen(HW.FindFrequencySweep.fidLog), logFilePath)
      % keep the file open during a session
      % (i.e. until the hook in HW is deleted)
      fidLog = fopen(logFilePath, 'a');
      HW.FindFrequencySweep.fidLogHook = onCleanup(@() fclose(fidLog));
      HW.FindFrequencySweep.fidLog = fidLog;
      fprintf(fidLog, 'time [Matlab], B0 [T], fLarmor [Hz], std. uncertainty [Hz]\n');
    end

    fprintf(HW.FindFrequencySweep.fidLog, '%.13f, %.6f, %.1f, %.3f\n', ...
      now, HW.B0, Seq.fLarmor, SeqOut.fLarmorStdev);

  else
    HW.B0 = Seq.fLarmor * 2*pi/HW.FindFrequencyGamma;  % calculate the new magnet field
    mySave.HW.B0 = Seq.fLarmor * 2*pi/HW.FindFrequencyGamma;  % save B0 to mySave
    mySave.HW.fLarmor = mySave.HW.B0 * HW.GammaDef/2/pi;  % save fLarmor of GammaDef (actual nuclei) to mySave
    fprintf('Failed. Old frequency: %12.3f Hz\n', mySave.HW.fLarmor);
    if newfLarmor < 0
      warning('Find_Frequency_Sweep:negativeFrequency', ...
        ['Please check the sample.\nPlease check HW.fLarmor.\n', ...
        'Detected Larmor frequency looks unreasonable:\n', ...
        '%.0f Hz'], ...
        newfLarmor);
    else
      warning('Find_Frequency_Sweep:fOffsetFIDsStdMaxValue', ...
        ['Please check the sample.\nPlease check HW.fLarmor.\n', ...
        'Standard deviation of offset frequency too high:\n', ...
        '%.3f Hz > HW.FindFrequencySweep.fOffsetFIDsStdMaxValue'], ...
        SeqOut.fLarmorStdev);
    end
  end

  if isemptyfield(Seq, 'TimeToNextSequence')
    % wait for HW.FindFrequencyPause; useful if using samples with a high T1 time (e.g. water)
    % gives the sample time to relax before starting the next measurement
    if HW.FindFrequencyPause >= 1
      fprintf('Waiting %2d seconds after finding Frequency... ', round(HW.FindFrequencyPause));
    end
    if HW.FindFrequencyPause > 0
      sleep(HW.FindFrequencyPause);
    end
    if HW.FindFrequencyPause >= 1
      fprintf('done\n');
    end
  end
end


%% save the new frequency to HW parameters
HW.B0 = mySave.HW.B0;
HW.fLarmor = mySave.HW.B0 * HW.GammaDef/2/pi;
HW.Grad(iDevice).AmpOffset(4) = mySave.HW.B0shift;


end
