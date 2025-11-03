function [HW, mySave, SeqOut] = Find_Frequency_Sweep(HW, mySave, minTime, span, doPlot, tPulse90, nMeasurements, nSamples)
%% Search the Larmor frequency and store it in HW
%
%   [HW, mySave, SeqOut] = Find_Frequency_Sweep(HW, mySave, minTime, span, doPlot, tPulse90, nMeasurements, nSamples)
%
% This function searches the Larmor frequency around the currently set Larmor
% frequency (HW.B0*HW.FindFrequencyGamma/2/pi) and saves it in HW and mySave.
% It uses a train of several TX and AQ windows with increasing mixer frequency:
% The offset frequency to the mixer frequency is calculated at each point of the
% acquired signal from the gradient of its phase. The Larmor frequency is
% calculated from the AQ window with the lowest standard deviation for the
% offset frequency.
% This frequency is used to calculate the magnitude of the B0 field (HW.B0). The
% frequency for the default nucleus (HW.GammaDef) is saved in HW.fLarmor.
%
% After the frequency is determined, the execution is paused for
% HW.FindFrequencyPause seconds.
%
% INPUT:
% All input parameters but HW are optional. If they are omitted or empty, the
% default values as set in the structure HW.FindFrequencySweep are used.
%
%   HW        HW structure or object
%
%   mySave    mySave structure (necessary for minTime)
%
%   minTime   Minimum time in seconds since the last time Find_Frequency_Sweep
%             was executed before new value for the Larmor frequency is searched
%             again.
%
%   span      searching span of frequency around fLarmor
%             (HW.B0*HW.FindFrequencyGamma/2/pi)
%
%   doPlot    Plot data (bool) or handle to the parent figure or uipanel in
%             which the data is displayed.
%
%   tPulse90  duration of the used (90 degrees) pulse in seconds
%
%   nMeasurements
%             number of measurements (TX and AQ trains)
%
%   nSamples  number of samples in each measurement
%
% OUTPUT:
%
%   HW        HW structure or object
%
%   mySave    mySave structure
%
%   SeqOut    structure with sequence parameters and results
%
% ------------------------------------------------------------------------
% (C) Copyright 2012-2021 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------

%% default input
if (nargin < 2),                     mySave  = []; end
if (nargin < 3) || isempty(minTime), minTime = HW.FindFrequencySweep.maxTime; end
% search span of frequency around fLarmor (HW.B0*HW.FindFrequencyGamma/2/pi)
if (nargin < 4) || isempty(span),    span    = HW.FindFrequencySweep.span; end
if (nargin < 5) || isempty(doPlot),  doPlot  = HW.FindFrequencySweep.doPlot; end
if (nargin < 6) || isempty(tPulse90)
  if isempty(HW.FindFrequencySweep.tPulse90)
    tPulse90 = HW.tFlip90Def;
  else
    tPulse90 = HW.FindFrequencySweep.tPulse90;
  end
end
if (nargin < 7) || isempty(nMeasurements), nMeasurements = HW.FindFrequencySweep.nMeasurements; end
if (nargin < 8) || isempty(nSamples), nSamples = HW.FindFrequencySweep.nSamples; end

%% initialization
% first time clear mySave
if isempty(mySave), mySave = struct(); end
if isemptyfield(mySave, 'lastTime'), mySave.lastTime = 0; end
if isemptyfield(mySave, 'HW'), mySave.HW = struct(); end
if isemptyfield(mySave.HW, 'B0'), mySave.HW.B0 = HW.B0; end

Dummy = false;
if ~isemptyfield(mySave, 'DummySerial') && mySave.DummySerial > 0
  Dummy = true;
  mySave.HW.B0 = HW.B0;
end

SeqOut = struct();

% FIXME: Handle measurements on multiple devices
iDevice = 1;

%% actual measurement
% If time since last frequency search has exceeded "minTime", search again.
if ~Dummy && (now*24*3600-mySave.lastTime > minTime)

  fprintf('Searching frequency...\n');

  % define parameters
  Seq.fLarmor     = HW.B0 * HW.FindFrequencyGamma/2/pi;  % center frequency for fLarmor sweep
  Seq.nMeasurements = nMeasurements;                    % number of measurements
  Seq.nSamples    = nSamples;                           % number of samples per AQ windows
  Seq.tPuls       = tPulse90/HW.FindFrequencySweep.tPulseDivider;  % use small flip angle depending on "tPulse90"

  % Generate a sequence in the following order:
  % TX AQ TX AQ TX AQ ...
  % the transmitting and the mixing frequency is increased step by step

  % AQ
  AQ.fSample    = HW.FindFrequencySweep.fSample;        % sampling rate (should be twice as high as the frequency stepsize)
  AQ.Start      = cumsum( (Seq.tPuls+get_DeadTimeTX2RX(HW,AQ.fSample(1)) + ...
                           Seq.nSamples/AQ.fSample(1) + ...
                           get_DeadTimeRX2TX(HW, AQ.fSample(1))) * ones(Seq.nMeasurements, iDevice)) - ...
                  (Seq.nSamples/AQ.fSample(1) + get_DeadTimeTX2RX(HW, AQ.fSample(1), iDevice));
  AQ.nSamples   = max(5, Seq.nSamples);                 % each window has the same number of samples
  if Seq.nMeasurements == 1, span = 0; end
  % mixing Frequency is increased step by step (to match rf pulse frequency)
  AQ.Frequency  = linspace(-span/2, span/2, Seq.nMeasurements)' + Seq.fLarmor;
  AQ.Device     = iDevice;

  % TX
  TX(1).Duration  = Seq.tPuls * ones(Seq.nMeasurements, 1);
  TX(1).Start     = AQ.Start - get_DeadTimeTX2RX(HW, AQ.fSample(1), iDevice) - TX(1).Duration;
  TX(1).Amplitude = HW.TX.AmpDef;                       % default TX amplitude
  TX(1).Frequency = AQ.Frequency;                       % rf pulse frequency matches mixing frequency
  TX(1).Device    = iDevice;

  % Seq
  Seq.tRep = AQ.Start(end)+(AQ.nSamples+10)/AQ.fSample+0.001;

  % plot pulse program
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
  % frequency offset to mixing frequency (derived from difference of angle between datapoints)
  fOffsetFIDs = diff(unwrap(angle(data(1).data(1:SeqOut.AQ(1).nSamples(1),:,1)))) * ...
    SeqOut.AQ(1).fSample(1)/2/pi;

  [maxA, maxI] = max(abs(data(1).data), [], 1);
  maxA(maxA<max(maxA)/10) = NaN;
  highAmp = double(bsxfun(@le, maxA/4, abs(data(1).data)));
  tt = 1:numel(maxA);
  firstI = NaN(size(maxA));
  secondI = NaN(size(maxA));
  nhighAmp = NaN(size(maxA));
  fOffsetFIDsStd = NaN(size(maxA));
  for t = tt(~isnan(maxA));
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
    if nhighAmp(t) > 3
      % standard deviation of frequency offset of each AQ window
      fOffsetFIDsStd(t) = std(fOffsetFIDs(firstI(t):(secondI(t)-1),t), 0, 1);
    end
  end

  % find the frequency
  % find index of the window with the lowest standard deviation
  [fOffsetFIDsStdMinValue, fOffsetFIDsStdIndx] = min(fOffsetFIDsStd, [], 2);
  SeqOut.fLarmorStdev = fOffsetFIDsStdMinValue/nhighAmp(fOffsetFIDsStdIndx)^0.5;  % standard uncertainty of measurement
  if fOffsetFIDsStdMinValue/nhighAmp(fOffsetFIDsStdIndx)^0.5 < HW.FindFrequencySweep.fOffsetFIDsStdMaxValue
    % mean of offset frequency of each AQ window
    fOffsetFID = mean(fOffsetFIDs, 1) - SeqOut.AQ(1).Frequency' + Seq.fLarmor;
    Seq.fLarmor = Seq.fLarmor - double(fOffsetFID(fOffsetFIDsStdIndx)); % calculate new Larmor frequency
    % save the new frequency to mySave
    HW.B0 = Seq.fLarmor * 2*pi/HW.FindFrequencyGamma;  % calculate the new magnetic field strength
    mySave.HW.B0 = Seq.fLarmor * 2*pi/HW.FindFrequencyGamma;  % save B0 to mySave
    mySave.HW.fLarmor = mySave.HW.B0 * HW.GammaDef/2/pi;  % save fLarmor of GammaDef (actual nuclei) to mySave
    mySave.lastTime = now*24*3600;  % set lastTime to now in seconds

    fprintf('New frequency: %12.3f Hz (SD %6.3f Hz)\n', ...
      mySave.HW.fLarmor, SeqOut.fLarmorStdev);

    % add entry to log file in temp folder
    logFileDir = fullfile(tempdir(), 'Pure Devices', 'logs');
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
    HW.B0 = Seq.fLarmor * 2*pi/HW.FindFrequencyGamma;         % calculate the new magnet field
    mySave.HW.B0 = Seq.fLarmor * 2*pi/HW.FindFrequencyGamma;  % save B0 to mySave
    mySave.HW.fLarmor = mySave.HW.B0 * HW.GammaDef/2/pi;    % save fLarmor of GammaDef (actual nuclei) to mySave
    fprintf('Failed. Old frequency: %12.3f Hz\n', mySave.HW.fLarmor);
    warning('Find_Frequency_Sweep:fOffsetFIDsStdMaxValue', ...
      ['Please check the sample.\nPlease check HW.fLarmor.\n', ...
      'Standard deviation of offset frequency too high:\n', ...
      '%.3f Hz > HW.FindFrequencySweep.fOffsetFIDsStdMaxValue'], ...
      SeqOut.fLarmorStdev);
  end

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

%% save the new frequency to HW parameters
HW.B0 = mySave.HW.B0;
HW.fLarmor = mySave.HW.B0 * HW.GammaDef/2/pi;

end
