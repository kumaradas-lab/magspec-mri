function [t2, T2, data, SeqOut] = get_T2(HW, SeqIn)
%% Measure T2 spectrum of sample using an Echo train
%
%       [t2, T2, data, SeqOut] = get_T2(HW, Seq)
%
% This function acquires an Echo train (CPMG) and evaluates the T2 time of a
% sample.
%
% INPUT:
%   HW      HW structure or object.
%   Seq     A structure with the following optional fields. If the fields are
%           omitted or empty, default values are used:
%     T2Estimated
%             Estimated T2 value of the sample in seconds (default: 0.1)
%     T2EstimatedMin
%             Shortest estimated T2 value of a sample with a T2-spectrum in
%             seconds (default: Seq.T2Estimated)
%     T2EstimatedMax
%             Longest estimated T2 value of a sample with a T2-spectrum in
%             seconds (default: Seq.T2Estimated)
%     tEcho   Echo time in seconds (default: Seq.T2EstimatedMin/10 but not
%             shorter than 0.5e-3 and not larger than  4e-3).
%     tEchoTrain
%             Duration of echo train in seconds (default: Seq.T2EstimatedMax*5)
%     preparationPulse
%             Handle to a pulse shape function for the 90 degrees flip pulse
%             (default: @Pulse_Rect).
%     inversionPulse
%             Handle to a pulse shape function for the 180 degrees inversion
%             pulse (default: @Pulse_Rect).
%     plot    Boolean value that enables plotting the results (default: false).
%     ConsoleOut
%             Boolean value. If true, the result of the T2 measurement is
%             displayed in the "Command Window" (default: false).
%     tAQEcho Duration of the acquisition window at the Echo times after the
%             inversion pulses in seconds where 0 corresponds to a single
%             acquisition sample and negative values correspond to no
%             acquisition
%             (default: max(1e-6, min(100e-6, Seq.tEcho-20e-6-HW.tFlip180Def)) )
%     fSample Sampling frequency of the acquisition window in Hz
%             (default: min(1000e3, max(20e3, 50/Seq.tAQEcho)) )
%     FitT2   Boolean value. If true, a (double) exponential fit to the data is
%             performed (default: true).
%     RingFilter
%             Boolean value. If true, the amplitude of two subsequent Echoes is
%             averaged to a "virtual" Echo at the mean time (default: true).
%     fitExp  Structure with the settings for fit_exp. For default values read
%             the help for that function.
%     RawData (deprecated)
%             Boolean value. (default: ~Seq.plot).
%     CalculateFFTOfData
%             Boolean value. If false, the FFT of the acquired data is NOT
%             calculated. That means that the CIC correction is also NOT
%             performed. On the other hand, the processing is faster in this
%             mode. (Default: ~Seq.RawData)
%     SeqAverage
%             Structure with the following optional fields. If the fields are
%             omitted or empty, default values are used:
%       average Number of averages (default: 2).
%       RandomTXRXPhaseOffset
%               Boolean value. If true, a random phase is added to all TX pulses
%               and acquisition windows (default: true).
%       TXPhaseInvertIncrement
%               Phase angle in degrees that is added to the invertion pulses at
%               each increment. The default value of 180 is useful to suppress
%               the FID of the inversion pulse.
%       SaveSeqAverageData
%               Boolean value. If true, the data of each average is saved in
%               data.SeqAverage.data (default: true).
%     SubtractEmptyMeasurement
%             Boolean value. If true, the measurement is repeated without the
%             sample to subtract the signal of the probe head (default: false).
%             This setting is handy when measuring samples with T2-values that
%             compare to the T2-values of the probe head materials (approx.
%             below a few milliseconds).
%     tAQFID  Duration of the acquisition window after the excitation pulse in
%             seconds (default: Seq.tEcho/8). If this value is larger than
%             Seq.tEcho/4, a separate measurement is performed to acquire just
%             the FID (followed by the Echo train as defined above). If the time
%             is set to 0, a single sample of the FID is acquired. For negative
%             values, the acquisition of the FID is omitted.
%     AQFIDMean
%             Boolean value. If true, the average of the acquired FID is used
%             for the T2 fit. Otherwise, each sample of the FID is used
%             (default: false).
%
% OUTPUT:
%   t2      Single exponential T2 time (if Seq.FitT2 is true. NaN otherwise.)
%   T2      Structure with results of the exponential fit. See: fit_exp. The
%           structure with the results for a potential separate FID sequence
%           (see tAQFID above) may be stored in T2.FID,
%   data    2x1 cell of structures with measured data with - among others - the
%           following fields (all as after the corrections for the exponential
%           fit, see fit_exp):
%     data{1}
%       DataAmplitudeEchoTrain
%               The complex amplitude data of the Echo train.
%       DataTimeEchoTrain
%               The time of the Echo train.
%       DataAmplitude
%               The complex amplitude of the Echo train potentially interleaved
%               with the data from the separate FID measurement.
%       DataTime
%               The time of the Echo train potentially interleaved with the data
%               from the separate FID measurement.
%     data{2}:
%             Structure with data of the potential separate FID measurement.
%   SeqOut  2x1 cell of structure with actually used measurement parameters (for
%           the Echo train and potentially for the separate FID measurement).
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% Default parameters
if nargin == 1, SeqIn.plot=0; end
if iscell(SeqIn), Seq = SeqIn{1}; else Seq = SeqIn; end
Seq = set_EmptyField(Seq, 'T2Estimated',      0.1);             % estimated T2, e.g. 100e-3
Seq = set_EmptyField(Seq, 'T2EstimatedMin',   Seq.T2Estimated); % shortest estimated T2, e.g. 100e-3
Seq = set_EmptyField(Seq, 'T2EstimatedMax',   Seq.T2Estimated); % longest estimated T2, e.g. 100e-3
Seq = set_EmptyField(Seq, 'tEcho',            min(max(0.5e-3, Seq.T2EstimatedMin/10), 4e-3));  % echo time, e.g. 10e-3           (T2/10)
Seq = set_EmptyField(Seq, 'tEchoTrain',       Seq.T2EstimatedMax*5);  % duration of echo train, e.g. 1  (T2*10)
Seq = set_EmptyField(Seq, 'preparationPulse', @Pulse_Rect);     % pulse shape used for flip pulse
Seq.FlipPulse = Seq.preparationPulse;  % FIXME: Remove this line when sequence_Recovery can be used.
Seq = set_EmptyField(Seq, 'inversionPulse',   @Pulse_Rect);     % pulse shape used for inversion
Seq.InvertPulse = Seq.inversionPulse;  % FIXME: Remove this line when sequence_Recovery can be used.
Seq = set_EmptyField(Seq, 'plot',             0);               % plot data
Seq = set_EmptyField(Seq, 'ConsoleOut',       0);               % display results in console
Seq = set_EmptyField(Seq, 'tAQEcho', max(1e-6, min(100e-6, Seq.tEcho-20e-6-HW.tFlip180Def)));  % duration of the acquisition between the inversions pulses [0...x] s (0 => only one sample, negative no sample)

Seq = set_EmptyField(Seq, 'fSample', min(1000e3, max(20e3, 50/Seq.tAQEcho)));  % Set the AQ sampling frequency
Seq = set_EmptyField(Seq, 'tAQFID',           -1);              % length of acquisition after excitation
Seq = set_EmptyField(Seq, 'AQFIDMean',        0);

Seq = set_EmptyField(Seq, 'FitT2',            1);               % FitT2 0 or 1
Seq = set_EmptyField(Seq, 'RingFilter',       1);               % RingFilter 0 or 1
if ~isfield(Seq, 'fitExp'),           Seq.fitExp      = []; end
Seq.fitExp = set_EmptyField(Seq.fitExp, 'hParent', Seq.plot);
Seq.fitExp = set_EmptyField(Seq.fitExp, 'CorrectFrequencyOffset', 0);
Seq.fitExp = set_EmptyField(Seq.fitExp, 'CorrectPhaseOffset', 1);
Seq.fitExp = set_EmptyField(Seq.fitExp, 'EndOffset', 0);
Seq.fitExp = set_EmptyField(Seq.fitExp, 'RingFilter', Seq.RingFilter);
Seq.fitExp = set_EmptyField(Seq.fitExp, 'hasFID', Seq.tAQFID >= 0);
Seq.fitExp = set_EmptyField(Seq.fitExp, 'FIDMean', 0);          % Use average of FID. If false, each sample of the FID is used for the fit.

Seq = set_EmptyField(Seq, 'RawData',            ~Seq.plot);     % faster processing
Seq = set_EmptyField(Seq, 'CalculateFFTOfData', ~Seq.RawData);  % faster processing
if ~isfield(Seq, 'SeqAverage'), Seq.SeqAverage = []; end        % use average
Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'average',                2);   % average number
Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'RandomTXRXPhaseOffset',  1);   %
Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'TXPhaseInvertIncrement', 180); % Invertpulse Phase increment average to average
Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'SaveSeqAverageData',     1);   % save average data

Seq = set_EmptyField(Seq, 'SubtractEmptyMeasurement', 0);       % Run the same measurement without sample to subtract signal caused by probe head.


if iscell(SeqIn)
  if numel(SeqIn) < 2
    SeqFID = SeqIn{1};
  else
    SeqFID = SeqIn{2};
  end
else
  SeqFID = SeqIn;
end


% relative part of the time tEcho/2 to be acquired between the first pulse and
% the inversion pulse [0...1[ (0 => only one sample, negative no sample)
Seq.AQFID = Seq.tAQFID / (Seq.tEcho/2);
if Seq.AQFID > 0.5
  % use two experiments for acquiring the FID and the Echo train
  Seq.separateFID = true;
  SeqFID.tEcho = 4*Seq.tAQFID;
  SeqFID.AQFID = 0.5;
  % pre-shots
  [~, ~] = sequence_EchoStandard(HW, SeqFID);
  [~, ~] = sequence_EchoStandard(HW, SeqFID);
  % FID and sparse Echoes
  [data{2}, SeqOut{2}] = sequence_EchoStandard(HW, SeqFID);
  Seq.AQFID = 0.5;
else
  % use one single experiment for acquiring FID and Echo train
  Seq.separateFID = false;
  data{2}.data = [];
  data{2}.time_all = [];
end
[data{1}, SeqOut{1}] = sequence_EchoStandard(HW, Seq);
% data.dataFID = dataFID.data;
% data.time_allFID = dataFID.time_all;
if SeqOut{1}.SubtractEmptyMeasurement
  disp('Please remove sample tube and press enter.')
  beep();
  pause();
  if Seq.separateFID
    [dataEmpty{2}] = sequence_EchoStandard(HW, SeqFID);
  else
    dataEmpty{2}.data = [];
  end
  [dataEmpty{1}] = sequence_EchoStandard(HW, Seq);
  % data with desired spacing between Echoes
  data{1}.dataFull = data{1}.data;
  data{1}.dataEmpty = dataEmpty{1}.data;
  data{1}.data = data{1}.data - dataEmpty{1}.data;
  % data for separate FID (optional)
  data{2}.dataFull = data{2}.data;
  data{2}.dataEmpty = dataEmpty{2}.data;
  data{2}.data = data{2}.dataFull - data{2}.dataEmpty;
end

%% (Double) exponential fit
t2 = NaN; T2 = NaN;
if SeqOut{1}.FitT2
  if SeqOut{1}.nEchos <= 2 % check for sufficient number of echoes
    warning('get_T2: Number of Echoes too low for accurate fitting. Must be > 2.');
     return
  end
  if Seq.separateFID
    % FIXME: Instead of fitting, the corrected FID should be "attached" to the
    % corrected data of the Echo train.
    SeqOut{2}.fitExp.hParent = SeqOut{2}.plot;
    [data{2}, SeqOut{2}, T2FID] = fitT2(data{2}, SeqOut{2});
    if isemptyfield(SeqOut{1}.fitExp, 'hParent') || ...
        SeqOut{1}.fitExp.hParent == SeqOut{2}.fitExp.hParent || ...
        (isnumeric(SeqOut{1}.fitExp.hParent) && SeqOut{1}.fitExp.hParent == 1)
      SeqOut{1}.fitExp.hParent = get(SeqOut{2}.fitExp.hParent, 'Number') + 1;
    end
  else
    data{2}.DataAmplitude = [];
    data{2}.DataTime = [];
    T2FID = NaN;
  end
  iAQ = SeqOut{1}.AQ(1).Device;
  [data{1}, SeqOut{1}, T2] = fitT2(data{1}(iAQ), SeqOut{1});
  T2.FID = T2FID;
  data{1}(iAQ).DataAmplitudeEchoTrain = data{1}(iAQ).DataAmplitude;
  data{1}(iAQ).DataTimeEchoTrain = data{1}(iAQ).DataTime;
  data{1}(iAQ).DataAmplitude = [data{2}(iAQ).DataAmplitude; data{1}(iAQ).DataAmplitudeEchoTrain];
  data{1}(iAQ).DataTime = [data{2}(iAQ).DataTime; data{1}(iAQ).DataTimeEchoTrain];
  [data{1}(iAQ).DataTime, iSort] = sort(data{1}(iAQ).DataTime);
  data{1}(iAQ).DataAmplitude = data{1}(iAQ).DataAmplitude(iSort);
  t2 = T2.tau;
end

end


function [data, SeqOut, T2] = fitT2(data, SeqOut)
%% Fit T2 with single and double exp functions

% fit exponential functions to the data
[T2, SeqOut.fitExp] = fit_exp(data.data, ...
  SeqOut.tEcho - SeqOut.tEchoFirst + data.time_all, ...
  SeqOut.fitExp);
if SeqOut.ConsoleOut
  fprintf('T2   = %10.1f ms\n', T2.tau*1000);
  if SeqOut.fitExp.DoubleExp
    fprintf('T2_1 = %10.1f ms; weighting factor = %5.1f%%\n', T2.xminDouble(3)*1000, T2.tau1w*100);
    fprintf('T2_2 = %10.1f ms; weighting factor = %5.1f%%\n', T2.xminDouble(5)*1000, T2.tau2w*100);
  end
end
data.DataAmplitude = T2.dataPhaseCorrected;
data.DataTime = T2.timeCorrected;
SeqOut.iLaplace1D.Problem = 'Decay';

end
