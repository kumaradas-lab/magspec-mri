function [HW, mySave] = Find_PulseDurationFID(HW, mySave, minTime, doPlot, iterations, tPulse90, T1)
%% Find 90 degrees pulse duration and store corresponding calibration value in file
%
%   [HW, mySave] = Find_PulseDuration(HW, mySave, minTime, doPlot, iterations, tPulse90, T1)
%
% This function searches the correct pulse length for a 90 degrees pulse (i.e.
% the coil efficiency) and stores it in HW and mySave. Additionally, it appends
% the corresponding factor for output voltage on the TX channel to B1 field
% strength to the PaUout2AmplitudeCal.m file in the User folder.
% For this, it iteratively measures the amplitude of an FID after a pulse with
% increasing pulse length.
%
% INPUT:
%
% All input parameters but HW are optional. If they are omitted or empty,
% default values are used.
% Do not use talker as first argument.
%
%   HW        HW structure or object
%
%   mySave    mySave structure (necessary for minTime)
%
%   minTime   Minimum time in seconds since the last time Find_PulseDurationFID
%             was executed before a new value for the coil efficiency is
%             determined again (default: 1000).
%
%   doPlot    Plot sequence and data (figure number). If 0, don't plot.
%             (Default: 320)
%
%   iterations
%             Number of iterations used for frequency determination.
%             (Default: 10)
%
%   tPulse90  Duration of the used (90 degrees) pulse in seconds (default:
%             HW.tFlip90Def)
%
%   T1        Approximate T1 time in seconds of the sample (default: 0.2)
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2025 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------


%% input check
if nargin < 1
  error('PD:Find_PulseDurationFID:NoInput', 'At least input argument HW is mandatory.');
end
if nargin < 2, mySave = []; end
if nargin < 3 || isempty(minTime), minTime = 1000; end
if nargin < 4 || isempty(doPlot), doPlot = 320; end
if nargin < 5 || isempty(iterations), iterations = 10; end
if nargin < 6 || isempty(tPulse90), tPulse90 = HW.tFlip90Def; end
if nargin < 7 || isempty(T1), T1 = 0.2; end


%% measurement and evaluation
if iterations > 0
  [HW, mySave] = Find_Frequency_Sweep(HW, mySave, minTime, [], 1);  % find magnet frequency
  t = 0;

  Seq.RepetitionTime = max(300e-3, 3*T1);

  minP90 = tPulse90/3;
  maxP90 = tPulse90*2;
  nStepsP90 = iterations;
  nSamples = 5;
  flipAngle = NaN(nStepsP90, 1);
  maxFID = NaN(nStepsP90, 1);
  dataFID = NaN(1, 1);
  timeFID = NaN(1, 1);

  if doPlot
    hf = figure(doPlot);
    clf(hf, 'reset')

    ax(1) = subplot(2,1,1, 'Parent', hf);
    plot(ax(1), ...
      timeFID*1e3, abs(dataFID)*1e6, ...
      timeFID*1e3, real(dataFID)*1e6, ...
      timeFID*1e3, imag(dataFID)*1e6);
    hold(ax(1), 'on');
    title(ax(1), 'Acquired signal');
    ylabel(ax(1), sprintf('Amplitude in %cV', 181));
    xlabel(ax(1), 'Time in ms');
    ylim(ax(1), [0, Inf]);
    grid(ax(1), 'on');

    ax(2) = subplot(2,1,2, 'Parent', hf);
    plot(ax(2), flipAngle*1e6, maxFID*1e6, 'x');
    hold(ax(2), 'on');
    ylabel(ax(2), sprintf('Amplitude in %cV', 181));
    xlabel(ax(2), sprintf('Pulse duration in %cs', 181));
    ylim(ax(2), [0, Inf]);
    grid(ax(2), 'on');
  end

  iDevice = 1;  % FIXME: Support multiple MMRT devices

  for p90 = linspace(minP90, maxP90, nStepsP90)
    t = t+1;
    % Parameters used for timing calculations
    Seq.p90         = p90;                          % duration of 1st TX pulse
    Seq.plotSeq     = [];                           % plot sequence off

    % Sequence parameters
    Seq.tRep        = 50e-3;                        % repetition time
    Seq.tOffset     = maxP90;                       % offset of tRep

    % RF transmission parameters
    TX.Start        = -Seq.p90/2;                   % start time of rf-pulse
    TX.Duration     = Seq.p90;                      % duration of rf-pulse
    TX.Frequency    = HW.fLarmor;                   % frequency of rf-pulse
    TX.Phase        = 0;                            % phase of rf-pulse
    % TX.Amplitude    = 9.9 * HW.TX.PaUout2Amplitude(2);

    % Acquisition parameters
    AQ.fSample      = 20e3;                         % sample rate of AQ window
    % AQ.Start        = 10e-6;                        % acquisition start time
    AQ.Start        = maxP90/2 + get_DeadTimeTX2RX(HW, AQ.fSample(1)) + 1/AQ.fSample + 0e-6;  % acquisition start time
    AQ.nSamples     = 500;                          % number of samples in AQ window
    AQ.Frequency    = HW.fLarmor;                   % frequency of AQ window
    AQ.Phase        = 0;                            % phase of AQ window
    % AQ.Gain         = HW.RX(iDevice).Amplitude2Uin / 20e-3;  % maximum input voltage
    % AQ.Gain         = HW.RX(iDevice).Amplitude2Uin / 500e-3;  % maximum input voltage

    [Grad(1:HW.Grad(iDevice).n).Time] = deal([]);   % time the values in Grad(t).Amp with the corresponding index are set
    [Grad(1:HW.Grad(iDevice).n).Amp] = deal([]);    % amplitude of the gradient at the time Grad(t).Time (linearly interpolated)
    [Grad(1:HW.Grad(iDevice).n).Shim] = deal([]);   % additional shim; magnet shim is already considered in HW.MagnetShim. Caution: Using high values over a long time will damage the gradient coils and amplifiers!
    [Grad(1:HW.Grad(iDevice).n).Repeat] = deal([]); % if there is a one in the array the gradients of the prior TR are used. (reduces data traffic between the PC and the MRI device)% Start measurement

    Seq.TimeToNextSequence = Seq.RepetitionTime;

    [~, SeqOut, data, ~] = set_sequence(HW, Seq, AQ, TX, Grad);

    Seq.TimeToNextSequence = Seq.RepetitionTime - SeqOut.SequenceTime;
    Seq.TimeFromLastSequence = SeqOut.TimeToNextSequence;
    Seq.Reinitialize = 0;

    Channel = 1;  % FIXME: Add support for multiple acquisition channels?
    iAQ = find([SeqOut.AQ(:).Channel] == Channel & [SeqOut.AQ(:).Device] == iDevice, 1, 'first');

    dataFID = data(iAQ).data * data(iAQ).Amplitude2Uin(1)/HW.RX(iDevice).LNAGain;
    maxFID(t) = mean(abs(dataFID(1:nSamples)));
    flipAngle(t) = p90;
    if doPlot
      % Plot results
      % plot_data_1D(HW,data_1D);

      % figure(10);
      % plot(data(iAQ).time_of_tRep,angle((data(iAQ).data)*data(iAQ).Amplitude2Uin(1)));
      % title('Phase of acquired signal');
      % ylabel('phase in rad');
      % xlabel('Time in s');

      % figure(11);
      % plot(data(iAQ).f_fft1_data,(abs(data(iAQ).fft1_data)./data(iAQ).cic_corr.*data(iAQ).Amplitude2Uin(1)));  %without filter correction
      % % plot(data(iAQ).f_fft1_data,abs(data(iAQ).fft1_data).*data(iAQ).Amplitude2Uin(1));               %with filter correction
      % xlim([data(iAQ).f_fft1_data(1) data(iAQ).f_fft1_data(end)])
      % title('FFT of acquired signal');
      % xlabel('Frequency in Hz');

      % figure(12);
      % % plot(data(iAQ).f_fft1_data-SeqOut.AQ.Frequency(1),abs(data(iAQ).fft1_data).*data(iAQ).Amplitude2Uin(1)); xlabel('Frequency in Hz');  % offset frequency
      % plot((data(iAQ).f_fft1_data/SeqOut.AQ.Frequency(1)-1)*1e6,abs(data(iAQ).fft1_data).*data(iAQ).Amplitude2Uin(1)); xlabel('Frequency in ppm'); % offset ppm
      % title('FFT of acquired signal');

      timeFID = data.time_of_tRep;
      plot(ax(1), ...
        timeFID*1e3, abs(dataFID)*1e6, '-', ...
        timeFID*1e3, real(dataFID)*1e6, ':', ...
        timeFID*1e3, imag(dataFID)*1e6, '-.');

      plot(ax(2), flipAngle(t)*1e6, maxFID(t)*1e6, 'x');
    end
  end
  MaxFIDinterp = interpn(maxFID, 4, 'cubic');
  flipAngleinterp = interpn(flipAngle, 4, 'cubic');
  [maxFIDAmp, maxi] = max(MaxFIDinterp);
  tFlip90 = flipAngleinterp(maxi);
  if doPlot
    hold(ax(1), 'off');
    plot(ax(2), flipAngleinterp*1e6, MaxFIDinterp*1e6, '-');
    plot(ax(2), repmat(tFlip90*1e6,2,1), [0,maxFIDAmp*1e6], '-o');
    ylim(ax(2), [0, ceil(maxFIDAmp*1e6/10)*10]);
    hold(ax(2), 'off');

    title(ax(2), sprintf('p90 = %4f %cs', tFlip90*1e6, 181));
  end
else
  tFlip90 = tPulse90;
end
% calculate magnetic flux density of the transmit coil at HW.TX.AmpDef amplitude
B1 = (1/(tFlip90*4)) / (HW.GammaDef/2/pi);


%% Store calculated coil efficiency in file or display it on screen
newPaUout2Amplitude = HW.TX(iDevice).PaUout2Amplitude;
TXVoltage = HW.TX(iDevice).AmpDef / HW.TX(iDevice).PaUout2Amplitude(HW.TX(iDevice).ChannelDef);
newPaUout2Amplitude(HW.TX(iDevice).ChannelDef) = B1 / TXVoltage;
comment = sprintf('%s (tFlip90 = %.3f us @ %.3f V) from bulk FID by %s', ...
  datestr(now, 'yyyy-mm-ddTHH:MM:SS'), tFlip90*1e6, TXVoltage, mfilename());

newCalLine = sprintf('HW.TX(%d).PaUout2Amplitude = [%.6f, %.6f]*1e-6;', ...
  iDevice, newPaUout2Amplitude*1e6);

if ~isempty(HW.TX(iDevice).CoilName)
  newCalLine = sprintf('if strcmp(HW.TX(%d).CoilName, ''%s''),  %s  end', ...
    iDevice, HW.TX(iDevice).CoilName, newCalLine);
end

if iDevice > 1
  newCalLine = sprintf('if numel(HW.TX) >= %d,  %s  end', ...
    iDevice, newCalLine);
end

newCalLine = sprintf('%s  %% %s\n', newCalLine, comment);
% FIXME: Add consistency checks
savePulseFile = true;
if savePulseFile
  if ~isempty(HW.TX(iDevice).PaUout2AmplitudePath) ...
      && (isemptyfield(mySave, 'DummySerial') ...
          || mySave.DummySerial(min(iDevice, numel(mySave.DummySerial))) <= 0)
    if ~exist(HW.TX(iDevice).PaUout2AmplitudePath, 'file')
      newCalLine = ['% factor from voltage amplitude at the coil input to B1+ field strength in T/V', ...
        sprintf('\n'), newCalLine];
    end
    if ~exist(fileparts(HW.TX(iDevice).PaUout2AmplitudePath), 'dir')
      mkdir(fileparts(HW.TX(iDevice).PaUout2AmplitudePath));
    end
    fid = fopen(HW.TX(iDevice).PaUout2AmplitudePath, 'a+');
    fwrite(fid, newCalLine);
    [~, name, ~] = fileparts(fopen(fid));
    fclose(fid);
    clear(name); clear('name');
    fprintf('\nA new line was added to the following file:\n%s\n%s\n', ...
      HW.TX(iDevice).PaUout2AmplitudePath, newCalLine);
  else
    fprintf('\n\nAdd the following line to your LoadMySystem.m file:\n%s\n', ...
      newCalLine);
  end
else
  fprintf('\n');
  warnStr = ['Determination of pulse length unsuccessful!\n', ...
    'If you want to use the uncertain best guess value anyway, ', ...
    'please manually add the following line to your LoadMySystem.m file:\n%s\n'];
  warning('PD:sequence_PulseDuration', warnStr, newCalLine);
end

end
