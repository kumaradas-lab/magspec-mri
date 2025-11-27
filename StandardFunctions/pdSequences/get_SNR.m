function [AmpSnrHz, FdB, F, AllOutputs] = get_SNR(HW, showPlot, tEcho, fSampleDivider, nEchoes, BW, nSmoothing)
%% Measure SNR
%
%  [AmpSnrHz, FdB, F, AllOutputs] = get_SNR(HW, showPlot, tEcho, fSampleDivider, nEchoes, BW, nSmoothing)
%
% Two similar Spin Echo measurements with equivalent acquisition windows at
% their respective echo time are performed where the first is *not* preceeded by
% any HF pulses while the second one is. From these measurements the SNR is
% calculated and the noise spectrum is displayed.
%
% For this to work correctly, a suitable sample should be inside the probe head.
% The delay to any prior measurement should be large enough to account for the
% T1 relaxation time of the sample.
%
%
% INPUT:
%
%   HW
%         HW object or structure
%         The structure HW.DefSeqValues.get_SNR can have the following fields
%         that influence the measurement:
%     SeqN
%           Structure with settings for the noise measurement. See documentation
%           of sequence_EchoStandard for available settings.
%     SeqS
%           Structure with settings for the signal measurement. See
%           documentation of sequence_EchoStandard for available settings.
%     iDevice
%           MMRT device to be used if multiple devices are connected.
%     iAQChannel
%           Acquisition channel that is used for the SNR measurement in case the
%           device supports multiple AQ channels.
%
%   showPlot
%         Boolean, figure handle, or (integer valued) double. If it is set to
%         false (the default), the SNR is measured without displaying a result
%         plot. If it is a figure handle or a positive (integer valued) double,
%         the results are shown in a plot which gets the focus. If it is a
%         negative (integer valued) double, the corresponding figure does not
%         get the focus if it already exists.
%
%   tEcho
%         Echo time in seconds. (Default: 2e-3)
%
%   fSampleDivider
%         Divider for the system frequency HW.RX.fSample (usually 125 MHz) to
%         get the sampling frequency of the acquisition window. (Default: 1250)
%
%   nEchoes
%         Number of echoes in the echo train (see Seq.nEchos for
%         sequence_EchoStandard). (Default: 2)
%
%   BW
%         Bandwidth of measured noise spectrum in Hertz. (Default: 5000)
%
%   nSmoothing
%         Size of the smoothing window to calculate the noise spectrum (odd
%         integer). (Default: 51)
%
%
% OUTPUT:
%
%   AmpSnrHz
%         Amplitude of SNR per sqrt(Hz).
%
%   FdB
%         Noise figure in dB.
%
%   F
%         Noise figure.
%
%   AllOutputs
%         Structure with the following fields:
%
%     AmpSnrHz
%           Amplitude of SNR per sqrt(Hz).
%     FdB
%           Noise figure in dB.
%     F
%           Noise figure.
%     PnoisePeakdBm
%           Peak noise power in dBm.
%     PsignaldBm
%           Power of signal in dBm.
%     Usignal
%           Amplitude of the signal in Volts.
%     SNRHz
%           Power SNR per Hz
%     UnoiseDensityMean
%           Mean noise density in Volts.
%     dfs
%           Frequency resolution of signal in Hertz.
%     dfn
%           Frequency resolution of noise in Hertz.
%     showPlot
%           See input argument.
%     tEcho
%           Used echo time in seconds.
%     nEchoes
%           Number of echoes (see input).
%     BW
%           Used bandwidth of measured noise in Hertz.
%     dataN
%           Structure with data of noise measurement (see "get_data").
%     SeqOutN
%           Structure with settings of noise measurement (see
%           "sequence_EchoStandard").
%     dataS
%           Structure with data of signal measurement (see "get_data").
%     SeqOutS
%           Structure with settings of signal measurement (see
%           "sequence_EchoStandard").
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%%

if nargin < 2 || isempty(showPlot)
  showPlot = 0;
end
if nargin < 3 || isempty(tEcho)
  tEcho = 2e-3;
end
fSampleDividerSet = true;
if nargin < 4 || isempty(fSampleDivider)
  fSampleDivider = 1250;
  fSampleDividerSet = false;
end
if nargin < 5 || isempty(nEchoes)
  nEchoes = 2;
end
if nargin < 6 || isempty(BW)
  BW = 5000;
end
if nargin < 7 || isempty(nSmoothing)
  nSmoothing = 51;
end


if isnumeric(showPlot) && showPlot < 0
  stealFocus = false;
  showPlot = abs(showPlot);
else
  stealFocus = true;
end
if showPlot == 1
  showPlot = 7134;
end


if isemptyfield(HW.DefSeqValues, {'get_SNR', 'iDevice'})
  iDevice = 1;
else
  iDevice = HW.DefSeqValues.get_SNR.iDevice;
end


%% noise measurement
if isemptyfield(HW.DefSeqValues, {'get_SNR', 'SeqN'})
  SeqN = struct();
else
  SeqN = HW.DefSeqValues.get_SNR.SeqN;
end


% FIXME: Should the noise measurement be equivalent to the signal measurement
%        (tEcho, AQEcho, ...)?
if isemptyfield(SeqN, 'average')
  SeqN.average = 1;
end
if ~isfield(SeqN, 'plotSeq')
  SeqN.plotSeq = [];  % plot pulse program; [] for no plot, 0 for TX and RX only
end
if isemptyfield(SeqN, 'tEcho')
  SeqN.tEcho = max(200e-3, tEcho*8);  % echo time in seconds
end
if isemptyfield(SeqN, 'plot')
  SeqN.plot = 0;  % plot data
end
if isemptyfield(SeqN, 'p90')
  % duration of excitation pulse in seconds
  SeqN.p90 = HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef/2;
end
if isemptyfield(SeqN, 'p180')
  % duration of inversion pulse in seconds
  SeqN.p180 = HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef;
end
if isemptyfield(SeqN, 'nEchos')
  SeqN.nEchos = 1;  % number of echoes in echo train
end
if isemptyfield(SeqN, 'fSample')
  SeqN.fSample = HW.RX(iDevice).fSample / min(8191, fSampleDivider) / 1;  % HW.RX.fSample/ 4,5,6....8191
else
  if fSampleDividerSet && (fSampleDivider ~= HW.RX(iDevice).fSample / SeqN.fSample)
    warning('PD:get_SNR:fSampleDividerChanged', ...
      'fSampleDivider adjusted to match SeqN.fSample');
  end
  fSampleDivider = HW.RX(iDevice).fSample / SeqN.fSample;
end
if isemptyfield(SeqN, 'AQEcho')
  SeqN.AQEcho = 0.5;  % relative part of echo time with data acquisiton (around echo time) ]0...1[
end
if isemptyfield(SeqN, 'FlipPulse')
  SeqN.FlipPulse = @Pulse_Rect_DummyNAN;  % dummy 90 degrees pulse function handle
end
if isemptyfield(SeqN, 'InvertPulse')
  SeqN.InvertPulse = @Pulse_Rect_DummyNAN;  % dummy 180 degrees pulse function handle
end

[dataN, SeqOutN] = sequence_EchoStandard(HW, SeqN);


if isemptyfield(HW.DefSeqValues, {'get_SNR', 'iChannel'})
  iChannel = 1;
else
  iChannel = HW.DefSeqValues.get_SNR.iAQChannel;
end
iAQ = find([SeqOutN.AQ(:).Channel] == iChannel & [SeqOutN.AQ(:).Device] == iDevice, 1, 'first');
dataN = dataN(iAQ);


%% signal measurement
if isemptyfield(HW.DefSeqValues, {'get_SNR', 'SeqS'})
  Seq = struct();
else
  Seq = HW.DefSeqValues.get_SNR.SeqS;
end


if ~isfield(Seq, 'plotSeq')
  Seq.plotSeq = [];  % plot pulse program; [] for no plot, 0 for TX and RX only
end
if isemptyfield(Seq, 'tEcho')
  Seq.tEcho = tEcho;  % echo time in seconds
end
if isemptyfield(Seq, 'plot')
  Seq.plot = 0;  % plot data
end
if isemptyfield(Seq, 'p90')
  % duration of excitation pulse in seconds
  Seq.p90 = HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef/2;
end
if isemptyfield(Seq, 'p180')
  % duration of inversion pulse in seconds
  Seq.p180 = HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef;
end
if isemptyfield(Seq, 'nEchos')
  Seq.nEchos = nEchoes;  % number of echoes in echo train
end
if isemptyfield(Seq, 'fSample')
  Seq.fSample = HW.RX(iDevice).fSample / min(8191, fSampleDivider) / 1;  % HW.RX.fSample/ 4,5,6....8191
end
if Seq.fSample ~= SeqN.fSample
  error('PD:get_SNR:fSampleMismatch', ...
    'fSample for noise and signal measurement must match');
end
if isemptyfield(Seq, 'AQEcho')
  % relative part of echo time with data acquisiton (around echo time) ]0...1[
  % Seq.AQEcho = (4*floor(Seq.tEcho*Seq.fSample*0.1/4)) / Seq.fSample/Seq.tEcho;
  Seq.AQEcho = (4*floor(Seq.tEcho*Seq.fSample*0.5/4)) / Seq.fSample/Seq.tEcho;
end
if Seq.nEchos == 0 && isemptyfield(Seq, 'AQFID')
  Seq.AQFID = 1;
end

[dataS, SeqOutS] = sequence_EchoStandard(HW, Seq);

iAQ = find([SeqOutS.AQ(:).Channel] == iChannel & [SeqOutS.AQ(:).Device] == iDevice, 1, 'first');
dataS = dataS(iAQ);


%% SNR calculation
PnoiseRinDensity = HW.Constant.Boltzmann * HW.TempR50 * 1;
UnoiseRinDensity = sqrt(PnoiseRinDensity * HW.RX(iDevice).Rin);

dfn = dataN.f_fft1_data(2,1,end) - dataN.f_fft1_data(1,1,end);

relFSD = fSampleDivider / min(8191, fSampleDivider);

INmin = ceil(size(dataN.fft1_data, 1) * (1/2 - 1/4/relFSD)) + 1;
INmax = ceil(size(dataN.fft1_data, 1) * (1/2 + 1/4/relFSD));
INroi = INmin:INmax;

ZFF = 1;
nSmoothing = floor(nSmoothing*ZFF/2)*2 + 1;  % only odd nSmoothing

smoothWindow = ones(nSmoothing, 1) / nSmoothing;
% smoothWindow=RaisedCosine(nSmoothing+2).^2/sum(RaisedCosine(nSmoothing+2).^2);
% dataN.fft1_dataZ=zeroFill_image(dataN.fft1_data(:,1,end),[size(dataN.fft1_data(:,1,end),1)*ZFF,1],1);
PnoiseDensityPerBin = abs(dataN.fft1_data(:,1,end).*(dataN.Amplitude2Uin(1)/HW.RX(iDevice).LNAGain./2^0.5)).^2 / HW.RX(iDevice).Rin / dfn * SeqN.average; % P = Ueff^2 / Rin;
PnoiseDensityMean = mean(PnoiseDensityPerBin(INroi));
PnoiseDensityPerBinSmooth = conv(PnoiseDensityPerBin, smoothWindow, 'same');
PnoisePeak = max(PnoiseDensityPerBinSmooth(INroi)*dfn/smoothWindow(round(end/2)));
PnoisePeakdBm = 10 * log10(PnoisePeak/1e-3);
UnoiseDensityMean = sqrt(PnoiseDensityMean*HW.RX(iDevice).Rin);

dfs = dataS.f_fft1_data(2,1,end) - dataS.f_fft1_data(1,1,end);

ISmin = ceil(size(dataS.fft1_data,1)*(1/2-1/4/relFSD)) + 1;
ISmax = ceil(size(dataS.fft1_data,1)*(1/2+1/4/relFSD));
ISroi = ISmin:ISmax;

s1m = floor(size(dataS.fft1_data,1)/2) + 1;
s1BW = ceil((BW/dfs)/2);
s1s = max(s1m-s1BW, ISmin);
s1e = min(s1m+s1BW, ISmax);

% Usignal= sqrt((  sum(abs(dataS.fft1_data(s1s :s1e ,1,end)).^2) ... %   center BW
%                 -(s1e-s1s+1)*mean(abs(dataS.fft1_data([ISmin:s1s-1,s1e+1:ISmax],1,end)).^2) ...  % - noise power
%               ))*(dataS.Amplitude2Uin(end)/HW.RX(iDevice).LNAGain/2^0.5);
if Seq.nEchos == 0
  %%
  nSamplesGrid=0;
  myTime = dataS.time_all(:,1,end);
  FIDs = abs(dataS.data(:,1,end))*(dataS.Amplitude2Uin(end)/HW.RX(iDevice).LNAGain/2^0.5);
  % StopI=min([max([find(FIDs>FIDs(1)*3/4,1,'last'),2]),numel(FIDs)]);
  StopI = min([max([find(myTime<myTime(1)*4,1,'last'),2]),numel(FIDs)]);

  [p, S, mu] = polyfit([-flipud(myTime(nSamplesGrid+(1:StopI)));myTime(nSamplesGrid+(1:StopI))], ...
                       [ flipud(abs(FIDs(nSamplesGrid+(1:StopI))));abs(FIDs(nSamplesGrid+(1:StopI)))],8);  % extrapolate to horizontal
  Usignal = polyval(p, 0, S, mu);
  tFit = linspace(0, myTime(StopI), 1000).';
  UsignalFit = polyval(p, tFit, S, mu);
else
  Usignal= sqrt(  abs(sum(dataS.fft1_data(s1s :s1e ,1,end))).^2 ... %   center BW
                 -(s1e-s1s+1)*mean(abs(dataS.fft1_data([ISmin:s1s-1,s1e+1:ISmax],1,end)).^2) ...  % - noise power
                ) * (dataS.Amplitude2Uin(end)/HW.RX(iDevice).LNAGain/2^0.5);
end
if ~isreal(Usignal)
  Usignal = 0;
end
Psignal = Usignal^2 / HW.RX(iDevice).Rin;
PsignaldBm = 10*log10(Psignal/1e-3);

SNRHz = (Psignal) / PnoiseDensityMean;
AmpSnrHz = SNRHz^0.5;
F = UnoiseDensityMean.^2 / UnoiseRinDensity.^2;
FdB = 10 * log10(F);


%% plot results
if ishghandle(showPlot, 'figure') || showPlot
  if ~ishghandle(showPlot, 'figure')
    hf = figure(showPlot);
  else
    hf = showPlot;
  end
  hf = clf(hf);
  if stealFocus
    hf = figure(hf);
  end
  ax1 = subplot(3,1,1, 'Parent', hf);
  plot(ax1, dataN.f_fft1_data(INroi,1,end), sqrt(PnoiseDensityPerBinSmooth(INroi)*HW.RX(iDevice).Rin), 'b'...
          , dataN.f_fft1_data(INroi,1,end), sqrt(PnoiseDensityMean*HW.RX(iDevice).Rin) * ones(size((dataN.f_fft1_data(INroi,1,end)))), '--b'...
          , dataN.f_fft1_data(INroi,1,end), UnoiseRinDensity * ones(size((dataN.f_fft1_data(INroi,1,end)))), '--g');
  title(ax1, 'Noise Vrms/\surd{Hz}');
  xlabel(ax1, 'frequency in Hz');
  ylabel(ax1, 'Vrms/\surd{Hz}');
  ylim(ax1, [0, Inf]);
  grid(ax1, 'on');
  legend(ax1, {'noise', 'mean noise', sprintf('%.0f Ohm noise', HW.RX(iDevice).Rin)}, 'Location', 'SouthEast');

  ax2 = subplot(3,1,2, 'Parent', hf);
  plot(ax2, dataS.f_fft1_data(ISroi,1,end), (dataS.Amplitude2Uin(1)/HW.RX(iDevice).LNAGain./2^0.5)*abs(dataS.fft1_data(ISroi,1,end)), 'c'...
          , dataS.f_fft1_data(s1s:s1e,1,end), (dataS.Amplitude2Uin(1)/HW.RX(iDevice).LNAGain./2^0.5)*abs(dataS.fft1_data(s1s:s1e,1,end)), 'b'...
          , dataS.f_fft1_data([ISmin;s1m;ISmax],1,end), [Usignal;Usignal;Usignal], '-xg'...
          , dataS.f_fft1_data([ISmin;ISmax],1,end), [sqrt(PnoiseDensityMean*HW.RX(iDevice).Rin);sqrt(PnoiseDensityMean*HW.RX(iDevice).Rin)]*sqrt(dfs), '--b'...
          , dataS.f_fft1_data([ISmin;ISmax],1,end), [UnoiseRinDensity;UnoiseRinDensity]*sqrt(dfs), '--g');
  xlabel(ax2, 'frequency in Hz');
  ylabel(ax2, 'Vrms');
  ylim(ax2, [0, Inf]);
  grid(ax2, 'on');
  legend(ax2, {'signal', 'signal BW', 'max signal', 'AQ noise', sprintf('%.0f Ohm noise', HW.RX(iDevice).Rin)}, 'Location', 'NorthEast');
  linkaxes([ax1, ax2], 'x');
  xlim(ax2, [dataN.f_fft1_data(INmin,1,end), dataN.f_fft1_data(INmax,1,end)]);

  ax3 = subplot(3,1,3, 'Parent', hf);
  if Seq.nEchos == 0
    plot(ax3, dataS.time_all(:,1,end), abs(dataS.data(:,1,end))*(dataS.Amplitude2Uin(end)/HW.RX(iDevice).LNAGain/2^0.5), 'c'...
            , tFit, UsignalFit, '-.b'...
            , [0;dataS.time_all(end,1,end)], abs([Usignal,Usignal]), 'g'...
            , dataS.time_all([1,end],1,end), [sqrt(PnoiseDensityMean*HW.RX(iDevice).Rin);sqrt(PnoiseDensityMean*HW.RX(iDevice).Rin)]*sqrt(SeqOutS.fSample)*mean((1./dataS.cic_corr(:,1,end)).^2).^(0.5), '--b'...
            , dataS.time_all([1,end],1,end), [UnoiseRinDensity;UnoiseRinDensity]*sqrt(SeqOutS.fSample)*mean((1./dataS.cic_corr(:,1,end)).^2).^(0.5), '--g');
    legend(ax3, {'signal', 'extrapolation to t=0', 't=0 signal', 'AQ noise', sprintf('%.0f Ohm noise', HW.RX(iDevice).Rin)}, 'Location', 'NorthEast');
  else
    fftsm = dataS.fft1_data(:,1,end)*0;
    fftsm(s1s:s1e,1,end) = dataS.fft1_data(s1s:s1e,1,end);
    plot(ax3, dataS.time_all(:,1,end), abs(dataS.data(:,1,end))*(dataS.Amplitude2Uin(end)/HW.RX(iDevice).LNAGain/2^0.5), 'c'...
            , dataS.time_all(:,1,end), abs(ifftshift(fft(fftshift(fftsm))))*(dataS.Amplitude2Uin(end)/HW.RX(iDevice).LNAGain/2^0.5), 'b'...
            , dataS.time_all([1,end],1,end), abs([Usignal,Usignal]), 'g'...
            , dataS.time_all([1,end],1,end), [sqrt(PnoiseDensityMean*HW.RX(iDevice).Rin);sqrt(PnoiseDensityMean*HW.RX(iDevice).Rin)]*sqrt(SeqOutS.fSample)*mean((1./dataS.cic_corr(:,1,end)).^2).^(0.5), '--b'...
            , dataS.time_all([1,end],1,end), [UnoiseRinDensity;UnoiseRinDensity]*sqrt(SeqOutS.fSample)*mean((1./dataS.cic_corr(:,1,end)).^2).^(0.5), '--g');
    legend(ax3, {'signal', 'signal BW', 'max signal', 'AQ noise', sprintf('%.0f Ohm noise', HW.RX(iDevice).Rin)}, 'Location', 'NorthEast');
  end

  xlabel(ax3, 'time in s');
  ylabel(ax3, 'Vrms');
  title(ax1, ['noise peak power = ' num2str(PnoisePeakdBm,'%.2f') ' dBm, noise figure = ' num2str(F,'%.2f') ' (' num2str(FdB,'%.1f') ' dB)']);
  title(ax2, ['signal power = ' num2str(PsignaldBm,'%.2f') ' dBm, amplitude SNR = ' num2str(AmpSnrHz,'%.0f') ' \surd{Hz}']);
  title(ax3, sprintf('max signal = %.3f %cVrms, noise = %.3f nVrms/\\surd{Hz}', Usignal*1e6, 181, UnoiseDensityMean*1e9));
  grid(ax3, 'on');
  drawnow();
end


if nargout > 3
  AllOutputs.AmpSnrHz = AmpSnrHz;
  AllOutputs.FdB = FdB;
  AllOutputs.F = F;
  AllOutputs.PnoisePeakdBm = PnoisePeakdBm;
  AllOutputs.PsignaldBm = PsignaldBm;
  AllOutputs.Usignal = Usignal;
  AllOutputs.SNRHz = SNRHz;
  AllOutputs.UnoiseDensityMean = UnoiseDensityMean;
  AllOutputs.dfs = dfs;
  AllOutputs.dfn = dfn;
  AllOutputs.showPlot = showPlot;
  AllOutputs.tEcho = tEcho;
  AllOutputs.fSampleDivider = fSampleDivider;
  AllOutputs.nEchoes = Seq.nEchos;
  AllOutputs.BW = BW;
  AllOutputs.dataN = dataN;
  AllOutputs.SeqOutN = SeqOutN;
  AllOutputs.dataS = dataS;
  AllOutputs.SeqOutS = SeqOutS;
end

%%
% disp(['Rauschzahl = ' num2str(F) '   Rauschzahl in dB = ' num2str(FdB) ' dB'])

% 2dB von Grad

% PnoiseHz150=PnoiseHz
% PnoiseHz100=PnoiseHz
% PnoiseHz50=PnoiseHz
% PnoiseHz30=PnoiseHz
% PnoiseGemessen=[2.129672935407666e-20,2.157812871605913e-20,2.263034825142230e-20,2.301434980606803e-20]
% Temp=[30,50,100,150]+273
% hold on
% figure(30);plot(Temp,PnoiseGemessen)
% PnoiseGemessen=[PnoiseHz30,PnoiseHz50,PnoiseHz100,PnoiseHz150]
% Steigung=mean(diff(PnoiseGemessen)./diff(Temp))
% figure(30);plot(0:500,HW.Constant.Boltzmann*(0:500))
% figure(30);plot(0:500,Steigung*(0:500))
% figure(30);plot(0:2000,HW.Constant.Boltzmann*(0:2000))
% Offset=mean(PnoiseGemessen)-mean(Temp)*Steigung
% figure(30);plot(0:500,Steigung*(0:500)+Offset)
% TnoisOffset=Offset/HW.Constant.Boltzmann

end
