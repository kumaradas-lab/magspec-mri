%% Demo Noise Figure from measured data
% Sample data is acquired. The noise figure is calculated in time and frequency
% domain.
% See also: get_SNR(HW, 1);


%% Preparations
LoadSystem;                                     % load system parameters

%  select frequency, there may be some extra noise from the usb connection at
%  e.g. 24, 48, 53, 29, 5, 19, 43 MHz harmonics of 24 MHz acquired with 125 MHz sampling rate
%  HW.fLarmor = 24.1e6;     % set Larmor frequency
%  HW.fLarmor = 6e6;        % set Larmor frequency

%  HW.RX.LNAGain = HW.RX.LNAGain*1.28; gain correction
%  [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 10, [], 1);  % find magnet frequency
%  HW.Grad.PowerDown  = 1;                             % turn off gradient amplifier

% Parameters used for timing calculations
Seq.plotSeq     = [];                               % plot sequence off

% Sequence parameters
Seq.average     = 1;                                % number of averages 

% RF transmission parameters
TX.Start        = NaN;                              % no TX

% Acquisition parameters
RawData = 0;                                        % whole spectrum or around HW.fLarmor
if RawData
  AQ.fSample      = HW.RX(1).fSample;               % sampling rate of AQ window
  AQ.nSamples     = 1024*64;                        % number of samples in AQ window
  AQ.Start        = 0;                              % acquisition start time
else
  AQ.fSample      = 100e3;                          % sampling rate of AQ window (1/BW)
  AQ.nSamples     = 1*AQ.fSample;                   % number of samples in AQ window
  AQ.Start        = -0.5/AQ.fSample;                % acquisition start time
end

AQ.Frequency    = HW.fLarmor;                       % frequency of AQ window
AQ.Phase        = 0;                                % phase of AQ window
% AQ.Gain         = HW.RX.Amplitude2LnaUin / 20e-3; % maximum input voltage

Seq.tRep        = (AQ.nSamples+10)/AQ.fSample+20e-3;     % repetition time
Seq.tOffset     = 10e-3;

% Start measurement
[Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);

% Plot results
HW.RX(1).AmplitudeUnitScale = 1e-6./SeqOut.AQ(1).Amplitude2Uin.*SeqOut.HW.RX(1).LNAGain;
HW.RX(1).AmplitudeUnit = sprintf('%cV', char(181));
% HW.RX(1).AmplitudeUnitScale = 1/SeqOut.AQ(1).Amplitude2Norm; HW.RX(1).AmplitudeUnit = 'FS';

plot_data_1D(HW, data_1D);

%% Load Gain Correction form CalibrationRxNoise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('CalibrationRxNoise.mat','file')
  load('CalibrationRxNoise.mat');
  data(1).fft1_data_GainCorrection = ...
    interp1(CalibrationRx.Frequency, CalibrationRx.Gain, data(1).f_fft1_data(:,1,end));
  data(1).data_GainCorrection = ...
    interp1(CalibrationRx.Frequency, CalibrationRx.Gain, SeqOut.HW.fLarmor);
else
  clear CalibrationRx
  data(1).fft1_data_GainCorrection = 1;
  data(1).data_GainCorrection = 1;
end

%% Noise Figure estimated from time domain data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nTime domain:\n');

% thermal noise at HW.TempR50
ThermalNoisePowerRinDensity = SeqOut.HW.Constant.Boltzmann * SeqOut.HW.TempR50 * 1;
ThermalNoiseUrmsRinDensity = sqrt(ThermalNoisePowerRinDensity * HW.RX(1).Rin);
BwNoise = SeqOut.AQ(1).fSample;

%% Power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             amplitude              to Uin                     LNA gain                   gain cal                       average^0.5
dataUpeak = data(1).data(:,1,end) .* data(1).Amplitude2Uin(1) / SeqOut.HW.RX(1).LNAGain ./ data(1).data_GainCorrection .* SeqOut.average^0.5;
% P = Ueff^2 / Rin;
% PowerDensity                   Ueff           ^2  /     Rin            /   BW
DataPowerDensityCic = mean(abs(dataUpeak/2^0.5).^2 / SeqOut.HW.RX(1).Rin / BwNoise );

% The CIC filter reduces the amplitude of off-center frequencies.
% Correct that by the mean reduction factor.
CicPowerFactor = mean((1./data(1).cic_corr).^2);
DataPowerDensity = DataPowerDensityCic ./ CicPowerFactor;
fprintf('Noise power density:       %10.3e W/Hz\n', DataPowerDensity);

NoiseFigure = DataPowerDensity / ThermalNoisePowerRinDensity;
fprintf('Noise figure (power):      %4.1f (%.2f dB)\n', NoiseFigure, 10*log10(NoiseFigure));

% Voltage
DataUrmsDensity = sqrt(DataPowerDensity * SeqOut.HW.RX(1).Rin);
fprintf('RMS Voltage noise density: %6.3f nVrms/sqrt(Hz)\n', DataUrmsDensity*1e9);
% DataUrmsNoise = DataUrmsDensity.*SeqOut.AQ.fSample.^0.5;
DataUrmsCic = mean(abs(dataUpeak/2^0.5).^2).^0.5;
% DataUrmsCic = 92e-9
DataUrms = DataUrmsCic ./ (CicPowerFactor.^0.5);
fprintf('RMS Voltage noise:         %6.3f %cVrms\n', DataUrms*1e6, char(181));
ThermalNoiseUrms = ThermalNoiseUrmsRinDensity .* BwNoise.^0.5;
NoiseFigure = (DataUrms/ThermalNoiseUrms).^2;
fprintf('Noise figure (voltage):    %4.1f (%.2f dB)\n', NoiseFigure, 10*log10(NoiseFigure));


%% Noise Figure estimated from corrected FFT data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nFrequency domain:\n');

FrequencyPerBin = SeqOut.AQ(1).fSample(end) / SeqOut.AQ(1).nSamples(end);

% Select the center part of the FFT (bandwidth of interest)
if SeqOut.AQ(1).fSample < SeqOut.HW.RX(1).fSample
  BwOI = FrequencyPerBin * SeqOut.AQ(1).nSamples/2;
  IBwOI = (data(1).f_fft1_data(:,1,end) - SeqOut.AQ(1).Frequency(end)<BwOI/2) & ...
    (data(1).f_fft1_data(:,1,end) - SeqOut.AQ(1).Frequency(end)>-BwOI/2);
else
  BwOI = SeqOut.HW.RX(1).fSample*23/25/2;
  % BwOI = SeqOut.HW.RX(1).fSample/2;
  IBwOI = (data(1).f_fft1_data(:,1,end) - SeqOut.HW.RX(1).fSample/4<BwOI/2) & ...
    (data(1).f_fft1_data(:,1,end) - SeqOut.HW.RX(1).fSample/4>-BwOI/2);
end

% Note: data.fft1_data is already CIC corrected
%                 amplitude                   to Uin                     LNA Gain                   gain cal                            average^0.5     
fft_dataUpeak = data(1).fft1_data(:,1,end) .* data(1).Amplitude2Uin(1) / SeqOut.HW.RX(1).LNAGain ./ data(1).fft1_data_GainCorrection .* SeqOut.average^0.5;
% P = Ueff^2 / Rin;
% PowerDensity                            Ueff     ^2 /     Rin             /   BW
DataPowerDensityPerBin = abs(fft_dataUpeak./2^0.5).^2 / SeqOut.HW.RX(1).Rin / FrequencyPerBin ;  
DataUrmsDensity = sqrt(DataPowerDensity * SeqOut.HW.RX(1).Rin);
DataPowerDensity = mean(DataPowerDensityPerBin(IBwOI));
fprintf('Noise power density:       %10.3e Vrms/Hz\n', DataPowerDensity);
DataPower = DataPowerDensity * BwOI;
DataPowerdBm = 10*log10(DataPower/1e-3);
NoiseFigure = DataPowerDensity / ThermalNoisePowerRinDensity;
NoiseFiguredB = 10*log10(NoiseFigure);
fprintf('Noise figure (power):      %4.1f (%.2f dB)\n', NoiseFigure, NoiseFiguredB);


%% plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSmoothing = SeqOut.AQ(1).nSamples/64;  % smoothing over FFT bins
nSmoothing = floor(nSmoothing/2)*2 + 1;  % only odd nSmoothing
if nSmoothing > sum(~IBwOI)
  warning('"nSmoothing" too high. Increase number of acquired samples or reduce "nSmoothing".');
end
smoothWindow = ones(nSmoothing, 1) / nSmoothing;

DataPowerDensityPerBinSmooth = conv(DataPowerDensityPerBin, smoothWindow, 'same');
DataUrmsDensityPerBinSmooth = sqrt(DataPowerDensityPerBinSmooth * SeqOut.HW.RX(1).Rin);
DataPowerPeak = max(DataPowerDensityPerBinSmooth(IBwOI) * FrequencyPerBin / smoothWindow(ceil(end/2)));
DataPowerPeakdBm = 10*log10(DataPowerPeak/1e-3);
fprintf('Peak power:                 %4.2e W (%.2f dBm)\n', DataPowerPeak, DataPowerPeakdBm);
NoiseFigurePerBin = DataPowerDensityPerBinSmooth / ThermalNoisePowerRinDensity;


hf = figure(11);
clf(hf);
hax = axes(hf);
yyaxis(hax,'right');
plot(hax, data(1).f_fft1_data(IBwOI), DataPowerDensityPerBinSmooth(IBwOI) , ...
  data(1).f_fft1_data(IBwOI) , repmat(ThermalNoisePowerRinDensity, sum(IBwOI), 1));  % with CIC filter correction
ylabel(hax, 'Power noise density in W/Hz');
ylim(hax, [0, Inf]);
yyaxis(hax, 'left');
plot(hax, data(1).f_fft1_data(IBwOI), NoiseFigurePerBin(IBwOI) );  % with CIC filter correction
ylabel(hax, 'Noise figure F');
ylim(hax, [0, Inf]);
title(hax, {'Power noise density of acquired signal with CIC correction.', ...
  ['Noise figure F = ' num2str(NoiseFigure,2) ' (' num2str(NoiseFiguredB,3) ' dB)']});
xlabel(hax, 'Frequency in Hz');
grid(hax, 'on');


hf = figure(12);
clf(hf);
hax = axes(hf);
plot(hax, data(1).f_fft1_data(IBwOI), DataUrmsDensityPerBinSmooth(IBwOI) , ...
  data(1).f_fft1_data(IBwOI) , repmat(ThermalNoiseUrmsRinDensity, sum(IBwOI), 1), '--');  % with CIC filter correction
ylabel(hax, 'RMS Voltage noise density in Vrms/sqrt(Hz)');
ylim(hax, [0, Inf]);
title(hax, {'RMS Voltage noise density of acquired signal with CIC correction.', ...
  ['Noise figure F = ' num2str(NoiseFigure,2) ' (' num2str(NoiseFiguredB,3) ' dB)']});
xlabel(hax, 'Frequency in Hz');
grid(hax, 'on');

%%
hf = figure(13);
clf(hf);
hax = subplot(2,1,1, 'Parent', hf);
plot(hax, data(1).f_fft1_data, ...
  abs(data(1).fft1_data)./data(1).cic_corr .* data(1).Amplitude2Uin(1) / ...
  SeqOut.HW.RX(1).LNAGain*1e6 ./ data(1).fft1_data_GainCorrection);  % without CIC filter correction
% xlim(hax, [data.f_fft1_data(1), data.f_fft1_data(end)])
title(hax, 'FFT of acquired signal without CIC correction');
xlabel(hax, 'Frequency in Hz');
ylabel(hax, sprintf('Voltage in %cV', char(181)));
grid(hax, 'on');

hax = subplot(2,1,2, 'Parent', hf);
plot(hax, data(1).f_fft1_data, ...
  abs(data(1).fft1_data) .* data(1).Amplitude2Uin(1) / ...
  SeqOut.HW.RX(1).LNAGain * 1e6 ./ data(1).fft1_data_GainCorrection);
title(hax, 'FFT of acquired signal with CIC correction');
xlabel(hax, 'Frequency in Hz');  % offset ppm
ylabel(hax, sprintf('Voltage in %cV', char(181)));
grid(hax, 'on');


%% -----------------------------------------------------------------------------
% (C) Copyright 2020-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
