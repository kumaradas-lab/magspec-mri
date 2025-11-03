function [data, SeqOut] = get_CorrectedPhase(data, SeqOut)
%% Correct phase in k-space due to changing magnet frequency during measurement
%
%   [data, SeqOut] = get_CorrectedPhase(data, SeqOut)
%
% This function is used in sequence_Flash if Seq.CorrectPhase is true.
% In this case, the FID immediately after the excitation pulse (before the
% dephasing gradient) is acquired additionally to the actual gradient echo. The
% phase of the gradient echo is corrected by the frequency offset determined
% from the acquired (non-encoded) FID.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% return early if nothing to correct
if ~SeqOut.CorrectPhase
  return;
end

%% default settings
if isemptyfield(SeqOut.AQSlice(1), 'iDevice'), SeqOut.AQSlice(1).iDevice = 1; end
if isemptyfield(SeqOut.AQSlice(2), 'iDevice'), SeqOut.AQSlice(2).iDevice = 1; end

Channel = 1;  % FIXME: Is this always 1?

%%
% FIXME: Handling data from different devices is not yet completely implemented.
% acquisition channel with the gradient echo
iAQ(1) = find(([SeqOut.AQ(:).Device]==SeqOut.AQSlice(1).iDevice) & ([SeqOut.AQ(:).Channel]==Channel), 1, 'first');
% acquisition channel with the FID
iAQ(2) = find(([SeqOut.AQ(:).Device]==SeqOut.AQSlice(2).iDevice) & ([SeqOut.AQ(:).Channel]==Channel), 1, 'first');

% SeqOut.SteadyState_PreShots;
% SeqOut.SteadyState_PostShots;
SeqOut.Correct_nSamples = SeqOut.AQ(iAQ(2)).nSamples(1);
% SeqOut.AQ.fSample;
% SeqOut.AQ.nSamples;
% SeqOut.tRep;
% SeqOut.AQ.Start;
SeqOut.tSlice = SeqOut.Slice(1).CenterOfPulse;
SeqOut.Correct_imageAQWindow = SeqOut.AQSlice(1).UseAQWindow(1,1);
SeqOut.Correct_FidAQWindow = SeqOut.AQSlice(2).UseAQWindow;

if isemptyfield(SeqOut, 'CorrectPlot'), SeqOut.CorrectPlot = 0; end
if isemptyfield(SeqOut, 'CorrectPlotFrequency'), SeqOut.CorrectPlotFrequency = 0; end
if isemptyfield(SeqOut, 'LoopPlot'), SeqOut.LoopPlot = 1; end

% create weights for sliding window average
SeqOut.Correct_mywin = cos(linspace(-0.5*pi, 0.5*pi, (2*min(SeqOut.SteadyState_PreShots, SeqOut.SteadyState_PostShots-1)+3))).'.^2;
SeqOut.Correct_mywin = SeqOut.Correct_mywin(2:end-1);  % remove zeros
SeqOut.Correct_mywin = SeqOut.Correct_mywin/sum(SeqOut.Correct_mywin(:));
% sliding window average of the mean frequency derived from the FID
SeqOut.Correct_foffset = conv(squeeze(mean(diff(unwrap(angle(data(iAQ(2)).data(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep)))))) * ...
                              SeqOut.AQ(iAQ(2)).fSample(1) / 2/pi, SeqOut.Correct_mywin, 'same').';
% mean time of AQ window of FID
tRepStart = cumsum([0, SeqOut.tRep(1:end-1)]);
SeqOut.Correct_tfoffset = SeqOut.AQ(iAQ(2)).Start(1,SeqOut.AQSlice(2).UsetRep) + ...
  SeqOut.AQ(iAQ(2)).nSamples(1,SeqOut.AQSlice(2).UsetRep) / 2 ./ SeqOut.AQ(iAQ(2)).fSample(1,SeqOut.AQSlice(2).UsetRep) + ...
  tRepStart(SeqOut.AQSlice(2).UsetRep);

% correct average for missing data in sliding window
% FIXME: It is unclear what this block of code should do.
% dataAvailable = squeeze(double(any(abs(data(iAQ(2)).data(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep)) > 0, 1)));
% if any(dataAvailable(:) == 0)
%   foffsetWeight = conv2(dataAvailable, SeqOut.Correct_mywin, 'same').';
%   SeqOut.Correct_foffset(foffsetWeight>0) = SeqOut.Correct_foffset(foffsetWeight>0) ./ foffsetWeight(foffsetWeight>0);
% end

% time of excitation pulse
% Assumption: The center of the pulse is at t=0 (tRep time). If that should
% change, this must be corrected here.
SeqOut.Correct_tPulse = tRepStart(SeqOut.Slice(1).UseAtRepetitionTime);

% (interpolated) frequency offset at excitation pulse
if numel(SeqOut.Correct_foffset) == 1
  SeqOut.Correct_foffsetPulse = SeqOut.Correct_foffset;
else
  SeqOut.Correct_foffsetPulse = interp1(SeqOut.Correct_tfoffset, SeqOut.Correct_foffset, SeqOut.Correct_tPulse, 'linear', 'extrap');
end
% phase correction (due to momentary frequency at excitation pulse)
% Use mean frequency between pulse and AQ center to determine the magnitude of
% the dephasing at the recorded FID due to the offset frequency.
SeqOut.Correct_poffsetPulse = ...
  (SeqOut.Correct_foffsetPulse + SeqOut.Correct_foffset)./2 .* ...
  (2*pi) .* (SeqOut.Correct_tfoffset - SeqOut.Correct_tPulse);
% select pulses with image acquisition.
SeqOut.Correct_poffsetPulse = ...
  SeqOut.Correct_poffsetPulse(SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots);


% data in ordering as acquired (not necessarily the same order as in the image)
imageUsetRep = sort(SeqOut.AQSlice(1).UsetRep(:));
% time of samples in AQ windows (for image)
AQnSamples = SeqOut.AQ(iAQ(1)).nSamples(SeqOut.Correct_imageAQWindow,imageUsetRep(1));
SeqOut.Correct_tSamples = squeeze(data(iAQ(1)).time_all(1:AQnSamples,SeqOut.Correct_imageAQWindow,imageUsetRep));

% (interpolated) frequency offset for each AQ sample
SeqOut.Correct_foffsetSamples = zeros(size(SeqOut.Correct_tSamples));
if numel(SeqOut.Correct_foffset) == 1
  SeqOut.Correct_foffsetSamples(:) = SeqOut.Correct_foffset;
else
  SeqOut.Correct_foffsetSamples(:) = interp1(SeqOut.Correct_tfoffset, SeqOut.Correct_foffset, SeqOut.Correct_tSamples(:), 'linear', 'extrap');
end

% select valid data (i.e. the data at the image acquisition) after interpolation
% is complete
SeqOut.Correct_foffset = SeqOut.Correct_foffset(SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots);
SeqOut.Correct_tfoffset = SeqOut.Correct_tfoffset(SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots);

% number of echoes in EPI segments
numEchoes = SeqOut.AQSlice(1).nImages * (SeqOut.AQSlice(1).oddEvenEchoes+1) * SeqOut.AQSlice(1).EPIFactor;

% phase correction (at each sample)
% Use mean frequency between FID center and each sample in the AQ window for the
% image to determine the magnitude of the dephasing at the samples.
SeqOut.Correct_poffsetSamples = ...
  (repelem(SeqOut.Correct_foffset(:).', size(SeqOut.Correct_foffsetSamples,1), numEchoes) + SeqOut.Correct_foffsetSamples)./2 .* ...  % mean frequency offset
  (2*pi) .* (SeqOut.Correct_tSamples - repelem(SeqOut.Correct_tfoffset(:).', size(SeqOut.Correct_foffsetSamples,1), numEchoes));

% combined phase correction (i.e. the dephasing at each sample in the AQ window
% for the image with respect to the respective pulse center)
SeqOut.Correct_poffsetBoth = SeqOut.Correct_poffsetSamples ...
  + repelem(SeqOut.Correct_poffsetPulse, size(SeqOut.Correct_foffsetSamples,1), numEchoes);


if SeqOut.CorrectPlot
  % plot FID data
  hFigure = figure(29);
  hax(1) = subplot(4,1,1, 'Parent', hFigure);
  dt = data(iAQ(2)).data(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep);
  time = data(iAQ(2)).time_of_tRep(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep);
  plot(hax(1), squeeze(time), squeeze(abs(dt)));
  ylabel(hax(1), 'Amp in T');
  title(hax(1), 'Phase Correction FIDs');
  grid(hax(1), 'on');
  hax(2) = subplot(4,1,2, 'Parent', hFigure);
  plot(hax(2), squeeze(time), squeeze(unwrap(angle(dt))));
  ylabel(hax(2), '\phi in rad');
  grid(hax(2), 'on');
  hax(3) = subplot(4,1,3, 'Parent', hFigure);
  plot(hax(3), (time(1:end-1,:) + time(2:end,:))/2, ...
    squeeze(diff(unwrap(angle(dt)))) * SeqOut.AQ.fSample(1) / 2/pi);
  ylabel(hax(3), 'f_{Offset} in Hz');
  xlabel(hax(3), 'time of tRep in s');
  grid(hax(3), 'on');
  linkaxes(hax, 'x');

  hax = subplot(4,1,4, 'Parent', hFigure);
  t_all = data(iAQ(2)).time_all(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep);
  plot(hax, squeeze(mean(t_all, 1)), ...
    squeeze(mean(diff(unwrap(angle(dt))))) * SeqOut.AQ.fSample(1) / 2/pi);
  hold(hax, 'all');
  f_offset = conv2(squeeze(mean(diff(unwrap(angle(data(iAQ(2)).data(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep)))))) * SeqOut.AQ.fSample(1)/2/pi, ...
                   SeqOut.Correct_mywin, 'same');
  plot(hax, squeeze(mean(t_all, 1)), ...
    f_offset);
  hold(hax, 'off');
  ylabel(hax, 'f_{Offset} in Hz');
  xlabel(hax, 'time all in s');
  grid(hax, 'on');

elseif SeqOut.LoopPlot && SeqOut.CorrectPlotFrequency
  hFigure = figure(29);
  clf(hFigure, 'reset');
  hax = axes(hFigure);
  time = tRepStart(SeqOut.AQSlice(2).UsetRep);
  dt = data(iAQ(2)).data(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep);
  plot(hax, ...
    cumsum(time(SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)).', ...
    squeeze(mean(diff(unwrap(angle(dt(:,:,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)))))) * SeqOut.AQ.fSample(1)/2/pi);
  hold(hax, 'all');
  plot(hax, ...
    cumsum(time(SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)).', ...
    conv2(squeeze(mean(diff(unwrap(angle(dt))))) * SeqOut.AQ.fSample(1)/2/pi, ...
          SeqOut.Correct_mywin, 'valid'), ...
    'r');
  hold(hax, 'off');
  legend(hax, {'raw','filtered'});
  ylabel(hax, 'frequency in Hz');
  xlabel(hax, 'time in s');
  title(hax, 'Offset Frequency at FID');
end


if SeqOut.CorrectPlot
  % plot echo data
  dt = data(iAQ(1)).data(1:AQnSamples,SeqOut.Correct_imageAQWindow,imageUsetRep);
  time = data(iAQ(1)).time_of_tRep(1:AQnSamples,SeqOut.Correct_imageAQWindow,imageUsetRep);
  hFigure = figure(30);
  clear hax
  hax(1) = subplot(4,1,1, 'Parent', hFigure);
  plot(hax(1), squeeze(time), squeeze(abs(dt)));
  ylabel(hax(1), 'Amp in T');
  title(hax(1), 'Phase Analysis Echoes (image)');
  grid(hax(1), 'on');
  hax(2) = subplot(4,1,2, 'Parent', hFigure);
  plot(hax(2), squeeze(time), unwrap2Dmiddle(squeeze(angle(dt))));
  ylabel(hax(2), '\phi in rad');
  xlabel(hax(2), 'time of tRep in s');
  grid(hax(2), 'on');
  linkaxes(hax, 'x');

  clear hax
  hax(1) = subplot(4,1,3, 'Parent', hFigure);
  t_all = data(iAQ(1)).time_all(1:AQnSamples,SeqOut.Correct_imageAQWindow,imageUsetRep);
  plot(hax(1), squeeze(mean(t_all, 1)), unwrap2Dmiddle(squeeze(angle(dt(ceil(end/2)+1,:,:)))));
  hold(hax(1), 'on');
  plot(hax(1), squeeze(mean(t_all, 1)), ...
    unwrap2Dmiddle(angle(reshape((dt(floor(end/2)+1,:,:)), [], 1) .* ...
                         reshape(exp(-1i*SeqOut.Correct_poffsetBoth(ceil(end/2)+1,:)), [], 1) )));
  hold(hax(1), 'off');
  ylabel(hax(1), '\phi at center sample in rad');
  grid(hax(1), 'on');
  legend(hax(1), {'measured', 'corrected'});

  hax(2) = subplot(4,1,4, 'Parent', hFigure);
  mean_phase_diff = unwrap(get_MeanPhaseDiffWeighted(dt, 1), [], 3);
  plot(hax(2), squeeze(mean(t_all, 1)), ...
    squeeze(mean_phase_diff) * ...
    SeqOut.AQ(iAQ(1)).fSample(SeqOut.Correct_imageAQWindow(1),imageUsetRep(1)) / 2/pi);
  hold(hax(2), 'all');
  mean_phase_diff = unwrap(get_MeanPhaseDiffWeighted(squeeze(dt) .* exp(-1i*SeqOut.Correct_poffsetBoth), 1), [], 3);
  plot(hax(2), squeeze(mean(t_all, 1)), ...
    mean_phase_diff * ...
    SeqOut.AQ(iAQ(1)).fSample(SeqOut.Correct_imageAQWindow(1),imageUsetRep(1)) / 2/pi);
  t_all = data(iAQ(2)).time_all(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep);
  t_all = t_all(:,:,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots);
  plot(hax(2), squeeze(mean(t_all, 1)), SeqOut.Correct_foffset);
  hold(hax(2), 'off');
  ylabel(hax(2), 'f_{Offset} in Hz');
  xlabel(hax(2), 'time all in s');
  grid(hax(2), 'on');
  legend(hax(2), {'mean', 'mean corrected', 'measured offset'});
  linkaxes(hax, 'x');
end

% correct data for frequency offset
data(iAQ(1)).data(1:AQnSamples,SeqOut.Correct_imageAQWindow,imageUsetRep) = ...
  data(iAQ(1)).data(1:AQnSamples,SeqOut.Correct_imageAQWindow,imageUsetRep) .* ...
  reshape(exp(-1i*SeqOut.Correct_poffsetBoth), [AQnSamples, 1, numel(imageUsetRep)]);

end
