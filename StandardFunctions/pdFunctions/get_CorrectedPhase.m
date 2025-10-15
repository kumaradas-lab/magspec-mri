function [data, SeqOut] = get_CorrectedPhase(data, SeqOut)
%% Correct phase in k-space due to changing magnet frequency during measurement
%
%   [data, SeqOut] = get_CorrectedPhase(data, SeqOut)
%
% This function is used in sequence_Flash and in sequence_Spin_echo if
% Seq.CorrectPhase is true.
% In this case, the FID immediately after the excitation pulse (before the
% dephasing gradient or inversion pulse) is acquired additionally to the actual
% gradient or spin echo.
% Similarly (but only with sequence_Flash), additional excitations with
% non-encoded FIDs can be acquired if Seq.CorrectPhaseSeparate is also true.
% The phase of the gradient or spin echoes is corrected by the frequency offset
% determined from the acquired (non-encoded) FIDs.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2024 Pure Devices GmbH, Wuerzburg, Germany
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


%% spin echo or gradient echo
if isemptyfield(SeqOut, 'SteadyState_PreShots90')
  % gradient echo (FLASH)
  isGradEcho = true;
  steadyState_PreShots = SeqOut.SteadyState_PreShots;
  steadyState_PostShots = SeqOut.SteadyState_PostShots;
  correctPhaseSeparate = SeqOut.CorrectPhaseSeparate;
  echoFactor = SeqOut.AQSlice(1).EPIFactor;
  tRep_excitation = SeqOut.Slice(1).UseAtRepetitionTime;  % FIXME: SingleTrep!
else
  % spin echo
  isGradEcho = false;
  steadyState_PreShots = SeqOut.SteadyState_PreShots90;
  steadyState_PostShots = SeqOut.SteadyState_PostShots90;
  correctPhaseSeparate = false;
  echoFactor = SeqOut.AQSlice(1).TurboFactor;
  tRep_excitation = SeqOut.P90tReps;
end


%%
% FIXME: Handling data from different devices is not yet completely implemented.
% acquisition channel with the (gradient) echo
iAQ(1) = find(([SeqOut.AQ(:).Device]==SeqOut.AQSlice(1).iDevice) & ([SeqOut.AQ(:).Channel]==Channel), 1, 'first');
% acquisition channel with the FID
iAQ(2) = find(([SeqOut.AQ(:).Device]==SeqOut.AQSlice(2).iDevice) & ([SeqOut.AQ(:).Channel]==Channel), 1, 'first');
iData = iAQ;

% assume that the first acquisition on this device is a correction window
tRepAQCorr = find(~isnan(SeqOut.AQ(iAQ(2)).Start), 1, 'first');

% support frequency tracking at secondary frequency to correct image at primary frequency
echoFrequency = SeqOut.AQ(iAQ(1)).Frequency(SeqOut.AQSlice(1).UseAQWindow(1),SeqOut.AQSlice(1).UsetRep(1));
fidFrequency = SeqOut.AQ(iAQ(2)).Frequency(tRepAQCorr);
if abs(echoFrequency - fidFrequency) > 5*eps(echoFrequency)
  % use x-frequency FID for correction
  iData(2) = iData(2) + numel(SeqOut.AQ);
  % disp('Correction with X-frequency FID');
end

% SeqOut.SteadyState_PreShots;
% SeqOut.SteadyState_PostShots;
SeqOut.Correct_nSamples = SeqOut.AQ(iAQ(2)).nSamples(tRepAQCorr);
% SeqOut.Correct_nSamples = SeqOut.AQSlice(2).nRead * SeqOut.AQSlice(2).ReadOS;
% SeqOut.AQ.fSample;
% SeqOut.AQ.nSamples;
% SeqOut.tRep;
% SeqOut.AQ.Start;
SeqOut.tSlice = SeqOut.Slice(1).CenterOfPulse;
SeqOut.Correct_imageAQWindow = SeqOut.AQSlice(1).UseAQWindow;
SeqOut.Correct_FidAQWindow = SeqOut.AQSlice(2).UseAQWindow;

if isemptyfield(SeqOut, 'CorrectPlot'), SeqOut.CorrectPlot = 0; end
if isemptyfield(SeqOut, 'CorrectPlotFrequency'), SeqOut.CorrectPlotFrequency = 0; end
if isemptyfield(SeqOut, 'LoopPlot'), SeqOut.LoopPlot = 1; end

% create weights for sliding window average
windowSize = (2*max(0, min(steadyState_PreShots, steadyState_PostShots-1))+1);
SeqOut.Correct_mywin = cos(linspace(-0.5*pi, 0.5*pi, windowSize+2)).^2;
SeqOut.Correct_mywin = SeqOut.Correct_mywin(2:end-1);  % remove zeros
SeqOut.Correct_mywin = SeqOut.Correct_mywin/sum(SeqOut.Correct_mywin(:));
% sliding window average of the mean frequency derived from the FID
[phaseDiffWeightedMean, phaseDiffWeightedStd] = get_MeanPhaseDiffWeighted(data(iData(2)).data(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep), 1, 'omitnan');
phaseDiffWeightedMean = reshape(phaseDiffWeightedMean, 1, []);
phaseDiffWeightedStd = reshape(phaseDiffWeightedStd, 1, []);
phaseDiffSmoothed = conv(phaseDiffWeightedMean, SeqOut.Correct_mywin, 'same');
% correct (start and end) of smoothed data for missing data in sliding window
w = conv(ones(1, numel(phaseDiffWeightedMean)), SeqOut.Correct_mywin, 'same');
SeqOut.Correct_foffset = phaseDiffSmoothed./w * SeqOut.AQ(iAQ(2)).fSample(tRepAQCorr) / 2/pi ...
  * echoFrequency / fidFrequency;
if any(isnan(phaseDiffWeightedStd))
  SeqOut.Correct_foffsetStd = ...
    sqrt(conv((phaseDiffWeightedMean-phaseDiffSmoothed./w).^2, SeqOut.Correct_mywin, 'same')) ...
    * SeqOut.AQ(iAQ(2)).fSample(tRepAQCorr) / 2/pi ...
    * echoFrequency / fidFrequency;
else
  SeqOut.Correct_foffsetStd = ...
    sqrt(conv(phaseDiffWeightedStd.^2 + ...
              (phaseDiffWeightedMean-phaseDiffSmoothed./w).^2, SeqOut.Correct_mywin, 'same')) ...
    * SeqOut.AQ(iAQ(2)).fSample(tRepAQCorr) / 2/pi ...
    * echoFrequency / fidFrequency;
end


% mean time of AQ window of FID
tRepStart = cumsum([0, SeqOut.tRep(1:end-1)]);
SeqOut.Correct_tfoffset = SeqOut.AQ(iAQ(2)).Start(1,SeqOut.AQSlice(2).UsetRep) + ...
  SeqOut.AQ(iAQ(2)).nSamples(1,SeqOut.AQSlice(2).UsetRep) / 2 ./ SeqOut.AQ(iAQ(2)).fSample(1,SeqOut.AQSlice(2).UsetRep) + ...
  tRepStart(SeqOut.AQSlice(2).UsetRep);

if correctPhaseSeparate
  numFIDs = floor((numel(SeqOut.Slice(1).UseAtRepetitionTime) - steadyState_PreShots - steadyState_PostShots) ...
    / (SeqOut.CorrectPhaseBlockSize+1)) + steadyState_PreShots + steadyState_PostShots;
else
  numFIDs = numel(SeqOut.Slice(1).UseAtRepetitionTime);
  % time of excitation pulse
  % Assumption: The center of the pulse is at t=0 (tRep time). If that should
  % change, this must be corrected here.
  SeqOut.Correct_tPulse = tRepStart(tRep_excitation) + SeqOut.Slice(1).CenterOfPulse;
  % (interpolated) frequency offset at excitation pulse
  if isscalar(SeqOut.Correct_foffset)
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
    SeqOut.Correct_poffsetPulse(steadyState_PreShots+1:end-steadyState_PostShots);
end


% data in ordering as acquired (not necessarily the same order as in the image)
[imageUsetRep, b] = sort(SeqOut.AQSlice(1).UsetRep(:));
if numel(SeqOut.AQSlice(1).UseAQWindow) == 1
  imageAQWindow = SeqOut.Correct_imageAQWindow;
else
  imageAQWindow = reshape(SeqOut.Correct_imageAQWindow(b), size(b));
end
% time of samples in AQ windows (for image)
szData = [size(data(iData(1)).time_all, 1), size(data(iData(1)).time_all, 2), size(data(iData(1)).time_all, 3)];
AQnSamples = SeqOut.AQ(iAQ(1)).nSamples(SeqOut.Correct_imageAQWindow(1),imageUsetRep(1));
if numel(SeqOut.AQSlice(1).UseAQWindow) == 1
  SeqOut.Correct_tSamples = reshape(data(iData(1)).time_all(1:AQnSamples,imageAQWindow,imageUsetRep(:)), AQnSamples, []);
else
  SeqOut.Correct_tSamples = reshape(data(iData(1)).time_all(1:AQnSamples,sub2ind(szData(2:end), imageAQWindow(:), imageUsetRep(:))), AQnSamples, []);
end


% (interpolated) frequency offset for each AQ sample
SeqOut.Correct_foffsetSamples = zeros(size(SeqOut.Correct_tSamples));
if isscalar(SeqOut.Correct_foffset)
  SeqOut.Correct_foffsetSamples(:) = SeqOut.Correct_foffset;
else
  SeqOut.Correct_foffsetSamples(:) = interp1(SeqOut.Correct_tfoffset, SeqOut.Correct_foffset, SeqOut.Correct_tSamples(:), 'linear', 'extrap');
end

% number of echoes in EPI segments
numEchoes = SeqOut.AQSlice(1).nImages * (SeqOut.AQSlice(1).oddEvenEchoes+1) * echoFactor;

idxNoPrePost = (steadyState_PreShots+1):(numFIDs-steadyState_PostShots);

if correctPhaseSeparate
  % interpolate to rf pulses with encoding
  Correct_tfoffset = tRepStart(1,SeqOut.AQSlice(1).UsetRep);
  Correct_foffset = ...
    interp1(SeqOut.Correct_tfoffset, SeqOut.Correct_foffset, Correct_tfoffset);
else
  Correct_tfoffset = SeqOut.Correct_tfoffset(idxNoPrePost);
  Correct_foffset = SeqOut.Correct_foffset(idxNoPrePost);
end

% phase correction (at each sample)
% Use mean frequency between FID center and each sample in the AQ window for the
% image to determine the magnitude of the dephasing at the samples.
if isemptyfield(SeqOut, 'Correct_ZStorageDuration'), SeqOut.Correct_ZStorageDuration = 0; end
if isemptyfield(SeqOut, 'Correct_ZStorageAtInversionNumber'), SeqOut.Correct_ZStorageAtInversionNumber = 1; end
% if isemptyfield(SeqOut, 'Correct_ZStorageAtSpinEcho'), SeqOut.Correct_ZStorageAtSpinEcho = 0; end

SeqOut.Correct_poffsetSamples = ...
  (repelem(reshape(Correct_foffset,1,[]), size(SeqOut.Correct_foffsetSamples,1), numEchoes) + SeqOut.Correct_foffsetSamples)./2 .* ...  % mean frequency offset
  (2*pi) .* (SeqOut.Correct_tSamples-SeqOut.Correct_ZStorageDuration - repelem(reshape(Correct_tfoffset,1,[]), size(SeqOut.Correct_foffsetSamples,1), numEchoes));
if isGradEcho
  if correctPhaseSeparate
    SeqOut.Correct_poffsetBoth = SeqOut.Correct_poffsetSamples;
  else
    % combined phase correction (i.e. the dephasing at each sample in the AQ window
    % for the image with respect to the respective pulse center)
    SeqOut.Correct_poffsetBoth = SeqOut.Correct_poffsetSamples ...
      + repelem(SeqOut.Correct_poffsetPulse, size(SeqOut.Correct_foffsetSamples,1), numEchoes);
  end

else
  % center of inversion pulse
  if SeqOut.SingletRep
    t_invPulse = SeqOut.singletRepStart(SeqOut.Slice(2).UseAtRepetitionTime);
  else
    t_invPulse = tRepStart(SeqOut.Slice(2).UseAtRepetitionTime);
  end
  t_invPulse = reshape(t_invPulse, size(SeqOut.Slice(2).UseAtRepetitionTime));
  % (interpolated) frequency offset at inversion pulse
  if isscalar(SeqOut.Correct_foffset)
    fOffset_invPulse = SeqOut.Correct_foffset;
  else
    fOffset_invPulse = interp1(SeqOut.Correct_tfoffset, SeqOut.Correct_foffset, t_invPulse, 'linear', 'extrap');
  end

  % calculate phase offset from center of frequency tracking window to center of
  % inversion pulse
  t_invPulseOffset=t_invPulse*0;
  t_invPulseOffset(SeqOut.Correct_ZStorageAtInversionNumber+1:end,:)=t_invPulseOffset(SeqOut.Correct_ZStorageAtInversionNumber+1:end,:)+SeqOut.Correct_ZStorageDuration;

  pOffset_InvPulse = ...
    bsxfun(@plus, fOffset_invPulse, SeqOut.Correct_foffset)./2 .* ...
    (2*pi) .* bsxfun(@minus, t_invPulse-t_invPulseOffset, SeqOut.Correct_tfoffset);

  % select inversion pulses with image acquisition.
  pOffset_InvPulse = pOffset_InvPulse(:,steadyState_PreShots+1:end-steadyState_PostShots);


  % handle echo trains

  % calculate phase offset at each inversion pulse relative to center of
  % excitation pulse
  pOffset_InvPulse = bsxfun(@plus, pOffset_InvPulse, SeqOut.Correct_poffsetPulse);

  % calculate steps at each inversion pulse in echo train
  signs = ones(size(pOffset_InvPulse));
  signs(2:2:end,:) = -1;
  pOffset_InvPulse_step = cumsum(2*pOffset_InvPulse.*signs, 1) .* signs;

  % select only inversion pulses with acquisition window
  pOffset_InvPulse_step = pOffset_InvPulse_step(SeqOut.SteadyState_PreShots180+1:end-SeqOut.SteadyState_PostShots180,:);

  pOffset_InvPulse_step = bsxfun(@plus, -SeqOut.Correct_poffsetPulse, pOffset_InvPulse_step);

  % combined phase correction (i.e. the dephasing at each sample in the AQ
  % window for the image with respect to the respective excitation pulse center)
  SeqOut.Correct_poffsetBoth = SeqOut.Correct_poffsetSamples ...
    - repelem(reshape(pOffset_InvPulse_step,1,[]), size(SeqOut.Correct_poffsetSamples,1), 1);
end


if SeqOut.CorrectPlot
  % plot FID data
  hFigure = figure(29);
  clf(hFigure);
  hax(1) = subplot(4,1,1, 'Parent', hFigure);
  dt = reshape(data(iData(2)).data(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep), SeqOut.Correct_nSamples, []);
  time = reshape(data(iData(2)).time_of_tRep(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep), SeqOut.Correct_nSamples, []);
  plot(hax(1), time, abs(dt));
  ylabel(hax(1), 'Amp in T');
  title(hax(1), 'Frequency Tracking Windows');
  grid(hax(1), 'on');
  hax(2) = subplot(4,1,2, 'Parent', hFigure);
  plot(hax(2), time, unwrap(angle(dt)));
  ylabel(hax(2), '\phi in rad');
  grid(hax(2), 'on');
  hax(3) = subplot(4,1,3, 'Parent', hFigure);
  hl = plot(hax(3), [(time(1:end-1,:) + time(2:end,:))/2; NaN(1,size(time,2))], ...
    [diff(unwrap(angle(dt))); NaN(1,size(time,2))] * SeqOut.AQ(iAQ(2)).fSample(tRepAQCorr) / 2/pi);
  if size(time,1) == 2
    set(hl, 'Marker', 'x');
  end
  ylabel(hax(3), 'f_{Offset} in Hz');
  xlabel(hax(3), 'time of tRep in s');
  grid(hax(3), 'on');
  linkaxes(hax, 'x');

  hax = subplot(4,1,4, 'Parent', hFigure);
  t_all = data(iData(2)).time_all(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep);
  t_center_corr = reshape(mean(t_all, 1), 1, []);
  hl(1) = plot(hax, t_center_corr, ...
    phaseDiffWeightedMean * SeqOut.AQ(iAQ(2)).fSample(1) / 2/pi);
  if ~isempty(idxNoPrePost)
    hold(hax, 'all');
    hl(2) = plot(hax, t_center_corr(idxNoPrePost), SeqOut.Correct_foffset(idxNoPrePost));
  end
  plot(hax, t_center_corr, ...
    (phaseDiffWeightedMean+phaseDiffWeightedStd) * SeqOut.AQ(iAQ(2)).fSample(tRepAQCorr) / 2/pi, ...
    'Color', get(hl(1), 'Color'), 'LineStyle', '--');
  plot(hax, t_center_corr, ...
    (phaseDiffWeightedMean-phaseDiffWeightedStd) * SeqOut.AQ(iAQ(2)).fSample(tRepAQCorr) / 2/pi, ...
    'Color', get(hl(1), 'Color'), 'LineStyle', '--');
  if ~isempty(idxNoPrePost)
    plot(hax, t_center_corr(idxNoPrePost), ...
      SeqOut.Correct_foffset(idxNoPrePost)+SeqOut.Correct_foffsetStd(idxNoPrePost), ...
      'Color', get(hl(2), 'Color'), 'LineStyle', '--');
    plot(hax, t_center_corr(idxNoPrePost), ...
      SeqOut.Correct_foffset(idxNoPrePost)-SeqOut.Correct_foffsetStd(idxNoPrePost), ...
      'Color', get(hl(2), 'Color'), 'LineStyle', '--');
    hold(hax, 'off');
  end
  ylabel(hax, 'f_{Offset} in Hz');
  xlabel(hax, 'time all in s');
  grid(hax, 'on');

elseif SeqOut.LoopPlot && SeqOut.CorrectPlotFrequency
  hFigure = figure(29);
  clf(hFigure, 'reset');
  hax = axes(hFigure);
  t_all = data(iData(2)).time_all(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep);
  t_center_corr = reshape(mean(t_all, 1), 1, []);
  errorbar(hax, t_center_corr, ...
    reshape(phaseDiffWeightedMean,1,[]) * SeqOut.AQ(iAQ(2)).fSample(tRepAQCorr) / 2/pi, ...
    reshape(phaseDiffWeightedStd,1,[]) * SeqOut.AQ(iAQ(2)).fSample(tRepAQCorr) / 2/pi);
  if ~isempty(idxNoPrePost)
    hold(hax, 'all');
    errorbar(hax, squeeze(mean(t_all(:,idxNoPrePost), 1)), SeqOut.Correct_foffset(idxNoPrePost).', SeqOut.Correct_foffsetStd(idxNoPrePost).');
    hold(hax, 'off');
  end
  ylabel(hax, 'f_{Offset} in Hz');
  xlabel(hax, 'time all in s');
  grid(hax, 'on');
  legend(hax, {'raw','filtered'});
  title(hax, 'Offset Frequency at FID');
end


if SeqOut.CorrectPlot
  % plot echo data
  if isscalar(imageAQWindow)
    dt = reshape(data(iData(1)).data(1:AQnSamples,imageAQWindow,imageUsetRep), AQnSamples, []);
    time = reshape(data(iData(1)).time_of_tRep(1:AQnSamples,imageAQWindow,imageUsetRep), AQnSamples, []);
  else
    dt = reshape(data(iData(1)).data(1:AQnSamples,sub2ind(szData(2:end), imageAQWindow, imageUsetRep)), AQnSamples, []);
    time = reshape(data(iData(1)).time_of_tRep(1:AQnSamples,sub2ind(szData(2:end), imageAQWindow, imageUsetRep)), AQnSamples, []);
  end
  hFigure = figure(30);
  clf(hFigure);
  clear hax
  hax(1) = subplot(4,1,1, 'Parent', hFigure);
  plot(hax(1), time, abs(dt));
  ylabel(hax(1), 'Amp in T');
  title(hax(1), {'Phase Analysis Echoes (image)', '(Only interpretable without phase encoding!)'});
  grid(hax(1), 'on');
  hax(2) = subplot(4,1,2, 'Parent', hFigure);
  plot(hax(2), time, unwrap2Dmiddle(angle(dt)));
  ylabel(hax(2), '\phi in rad');
  xlabel(hax(2), 'time of tRep in s');
  grid(hax(2), 'on');
  linkaxes(hax, 'x');

  clear hax
  hax(1) = subplot(4,1,3, 'Parent', hFigure);
  if isscalar(imageAQWindow)
    t_all = data(iData(1)).time_all(1:AQnSamples,imageAQWindow,imageUsetRep);
  else
    t_all = data(iData(1)).time_all(1:AQnSamples,sub2ind(szData(2:end), imageAQWindow, imageUsetRep));
  end
  t_center_image = reshape(mean(t_all, 1), 1, []);
  plot(hax(1), t_center_image, angle(dt(floor(end/2)+1,:)), '-x');
  hold(hax(1), 'on');
  plot(hax(1), t_center_image, ...
    angle(dt(floor(end/2)+1,:) .* ...
          reshape(exp(-1i*SeqOut.Correct_poffsetBoth(floor(end/2)+1,:)), 1, []) ), '-x');
  hold(hax(1), 'off');
  ylabel(hax(1), '\phi at center sample in rad');
  grid(hax(1), 'on');
  legend(hax(1), {'measured', 'corrected'});

  hax(2) = subplot(4,1,4, 'Parent', hFigure);
  mean_phase_diff = unwrap(get_MeanPhaseDiffWeighted(dt, 1), [], 2);
  plot(hax(2), t_center_image, ...
    reshape(mean_phase_diff, 1, []) * ...
    SeqOut.AQ(iAQ(1)).fSample(SeqOut.Correct_imageAQWindow(1),imageUsetRep(1)) / 2/pi, '-x');
  hold(hax(2), 'all');
  mean_phase_diff = unwrap(get_MeanPhaseDiffWeighted(squeeze(dt) .* exp(-1i*SeqOut.Correct_poffsetBoth), 1), [], 2);
  plot(hax(2), t_center_image, ...
    mean_phase_diff * ...
    SeqOut.AQ(iAQ(1)).fSample(SeqOut.Correct_imageAQWindow(1),imageUsetRep(1)) / 2/pi, '-x');
  hold(hax(2), 'off');
  ylabel(hax(2), 'f_{offset} in Hz');
  xlabel(hax(2), 'time all in s');
  legend_foffset = {'mean', 'mean corrected'};
  grid(hax(2), 'on');
  if ~isempty(idxNoPrePost)
    yyaxis(hax(2), 'right');
    t_all = data(iData(2)).time_all(1:SeqOut.Correct_nSamples,SeqOut.Correct_FidAQWindow,SeqOut.AQSlice(2).UsetRep);
    t_center_corr = reshape(mean(t_all, 1), 1, []);
    plot(hax(2), t_center_corr(idxNoPrePost), SeqOut.Correct_foffset(idxNoPrePost), '-x');
    grid(hax(2), 'on');
    ylabel(hax(2), 'f_{correction} in Hz');
    legend_foffset = [legend_foffset, {'correction'}];
  end
  legend(hax(2), legend_foffset);
  linkaxes(hax, 'x');
end


% correct encoded data for frequency offset
if isscalar(imageAQWindow)
  data(iData(1)).data(1:AQnSamples,imageAQWindow,imageUsetRep(:)) = ...
    data(iData(1)).data(1:AQnSamples,imageAQWindow,imageUsetRep(:)) .* ...
    reshape(exp(-1i*SeqOut.Correct_poffsetBoth), [AQnSamples, numel(imageUsetRep)]);
else
  data(iData(1)).data(1:AQnSamples,sub2ind(szData(2:end), imageAQWindow, imageUsetRep)) = ...
    data(iData(1)).data(1:AQnSamples,sub2ind(szData(2:end), imageAQWindow, imageUsetRep)) .* ...
    exp(-1i*SeqOut.Correct_poffsetBoth);
end

end
