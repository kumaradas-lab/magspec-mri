function [T, fitExpSettings] = fit_exp(data, time, hParent, CorrectFrequencyOffset, CorrectPhaseOffset, EndOffset, RingFilter)
%% Perform single and double exponential fit on data
%
%   [T, fitExpSettings] = fit_exp(data, time, fitExpSettings)
% alternatively (not recommended):
%   [T, fitExpSettings] = fit_exp(data, time, hParent, CorrectFrequencyOffset, CorrectPhaseOffset, EndOffset, RingFilter)
%
% The function performs a single and (optionally) a double exponential fit to
% the real part of data. By default, some corrections are applied before the fit
% (see below).
% The objective function for the single exponential fit is:
%                 f(t) = A * exp(t/T) [+ B]
% For the optional double exponential fit, this objective function is used:
%                 f(t) = A_1 * exp(t/T_1) + A_2 * exp(t/T_2) [+ B]
% Free parameters of the objective functions are written using capital letters.
%
% INPUT:
%   data      Vector or array with complex data to be fitted. In case this is an
%             array, "time" must also be passed as a array where the samples
%             of the same acquisition window run along the first dimension and
%             the acquisitions run along the third (or second) dimension.
%   time      Vector or array of the same size as "data" with the corresponding
%             time.
%
% The following input arguments are optional. If they are omitted or empty,
% default values are used. If the third argument is a structure fitExpSettings,
% it can have the following fields. Alternatively (but not recommended), you can
% pass some of them as separate arguments:
%
%   hParent   Handle to a figure which is used to display the measured data and
%             the results of the fit. If 0, no results are displayed (default:
%             0).
%   CorrectFrequencyOffset
%             Boolean value that enables the correction of the signal phase by
%             the offset to the Larmor frequency (default: 1 if data is a 3d
%             array, 0 otherwise).
%   CorrectMeanFrequencyOffset
%             Boolean value that selects whether the average frequency offset of
%             all acquisitions will be corrected (true), or whether the
%             frequency offset is corrected independently for each acquisition
%             window (false). It is always treated as false if
%             CorrectFrequencyDrift is not false. (Default: true)
%   CorrectFrequencyDrift
%             Boolean value or array the same size as "data" that enables the
%             correction of the signal phase by the linear drift with respect to
%             the Larmor frequency (default: 0). If the real time of the
%             measurement differs from the values in "time", pass the values
%             that should be used for the linear regression in this switch.
%   CorrectFrequencyReferenceTime
%             Scalar or array with the reference time for each echo (or FID) in
%             the input data. If omitted or empty, the time of the sample with
%             the highest amplitude (on average) for each echo (or FID) is used
%             as the reference time.
%   CorrectPhaseOffset
%             Boolean value that enables a rotation of the measurement signal
%             such that the average of the imaginary part of the first or last
%             part of the data is zero. "fitExpSettings.CorrectPhaseOffsetAtEnd"
%             determines whether the initial eighth part of the data or the last
%             eighth part of the data is used to determine the phase angle for
%             the correction. (Default: 1)
%   EndOffset
%             Boolean value that indicates whether the expected exponential
%             function does not converge to zero (after long times), i.e. a
%             recovery experiment. If true, the offset "B" of the exponential
%             curve is a free fit parameter. If false, "B" is forced to be zero,
%             i.e. the data is from a decay experiment. (Default: 0)
%   CorrectPhaseOffsetAtEnd
%             Boolean value that indicates whether the phase is corrected using
%             the last eighth part of the data (use this for recovery
%             experiments). Otherwise, the initial eighth part of the data is
%             used (decay experiments). (Default: fitExpSettings.EndOffset)
%   omitFirstnEchoes
%             Omit the first n Echoes in the fit. (Default: 0)
%   Factors
%             1xN vector with the correction factors for the amplitude of the
%             first N Echoes. A NaN value means that the respective Echo is
%             omitted. (Default: nan(1, fitExpSettings.omitFirstnEchoes) )
%   RingFilter
%             Boolean value that enables the correction of the "ringing effect"
%             that leads to the signal amplitudes being alternately too high and
%             too low by applying a sliding average of length 2 on the real part
%             of "data". (Default: 1).
%   RingFilterAbsAmpMeanPhase
%             Boolean value that indicates how the ring filter is applied on
%             complex data. If false, the average value of pairs of complex data
%             is calculated in a mathematical sense. If true, the average of the
%             absolute value of these pairs is calculated and given the average
%             phase of these respective pairs.
%             (Default: false)
%   SingleExp Boolean value that enables a single exponential fit of the data
%             (default: true).
%   SingleExpFitType
%             Integer that selects the fitting algorithm for the single
%             exponential fit:
%             0:  Least squares search using Matlab's fminsearch (default).
%             1:  Weighted linearized fit for the single exponential fit.
%                 Faster than fminsearch but non-positive values are ignored.
%                 This might lead to deviations if noisy data at low amplitudes
%                 crosses the base line.
%             2:  Regression using numeric integration (c.f.
%                 https://math.stackexchange.com/a/1897000 )
%                 This is not strictly a least squares fit but it is faster than
%                 fminsearch while also including data below the base line.
%   SingleExpFixTau
%             Integer value with a fixed decay constant in seconds. If this is
%             non-empty, the decay constant isn't determined from the data, but
%             the given decay constant is used. All other parameters are still
%             optimized. (Default: [])
%   DoubleExp Boolean value that enables an additional double exponential fit of
%             the data (default: true). If true, a single exponential fit is
%             also performed.
%   GaussExp  Boolean value that enables an additional fit of the data using a
%             exponential function (for the longer tau) and a Gaussian function
%             (for the shorter tau). This can be useful, e.g., when fitting T2
%             data containing a solid echo signal. If true, a single exponential
%             fit is also performed. (Default: false)
%   SingleGauss
%             Boolean value that enables a fit of the data with a single normal
%             (Gaussian) distribution function.
%             (Default: false)
%   hasFID    Boolean value. If true, the data in the first tRep is treated as
%             an FID (exp. decaying) signal (default: false).
%   FIDMean   Boolean value. If true, the average value of the FID is used at
%             the average time. Otherwise, each sample of the FID is used in the
%             exponential fits. (Default: false).
%   FrequencyFilter
%             A structure with the following fields:
%     UseFilter
%               Boolean value. If true, a bandpass filter is applied to the
%               signal at each Echo (default: false). For this, data and time
%               must be matrices.
%     CenterOffset
%               The offset of the center frequency (after mixing) in Hz
%               (Default: 0).
%     Bandwidth
%               The bandwidth of the filter in Hz. (No default value.)
%     ZeroFillFactor
%               Zero fill factor (in the time domain) for the band pass filter.
%               This must be an odd number. A value of 1 means no zero filling
%               (default: 5).
%
% OUTPUT:
%   T         Structure with the following fields:
%     dataPhaseCorrection
%             The phase in radians that was used to correct the data (see input
%             parameters "CorrectFrequencyOffset" and "CorrectPhaseOffset").
%     dataPhaseCorrected
%             The "data" with the phase correction applied.
%     dataPhaseCorrectedReal
%             The real part of the corrected data that was used for the
%             exponential fit.
%     timeCorrected
%             The time values used for the fit (see input argument
%             "RingFilter").
%     data    Same as input "data".
%     time    Same as input "time".
%     fitExpSettings
%             Same as input "fitExpSettings" or a structure with settings
%             generated by the alternative input arguments.
%     hFigure Same as input "hParent" (deprecated and will be removed in a
%             future version).
%     CorrectFrequencyOffset
%             Same as input "CorrectFrequencyOffset" (deprecated and will be
%             removed in a future version).
%     CorrectPhaseOffset
%             Same as input "CorrectPhaseOffset" (deprecated and will be removed
%             in a future version).
%     EndOffset
%             Same as input "EndOffset" (deprecated and will be removed in a
%             future version).
%     RingFilter
%             Same as input "RingFilter" (deprecated and will be removed in a
%             future version).
%     xstart  Start values of the free fit parameters that were used for the
%             single exponential fit.
%     xminSingle
%             Best fit values of these parameters.
%     errorSingle
%             RMS deviation of the single exponential fit to the data.
%     tau     Best fit single exponential decay constant.
%     functionTime
%             Vector with times with linear spacing used for plotting the fit
%             curves.
%     functionAmpSingle
%             Vector with best fit amplitudes of the single exponential curve
%             that matches "functionTime".
%     functionTimeLog
%             Vector with times with logarithmic spacing used for plotting the
%             fit curves.
%     functionAmpSingleLog
%             Vector with best fit amplitudes of the single exponential curve
%             that matches "functionTimeLog".
% If DoubleExp is true, the results include the following fields:
%     xstart5 Start values of the free fit parameters that were used for the
%             double exponential fit.
%     xminDouble
%             Best fit values of these parameters.
%     errorDouble
%             RMS deviation of the double exponential fit to the data.
%     tau1, tau2
%             Best fit double exponential decay constants.
%     tau1w, tau2w
%             Amplitudes ("weights") of the two exponential curves normalized to
%             1.
%     functionAmpDouble
%             Vector with best fit amplitudes of the double exponential curve
%             that matches "functionTime".
%     functionAmpDoubleLog
%             Vector with best fit amplitudes of the double exponential curve
%             that matches "functionTimeLog".
% If GaussExp is true, the results include the following fields:
%     xstartGaussExp
%             Start values of the free fit parameters that were used for the
%             exponential+Gauss fit.
%     xminGauss
%             Best fit values of these parameters.
%     errorGauss
%             RMS deviation of the exponential+Gaussian fit to the data.
%     tauG1, tauG2
%             Best fit decay constants for exponential (tauG1) and Gaussian
%             (tauG2) part.
%     tauG1w, tauG2w
%             Amplitudes ("weights") of the exponential (tauG1w) and the
%             Gaussian (tauG2w) curves normalized to 1.
%     functionAmpGauss
%             Vector with best fit amplitudes of the exponential + Gaussian
%             curve that matches "functionTime".
%     functionAmpGaussLog
%             Vector with best fit amplitudes of the exponential + Gaussian
%             curve that matches "functionTimeLog".
% If SingleGauss is true, the results include the following fields:
%     xstartGaussExp
%             Start values of the free fit parameters that were used for the
%             exponential+Gauss fit.
%     xminSingleGauss
%             Best fit values of these parameters.
%     errorSingleGauss
%             RMS deviation of the exponential+Gaussian fit to the data.
%     tauSingleGauss
%             Best fit decay constants for Gaussian distribution function.
%     functionAmpSingleGauss
%             Vector with best fit amplitudes of the Gaussian curve that matches
%             "functionTime".
%     functionAmpSingleGaussLog
%             Vector with best fit amplitudes of the Gaussian curve that matches
%             "functionTimeLog".
%
%   fitExpSettings
%           Structure (corresponding to the input argument of same name) with
%           the actually used settings.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% input check and default parameters
if nargin < 3 || isempty(hParent), hParent = 0; end
if isstruct(hParent)
  fitExpSettings = hParent;
else
  if nargin < 4 || isempty(CorrectFrequencyOffset), CorrectFrequencyOffset = (ndims(data) == 3) && (size(data, 1) > 1); end
  if nargin < 5 || isempty(CorrectPhaseOffset),     CorrectPhaseOffset = 1; end
  if nargin < 6 || isempty(EndOffset),              EndOffset = 0; end
  if nargin < 7 || isempty(RingFilter),             RingFilter = 1; end
  fitExpSettings.hParent = hParent;
  fitExpSettings.CorrectFrequencyOffset = CorrectFrequencyOffset;
  fitExpSettings.CorrectPhaseOffset = CorrectPhaseOffset;
  fitExpSettings.EndOffset = EndOffset;
  fitExpSettings.RingFilter = RingFilter;
end

fitExpSettings = set_EmptyField(fitExpSettings, 'hParent', 0);
if isemptyfield(fitExpSettings, 'CorrectFrequencyOffset')
  fitExpSettings.CorrectFrequencyOffset = ...
    (ndims(data) == 3) && (size(data, 1) > 1);
end
if isemptyfield(fitExpSettings, 'CorrectMeanFrequencyOffset')
  fitExpSettings.CorrectMeanFrequencyOffset = true;
end
if isemptyfield(fitExpSettings, 'CorrectFrequencyDrift')
  fitExpSettings.CorrectFrequencyDrift = 0;
end
if ~isfield(fitExpSettings, 'CorrectFrequencyReferenceTime')
  fitExpSettings.CorrectFrequencyReferenceTime = [];
end
fitExpSettings = set_EmptyField(fitExpSettings, 'CorrectPhaseOffset', 1);
fitExpSettings = set_EmptyField(fitExpSettings, 'EndOffset', 0);
fitExpSettings = set_EmptyField(fitExpSettings, 'CorrectPhaseOffsetAtEnd', fitExpSettings.EndOffset);
fitExpSettings = set_EmptyField(fitExpSettings, 'RingFilter', 1);
fitExpSettings = set_EmptyField(fitExpSettings, 'RingFilterAbsAmpMeanPhase', 0);
fitExpSettings = set_EmptyField(fitExpSettings, 'SingleExp', 1);
fitExpSettings = set_EmptyField(fitExpSettings, 'SingleExpFitType', 0);
if ~isfield(fitExpSettings, 'SingleExpFixTau')
  fitExpSettings.SingleExpFixTau = [];
end
fitExpSettings = set_EmptyField(fitExpSettings, 'DoubleExp', 1);
fitExpSettings = set_EmptyField(fitExpSettings, 'GaussExp', false);
fitExpSettings = set_EmptyField(fitExpSettings, 'SingleGauss', false);
fitExpSettings = set_EmptyField(fitExpSettings, 'hasFID', 0);
fitExpSettings = set_EmptyField(fitExpSettings, 'FIDMean', 0);
fitExpSettings = set_EmptyField(fitExpSettings, 'omitFirstnEchoes', 0);
fitExpSettings = set_EmptyField(fitExpSettings, 'Factors', []);
fitExpSettings = set_EmptyField(fitExpSettings, 'PhaseOffset', 0);

fitExpSettings = set_EmptyField(fitExpSettings, 'FrequencyFilter', struct());
fitExpSettings.FrequencyFilter = set_EmptyField(fitExpSettings.FrequencyFilter, 'UseFilter', false);
if fitExpSettings.FrequencyFilter.UseFilter
  fitExpSettings.FrequencyFilter = set_EmptyField(fitExpSettings.FrequencyFilter, 'CenterOffset', 0);
  if isemptyfield(fitExpSettings.FrequencyFilter, 'Bandwidth')
    error('PD:fit_exp:UseFilterParameters', ['If FrequencyFilter.UseFilter is true, ' ...
      'FrequencyFilter.Bandwidth must be set.']);
  end
  fitExpSettings.FrequencyFilter = set_EmptyField(fitExpSettings.FrequencyFilter, 'ZeroFillFactor', 5);
  if (fitExpSettings.FrequencyFilter.ZeroFillFactor-1)/2 ...
      ~= round((fitExpSettings.FrequencyFilter.ZeroFillFactor-1)/2)
    error('PD:fit_exp:UseFilterZeroFillFactorNotOdd', ...
      'FrequencyFilter.ZeroFillFactor must be an odd number.');
  end
end

if fitExpSettings.hParent == 1, fitExpSettings.hParent = 41; end


% reshape 2d matrices to 3d arrays
% (assuming that the 2d matrix is samples x AQ windows)
% FIXME: Is this safe?
if ismatrix(data) && size(data, 2) > 1
  data = reshape(data, [size(data, 1), 1, size(data, 2)]);
end
if ismatrix(time) && size(time, 2) > 1
  time = reshape(time, [size(time, 1), 1, size(time, 2)]);
end
if ismatrix(fitExpSettings.CorrectFrequencyDrift) && size(fitExpSettings.CorrectFrequencyDrift, 2) > 1
  fitExpSettings.CorrectFrequencyDrift = ...
    reshape(fitExpSettings.CorrectFrequencyDrift, ...
    [size(fitExpSettings.CorrectFrequencyDrift, 1), 1, size(fitExpSettings.CorrectFrequencyDrift, 2)]);
end


% remove tReps without acquisition windows
noData = all(isnan(data), 1);
noTime = all(isnan(time), 1);
data(:,:,noData|noTime) = [];
time(:,:,noData|noTime) = [];


%% data corrections
% Correct frequency offset (phase drift)
if fitExpSettings.CorrectFrequencyOffset
  if ndims(data) == 3 && size(data, 1) < 2
    % frequency drift over (mean value of) several echoes
    if any(abs(diff(time) - diff(time(1:2))) > 8e-9)
      % FIXME: Only works if time raster is linear.
      warning('PD:fit_exp:FreqDriftNonLinear', ...
        ['Frequency offset (and drift) can only be corrected if echoes are ', ...
        'spaced linearly.']);
      fitExpSettings.CorrectFrequencyOffset = 0;
      fitExpSettings.CorrectFrequencyDrift = 0;
      dataFrequencyCorrected = data;
    else
      phasePerSec = get_MeanPhaseDiffWeighted(data, 3, 'omitnan') / diff(time(1:2));
      dataFrequencyCorrected = data .* exp(-1i * phasePerSec * time);
      if numel(fitExpSettings.CorrectFrequencyDrift) > 1 || fitExpSettings.CorrectFrequencyDrift
        warning('PD:fit_exp:FreqDriftNotCorrectable', ...
          ['Frequency drift can only be corrected if "data" and "time" ' ...
          'are 3d matrices that are non-singleton in the first dimension.']);
        fitExpSettings.CorrectFrequencyDrift = 0;
      end
    end
  else
    % determine frequency drift within each echo
    % (and use that to calculate the drift in the train)

    % get weighted average of phase drift in each echo
    [phasePerSec, phasePerSecStd] = get_MeanPhaseDiffWeighted(data, 1, 'omitnan');
    phasePerSec = bsxfun(@rdivide, phasePerSec, diff(time(1:2,1,:), 1, 1));
    phasePerSecStd = bsxfun(@rdivide, phasePerSecStd, diff(time(1:2,1,:), 1, 1));
    if isempty(fitExpSettings.CorrectFrequencyReferenceTime)
      % Use sample with highest amplitude in most tReps as the reference time.
      % FIXME: This should either be the echo time (for echoes) or the center of
      % the excitation pulse (for FIDs).
      isMaxDataCount = sum(bsxfun(@eq, abs(data), max(abs(data), [], 1)), 3);
      refFreqTime = time(find(isMaxDataCount == max(isMaxDataCount), 1, 'first'),:,:);
      if fitExpSettings.hasFID
        if fitExpSettings.FIDMean
          refFreqTime(1) = mean(time(:,:,1), 1, 'omitnan');
        else
          refFreqTime(1) = time(1);
        end
      end
    elseif isscalar(fitExpSettings.CorrectFrequencyReferenceTime)
      refFreqTime = ...
        fitExpSettings.CorrectFrequencyReferenceTime ...
        .* ones(size(phasePerSec));
    elseif all(size(fitExpSettings.CorrectFrequencyReferenceTime) ...
        == size(phasePerSec))
      refFreqTime = fitExpSettings.CorrectFrequencyReferenceTime;
    else
      error('PD:fit_exp:CorrectFrequencyReferenceTimeInvalidSize', ...
        ['CorrectFrequencyReferenceTime must be a scalar or ', ...
        'must match the number of echoes/FIDs in all trains.']);
    end
    % Use mean amplitude and standard error of mean phase drift as weights.
    weights = (double(mean(abs(data), 1, 'omitnan'))./phasePerSecStd).^2;
    validData = ~isnan(phasePerSec);
    validData(1,1,fitExpSettings.hasFID + (1:fitExpSettings.omitFirstnEchoes)) = false;
    weights = weights(validData);
    if numel(fitExpSettings.CorrectFrequencyDrift) > 1 || fitExpSettings.CorrectFrequencyDrift
      % Correct linear frequency drift as well.
      if numel(fitExpSettings.CorrectFrequencyDrift) > 1
        refFreqDriftTime = ...
          fitExpSettings.CorrectFrequencyDrift(find(isMaxDataCount == max(isMaxDataCount), 1, 'first'),:,:);
      else
        refFreqDriftTime = refFreqTime;
      end
      P = (diag(weights(:)) * [squeeze(refFreqDriftTime(validData)), ones(sum(validData), 1)]) \ ...
        (weights(:) .* squeeze(phasePerSec(validData)));
      % P = polyfit(refFreqDriftTime(:), phasePerSec(:), 1);
      dataFrequencyCorrected = bsxfun(@times, data, ...
        exp(-1i * bsxfun(@times, polyval(P, refFreqDriftTime), bsxfun(@minus, time, refFreqTime))));
    else
      if fitExpSettings.CorrectMeanFrequencyOffset
        % Correct all acquisition windows with the (same) average frequency
        % offset.
        meanPhasePerSec = sum(phasePerSec(validData).*weights) ./ sum(weights);
      else
        % Use frequency offset determined independently for each acquisition
        % window.
        meanPhasePerSec = phasePerSec;
      end
      dataFrequencyCorrected = data ...
        .* exp(-1i * bsxfun(@times, meanPhasePerSec, bsxfun(@minus, time, refFreqTime)));
    end
  end
else
  dataFrequencyCorrected = data;
end


% band pass filter
if fitExpSettings.FrequencyFilter.UseFilter
  if ndims(data) ~= 3 || ndims(time) ~= 3
      % (ndims(time) ~= 3 && ndims(fitExpSettings.CorrectFrequencyDrift) ~= 3)
    warning('PD:fit_exp:BandpassFilterNotCorrectable', ...
      ['Frequency filter can only be applied if "data" and "time" are 3d matrices ' ...
      'that are non-singleton in the first dimension.']);
    dataBandPass = dataFrequencyCorrected;
  else
    % FIXME: It would also be possible to use
    % fitExpSettings.CorrectFrequencyDrift as the time matrix. This isn't
    % implemented yet.
    dataFFT = fftshift(fft(cat(1, ...
      zeros(size(dataFrequencyCorrected, 1)*(fitExpSettings.FrequencyFilter.ZeroFillFactor-1)/2, size(dataFrequencyCorrected, 2), size(dataFrequencyCorrected, 3)), ...
      dataFrequencyCorrected, ...
      zeros(size(dataFrequencyCorrected, 1)*(fitExpSettings.FrequencyFilter.ZeroFillFactor-1)/2, size(dataFrequencyCorrected, 2), size(dataFrequencyCorrected, 3))), ...
      [], 1));
    freqFFT = get_FFTGrid(1/diff(time(1:2,1,2)), size(time, 1)*fitExpSettings.FrequencyFilter.ZeroFillFactor);
    dataFFT(abs(freqFFT-fitExpSettings.FrequencyFilter.CenterOffset) > fitExpSettings.FrequencyFilter.Bandwidth,:,:) = 0;
    dataBandPassZ = ifft(ifftshift(dataFFT), [], 1);
    k = size(dataBandPassZ,1)/fitExpSettings.FrequencyFilter.ZeroFillFactor/2;
    dataBandPass = dataBandPassZ(ceil(end/2+(-k:(k-1))),:,:);
  end
else
  dataBandPass = dataFrequencyCorrected;
end


numSamplesFID = 1;
if ndims(dataBandPass) == 3
  if fitExpSettings.hasFID && ~fitExpSettings.FIDMean
    numSamplesFID = sum(~isnan(dataBandPass(:,1,1)));
    dataFrequencyCorrectedMean = [dataBandPass(~isnan(dataBandPass(:,1,1)),1,1); ...
      reshape(mean(dataBandPass(:,:,2:end), 1, 'omitnan'), [], 1)];
  else
    dataFrequencyCorrectedMean = reshape(mean(dataBandPass, 1, 'omitnan'), [], 1);
  end
else
  dataFrequencyCorrectedMean = dataBandPass;
end
if ndims(time) == 3
  if fitExpSettings.hasFID && ~fitExpSettings.FIDMean
    timeUsed = [time(~isnan(time(:,1,1)),1,1); ...
      reshape(mean(time(:,:,2:end), 1, 'omitnan'), [], 1)];
  else
    timeUsed = reshape(mean(time, 1, 'omitnan'), [], 1);
  end
else
  timeUsed = time;
end

% (manually) correct amplitude
if ~isempty(fitExpSettings.Factors)
  nFactors = min(numel(fitExpSettings.Factors), size(dataFrequencyCorrectedMean, 1));
  dataFrequencyCorrectedMean((1:nFactors)+fitExpSettings.hasFID*numSamplesFID) = ...
    dataFrequencyCorrectedMean((1:nFactors)+fitExpSettings.hasFID*numSamplesFID) .* ...
    fitExpSettings.Factors(1:nFactors);
end

dataFrequencyCorrectedMean=dataFrequencyCorrectedMean .* exp(-1i*fitExpSettings.PhaseOffset);

% correct phase offset
if fitExpSettings.CorrectPhaseOffset
  % [~, maxI] = max(abs(dataFrequencyCorrectedMean));
  % if maxI < numel(dataFrequencyCorrectedMean)/2
  if fitExpSettings.CorrectPhaseOffsetAtEnd
    % inversion or saturation recovery approaching full recovery
    % take last eighth part
    selectDataIndex = max(ceil(numel(dataFrequencyCorrectedMean)*7/8), 2):numel(dataFrequencyCorrectedMean);
    selectDataIndex(1:min(find(abs(dataFrequencyCorrectedMean(selectDataIndex)) < abs(dataFrequencyCorrectedMean(end))/2, 1, 'last'), end-2)) = [];
    addPhase = 0;
  else
    % take initial eighth part
    selectDataIndex = [ones(fitExpSettings.hasFID, 1)*(1:numSamplesFID), ...
      (fitExpSettings.hasFID*numSamplesFID+fitExpSettings.omitFirstnEchoes+1):max(floor((numel(dataFrequencyCorrectedMean)-fitExpSettings.omitFirstnEchoes)/8+1), fitExpSettings.omitFirstnEchoes+2)];
    selectDataIndex(max(find(abs(dataFrequencyCorrectedMean(selectDataIndex)) < abs(dataFrequencyCorrectedMean(1))/4, 1, 'first'), 2):end) = [];
    if fitExpSettings.EndOffset
      % inversion recovery with only short recovery times compared to T1
      addPhase = pi;
    else
      % decay (T2)
      addPhase = 0;
    end
  end

  pOffset = mean(unwrap(angle(dataFrequencyCorrectedMean(selectDataIndex))))+addPhase;
else
  pOffset = 0;
end

T.dataPhaseCorrection = pOffset;
T.dataPhaseCorrected = dataFrequencyCorrectedMean.*exp(-1i*pOffset);

% remove NaN values
% FIXME: If "fit_exp" is called repeatedly, e.g., for parts of a T1-T2
%        measurement, removing NaN values only for some of them might lead to
%        inconsistent matrix dimensions.
%        Consider handling NaN values differently to keep consistent matrix
%        dimensions.
timeUsed(isnan(T.dataPhaseCorrected)) = [];
T.dataPhaseCorrected(isnan(T.dataPhaseCorrected)) = [];

if isempty(T.dataPhaseCorrected)
  % continue with dummy data
  % FIXME: Consider removing this work-around when the above FIXME note has been
  %        addressed.
  timeUsed = NaN(size(dataFrequencyCorrectedMean));
  T.dataPhaseCorrected = timeUsed;
end

% ring filter (reduce difference between odd and even Echoes)
if fitExpSettings.RingFilter
  if fitExpSettings.RingFilterAbsAmpMeanPhase
    if fitExpSettings.hasFID
      if fitExpSettings.FIDMean
        T.dataPhaseCorrected = [T.dataPhaseCorrected(1); conv(abs(T.dataPhaseCorrected(2:end)), [0.5,0.5], 'valid').*exp(1i*angle(conv(T.dataPhaseCorrected(2:end), [0.5,0.5], 'valid')))];
        timeCorrected = [timeUsed(1); conv(timeUsed(2:end), [0.5,0.5], 'valid')];
      else
        T.dataPhaseCorrected = [T.dataPhaseCorrected(1:numSamplesFID); conv(abs(T.dataPhaseCorrected(numSamplesFID+1:end)), [0.5,0.5], 'valid').*exp(1i*angle(conv(T.dataPhaseCorrected(numSamplesFID+1:end), [0.5,0.5], 'valid')))];
        timeCorrected = [timeUsed(1:numSamplesFID); conv(timeUsed(numSamplesFID+1:end), [0.5,0.5], 'valid')];
      end
    else
      T.dataPhaseCorrected = conv(abs(T.dataPhaseCorrected), [0.5,0.5], 'valid').*exp(1i*angle(conv(T.dataPhaseCorrected, [0.5,0.5], 'valid')));
      timeCorrected = conv(timeUsed, [0.5,0.5], 'valid');
    end
  else
    if fitExpSettings.hasFID
      if fitExpSettings.FIDMean
        T.dataPhaseCorrected = [T.dataPhaseCorrected(1); conv(T.dataPhaseCorrected(2:end), [0.5,0.5], 'valid')];
        timeCorrected = [timeUsed(1); conv(timeUsed(2:end), [0.5,0.5], 'valid')];
      else
        T.dataPhaseCorrected = [T.dataPhaseCorrected(1:numSamplesFID); conv(T.dataPhaseCorrected(numSamplesFID+1:end), [0.5,0.5], 'valid')];
        timeCorrected = [timeUsed(1:numSamplesFID); conv(timeUsed(numSamplesFID+1:end), [0.5,0.5], 'valid')];
      end
    else
      T.dataPhaseCorrected = conv(T.dataPhaseCorrected, [0.5,0.5], 'valid');
      timeCorrected = conv(timeUsed, [0.5,0.5], 'valid');
    end
  end
else
  timeCorrected = timeUsed;
end
T.dataPhaseCorrectedReal = real(T.dataPhaseCorrected);


%% fit (and plot) corrected data
% remove omitted Echoes from data for fit
dataPhaseCorrectedFit = T.dataPhaseCorrected(fitExpSettings.omitFirstnEchoes+1:end);
timeCorrectedFit = timeCorrected(fitExpSettings.omitFirstnEchoes+1:end);


% plot data used for fit
if ishghandle(fitExpSettings.hParent, 'figure') || ...
    ishghandle(fitExpSettings.hParent, 'uipanel') || fitExpSettings.hParent
  if ishghandle(fitExpSettings.hParent, 'uipanel')
    delete(get(fitExpSettings.hParent, 'Children'));
  else
    fitExpSettings.hParent = figure(fitExpSettings.hParent);
    clf(fitExpSettings.hParent)
  end
  ax(1) = subplot(2,1,1, 'Parent', fitExpSettings.hParent);
  if size(time, 1) > 1
    timeNAN = cat(1, time, nan(1,size(time,2),size(time,3),size(time,4)));
    dataNAN = cat(1, data, nan(1,size(data,2),size(data,3),size(data,4)));
    dataFrequencyCorrectedNAN = cat(1, dataBandPass, nan(1,size(dataBandPass,2),size(dataBandPass,3),size(dataBandPass,4)));
  else
    timeNAN = time;
    dataNAN = data;
    dataFrequencyCorrectedNAN = dataBandPass;
  end
  hLines1 = plot(ax(1), ...
    timeNAN(:), abs(dataFrequencyCorrectedNAN(:)), ...
    timeNAN(:), real(dataFrequencyCorrectedNAN(:)*exp(-1i*pOffset)), ...
    timeNAN(:), imag(dataFrequencyCorrectedNAN(:)*exp(-1i*pOffset)), ...
    timeCorrectedFit, real(dataPhaseCorrectedFit));
  set(hLines1(4), 'LineStyle', 'none', 'Marker', 'x');
  hold(ax(1), 'on');
  set(ax(1), 'ColorOrderIndex', 1);
  %   plot(ax(1), timeNAN(round(end/2),:), abs(dataFrequencyCorrectedNAN(round(end/2),:)), '.', ...
  %               timeNAN(round(end/2),:), real(dataFrequencyCorrectedNAN(round(end/2),:)*exp(-1i*pOffset)), '.', ...
  %               timeNAN(round(end/2),:), imag(dataFrequencyCorrectedNAN(round(end/2),:)*exp(-1i*pOffset)), '.');
  tCenter=reshape(mean(timeNAN,1,'omitnan'),1,[]);
  dMean=reshape(mean(dataFrequencyCorrectedNAN,1,'omitnan'),1,[])*exp(-1i*pOffset);
  plot(ax(1), ...
    tCenter,  abs(dMean), '.', ...
    tCenter, real(dMean), '.', ...
    tCenter, imag(dMean), '.');
  legendCell = {'abs', 'real', 'imag', 'amp for fit'};
  legend(hLines1, legendCell);
  xlabel(ax(1), 'time in s');
  ylabel(ax(1), 'amplitude in T');
  title_str = 'Corrected Measurement Data';


  ax(2) = subplot(2,1,2, 'Parent', fitExpSettings.hParent);
%   if fitExpSettings.CorrectFrequencyOffset
%     plot(ax(2), timeNAN(:), unwrap(angle(dataNAN(:))), 'r', ...
%       timeNAN(:), unwrap(angle(dataFrequencyCorrectedNAN(:))), 'b', ...
%       timeCorrectedFit, unwrap(angle(dataPhaseCorrectedFit)), 'kx-');
%     legend(ax(2), {'phase of data', 'phase frequency corrected', 'phase for fit'});
%   else
%     plot(ax(2), timeNAN(:), unwrap(angle(dataNAN(:))), 'r', ...
%       timeCorrectedFit, unwrap(angle(dataPhaseCorrectedFit)), 'kx-');
%     legend(ax(2), {'phase of data', 'phase for fit'});
%   end
  ax(2) = subplot(2,1,2, 'Parent', fitExpSettings.hParent);
  if fitExpSettings.CorrectFrequencyOffset
    plot(ax(2), ...
      timeNAN(:), angle(dataNAN(:)), 'r', ...
      timeNAN(:), angle(dataFrequencyCorrectedNAN(:)*exp(-1i*pOffset)), 'b', ...
      tCenter, angle(dMean), '.b', ...
      timeCorrectedFit, (angle(dataPhaseCorrectedFit)), 'kx-');
    legend(ax(2), {'phase of data', 'phase frequency corrected', 'phase for fit'});
  else
    plot(ax(2), timeNAN(:), (angle(dataNAN(:))), 'r', ...
      timeCorrectedFit, (angle(dataPhaseCorrectedFit)), 'kx-');
    legend(ax(2), {'phase of data', 'phase for fit'});
  end
  xlabel(ax(2), 'time in s')
  ylabel(ax(2), 'angle in rad')
  % linkaxes(ax,'x');
  % set(zoom(ax(1)),'Motion','both','Enable','on');
end


% Fit the corrected data
if length(T.dataPhaseCorrectedReal) - fitExpSettings.omitFirstnEchoes < 2+double(fitExpSettings.EndOffset)
  error('PD:fit_exp:NotEnoughData', 'Not enough data to apply the exponential fit.');
end


if fitExpSettings.EndOffset
  % offset
  T.xstart(1) = double(max(abs(T.dataPhaseCorrectedReal([ones(1, fitExpSettings.hasFID), (fitExpSettings.hasFID+fitExpSettings.omitFirstnEchoes+1):end])))); % offset
  % amplitude
  if T.dataPhaseCorrectedReal(fitExpSettings.hasFID+fitExpSettings.omitFirstnEchoes+1) < -T.xstart(1)/10
    T.xstart(2) = -2*T.xstart(1); % amplitude inversion
  else
    T.xstart(2) = -1*T.xstart(1); % amplitude saturation
  end
else
  % offset
  T.xstart(1) = 0;
  % amplitude
  if fitExpSettings.hasFID
    T.xstart(2) = T.dataPhaseCorrectedReal(1);
  else
    T.xstart(2) = T.dataPhaseCorrectedReal(fitExpSettings.omitFirstnEchoes+1);
  end
end

% decay constant
% Find index of first sample for which the amplitude decayed by 1/e.
timeIndex = find((T.dataPhaseCorrectedReal-T.xstart(1))<=(T.xstart(2)-T.xstart(1))/exp(1), 1, 'first');
if isempty(timeIndex)
  % If that point was not yet reached, estimate an approximate decay constant
  % from the amplitudes of the first and the last data point.
  T.xstart(3) = - (timeCorrected(1) - timeCorrected(end)) ...
    / log( (T.dataPhaseCorrectedReal(1) - T.xstart(1)) ...
           / (T.dataPhaseCorrectedReal(end) - T.xstart(1)));
  if ~isreal(T.xstart(3)) || T.xstart(3) < 0 || ~isfinite(T.xstart(3))
    % If that results in an unreasonable decay constant, just use the time of
    % the last data point as an estimate for the decay constant.
    T.xstart(3) = timeCorrected(end);
  end
else
  T.xstart(3) = timeCorrected(timeIndex);
end

AmpScale = T.xstart(3)/T.xstart(2);
dataPhaseCorrectedRealFit = double(real(dataPhaseCorrectedFit)*AmpScale);
T.xstart(1:2) = T.xstart(1:2)*AmpScale;


if fitExpSettings.SingleExp || fitExpSettings.DoubleExp || fitExpSettings.GaussExp
  maxerror = max(abs(dataPhaseCorrectedRealFit)/1e2)^2 * numel(dataPhaseCorrectedRealFit);
  % single exponential fit
  if fitExpSettings.EndOffset
    switch fitExpSettings.SingleExpFitType
      case 1  % Weighted Linearized Fit
        if isempty(fitExpSettings.SingleExpFixTau)
          % estimate offset C of the exponential decay: y = A*exp(-t/tau) + C
          maxAmp = max(dataPhaseCorrectedRealFit) - min(dataPhaseCorrectedRealFit);
          C_estimated = min(dataPhaseCorrectedRealFit) - maxAmp/100;
          % Solve v(t) = y(t) - C = A*e^(-t/tau) <=> log(v(t)) = log(A) - 1/tau * t
          X = [ones(length(timeCorrectedFit), 1), reshape(timeCorrectedFit, length(timeCorrectedFit), 1)];

          % search best offset C
          C = fminsearch(@(c) linearizedFit(timeCorrectedFit, dataPhaseCorrectedRealFit, X, c), ...
                         C_estimated, optimset('TolX', 1e-4*abs(maxAmp)));

          [T.errorSingle, coeff] = linearizedFit(timeCorrectedFit, dataPhaseCorrectedRealFit, X, C);
          T.xminSingle = [C, coeff];
        else
          % estimate offset C of the exponential decay: y = A*exp(-t/tau) + C
          maxAmp = max(dataPhaseCorrectedRealFit) - min(dataPhaseCorrectedRealFit);
          C_estimated = min(dataPhaseCorrectedRealFit) - maxAmp/100;
          % Solve v(t) = y(t) - C = A*e^(-t/tau) <=> log(v(t)) = log(A) - 1/tau * t
          X = ones(length(timeCorrectedFit), 1);

          % search best offset C
          C = fminsearch(@(c) linearizedFitFixDecay(timeCorrectedFit, dataPhaseCorrectedRealFit, X, c, ...
                                                    fitExpSettings.SingleExpFixTau), ...
                         C_estimated, optimset('TolX', 1e-4*abs(maxAmp)));

          [T.errorSingle, coeff] = ...
            linearizedFitFixDecay(timeCorrectedFit, dataPhaseCorrectedRealFit, X, ...
                                  C, fitExpSettings.SingleExpFixTau);
          T.xminSingle = [C, coeff];
        end

      case 2  % Numerical Integration
        [T.errorSingle, T.xminSingle] = ...
          integrNumExp(timeCorrectedFit(:).', dataPhaseCorrectedRealFit(:).', ...
                       true, fitExpSettings.SingleExpFixTau);

      otherwise  % least squares search
        if isempty(fitExpSettings.SingleExpFixTau)
          [T.xminSingle, T.errorSingle, T.exitflagSingle, T.outputSingle] = ...
            fminsearch(@(x) lsqFitExp(timeCorrectedFit, dataPhaseCorrectedRealFit, x), T.xstart, ...
                       optimset('MaxFunEvals', 1000, 'TolX', 1e-4*abs(T.xstart(2)), 'TolFun', maxerror, 'Display', 'off'));
        else
          [T.xminSingle, T.errorSingle, T.exitflagSingle, T.outputSingle] = ...
            fminsearch(@(x) lsqFitExpFixDecay(timeCorrectedFit, dataPhaseCorrectedRealFit, x, ...
                                              fitExpSettings.SingleExpFixTau, fitExpSettings.EndOffset), ...
                       T.xstart(1:2), ...
                       optimset('MaxFunEvals', 1000, 'TolX', 1e-4*abs(T.xstart(2)), 'TolFun', maxerror, 'Display', 'off'));
          T.xminSingle = [T.xminSingle, fitExpSettings.SingleExpFixTau];
        end
    end
  else
    switch fitExpSettings.SingleExpFitType
      case 1 % Weighted Linearized Fit
        if isempty(fitExpSettings.SingleExpFixTau)
          % linearize exponential function: y(t) = A*e^(-t/tau) <=> log(y) = log(A) - t/tau
          X = [ones(length(timeCorrectedFit), 1), reshape(timeCorrectedFit, length(timeCorrectedFit), 1)];
          [T.errorSingle, coeff] = linearizedFit(timeCorrectedFit, dataPhaseCorrectedRealFit, X, 0);
          T.xminSingle = [0, coeff];
        else
          % linearize exponential function: y(t) = A*e^(-t/tau) <=> log(y) = log(A) - t/tau
          X = ones(length(timeCorrectedFit), 1);
          [T.errorSingle, coeff] = ...
            linearizedFitFixDecay(timeCorrectedFit, dataPhaseCorrectedRealFit, X, ...
                                  0, fitExpSettings.SingleExpFixTau);
          T.xminSingle = [0, coeff];
        end

      case 2 % Numerical Integration
        [T.errorSingle, coeff] = ...
          integrNumExp(timeCorrectedFit(:).', dataPhaseCorrectedRealFit(:).', ...
                       false, fitExpSettings.SingleExpFixTau);
        T.xminSingle = [0, coeff];

      otherwise % least squares search
        if isempty(fitExpSettings.SingleExpFixTau)
          [T.xminSingle, T.errorSingle, T.exitflagSingle, T.outputSingle] = ...
            fminsearch(@(x) lsqFitExp(timeCorrectedFit, dataPhaseCorrectedRealFit, x), T.xstart(2:3), ...
                       optimset('MaxFunEvals', 1000, 'TolX', 1e-4*abs(T.xstart(2)), 'TolFun', maxerror, 'Display', 'off'));
          T.xminSingle = [0, T.xminSingle(1:2)];
        else
          [T.xminSingle, T.errorSingle, T.exitflagSingle, T.outputSingle] = ...
            fminsearch(@(x) lsqFitExpFixDecay(timeCorrectedFit, dataPhaseCorrectedRealFit, x, ...
                                              fitExpSettings.SingleExpFixTau, fitExpSettings.EndOffset), ...
                       T.xstart(2), ...
                       optimset('MaxFunEvals', 1000, 'TolX', 1e-4*abs(T.xstart(2)), 'TolFun', maxerror, 'Display', 'off'));
          T.xminSingle = [0, T.xminSingle(1), fitExpSettings.SingleExpFixTau];
        end
    end
  end

  % rescale to original units
  T.xminSingle(1:2) = T.xminSingle(1:2) / AmpScale;
  T.errorSingle = nan; %sqrt(T.errorSingle/AmpScale^2 / numel(timeCorrectedFit));

  % interpolate fitted function
  T.functionTime = linspace(timeCorrectedFit(1)*0, timeCorrectedFit(end)*2, max(2048,numel(timeCorrectedFit)*2));
  [~, T.functionAmpSingle] = lsqFitExp(T.functionTime, length(T.functionTime), T.xminSingle);
  T.functionTimeLog = logspace(log10(timeCorrectedFit(1)), log10(timeCorrectedFit(end)*2), max(2048,numel(timeCorrectedFit)*2));
  [~, T.functionAmpSingleLog] = lsqFitExp(T.functionTimeLog, length(T.functionTime), T.xminSingle);
  if ishghandle(fitExpSettings.hParent, 'figure') || ...
    ishghandle(fitExpSettings.hParent, 'uipanel') || fitExpSettings.hParent
    hLines1(end+1) = plot(ax(1), T.functionTime, T.functionAmpSingle, 'b--');
    if T.xminSingle(3) > 10
      title_str = sprintf('\\tau = %.3f s', T.xminSingle(3));
    elseif T.xminSingle(3) > 10e-3
      title_str = sprintf('\\tau = %.3f ms', T.xminSingle(3)*1e3);
    else
      title_str = sprintf('\\tau = %.3f %cs', T.xminSingle(3)*1e6, char(181));
    end
    legendCell = [legendCell, 'single expo.'];
  end
  T.tau = T.xminSingle(3);

end

if fitExpSettings.DoubleExp
  % double exponential fit
  T.xstart5(1) = T.xminSingle(1)*AmpScale;                  % offset
  T.xstart5(2) = T.xminSingle(2)/2*AmpScale;                % amp1
  T.xstart5(3) = T.xminSingle(3)*2;                         % tau1
  T.xstart5(4) = T.xminSingle(2)/2*AmpScale;                % amp2
  T.xstart5(5) = max(T.xminSingle(3)/2,2*timeCorrectedFit(1)); % tau2

  if fitExpSettings.EndOffset
    [T.xminDouble, T.errorDouble, T.exitflagDouble, T.outputDouble] = ...
      fminsearch(@(x) lsqFitExp(timeCorrectedFit, dataPhaseCorrectedRealFit, x), T.xstart5, ...
      optimset('MaxFunEvals', 2000, 'TolX', 1e-4*abs(T.xstart(2)), 'TolFun', maxerror, 'Display', 'off'));
  else
    [T.xminDouble, T.errorDouble, T.exitflagDouble, T.outputDouble] = ...
      fminsearch(@(x) lsqFitExp(timeCorrectedFit, dataPhaseCorrectedRealFit, x), T.xstart5(2:5), ...
      optimset('MaxFunEvals', 1000, 'TolX', 1e-4*abs(T.xstart(2)), 'TolFun', maxerror, 'Display', 'off'));
    T.xminDouble = [0 T.xminDouble(1:4)];
  end
  T.xminDouble([1:2,4])=T.xminDouble([1:2,4])/AmpScale;
  T.errorDouble = sqrt(T.errorDouble/AmpScale^2 / numel(timeCorrectedFit));

  % interpolate fitted function
  [~, T.functionAmpDouble] = lsqFitExp(T.functionTime, length(T.functionTime), T.xminDouble);
  [~, T.functionAmpDoubleLog] = lsqFitExp(T.functionTimeLog, length(T.functionTime), T.xminDouble);
  if abs(T.xminDouble(2)) < abs(T.xminDouble(4))
    temp1 = T.xminDouble(2:3);
    T.xminDouble(2:3) = T.xminDouble(4:5);
    T.xminDouble(4:5) = temp1;
  end
  % decay constants
  T.tau1 = T.xminDouble(3);
  T.tau2 = T.xminDouble(5);
  % amplitude "weight"
  T.tau1w = T.xminDouble(2)/(T.xminDouble(2)+T.xminDouble(4));
  T.tau2w = T.xminDouble(4)/(T.xminDouble(2)+T.xminDouble(4));

  if ishghandle(fitExpSettings.hParent, 'figure') || ...
      ishghandle(fitExpSettings.hParent, 'uipanel') || fitExpSettings.hParent
    hLines1(end+1) = plot(ax(1), T.functionTime, T.functionAmpDouble, 'b-.');
    if T.xminDouble(3) > 10
      title_str = [title_str '', ...
        ' || \tau_1 = ' num2str(T.xminDouble(3), '%.3f') ' s @ ' num2str(T.tau1w*100, '%.1f') '%'];
    elseif T.xminDouble(3) > 10e-3
      title_str = [title_str '', ...
        ' || \tau_1 = ' num2str(T.xminDouble(3)*1e3, '%.3f') ' ms @ ' num2str(T.tau1w*100, '%.1f') '%'];
    else
      title_str = [title_str '', ...
        ' || \tau_1 = ' num2str(T.xminDouble(3)*1e6, '%.3f') ' ' char(181) 's @ ' num2str(T.tau1w*100, '%.1f') '%'];
    end
    if T.xminDouble(5) > 10
      title_str = [title_str '', ...
        ', \tau_2 = ' num2str(T.xminDouble(5), '%.3f') ' s @ ' num2str(T.tau2w*100, '%.1f') '%'];
    elseif T.xminDouble(5) > 10e-3
      title_str = [title_str '', ...
        ', \tau_2 = ' num2str(T.xminDouble(5)*1e3, '%.3f') ' ms @ ' num2str(T.tau2w*100, '%.1f') '%'];
    else
      title_str = [title_str '', ...
        ', \tau_2 = ' num2str(T.xminDouble(5)*1e6, '%.3f') ' ' char(181) 's @ ' num2str(T.tau2w*100, '%.1f') '%'];
    end
    legendCell = [legendCell, 'double expo.'];
  end
end

if fitExpSettings.GaussExp
  % fit with exponential function + Gaussian
  T.xstartGaussExp(1) = T.xminSingle(1)*AmpScale;                  % offset
  T.xstartGaussExp(2) = T.xminSingle(2)/2*AmpScale;                % amp1 (exponential)
  T.xstartGaussExp(3) = T.xminSingle(3)*2;                         % tau1 (exponential)
  T.xstartGaussExp(4) = T.xminSingle(2)/2*AmpScale;                % amp2 (Gaussian)
  T.xstartGaussExp(5) = max(T.xminSingle(3)/2,2*timeCorrectedFit(1)); % tau2 (Gaussian)

  if fitExpSettings.EndOffset
    [T.xminGauss, T.errorGauss, T.exitflagGauss, T.outputGauss] = ...
      fminsearch(@(x) lsqFitGaussExp(timeCorrectedFit, dataPhaseCorrectedRealFit, x), T.xstartGaussExp, ...
      optimset('MaxFunEvals', 2000, 'TolX', 1e-4*abs(T.xstartGaussExp(2)), 'TolFun', maxerror, 'Display', 'off'));
  else
    [T.xminGauss, T.errorGauss, T.exitflagGauss, T.outputGauss] = ...
      fminsearch(@(x) lsqFitGaussExp(timeCorrectedFit, dataPhaseCorrectedRealFit, x), T.xstartGaussExp(2:5), ...
      optimset('MaxFunEvals', 1000, 'TolX', 1e-4*abs(T.xstartGaussExp(2)), 'TolFun', maxerror, 'Display', 'off'));
    T.xminGauss = [0 T.xminGauss(1:4)];
  end
  T.xminGauss([1:2,4]) = T.xminGauss([1:2,4])/AmpScale;
  T.errorGauss = sqrt(T.errorGauss/AmpScale^2 / numel(timeCorrectedFit));

  % interpolate fitted function
  [~, T.functionAmpGauss] = lsqFitGaussExp(T.functionTime, length(T.functionTime), T.xminGauss);
  [~, T.functionAmpGaussLog] = lsqFitGaussExp(T.functionTimeLog, length(T.functionTime), T.xminGauss);
  % decay constants
  T.tauG1 = T.xminGauss(3);
  T.tauG2 = T.xminGauss(5);
  % amplitude "weight"
  T.tauG1w = T.xminGauss(2)/(T.xminGauss(2)+T.xminGauss(4));
  T.tauG2w = T.xminGauss(4)/(T.xminGauss(2)+T.xminGauss(4));

  if ishghandle(fitExpSettings.hParent, 'figure') || ...
      ishghandle(fitExpSettings.hParent, 'uipanel') || fitExpSettings.hParent
    hLines1(end+1) = plot(ax(1), T.functionTime, T.functionAmpGauss, 'b-.');
    if T.xminDouble(3) > 10
      title_str = [title_str '', ...
        ' || \tau_exp = ' num2str(T.xminGauss(3), '%.3f') ' s @ ' num2str(T.tauG1w*100, '%.1f') '%'];
    elseif T.xminDouble(3) > 10e-3
      title_str = [title_str '', ...
        ' || \tau_exp = ' num2str(T.xminGauss(3)*1e3, '%.3f') ' ms @ ' num2str(T.tauG1w*100, '%.1f') '%'];
    else
      title_str = [title_str '', ...
        ' || \tau_exp = ' num2str(T.xminGauss(3)*1e6, '%.3f') ' ' char(181) 's @ ' num2str(T.tauG1w*100, '%.1f') '%'];
    end
    if T.xminDouble(5) > 10
      title_str = [title_str '', ...
        ' || \tau_gauss = ' num2str(T.xminGauss(5), '%.3f') ' s @ ' num2str(T.tauG2w*100, '%.1f') '%'];
    elseif T.xminDouble(5) > 10e-3
      title_str = [title_str '', ...
        ' || \tau_gauss = ' num2str(T.xminGauss(5)*1e3, '%.3f') ' ms @ ' num2str(T.tauG2w*100, '%.1f') '%'];
    else
      title_str = [title_str '', ...
        ' || \tau_gauss = ' num2str(T.xminGauss(5)*1e6, '%.3f') ' ' char(181) 's @ ' num2str(T.tauG2w*100, '%.1f') '%'];
    end
    legendCell = [legendCell, 'expo. + Gaussian'];
  end
end

if fitExpSettings.SingleGauss
  % fit with single Gaussian distribution function
  T.xstartGauss(1:2) = T.xstart(1:2);  % offset, amp
  T.xstartGauss(3) = 2*T.xstart(3)^2;  % tau

  % FIXME: Implement other fit types (e.g. linearized)?
  if fitExpSettings.EndOffset
    [T.xminSingleGauss, T.errorSingleGauss, T.exitflagSingleGauss, T.outputSingleGauss] = ...
      fminsearch(@(x) lsqFitGauss(timeCorrectedFit, dataPhaseCorrectedRealFit, x), T.xstartGauss, ...
                 optimset('MaxFunEvals', 1000, 'TolX', 1e-4*abs(T.xstartGauss(2)), 'TolFun', maxerror, 'Display', 'off'));
  else
    [T.xminSingleGauss, T.errorSingleGauss, T.exitflagSingleGauss, T.outputSingleGauss] = ...
      fminsearch(@(x) lsqFitGauss(timeCorrectedFit, dataPhaseCorrectedRealFit, x), T.xstartGauss(2:3), ...
                 optimset('MaxFunEvals', 1000, 'TolX', 1e-4*abs(T.xstartGauss(2)), 'TolFun', maxerror, 'Display', 'off'));
    T.xminSingleGauss = [0, T.xminSingleGauss(1:2)];
  end

  % rescale to original units
  T.xminSingleGauss(1:2) = T.xminSingleGauss(1:2) / AmpScale;
  T.tauSingleGauss = sqrt(T.xminSingleGauss(3)/2);

  % interpolate fitted function
  if isemptyfield(T, 'functionTime')
    T.functionTime = linspace(timeCorrectedFit(1)*0, timeCorrectedFit(end)*2, max(2048,numel(timeCorrectedFit)*2));
  end
  [~, T.functionAmpSingleGauss] = lsqFitGauss(T.functionTime, length(T.functionTime), T.xminSingleGauss);
  if isemptyfield(T, 'functionTime')
    T.functionTimeLog = logspace(log10(timeCorrectedFit(1)), log10(timeCorrectedFit(end)*2), max(2048,numel(timeCorrectedFit)*2));
  end
  [~, T.functionAmpSingleGaussLog] = lsqFitGauss(T.functionTimeLog, length(T.functionTimeLog), T.xminSingleGauss);
  if ishghandle(fitExpSettings.hParent, 'figure') || ...
    ishghandle(fitExpSettings.hParent, 'uipanel') || fitExpSettings.hParent
    hLines1(end+1) = plot(ax(1), T.functionTime, T.functionAmpSingleGauss, 'b:');
    if exist('title_str', 'var')
      title_str = [title_str, ' || '];
    else
      title_str = '';
    end
    if T.tauSingleGauss > 10
      title_str = [title_str, sprintf('\\tau_{Gauss} = %.3f s', T.tauSingleGauss)];
    elseif T.tauSingleGauss > 10e-3
      title_str = [title_str, sprintf('\\tau_{Gauss} = %.3f ms', T.tauSingleGauss*1e3)];
    else
      title_str = [title_str, sprintf('\\tau_{Gauss} = %.3f %cs', T.tauSingleGauss*1e6, char(181))];
    end
    legendCell = [legendCell, 'single Gauss'];
  end

end

if ishghandle(fitExpSettings.hParent, 'figure') || ...
    ishghandle(fitExpSettings.hParent, 'uipanel') || fitExpSettings.hParent
  % finalize plot decoration
  hold(ax(1), 'off')
  title(ax(1), title_str);
  legend(hLines1, legendCell);
  grid(ax(1), 'on');
  grid(ax(2), 'on');
  linkaxes(ax, 'x');
  set(zoom(ax(1)), 'Motion', 'both', 'Enable', 'on');
  drawnow expose
end


%% collect data for output
T.data = data;
T.time = time;
T.timeCorrected = timeCorrected;
T.fitExpSettings = fitExpSettings;
% FIXME: Remove these deprecated fields
T.hParent = fitExpSettings.hParent;
T.hFigure = fitExpSettings.hParent;
T.CorrectPhaseOffset = fitExpSettings.CorrectPhaseOffset;
T.EndOffset = fitExpSettings.EndOffset;
T.RingFilter = fitExpSettings.RingFilter;

end


function [dev, coeff] = linearizedFit(time, data, X, C)
%% Find solution of the exponential decay: y = A*exp(-t/tau) + C

% linearize exponential function: y(t) = A*e^(-t/tau) <=> log(y) = log(A) - t/tau

% remove offset C before linearization
dataOffset = data - C;

% dismiss non-positive values
isPositive = dataOffset > eps;

% use squared amplitudes as weights in linearized fit
b = (diag(dataOffset(isPositive))*X(isPositive,:)) ...
  \ (dataOffset(isPositive) .* log(dataOffset(isPositive)));

% calculate fitted data
y = exp(b(1))*exp(time*b(2));

% RMS deviation to data
dev = sum((dataOffset - y).^2);

coeff = [exp(b(1)), -1/b(2)];

end


function [dev, coeff] = linearizedFitFixDecay(time, data, X, C, tau)
%% Find solution of the exponential decay: y = A*exp(-t/tau) + C with fix tau

% linearize exponential function: y(t) = A*e^(-t/tau) <=> log(y) = log(A) - t/tau

% remove offset C before linearization
dataOffset = data - C;

% dismiss non-positive values
isPositive = dataOffset > eps;

% use squared amplitudes as weights in linearized fit
b = (diag(dataOffset(isPositive))*X(isPositive,:)) ...
  \ (dataOffset(isPositive) .* (log(dataOffset(isPositive))+time(isPositive)/tau));

% calculate fitted data
y = exp(b(1))*exp(-time/tau);

% RMS deviation to data
dev = sum((dataOffset - y).^2);

coeff = [exp(b(1)), tau];

end


function [dev, coeff] = integrNumExp(time, data, hasOffset, tau)
%% numerically integrate the solution
% c.f.: https://de.scribd.com/doc/14674814/Regressions-et-equations-integrales

% numerical approximation of integral
S = cumsum([0, 1/2 * (data(2:end) + data(1:end-1)) .* (time(2:end) - time(1:end-1))]);

if isempty(tau)
  if hasOffset
    % fit with offset
    % y = a + b * exp(c * time)

    % corresponding integral equation is:
    % y - (a + b * exp(c * time(1))) = - a*c*(time-time(1)) + c * integ(y, x1, x)
    % Hence the approx. linear equation to solve is:
    % y - y(1) ?= -a*c*(time-time(1)) + c * S
    % quadratic deviation to real solution:
    % sum( (-a*c*(time-time(1)) + c*S - (y-y(1))).^2 )

    % solve linear equation for (-a*c) and (c)
    off_diag = sum(S.*(time-time(1)));
    k = [sum(S.^2), off_diag; off_diag, sum((time-time(1)).^2)] \ ...
      [sum(S.*(data-data(1))); sum((time-time(1)).*(data-data(1)))];

    % determine b
    T = exp(k(1)*time);
    l = [sum(T.^2), sum(T); sum(T), numel(time)] \ [sum(T.*data); sum(data)];

    coeff = [-k(2)/k(1), l(1), -1/k(1)];  % a, b, -1/c

    % calculate fitted data
    y = coeff(1) + coeff(2) * exp(k(1)*time);

  else
    % fit without offset
    % y = b * exp(c * time)

    % corresponding integral equation is:
    % y - b * exp(c * time(1)) = c * integ(y, x1, x)
    % Hence the approx. linear equation to solve is:
    % y - y(1) ?= c * S
    % quadratic deviation to real solution:
    % sum( (c*S - (y-y(1))).^2 )

    % solve linear equation for (c)
    c = S.' \ (data - data(1)).';
    % c = (S.*(data - data(1))).' \ ((data - data(1)).^2).';
    % c = sum(S.*(data - data(1))) / sum(S.^2);

    % determine b by linear regression
    b = exp(c * time).' \ data.';

    coeff = [b, -1/c];

    % calculate fitted data
    y = b * exp(c*time);

  end
else
  if hasOffset
    % fit with offset and given c (-1/tau)
    % y = a + b * exp(c * time)

    % corresponding integral equation is:
    % y - (a + b * exp(c * time(1))) - c * integ(y, x1, x) = - a*c*(time-time(1))
    % Hence the approx. linear equation to solve is:
    % y - y(1) - c * S ?= -a*c*(time-time(1))

    % solve linear equation for (-a*c)
    mac = (time-time(1)).' \ (data - data(1) + 1/tau*S).';

    % determine b
    T = exp(-time/tau);
    l = [sum(T.^2), sum(T); sum(T), numel(time)] \ [sum(T.*data); sum(data)];

    coeff = [-mac*tau, l(1), tau];  % a, b, -1/c

    % calculate fitted data
    y = coeff(1) + coeff(2) * exp(-time/tau);

  else
    % fit without offset and given c (-1/tau)
    % y = b * exp(c * time)

    % corresponding integral equation is:
    % y - b * exp(c * time(1)) = c * integ(y, x1, x)
    % Hence the approx. linear equation to solve is:
    % y - y(1) - c * S ?= 0
    % quadratic deviation to real solution:
    % sum( (c*S - (y-y(1))).^2 )

    % This case is not really solved via numerical integration but with a linear
    % regression.

    % determine b by linear regression
    b = exp(-time/tau).' \ data.';

    coeff = [b, tau];

    % calculate fitted data
    y = b * exp(-time/tau);

  end
end

% squared deviation sum to (measured) data
dev = sum((data - y).^2);

end


function [dev, fitAmp] = lsqFitExp(time, data, A)
%% function for least squares fit of a (bi-)exponential decay with fminsearch

switch numel(A)
  case 5
    % bi-exponential with offset
    % y(t) = A(1) + A(2)*exp(-t/A(3)) + A(4)*exp(-t/A(5))
    fitAmp = A(1) + A(2)*exp(-time/A(3)) + A(4)*exp(-time/A(5));
    dev = sum(abs(data - fitAmp).^2);
    % penalty for exponential growth
    if any(A(2:end) < 0), dev = dev * 1000; end
    % penalty for improper decay time constants
    if A(3) < time(1)/2 || A(5) < time(1)/2 ...
        || A(3) > time(end)*2 || A(5) > time(end)*2
      dev = dev * 1000;
    end

  case 4
    % bi-exponential without offset
    % y(t) = A(1)*exp(-t/A(2)) + A(3)*exp(-t/A(4))
    fitAmp = A(1)*exp(-time/A(2)) + A(3)*exp(-time/A(4));
    dev = sum(abs(data - fitAmp).^2);
    % penalty for exponential growth
    if any(A < 0), dev = dev * 1000; end
    % penalty for improper decay time constants
    if A(2) < time(1)/2 || A(4) < time(1)/2 ...
        || A(2) > time(end)*2 || A(4) > time(end)*2
      dev = dev * 1000;
    end

  case 3
    % mono-exponential with offset
    % y(t) = A(1) + A(2)*exp(-t/A(3))
    fitAmp = A(1) + A(2)*exp(-time/A(3));
    dev = sum(abs( data - fitAmp ).^2);
    % penalty for exponential growth
    if any(A(2:end) < 0), dev = dev * 1000; end
    % penalty for improper decay time constant
    if A(3) < time(1)/2 || A(3) > time(end)*2
      dev = dev * 1000;
    end

  case 2
    % mono-exponential without offset
    % y(t) = A(1)*exp(-t/A(2))
    fitAmp = A(1)*exp(-time/A(2));
    dev = sum(abs(data - fitAmp).^2);
    % penalty for exponential growth
    if any(A < 0), dev = dev * 1000; end
    % penalty for improper decay time constant
    if A(2) < time(1)/2 || A(2) > time(end)*2
      dev = dev * 1000;
    end

end

end


function [dev, fitAmp] = lsqFitExpFixDecay(time, data, A, tau, hasOffset)
%% function for least squares fit of a mono-exponential decay with fixed decay constant

switch numel(A)
  case 3
    % bi-exponential with offset
    % y(t) = A(1) + A(2)*exp(-t/A(3)) + A(4)*exp(-t/A(5))
    fitAmp = A(1) + A(2)*exp(-time/tau(1)) + A(3)*exp(-time/tau(2));
    dev = sum(abs(data - fitAmp).^2);
    % penalty for exponential growth
    if any(A(2:end) < 0), dev = dev * 1000; end

  case 2
    if hasOffset
      % mono-exponential with offset
      % y(t) = A(1) + A(2)*exp(-t/A(3))
      fitAmp = A(1) + A(2)*exp(-time/tau);
      dev = sum(abs( data - fitAmp ).^2);
      % penalty for exponential growth
      if any(A(2:end) < 0), dev = dev * 1000; end
    else
      % bi-exponential without offset
      % y(t) = A(1)*exp(-t/A(2)) + A(3)*exp(-t/A(4))
      fitAmp = A(1)*exp(-time/tau(1)) + A(2)*exp(-time/tau(1));
      dev = sum(abs(data - fitAmp).^2);
      % penalty for exponential growth
      if any(A < 0), dev = dev * 1000; end
    end

  case 1
    % mono-exponential without offset
    % y(t) = A(1)*exp(-t/A(2))
    fitAmp = A(1)*exp(-time/tau);
    dev = sum(abs(data - fitAmp).^2);
    % penalty for exponential growth
    if any(A < 0), dev = dev * 1000; end

end

end


function [dev, fitAmp] = lsqFitGaussExp(time, data, A)
%% function for least squares fit of a exponential decay + Gaussian with fminsearch

switch numel(A)
  case 5
    % exponential decay + Gaussian with offset
    % y(t) = A(1) + A(2)*exp(-t/A(3)) + A(4)*exp(-t^2/(2*A(5)^2))
    fitAmp = A(1) + A(2)*exp(-time/A(3)) + A(4)*exp(-time.^2/(2*A(5)^2));
    dev = sum(abs(data - fitAmp).^2);
    % penalty for exponential growth
    if any(A(2:end) < 0), dev = dev * 1000; end
    % penalty for improper decay time constants
    if A(3) < time(1)/2 || A(5) < time(1)/2 ...
        || A(3) > time(end)*2 || A(5) > time(end)*2
      dev = dev * 1000;
    end

  case 4
    % exponential decay + Gaussian without offset
    % y(t) = A(1)*exp(-t/A(2)) + A(3)*exp(-t^2/(2*A(4)^2))
    fitAmp = A(1)*exp(-time/A(2)) + A(3)*exp(-time.^2/(2*A(4)^2));
    % weights
    numSamples = numel(fitAmp);
    nWeighted = 0;
    weights = [1000*ones(nWeighted, 1); ones(numSamples-nWeighted, 1)];
    dev = sum(abs(data - fitAmp).^2.*weights);
    % penalty for exponential growth
    if any(A < 0), dev = dev * 1000; end
    % penalty for improper decay time constants
    if A(2) < time(1)/2 || A(4) < time(1)/2 ...
        || A(2) > time(end)*2 || A(4) > time(end)*2
      dev = dev * 1000;
    end

end

end


function [dev, fitAmp] = lsqFitGauss(time, data, A)
%% function for least squares fit of a normal distribution function with fminsearch

switch numel(A)
  case 3
    % with offset
    % y(t) = A(1) + A(2)*exp(-t^2/A(3))
    fitAmp = A(1) + A(2)*exp(-time.^2/A(3));
    dev = sum(abs( data - fitAmp ).^2);
    % penalty for exponential growth
    if any(A(2:end) < 0), dev = dev * 1000; end
    % penalty for improper decay time constant
    if A(3) < time(1)/2 || A(3) > time(end)*2
      dev = dev * 1000;
    end

  case 2
    % without offset
    % y(t) = A(1)*exp(-t^2/A(2))
    fitAmp = A(1)*exp(-time.^2/A(2));
    dev = sum(abs(data - fitAmp).^2);
    % penalty for exponential growth
    if any(A < 0), dev = dev * 1000; end
    % penalty for improper decay time constant
    if A(2) < time(1)/2 || A(2) > time(end)*2
      dev = dev * 1000;
    end

end

end
