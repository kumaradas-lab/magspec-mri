function [C, T1, T2, Seq] = get_iLaplace2D(tau1, tau2, amplitude, Seq)
%% 2-d inverse Laplace transform
%
%     [data, Seq] = get_iLaplace2D(data, Seq)
%   or (alternatively):
%     [C, T1, T2, Seq] = get_iLaplace2D(tau1, tau2, amplitude, Seq);
%
% This function calculates a T1-T2-map (corresponding to tau1 and tau2) from the
% data given in amplitude.
%
%
% INPUT:
%
%   data
%       Structure with the following fields:
%
%     Tau1Time
%         Row vector with the tau1 times in seconds. (Corresponds to "tau1" in
%         the alternative syntax.)
%
%     EchoTime
%         Column vector with the tau2 times in seconds. Usually this is the CPMG
%         train. (Corresponds to "tau2" in the alternative syntax.)
%
%     MeanEchoTau1PhaseCorrected
%         Matrix with the (real valued) amplitude where the first dimension
%         corresponds to tau1 and the second dimension corresponds to tau2.
%         (Corresponds to "amplitude" in the alternative syntax.)
%
%   Seq
%       Structure with the following fields (default values are used if omitted
%       or empty):
%
%     iLaplace2D
%          Structure with the following fields (default values are used if
%          omitted or empty):
%
%       QualityFactor
%           Regularization factor that adds a penalty for high values in the
%           spectrum. A low value effectively smoothes the spectrum. High values
%           separate more peaks in the spectrum.
%           (Default: 20)
%
%       n_Lcomp
%           Number of T1 values in the T1-T2-map.
%           (Default: 101)
%
%       n_Tcomp
%           Number of T2 values in the T1-T2-map.
%           (Default: 101)
%
%       n_tau2Smooth
%           To decrease computation time, a 1d inverse Laplace transform can be
%           done on each individual tau2 train. This T2 decomposition can be
%           used to compute n_tau2Smooth logarithmically spaced amplitude values
%           which are used for the 2-d inverse Laplace transform. If 0, the data
%           is not reduced in tau2 direction (the 2-d inversion can take very
%           long). See also "LastEchoTrainCorrection".
%           (Default: 100)
%
%       T1Start
%           Shortest T1 in the spectrum in seconds (corresponding to tau1).
%           (Default: tau1(1) )
%
%       T1End
%           Longest T1 in the spectrum in seconds (corresponding to tau1).
%           (Default: tau1(end)*10 )
%
%       T2Start
%           Shortest T2 in the spectrum in seconds (corresponding to tau2).
%           (Default: tau2(1) )
%
%       T2End
%           Longest T2 in the spectrum in seconds (corresponding to tau2).
%           (Default: tau2(end)*10 )
%
%       Recovery
%           Type of recovery used for the tau1 dimension. Can be 'Saturation',
%           'Inversion', or 'Decay'.
%           (Default: 'Inversion')
%
%       IgnoreFirstEcho
%           Do not use the first N echoes in the evaluation.
%           (Default: 0)
%
%       LastEchoTrainCorrection
%           If true, subtract the last echo train from all others (effectively
%           converting into a Decay problem). For this to work properly the
%           longest tau1 should be at least approx. 2*T1 of the sample. This is
%           always set to 1 for inversion recovery if
%           Seq.iLaplace2D.n_tau2Smooth > 0.
%           (Default: fasle for 'Decay', true otherwise)
%
%       RingFilter
%           If true, calculate mean value of each two neighboring echoes to
%           reduce ringing.
%           (Default: true)
%
%       PenalizeStepHeight
%           The discrete inverse Laplace transform is ill-conditioned in
%           general. Two regularization approaches are implemented to improve
%           the condition of the inversion.
%           If true, the step height between neighboring bins in the output
%           spectrum is penalized as a means of regularization.
%           If false, the absolute spectrum amplitude is penalized.
%           (Default: false)
%
%       FullScaleAmplitude
%           Amplitude for scaling the results.
%           (Default: max(abs(amplitude(:))) )
%
%       skip_iLaplace
%            If true, the function returns before actually computing the inverse
%            Laplace transform. This can be useful because computing the inverse
%            Laplace transform can take a considerable amount of time. In some
%            situations it is still useful to call this function to get the
%            correctly formatted time-domain data (i.e., the input for the
%            inverse Laplace transform).
%            (Default: false)
%
%   If Seq.iLaplace2D.ntau2Smooth > 0, the following field is also used:
%     iLaplace1D
%           Structure as used for get_iLaplace1D. Please, see the help of that
%           function for available options.
%
%
% OUTPUT:
%
%   data
%       Same as input "data" with the following added fields in the
%       sub-structure iLaplace2D:
%
%     DataAmplitude
%         Matrix with data amplitudes as used for the decomposition.
%
%     SpectrumAmplitude
%         Amplitude of the T1-T2-map. (Corresponds to "C" in the alternative
%         syntax.)
%
%     resnorm, residual, exitflag, output, lambda
%         See help for Matlab function "lsqnonneg".
%
%     FitAmplitude
%         Amplitude of the reverse of the inverse Laplace transform, i.e., the
%         actual T1-T2-map.
%
%     T1
%         Vector with the T1 values of the map. (Corresponds to "T1" in the
%         alternative syntax.)
%
%     T2
%         Vector with the T2 values of the map. (Corresponds to "T2" in the
%         alternative syntax.)
%
%   Seq
%       Same as input "Seq" with the actually used settings.
%       Additionally, the following fields are added (among others):
%
%     iLaplace2D.FullScaleAmplitude
%         The data is normalized such that the (extrapolated, 1d-fitted) value
%         at tau2=0 equals 1 for the last echo train (for inversion recovery
%         with Seq.iLaplace2D.n_tau2Smooth>0). The data is normalized to the
%         maximum input amplitude otherwise. This value is the "actual" value
%         where the normalized amplitude is 1.
%
%
% EXAMPLE:
%
%   alternative syntax:
%     tau1 = data.Tau1Time(1,:);
%     tau2 = data.EchoTime(:,1);
%     amplitude = real(data.MeanEchoTau1PhaseCorrected);
%     [s, T1, T2] = get_iLaplace2D(tau1, tau2, amplitude, Seq);
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% default input

if nargin==2
  data = tau1;
  Seq = tau2;
  tau1 = data.Tau1Time(1,:);
  tau2 = data.EchoTime(1:end,1);
  amplitude = data.MeanEchoTau1PhaseCorrected;
end

if ~isfield(Seq, 'iLaplace2D')
  Seq.iLaplace2D = [];
end
if isemptyfield(Seq.iLaplace2D, 'QualityFactor')
  % regularization parameter
  Seq.iLaplace2D.QualityFactor = 20;
end
if isemptyfield(Seq.iLaplace2D, 'n_Lcomp')
  % Number of relevant singular values T1
  Seq.iLaplace2D.n_Lcomp = 101;
end
if isemptyfield(Seq.iLaplace2D, 'n_Tcomp')
  % Number of relevant singular values T2
  Seq.iLaplace2D.n_Tcomp = 101;
end
if isemptyfield(Seq.iLaplace2D, 'n_tau2Smooth')
  % Number of data amplitudes in tau2 direction after smooth T2
  Seq.iLaplace2D.n_tau2Smooth = 100;
end
if isemptyfield(Seq.iLaplace2D, 'T1Start')
  % T1 grid start
  Seq.iLaplace2D.T1Start = min(tau1);
end
if isemptyfield(Seq.iLaplace2D, 'T1End')
  % T1 grid stop
  Seq.iLaplace2D.T1End = max(tau1)*10;
end
if isemptyfield(Seq.iLaplace2D, 'T2Start')
  % T2 grid start
  Seq.iLaplace2D.T2Start = min(tau2);
end
if isemptyfield(Seq.iLaplace2D, 'T2End')
  % T2 grid stop
  Seq.iLaplace2D.T2End = max(tau2)*10;
end
if isemptyfield(Seq.iLaplace2D, 'Recovery')
  % 'Saturation' or 'Inversion'
  Seq.iLaplace2D.Recovery = 'Inversion';
end
if isemptyfield(Seq.iLaplace2D, 'IgnoreFirstEcho')
  % ignore first n echoes
  Seq.iLaplace2D.IgnoreFirstEcho = 0;
end
if isemptyfield(Seq.iLaplace2D, 'LastEchoTrainCorrection')
  % subtract last Echo train from all others
  Seq.iLaplace2D.LastEchoTrainCorrection = ~strcmp(Seq.iLaplace2D.Recovery, 'Decay');
end
if isemptyfield(Seq.iLaplace2D, 'RingFilter')
  % mean value of two Echoes
  Seq.iLaplace2D.RingFilter = false;
end
if isemptyfield(Seq.iLaplace2D, 'PenalizeStepHeight')
  Seq.iLaplace2D.PenalizeStepHeight = false;
end

if isemptyfield(Seq.iLaplace2D, 'skip_iLaplace')
  % return before actually calculating the inverse Laplace transform
  Seq.iLaplace2D.skip_iLaplace = false;
end

if Seq.iLaplace2D.IgnoreFirstEcho > 0
  amplitude = amplitude(Seq.iLaplace2D.IgnoreFirstEcho+1:end,:);
  tau2 = tau2(Seq.iLaplace2D.IgnoreFirstEcho+1:end);
end

if isemptyfield(Seq.iLaplace2D, 'FullScaleAmplitude')
  % amplitude for scaling results
  Seq.iLaplace2D.FullScaleAmplitude = max(abs(amplitude(:)));
end

% plot parents
if isemptyfield(Seq.iLaplace2D, 'plotAmp')
  % figure for measured amplitude
  Seq.iLaplace2D.plotAmp = 83;
end
if isemptyfield(Seq.iLaplace2D, 'plotFitAmp')
  % figure for fitted amplitude
  Seq.iLaplace2D.plotFitAmp = 85;
end
if isemptyfield(Seq.iLaplace2D, 'plotResAmp')
  % figure for residual amplitude
  Seq.iLaplace2D.plotResAmp = 86;
end
if isemptyfield(Seq.iLaplace2D, 'plotT1T2Map')
  % figure for T1-T2 map
  Seq.iLaplace2D.plotT1T2Map = 84;
end


%% Initialization
tau1 = reshape(tau1, 1, []);
tau2 = reshape(tau2, [], 1);
amplitude = amplitude ./ Seq.iLaplace2D.FullScaleAmplitude;

if Seq.iLaplace2D.RingFilter
  amplitude = convn(amplitude, [0.5;0.5], 'valid');
  tau2 = conv(tau2, [0.5;0.5], 'valid');
end
amplitude = real(amplitude);

if Seq.iLaplace2D.n_tau2Smooth && strcmp(Seq.iLaplace2D.Recovery, 'Inversion')
  % FIXME: Can we do something similar for saturation recovery measurements?
  Seq.iLaplace2D.LastEchoTrainCorrection = 1;
  Seq.iLaplace1D.FitTime = logspace(log10(tau2(1)), log10(tau2(end)), ...
    Seq.iLaplace2D.n_tau2Smooth).';
  amplitudeSmooth = zeros(numel(Seq.iLaplace1D.FitTime), numel(tau1));
  FS = ones(1, numel(tau1));
  if ~isfield(Seq, 'iLaplace1D')
    Seq.iLaplace1D = [];
  end
  if isemptyfield(Seq.iLaplace1D, 'Problem')
    % type of problem
    Seq.iLaplace1D.Problem = 'Decay';
  end
  if isemptyfield(Seq.iLaplace1D, 'Plot')
    Seq.iLaplace1D.Plot = 0;
  end
  if isemptyfield(Seq.iLaplace1D, 'nSpectrum')
    % number of output grid
    Seq.iLaplace1D.nSpectrum = Seq.iLaplace2D.n_Tcomp;
  end
  % if isemptyfield(Seq.iLaplace1D, 'SpectrumTimeStart')
  %   % T2 grid start
  %   Seq.iLaplace1D.SpectrumTimeStart = Seq.iLaplace2D.T2Start;
  % end
  % if isemptyfield(Seq.iLaplace1D, 'SpectrumTimeEnd')
  %   % T2 grid stop
  %   Seq.iLaplace1D.SpectrumTimeEnd = Seq.iLaplace2D.T2End;
  % end
  % if isemptyfield(Seq.iLaplace1D, 'SpectrumTimeStartCut')
  %   % T2 start cut
  %   Seq.iLaplace1D.SpectrumTimeStartCut = Seq.iLaplace1D.SpectrumTimeStart*1;
  % end
  % if isemptyfield(Seq.iLaplace1D, 'SpectrumTimeEndCut')
  %   % T2 grid stop
  %   Seq.iLaplace1D.SpectrumTimeEndCut = Seq.iLaplace1D.SpectrumTimeEnd;
  % end
  for t = numel(tau1):-1:1
    Seq.iLaplace1D.FullScaleAmplitude = 1;
    [~, ~, Seq] = get_iLaplace1D(tau2, amplitude(:,t), Seq);
    % revert scaling done by get_iLaplace1D
    amplitudeSmooth(:,t) = Seq.iLaplace1D.FitAmplitude * Seq.iLaplace1D.FullScaleAmplitude;
    % save scaling factors
    FS(1,t) = Seq.iLaplace1D.FullScaleAmplitude;

    if t == numel(tau1)
      % Subtract last Echo train (~>2*T1) to ensure decay (lsqnonneg)
      Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude = ...
        Seq.iLaplace1D.FitAmplitudeAtDataTime * Seq.iLaplace1D.FullScaleAmplitude;
      Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth = ...
        Seq.iLaplace1D.FitAmplitude * Seq.iLaplace1D.FullScaleAmplitude;
      amplitude = -bsxfun(@minus, amplitude, Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude);
      amplitudeSmooth(:,t) = 0;
    end
  end
  % Scale such that (extrapolated) value at tau2=0 equals 1 for the last echo
  % train.
  amplitude = -amplitudeSmooth ./ FS(end);
  Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude = ...
    Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude ./ FS(end);
  Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth = ...
    Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth ./ FS(end);
  Seq.iLaplace2D.FullScaleAmplitude = Seq.iLaplace2D.FullScaleAmplitude * FS(end);

  tau2 = Seq.iLaplace1D.FitTime.';
elseif Seq.iLaplace2D.LastEchoTrainCorrection
  iLastValid = find(~isnan(amplitude(1,:)), 1, 'last');
  Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude = amplitude(:,iLastValid);
  Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth = ...
    Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude;

  amplitude = bsxfun(@minus, amplitude, Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude);
end
tau2 = tau2(:).';

data.iLaplace2D.tau1 = tau1;
data.iLaplace2D.tau2 = tau2;

% if Seq.iLaplace2D.n_tau2Smooth
%   tau2Smooth = logspace(log10(tau2(1)), log10(tau2(end)), Seq.iLaplace2D.n_tau2Smooth);
%   IndexPolifitEnd = find(tau2>=tau2(end)/Seq.iLaplace2D.n_tau2Smooth, 1, 'first');
%   IndexPolifitSmoothEnd = find(tau2Smooth<tau2(end)/Seq.iLaplace2D.n_tau2Smooth, 1, 'last');
%
%   if Seq.iLaplace2D.LastEchoTrainCorrection
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth = zeros(1, numel(tau2Smooth));
%     [p, ~, mu] = polyfit(...
%       tau2(1:IndexPolifitEnd), ...
%       Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude(1,1:IndexPolifitEnd), ...
%       round(max(1,min(6,IndexPolifitEnd/2))));
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth(1,(1:IndexPolifitSmoothEnd)) = ...
%       polyval(p, tau2Smooth(1:IndexPolifitSmoothEnd), [], mu);
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude(1,(1:IndexPolifitEnd)) = ...
%       polyval(p, tau2(1:IndexPolifitEnd), [], mu);
%     Seq.iLaplace2D.AmplitudeScalingFactorError = polyval(p, 0, [], mu);
%
%     [p, ~, mu] = polyfit(...
%       tau2(IndexPolifitEnd+1:end), ...
%       Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude(1,IndexPolifitEnd+1:end), ...
%       round(max(1,min(16,(numel(tau2)-IndexPolifitEnd)/2))));
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth(1,IndexPolifitSmoothEnd+1:end) ...
%       = polyval(p, tau2Smooth(IndexPolifitSmoothEnd+1:end), [], mu);
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude(IndexPolifitEnd+1:end) = ...
%       polyval(p, tau2(IndexPolifitEnd+1:end), [], mu);
%
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth = ...
%       Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth ...
%       ./ Seq.iLaplace2D.AmplitudeScalingFactorError;
%   end
% end
% Seq.iLaplace2D.FullScaleAmplitude = ...
%   Seq.iLaplace2D.FullScaleAmplitude * Seq.iLaplace2D.AmplitudeScalingFactorError;
% amplitude = amplitude ./ Seq.iLaplace2D.AmplitudeScalingFactorError;
% if Seq.iLaplace2D.LastEchoTrainCorrection
%   Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude = ...
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude ...
%     ./ Seq.iLaplace2D.AmplitudeScalingFactorError;
%   amplitude = bsxfun(@minus, amplitude, Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude);
%   if Seq.iLaplace2D.n_tau2Smooth
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude = ...
%       Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth;
%   end
% end

% if Seq.iLaplace2D.n_tau2Smooth
%   amplitudeSmooth = zeros(numel(tau1), numel(tau2Smooth));
%   % amplitudeSmoothError = zeros(numel(tau1), numel(tau2));
%   for t=1:numel(tau1)
%     [p, ~ ,mu] = polyfit(...
%       tau2(1:IndexPolifitEnd), ...
%       amplitude(t,1:IndexPolifitEnd), ...
%       round(max(1,min(6,IndexPolifitEnd/2))));
%     amplitudeSmooth(t,(1:IndexPolifitSmoothEnd)) = ...
%       polyval(p, tau2Smooth(1:IndexPolifitSmoothEnd), [], mu);
%     % amplitudeSmoothError(t,1:IndexPolifitEnd) = ...
%     %   polyval(p, tau2(1:IndexPolifitEnd), [], mu) - amplitude(t,1:IndexPolifitEnd);
%     [p, ~, mu] = polyfit(...
%       tau2(IndexPolifitEnd+1:end), ...
%       amplitude(t,IndexPolifitEnd+1:end), ...
%       round(max(1,min(16,(numel(tau2)-IndexPolifitEnd)/2))));
%     amplitudeSmooth(t,IndexPolifitSmoothEnd+1:end) = ...
%       polyval(p, tau2Smooth(IndexPolifitSmoothEnd+1:end), [], mu);
%     % amplitudeSmoothError(t,IndexPolifitEnd+1:end) = ...
%     %   polyval(p, tau2(IndexPolifitEnd+1:end), [], mu) - amplitude(t,IndexPolifitEnd+1:end);
%     % figure(1)
%     % subplot(2,1,1)
%     % plot(tau2, amplitude(t,:),'-',...
%     %      tau2Smooth,amplitudeSmooth(t,:),'-x')
%     % subplot(2,1,2)
%     % plot(  tau2, amplitudeSmoothError(t,:),'-d')
%     % pause
%   end
%
%   amplitude = amplitudeSmooth;
%   tau2 = tau2Smooth;
% end


%% generate matrix with objective functions
% span of T1 values (logarithmic)
T1 = logspace(log10(Seq.iLaplace2D.T1Start), log10(Seq.iLaplace2D.T1End), Seq.iLaplace2D.n_Lcomp);
% span of T2 values (logarithmic)
T2 = logspace(log10(Seq.iLaplace2D.T2Start), log10(Seq.iLaplace2D.T2End), Seq.iLaplace2D.n_Tcomp);

% T1 -Measurement ---------------------------------------------------------
if strcmp(Seq.iLaplace2D.Recovery, 'Inversion')
  % Inversion Recovery
  if Seq.iLaplace2D.LastEchoTrainCorrection
    % Effectively, a decay problem
    L = (-2)*exp(-tau1.' * (1./T1));
  else
    L = 1-2*exp(-tau1.' * (1./T1));
  end
elseif strcmp(Seq.iLaplace2D.Recovery, 'Saturation')
  % Saturation Recovery
  if Seq.iLaplace2D.LastEchoTrainCorrection
    % Effectively, a decay problem
    L = (-1)*exp(-tau1.' * (1./T1));
  else
    L = 1-1*exp(-tau1.' * (1./T1));
  end
elseif strcmp(Seq.iLaplace2D.Recovery, 'Decay')
  % Decay
  L = exp(-tau1.' * (1./T1));
else
  error('PD:get_iLaplace2D:UnknownRecovery', ...
    'Seq.iLaplace2D.Recovery must be set to ''Inversion'', ''Saturation'', or ''Decay''.');
end

T = exp(-tau2.' * (1./T2));  % T2 (decay) measurement
T = permute(T, [1,3,2]);
data.iLaplace2D.DataAmplitude = reshape(real(amplitude), [], 1);

if Seq.iLaplace2D.skip_iLaplace
  if nargin == 2
    data.iLaplace2D.DataAmplitude = reshape(data.iLaplace2D.DataAmplitude, size(T,1), size(L,1));

    if Seq.iLaplace2D.LastEchoTrainCorrection
      data.iLaplace2D.DataAmplitude = bsxfun(@plus, ...
        data.iLaplace2D.DataAmplitude, ...
        Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth);
    end
    C = data;
    T1 = Seq;
  end
  return;
end

% dimensions: tau2 x tau1 x T2 x T1
TL = zeros(size(T,1), size(L,1), size(T,3), size(L,2));
for t = 1:size(L, 2)
  TL(:,:,:,t) = bsxfun(@times, L(:,t).', T);
  % TLR(:,:,:,t) = bsxfun(@times, LR(:,t).', T);
end
A = reshape(TL, size(T,1)*size(L,1), size(T,3)*size(L,2));
% AR = reshape(TLR, size(T,1)*size(L,1), size(T,3)*size(L,2));


%% regularize LSE
if Seq.iLaplace2D.QualityFactor<1e6
  if Seq.iLaplace2D.PenalizeStepHeight
    % penalize large steps between neighboring bins in output spectrum

    % neighbors in T2 direction
    % The more different the model function for neighboring bins are, the lower
    % is the penalty.
    % Additionally, rescale with amplitude of model functions to have similar
    % effect of QualityFactor to other branch.
    diffT =  1 ...
      ./ sqrt(mean(mean(diff(abs(TL), [], 3).^2, 1), 2)) ...
      .* sqrt(mean(mean(((TL(:,:,1:end-1,:)+TL(:,:,2:end,:))/2).^2, 1), 2));
    neighborPenaltyT = (1/Seq.iLaplace2D.QualityFactor^4 * diffT).^0.5;
    % Allow higher steps if neighboring bins in spectrum are further apart.
    T2Step = diff(T2);
    meanT2Step = mean(T2Step);
    neighborPenaltyT = cat(3, ...
      bsxfun(@rdivide,  neighborPenaltyT, ...
             sqrt(reshape(T2Step, [1, 1, numel(T2Step), 1]) / meanT2Step)) ...
      / sqrt(numel(T2Step) / (numel(tau2))), ...
      zeros(1, 1, 1, size(TL, 4)));
    neighborPenaltyArrayT = diag(neighborPenaltyT(:)) ...
      - circshift(diag(neighborPenaltyT(:)), 1, 2);
    neighborPenaltyArrayT(numel(T2):numel(T2):end,:) = [];

    % neighbors in T1 direction
    diffL =  1 ...
      ./ sqrt(mean(mean(diff(abs(TL), [], 4).^2, 1), 2)) ...
      .* sqrt(mean(mean(((TL(:,:,:,1:end-1)+TL(:,:,:,2:end))/2).^2, 1), 2));
    neighborPenaltyL = (1/Seq.iLaplace2D.QualityFactor^4 * diffL).^0.5;
    % Allow higher steps if neighboring bins in spectrum are further apart.
    T1Step = diff(T1);
    meanT1Step = mean(T1Step);
    neighborPenaltyL = cat(4, ...
      bsxfun(@rdivide,  neighborPenaltyL, ...
             sqrt(reshape(T1Step, [1, 1, 1, numel(T1Step)]) / meanT1Step)) ...
      / sqrt(numel(T1Step) / (numel(tau1))), ...
      zeros(1, 1, size(TL, 3), 1));
    neighborPenaltyArrayL = diag(neighborPenaltyL(:)) ...
      - circshift(diag(neighborPenaltyL(:)), size(TL, 3), 2);
    neighborPenaltyArrayL(end-size(TL, 3)+1:end,:) = [];

    % create objective matrix and data with penalties
    % include regularization in design matrix
    G = [A; neighborPenaltyArrayT; neighborPenaltyArrayL];
    % include zeros in signal-vector
    d = [data.iLaplace2D.DataAmplitude; ...
      zeros(size(neighborPenaltyArrayT,1)+size(neighborPenaltyArrayL,1),1)];

    clear('neighborPenaltyArrayT', 'neighborPenaltyArrayL');

  else
    % penalize large amplitudes of bins in output spectrum
    % include regularization in design matrix
    G = [A; diag((1/Seq.iLaplace2D.QualityFactor^4 * sum(abs(A), 1)).^0.5)];
    % include zeros in signal-vector
    d = [data.iLaplace2D.DataAmplitude; zeros(size(A,2), 1)];

  end

else
  % no regularization in design matrix for very high QualityFactor
  G = A;
  d = data.iLaplace2D.DataAmplitude;

end


%% inversion (with non-negative factors)
options.TolX = 1*eps*norm(G,1)*length(G);
% do regularized non-negative least squares
[data.iLaplace2D.SpectrumAmplitude, data.iLaplace2D.resnorm, ...
  data.iLaplace2D.FitResidual, data.iLaplace2D.exitflag, ...
  data.iLaplace2D.output, data.iLaplace2D.lambda] = lsqnonneg(G, double(d), options);


%% process results
% data.iLaplace2D.FitAmplitude = AR*data.iLaplace2D.SpectrumAmplitude;
data.iLaplace2D.FitAmplitude = ...
  G(1:length(data.iLaplace2D.DataAmplitude),:) * data.iLaplace2D.SpectrumAmplitude;
% data.iLaplace2D.FitAmplitude = data.iLaplace2D.FitAmplitude(1:length(data.iLaplace2D.DataAmplitude));

data.iLaplace2D.FitAmplitude = reshape(data.iLaplace2D.FitAmplitude, size(T,1), size(L,1));
data.iLaplace2D.DataAmplitude = reshape(data.iLaplace2D.DataAmplitude, size(T,1), size(L,1));
data.iLaplace2D.SpectrumAmplitude = reshape(data.iLaplace2D.SpectrumAmplitude, size(T,3), size(L,2));

if Seq.iLaplace2D.LastEchoTrainCorrection
  data.iLaplace2D.FitAmplitude = bsxfun(@plus, ...
    data.iLaplace2D.FitAmplitude, ...
    Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth);
  data.iLaplace2D.DataAmplitude = bsxfun(@plus, ...
    data.iLaplace2D.DataAmplitude, ...
    Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth);
end


data.iLaplace2D.T1 = T1;
data.iLaplace2D.T2 = T2;
Seq = plot_iLaplace2D(data, Seq);

% residual
% Seq.iLaplace2D.res = norm(d - G*s(:),2)^2;

if nargin == 2
  C = data;
  T1 = Seq;
else
  C = data.iLaplace2D.SpectrumAmplitude;
end

end
