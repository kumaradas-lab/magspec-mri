function [SpectrumTime, SpectrumAmplitudeOut, Seq] = get_iLaplace1D(DataTime, DataAmplitude, Seq)
%% Calculate discrete inverse Laplace transform (ILT) of discrete data
%
%   [data, Seq] = get_iLaplace1D(data, Seq)
%
% alternative (deprecated) syntax:
%   [SpectrumTime, SpectrumAmplitude, Seq] = get_iLaplace1D(DataTime, DataAmplitude, Seq)
%
%
% INPUT:
%
%   DataTime
%       Row vector with the independent variable (time).
%   DataAmplitude
%       Row vector (or 2D matrix of row vectors) with the dependent variable
%       (amplitude).
%   data
%       Structure with the fields "DataTime" and "DataAmplitude" (see above).
%   Seq
%       Structure with the field "iLaplace1D" containing the following fields:
%     QualityFactor
%         Quality factor (Q-Factor) determining width of peaks in ILT spectrum
%         (default: 50)
%     nSpectrum
%         number of fit functions for the spectrum (default: 1000)
%     SpectrumTimeStart
%         lowest decay constant in the fit functions for the spectrum
%         (default: DataTime(1)/2)
%     SpectrumTimeEnd
%         highest decay constant in the fit functions for the spectrum
%         (default: max(DataTime)*10)
%     SpectrumTime
%         vector with decay constants in the fit functions for the spectrum
%         (default: logspace(log10(Seq.iLaplace1D.SpectrumTimeStart), ...
%                            log10(Seq.iLaplace1D.SpectrumTimeEnd), ...
%                            Seq.iLaplace1D.nSpectrum);)
%     SpectrumWeight
%         vector with weights for the model function for each decay constant
%         (default: ones(numel(iLaplace1D.SpectrumTime), 1);)
%     SpectrumTimeStartCut
%         Cut (ignore) decay constants below this value in the resulting
%         spectrum (default: Seq.iLaplace1D.SpectrumTimeStart)
%     SpectrumTimeEndCut
%         Cut (ignore) decay constants above this value in the resulting
%         spectrum (default: Seq.iLaplace1D.SpectrumTimeEnd)
%     Problem
%         Type of fit functions: 'Saturation' or 'Inversion' or 'Decay' or
%         'ProblemFunctionHandle'.
%         (Default: 'Decay')
%     ProblemFunctionHandle
%         Handle to a function that is used to generate the base matrix for the
%         inversion if Seq.iLaplace1D.Problem is set to 'ProblemFunctionHandle'.
%         The function handle must have the following signature
%             A = @(DataTime, SpectrumTime, Seq)
%         where "DataTime" is a vector with the time for each sample of the
%         acquired signal, "SpectrumTime" is a vector with the times in the
%         inverted spectrum, "Seq" is the structure as it is passed to
%         get_iLaplace1D, and "A" is the base matrix for the inversion. "A" must
%         be of the size (numel(DataTime)) x (numel(SpectrumTime)) containing
%         the amplitudes of the bas functions for the generalized discrete
%         inverse Laplace transform. (Default: [])
%     DataWeight
%         Vector of the same size as DataTime with weights for a weighted
%         (non-negative) least squares fit used in the (generalized) inverse
%         Laplace transform. (Default: 1)
%     Plot
%         Boolean (default: true)
%     hParent
%         Handle to a parent (figure or uipanel) for the plot with the Laplace
%         spectrum. If empty, figure 222 is used.
%     raiseFigure
%         If false and the figure already exists, the figure is not raised above
%         other windows and doesn't steal the focus when plotted.
%         (default: true)
%     IgnoreFirstEcho
%         Ignore n first datapoints in DataTime and DataAmplitudes (default: 0)
%     AllowDataOffset
%         Allow a global offset in DataAmplitude that is removed by the
%         inversion algorithm. (default: true)
%     QFactorPreFitThreshold
%         For smaller Q-Factors than this threshold, fit without regularization
%         first. Then re-grid the data on a looser logarithmic grid for a second
%         fit with the actual Q-Factor. (default: 50)
%
%
% OUTPUT:
%
%   SpectrumTime
%       vector with decay constants
%   SpectrumAmplitude
%       amplitudes for Laplace spectrum
%
%   If called with 2 input arguments, the following fields are added to
%   data.iLaplace. If called with 3 input arguments, they are stored in
%   Seq.iLaplace1D:
%     DataTime
%         see input
%     DataAmplitude
%         see input
%     FullScaleAmplitude
%         full scale amplitude
%     SpectrumTime
%         vector with decay constants
%     SpectrumAmplitude
%         vector with amplitude of Laplace spectrum in percent
%     SpectrumAmplitudeOffset
%         amplitude offset for inversion measurement in percent
%     FitTime
%         decay constants with tighter grid at small decay constants
%     FitAmplitudeAtDataTime
%         Fit amplitude in percent at DataTime
%     FitAmplitude
%         Fit amplitude in percent at FitTime
%     FitResidual
%         Fit residuals at DataTime
%     axHandle
%         handles to the axes in hParent with the results.
%         (In data.axHandle if called with 2 input arguments.)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% 11.12.2015

%% convert deprecated syntax to new syntax

if nargin == 2
  data = DataTime;
  Seq = DataAmplitude;
  DataTime = data.DataTime;
  DataAmplitude = real(data.DataAmplitude);
end


%% default parameters
if ~isfield(Seq, 'iLaplace1D')
  Seq.iLaplace1D = [];
end
iLaplace1D = Seq.iLaplace1D;

% Quality factor (Q-Factor) determining width of peaks in ILT spectrum
if isemptyfield(iLaplace1D, 'QualityFactor')
  iLaplace1D.QualityFactor = 50;
end
if ~isemptyfield(iLaplace1D, 'SpectrumTime')
  % Make sure properties are consistent and give precedence to SpectrumTime over
  % other related properties
  iLaplace1D.nSpectrum = numel(iLaplace1D.SpectrumTime);
  iLaplace1D.SpectrumTimeStart = min(iLaplace1D.SpectrumTime);
  iLaplace1D.SpectrumTimeEnd = max(iLaplace1D.SpectrumTime);
end
% number of fit functions for the spectrum
if isemptyfield(iLaplace1D, 'nSpectrum')
  iLaplace1D.nSpectrum = 1000;
end
% lowest decay constant in the fit functions for the spectrum
if isemptyfield(iLaplace1D, 'SpectrumTimeStart')
  iLaplace1D.SpectrumTimeStart = DataTime(1)/2;
end
% highest decay constant in the fit functions for the spectrum
if isemptyfield(iLaplace1D, 'SpectrumTimeEnd')
  iLaplace1D.SpectrumTimeEnd = max(DataTime)*10;
end
% vector with decay constants in the fit functions for the spectrum
if isemptyfield(iLaplace1D, 'SpectrumTime')
  % span of T2 values (logarithmic)
  iLaplace1D.SpectrumTime = ...
    logspace(log10(iLaplace1D.SpectrumTimeStart), ...
             log10(iLaplace1D.SpectrumTimeEnd), ...
             iLaplace1D.nSpectrum);
end
% vector with weights for the model function for each decay constant
if isemptyfield(iLaplace1D, 'SpectrumWeight')
  iLaplace1D.SpectrumWeight = ones(numel(iLaplace1D.SpectrumTime), 1);
end

iLaplace1D = set_EmptyField(iLaplace1D, 'SpectrumTimeStartCut',   iLaplace1D.SpectrumTimeStart*1);    % T2 start cut
iLaplace1D = set_EmptyField(iLaplace1D, 'SpectrumTimeEndCut',     iLaplace1D.SpectrumTimeEnd*1);      % T2 grid stop
iLaplace1D = set_EmptyField(iLaplace1D, 'Problem',                'Decay');                           % 'Saturation' or 'Inversion' or 'Decay' or 'ProblemFunctionHandle'
iLaplace1D = set_EmptyField(iLaplace1D, 'ProblemFunctionHandle',  []);                                % ProblemFunctionHandle for A and AFit
iLaplace1D = set_EmptyField(iLaplace1D, 'Plot',                   1);                                 % Plot
iLaplace1D = set_EmptyField(iLaplace1D, 'hParent',                222);                               % handle to parent
iLaplace1D = set_EmptyField(iLaplace1D, 'raiseFigure',            true);                              % raise figure each time function is called
iLaplace1D = set_EmptyField(iLaplace1D, 'IgnoreFirstEcho',        0);                                 % Ignore First x Echoes
iLaplace1D = set_EmptyField(iLaplace1D, 'DataWeight',             1);                                 % weighing factors for fit
% allow an offset in the input data
if isemptyfield(iLaplace1D, 'AllowDataOffset')
  iLaplace1D.AllowDataOffset = true;
end

iLaplace1D = set_EmptyField(iLaplace1D, 'QFactorPreFitThreshold', 50);                                % For smaller Q-Factors: Fit without Q-Factor first and re-grid the data for second fit with Q-Factor
Seq.iLaplace1D = iLaplace1D;

if Seq.iLaplace1D.IgnoreFirstEcho > 0
  DataAmplitude = DataAmplitude(Seq.iLaplace1D.IgnoreFirstEcho+1:end,:,:,:,:,:);
  DataTime = DataTime(Seq.iLaplace1D.IgnoreFirstEcho+1:end);
end
FullScaleAmplitude = max(abs(DataAmplitude), [], 1);


%% open parent figure if necessary
if Seq.iLaplace1D.Plot
  hParent = Seq.iLaplace1D.hParent;
  if ishghandle(hParent, 'figure') || (isa(hParent, 'double') && mod(hParent, 1) == 0)
    if Seq.iLaplace1D.raiseFigure || ~ishghandle(hParent, 'figure')
      hFigure = figure(hParent);
    else
      hFigure = hParent;
    end
  elseif ishghandle(hParent, 'uipanel')
    hFigure = ancestor(hParent, 'figure');
  else
    error('PD:get_iLaplace1D:invalidParent', ...
      '"Seq.iLaplace1D.hParent" must be a valid handle to a figure or uipanel.');
  end
end


%% Initialization
DataTime = reshape(DataTime, [], 1);  % reshape to column vector
if isvector(DataAmplitude)
  DataAmplitude = reshape(DataAmplitude, [], 1) ./ FullScaleAmplitude;
else
  DataAmplitude = bsxfun(@times, DataAmplitude, 1 ./ FullScaleAmplitude);
  % DataAmplitude=reshape(DataAmplitude,size(DataAmplitude,1),[]);
end
if isscalar(Seq.iLaplace1D.DataWeight)
  Seq.iLaplace1D.DataWeight = repmat(Seq.iLaplace1D.DataWeight, numel(DataTime), 1);
end

if isemptyfield(Seq.iLaplace1D, 'FitTime')
  Seq.iLaplace1D.FitTime = DataTime;
  FitTime = reshape(Seq.iLaplace1D.FitTime, [], 1);
  nEchoFitFine = 30;
  gridMultiplier = 30; % must be an even number
  FitTime = [linspace(0, Seq.iLaplace1D.FitTime(min(nEchoFitFine, length(Seq.iLaplace1D.FitTime))), ...
                      (nEchoFitFine+0.5)*gridMultiplier + 1)'; FitTime];
  FitTime(1+((1.5*gridMultiplier):gridMultiplier:((nEchoFitFine+0.5)*gridMultiplier))) = [];
  Seq.iLaplace1D.FitTime = sort(FitTime);
end
FitTime = reshape(Seq.iLaplace1D.FitTime, [], 1);


% reshape to row vector
SpectrumTime = reshape(iLaplace1D.SpectrumTime, 1, []);
iLaplace1D.SpectrumWeight = reshape(iLaplace1D.SpectrumWeight, 1, []);


%% prepare "base" of exponential functions
if strcmp(Seq.iLaplace1D.Problem, 'Inversion')
  % Inversion Problem
  A = 1 - 2*exp(-(DataTime)*(1./SpectrumTime));
  AFit = 1 - 2*exp(-(FitTime)*(1./SpectrumTime));
elseif strcmp(Seq.iLaplace1D.Problem, 'Saturation')
  % Saturation Problem
  A = 1 - 1*exp(-(DataTime)*(1./SpectrumTime));
  AFit = 1 - 1*exp(-(FitTime)*(1./SpectrumTime));
elseif strcmp(Seq.iLaplace1D.Problem, 'Decay')
  % Decay Problem
  A = exp(-(DataTime)*(1./SpectrumTime));
  AFit = exp(-(FitTime)*(1./SpectrumTime));

  % figure(5)
  % minval=-5;
  % surf(Seq.dropletRadius_spectrum*2*1e6,(sqrt(DataTime/2)./Seq.Gamma).^2.*Seq.bGradNorm,max(minval,log(A)),'LineStyle', 'none')
  % xlabel(sprintf('Diameter in %cm', char(181)));
  % ylabel(sprintf('Diffusion encoding b in s/m^2'));
  % view(2)
  % colorbar
  % grid on
  % grid minor

elseif strcmp(Seq.iLaplace1D.Problem, 'ProblemFunctionHandle')
  % ProblemFunctionHandle Problem
  if isa(Seq.iLaplace1D.ProblemFunctionHandle, 'function_handle')
    A = Seq.iLaplace1D.ProblemFunctionHandle(DataTime, SpectrumTime, Seq);
    AFit = Seq.iLaplace1D.ProblemFunctionHandle(FitTime, SpectrumTime, Seq);
  else
    error('Seq.iLaplace1D.ProblemFunctionHandle is not a function_handle');
  end
else
  error('Seq.iLaplace1D.Problem must be ''Saturation'', ''Inversion'', ''Decay'' or ''ProblemFunctionHandle''.');
end
A = bsxfun(@times, A, Seq.iLaplace1D.DataWeight);
DataAmplitude = bsxfun(@times, DataAmplitude, Seq.iLaplace1D.DataWeight);

A = bsxfun(@times, A, iLaplace1D.SpectrumWeight);
AFit = bsxfun(@times, AFit, iLaplace1D.SpectrumWeight);
if iLaplace1D.AllowDataOffset
  % separate pos and neg offset for lsqnonneg
  A = [A, ones(size(A,1),1), -ones(size(A,1),1)];
end


%% Solve linear equations
sizeDataAmplitude = [size(DataAmplitude,1),size(DataAmplitude,2),size(DataAmplitude,3),size(DataAmplitude,4),size(DataAmplitude,5)];
SpectrumAmplitudeOut=       zeros([numel(SpectrumTime)-2*double(strcmp(Seq.iLaplace1D.Problem,'Inversion')),sizeDataAmplitude(2:end)]);
SpectrumAmplitudeOffsetOut= zeros([1,sizeDataAmplitude(2:end)]);
FitAmplitudeOut=            zeros(size(FitTime));
FitAmplitudeAtDataTimeOut=  zeros(size(DataAmplitude));
FitResidualOut=             zeros(size(DataAmplitude));
FSOut=                      zeros([1,sizeDataAmplitude(2:end)]);
DataAmplitudeOut=           DataAmplitude;

maxQFactor = 1e6;
if (Seq.iLaplace1D.QualityFactor < Seq.iLaplace1D.QFactorPreFitThreshold) && ...
    (length(DataTime) > Seq.iLaplace1D.nSpectrum)
  QualityFactor = [maxQFactor Seq.iLaplace1D.QualityFactor];
else
  QualityFactor = Seq.iLaplace1D.QualityFactor;
end

G = cell(1, numel(QualityFactor));

for t=1:numel(DataAmplitudeOut)/size(DataAmplitudeOut,1)
  DataAmplitude = DataAmplitudeOut(:,t);
  DataAmplitudeFit = double(DataAmplitude);
  for iQF = 1:numel(QualityFactor)
    if QualityFactor(iQF) < maxQFactor
      if t == 1
        % include regularization in design matrix
        if iLaplace1D.AllowDataOffset
          SpectrumWeightWithOffset = [iLaplace1D.SpectrumWeight, 1, 1];
        else
          SpectrumWeightWithOffset = iLaplace1D.SpectrumWeight;
        end
        penalty = (1/QualityFactor(iQF)^4 * (sum(abs(A),1) ./ SpectrumWeightWithOffset)).^0.5;
        if iLaplace1D.AllowDataOffset
          % reduce penalty for offset
          penalty(end-1:end) = penalty(end-1:end)*1e-6;
        end
        G{iQF} = [A; diag(penalty)];
      end
      d = [DataAmplitudeFit; zeros(size(A,2),1)]; % include zeros in signal-vector
    else
      if t==1, G{iQF} = A; end     % no regularization term in design matrix
      d = DataAmplitudeFit;  % no zeros in signal-vector
    end

    % do (regularized) non-negative least squares fit
    options = struct();
    options.TolX = eps() * norm(G{iQF}, 1) * length(G{iQF});
    [SpectrumAmplitude, data.iLaplace1D.resnorm, data.iLaplace1D.FitResidual, ...
      data.iLaplace1D.exitflag, data.iLaplace1D.output, data.iLaplace1D.lambda] = ...
      lsqnonneg(G{iQF}, double(d), options);

    if iQF == 1 && length(QualityFactor) > 1
      % re-grid DataAmplitude with fitted data on reduced logarithmic timegrid
      DataTimeReduced = logspace(log10(min(DataTime)), log10(max(DataTime)), Seq.iLaplace1D.nSpectrum).';
      if strcmp(Seq.iLaplace1D.Problem, 'Inversion')
        % Inversion Problem
        A = 1 - 2*exp(-(DataTimeReduced)*(1./SpectrumTime));
      elseif strcmp(Seq.iLaplace1D.Problem, 'Saturation')
        % Saturation Problem
        A = 1 - 1*exp(-(DataTimeReduced)*(1./SpectrumTime));
      elseif strcmp(Seq.iLaplace1D.Problem, 'Decay')
        % Decay Problem
        A = exp(-(DataTimeReduced)*(1./SpectrumTime));
      elseif strcmp(Seq.iLaplace1D.Problem, 'ProblemFunctionHandle')
        % ProblemFunctionHandle Problem
        A = Seq.iLaplace1D.ProblemFunctionHandle(DataTimeReduced,SpectrumTime,Seq);
      else
        error('Seq.iLaplace1D.Problem must be ''Saturation'', ''Inversion'', ''Decay'' or ''ProblemFunctionHandle''.');
      end
      A = bsxfun(@times, A, iLaplace1D.SpectrumWeight);
      if iLaplace1D.AllowDataOffset
        % separate pos and neg offset for lsqnonneg
        A = [A, ones(size(A,1),1), -ones(size(A,1),1)];  %#ok<AGROW>
      end
      DataAmplitudeFit = A * SpectrumAmplitude;
   end
  end

  if iLaplace1D.AllowDataOffset
    % correct offset
    SpectrumAmplitudeOffset = SpectrumAmplitude(end-1) - SpectrumAmplitude(end);
    SpectrumAmplitude(end-1:end) = [];
  else
    SpectrumAmplitudeOffset = 0;
  end
  % cut spectrum
  SpectrumAmplitude((SpectrumTime<Seq.iLaplace1D.SpectrumTimeStartCut-1e-9) ...
                    | (1e-9+Seq.iLaplace1D.SpectrumTimeEndCut<SpectrumTime)) = 0;

  FS = sum(SpectrumAmplitude(:).*iLaplace1D.SpectrumWeight(:)) + SpectrumAmplitudeOffset;
  SpectrumAmplitude = SpectrumAmplitude ./ FS;
  SpectrumAmplitudeOffset = SpectrumAmplitudeOffset ./ FS;
  DataAmplitude = DataAmplitude ./ FS ./ Seq.iLaplace1D.DataWeight;
  FullScaleAmplitude(t) = FullScaleAmplitude(t) .* FS;

  %% plot results
  FitAmplitudeAtDataTime = ...
    G{1}(1:length(DataAmplitude),...
         1:end-(2*Seq.iLaplace1D.AllowDataOffset)) * SpectrumAmplitude ./ Seq.iLaplace1D.DataWeight...
    + SpectrumAmplitudeOffset;
  FitAmplitude = AFit * SpectrumAmplitude + SpectrumAmplitudeOffset;
  FitResidual = DataAmplitude - FitAmplitudeAtDataTime;
  ax = [];
  if  Seq.iLaplace1D.Plot
    delete(get(hParent, 'Children')); % clear parent
    ax(1) = subplot(4,1,1, 'Parent', hParent);
    plot(DataTime, DataAmplitude*100, '.', ...
      FitTime, FitAmplitude*100, '-', ...
      DataTime, FitResidual*100, 'k:', ...
      DataTime, 1./Seq.iLaplace1D.DataWeight*100, 'k:', ...
      'Parent', ax(1));
    title(ax(1), 'Signal, Fit and Residuals');
    ylabel(ax(1), {'Signal [% FS]'; ['FS = ', num2str(FullScaleAmplitude(t))]});
    xlabel(ax(1), 'Time [s]');
    legend(ax(1), {'Data', 'Fit', 'Residuals', '1/Weight'});
    grid(ax(1), 'on');

    ax(2) = subplot(4,1,2, 'Parent', hParent);
    plot(DataTime, FitResidual*100, 'k', 'Parent', ax(2));
    title(ax(2), 'Residuals (Signals - Fit)');
    ylabel(ax(2), 'Error [% FS]');
    xlabel(ax(2), 'Time [s]');
    legend(ax(2), {'Residuals'});
    grid(ax(2), 'on');
    linkaxes(ax(1:2), 'x')

    ax(3) = subplot(4,1,3, 'Parent', hParent);
    semilogx(SpectrumTime, SpectrumAmplitude./sum(SpectrumAmplitude(:))*100, '-', 'LineWidth', 2, 'Parent', ax(3));
    title(ax(3), 'Spectrum');
    xlabel(ax(3),'Relaxation Time [s]');
    ylabel(ax(3), 'Amplitude [% FS]');
    legend(ax(3), 'Fit');
    grid(ax(3), 'on');

    ax(4) = subplot(4,1,4, 'Parent', hParent);
    semilogx(SpectrumTime, cumsum(SpectrumAmplitude)./sum(SpectrumAmplitude(:))*100, 'LineWidth', 2, 'Parent', ax(4));
    title(ax(4), 'Cumulative Spectrum');
    ylabel(ax(4), {'Amplitude [% FS]';['FS = ',num2str(FullScaleAmplitude(t))]});
    xlabel(ax(4), 'Relaxation Time [s]');
    legend(ax(4), 'Fit');
    grid(ax(4), 'on');
    linkaxes(ax(3:4), 'x')
    ylim(ax(4), [0,100])

    zoom(hFigure, 'on');
    % call subplot again to correctly align the decorated axes
    subplot(4,1,1, ax(1));
    subplot(4,1,2, ax(2));
    subplot(4,1,3, ax(3));
    subplot(4,1,4, ax(4));
    drawnow expose;
  end

  if strcmp(Seq.iLaplace1D.Problem, 'Inversion')
    if t == numel(DataAmplitudeOut)/size(DataAmplitudeOut,1)
      SpectrumTime = SpectrumTime(2:end-1);
    end
    SpectrumAmplitudeOut(:,t) = SpectrumAmplitude(2:end-1);
  else
    SpectrumAmplitudeOut(:,t) = SpectrumAmplitude;
  end
  DataAmplitudeOut(:,t) = DataAmplitude;
  FSOut(:,t )= FS;
  FitAmplitudeOut(:,t) = FitAmplitude;
  FitAmplitudeAtDataTimeOut(:,t) = FitAmplitudeAtDataTime;
  SpectrumAmplitudeOffsetOut(:,t) = SpectrumAmplitudeOffset;
  FitResidualOut(:,t) = FitResidual;
end


%% assign results to output variables
if nargin == 2
  temp = SpectrumTime.';
  clear SpectrumTime;
  SpectrumTime = data;
  SpectrumTime.iLaplace1D.DataTime = DataTime;
  SpectrumTime.iLaplace1D.DataAmplitude = DataAmplitudeOut;
  SpectrumTime.iLaplace1D.FullScaleAmplitude = FullScaleAmplitude;
  SpectrumTime.iLaplace1D.SpectrumAmplitude = SpectrumAmplitudeOut;
  SpectrumTime.iLaplace1D.SpectrumAmplitudeOffset = SpectrumAmplitudeOffsetOut;
  SpectrumTime.iLaplace1D.SpectrumTime = temp;
  SpectrumTime.iLaplace1D.FitTime = FitTime;
  SpectrumTime.iLaplace1D.FitAmplitudeAtDataTime = FitAmplitudeAtDataTimeOut;
  SpectrumTime.iLaplace1D.FitAmplitude = FitAmplitudeOut;
  SpectrumTime.iLaplace1D.FitResidual = FitResidualOut;
  SpectrumAmplitudeOut = Seq;
  SpectrumTime.axHandle = ax;
else
  Seq.iLaplace1D.DataTime = DataTime;
  Seq.iLaplace1D.DataAmplitude = DataAmplitudeOut;
  Seq.iLaplace1D.FullScaleAmplitude = FullScaleAmplitude;
  Seq.iLaplace1D.SpectrumAmplitudeOffset = SpectrumAmplitudeOffsetOut;
  Seq.iLaplace1D.SpectrumAmplitude = SpectrumAmplitudeOut;
  Seq.iLaplace1D.SpectrumTime = SpectrumTime.';
  Seq.iLaplace1D.FitTime = FitTime;
  Seq.iLaplace1D.FitAmplitudeAtDataTime = FitAmplitudeAtDataTimeOut;
  Seq.iLaplace1D.FitAmplitude = FitAmplitudeOut;
  Seq.iLaplace1D.FitResidual = FitResidualOut;
  Seq.iLaplace1D.axHandle = ax;
end

end
