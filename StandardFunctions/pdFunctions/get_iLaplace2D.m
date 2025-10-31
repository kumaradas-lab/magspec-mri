function [C, T1, T2, Seq] = get_iLaplace2D(tau1, tau2, amplitude, Seq)
%% 2d inverse Laplace transform
%
%     [data, Seq] = get_iLaplace2D(data, Seq)
%   or (alternatively):
%     [C, T1, T2, Seq] = get_iLaplace2D(tau1, tau2, amplitude, Seq);
%
% This function calculates a T1-T2-map (corresponding to tau1 and tau2) from the
% data given in amplitude.
%
% INPUT:
%   data      Structure with the following fields:
%     Tau1Time  Row vector with the tau1 times in seconds. (Corresponds to
%               "tau1" in the alternative syntax.)
%     EchoTime  Column vector with the tau2 times in seconds. Usually this is
%               the CPMG train. (Corresponds to "tau2" in the alternative
%               syntax.)
%     MeanEchoTau1PhaseCorrected
%               Matrix with the (real valued) amplitude where the first
%               dimension corresponds to tau1 and the second dimension
%               corresponds to tau2. (Corresponds to "amplitude" in the
%               alternative syntax.)
%
%   Seq       Structure with the following fields (default values are used if
%             omitted or empty):
%     iLaplace2D
%               Structure with the following fields (default values are used if
%               omitted or empty):
%       QualityFactor
%                 Regularization factor that adds a penalty for high values in
%                 the spectrum. A low value effectively smoothes the spectrum.
%                 High values separate more peaks in the spectrum. (Default: 20)
%       n_Lcomp   Number of T1 values in the T1-T2-map (default: 101).
%       n_Tcomp   Number of T2 values in the T1-T2-map (default: 101).
%       n_tau2Smooth
%                 To decrease computation time, a 1d inverse Laplace transform
%                 can be done on each individual tau2 train. This T2
%                 decomposition can be used to compute n_tau2Smooth
%                 logarithmically spaced amplitude values which are used for the
%                 2d inverse Laplace transform. If 0, the data is not reduced in
%                 tau2 direction (the 2d inversion can take very long). See also
%                 "LastEchoTrainCorrection". (Default: 100)
%       T1Start   Shortest T1 in the spectrum in seconds (corresponding to
%                 tau1). (Default: tau1(1) )
%       T1End     Longest T1 in the spectrum in seconds (corresponding to tau1).
%                 (Default: tau1(end)*10 )
%       T2Start   Shortest T2 in the spectrum in seconds (corresponding to
%                 tau2). (Default: tau2(1) )
%       T2End     Longest T2 in the spectrum in seconds (corresponding to tau2).
%                 (Default: tau2(end)*10 )
%       Recovery  Type of recovery used for the tau1 dimension. Can be
%                 'Saturation', 'Inversion', or 'Decay'. (Default: 'Inversion')
%       IgnoreFirstEcho
%                 Don't use the first N echoes in the evaluation. (default: 0)
%       LastEchoTrainCorrection
%                 If true, subtract the last Echo train from all others
%                 (effectively converting into a Decay problem). For this to
%                 work properly the longest tau1 should be approx. 2*T1 of the
%                 sample. This is always set to 1 for inversion recovery if
%                 Seq.iLaplace2D.n_tau2Smooth > 0.
%                 (Default: 0 for 'Decay', 1 otherwise)
%       RingFilter
%                 If true, calculate mean value of each two neighboring Echoes
%                 to reduce ringing (default: 1).
%       FullScaleAmplitude
%                 Amplitude for scaling the results.
%                 (Default: max(abs(amplitude(:))) )
%     If Seq.iLaplace2D.ntau2Smooth > 0, the following field is also used:
%     iLaplace1D
%               Structure as used for get_iLaplace1D. Please, see the help of
%               that function for available options.
%
% OUTPUT:
%   data      Same as input "data" with the following added fields in the
%             sub-structure iLaplace2D:
%     DataAmplitude
%               Matrix with data amplitudes as used for the decomposition.
%     SpectrumAmplitude
%               Amplitude of the T1-T2-map. (Corresponds to "C" in the
%               alternative syntax.)
%     resnorm, residual, exitflag, output, lambda
%               See help for Matlab function "lsqnonneg".
%     FitAmplitude
%               Amplitude of the reverse of the inverse Laplace transform, i.e.
%               the actual T1-T2-map.
%     T1        Vector with the T1 values of the map. (Corresponds to "T1" in
%               the alternative syntax.)
%     T2        Vector with the T2 values of the map. (Corresponds to "T2" in
%               the alternative syntax.)
%
%   Seq       Same as input "Seq" with the actually used settings.
%
%
% EXAMPLE:
%   alternative syntax:
%     tau1 = data.Tau1Time(1,:);
%     tau2 = data.EchoTime(:,1);
%     amplitude = real(data.MeanEchoTau1PhaseCorrected);
%     [s, T1, T2] = get_iLaplace2D(tau1, tau2, amplitude, Seq);
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% 04.11.2015

%%

if nargin==2
  data = tau1;
  Seq = tau2;
  tau1 = data.Tau1Time(1,:);
  tau2 = data.EchoTime(1:end,1);
  amplitude = data.MeanEchoTau1PhaseCorrected;
end

if ~isfield(Seq, 'iLaplace2D'),                   Seq.iLaplace2D = [];                end
if isemptyfield(Seq.iLaplace2D, 'QualityFactor'), Seq.iLaplace2D.QualityFactor = 20;  end % regularization
if isemptyfield(Seq.iLaplace2D, 'n_Lcomp'),       Seq.iLaplace2D.n_Lcomp = 101;       end % Number of relevant singular values T1
if isemptyfield(Seq.iLaplace2D, 'n_Tcomp'),       Seq.iLaplace2D.n_Tcomp = 101;       end % Number of relevant singular values T2
if isemptyfield(Seq.iLaplace2D, 'n_tau2Smooth'),  Seq.iLaplace2D.n_tau2Smooth = 100;  end % Number of data amplitudes in tau2 direction after smooth T2
if isemptyfield(Seq.iLaplace2D, 'T1Start'),       Seq.iLaplace2D.T1Start = tau1(1);   end % T1 grid start
if isemptyfield(Seq.iLaplace2D, 'T1End'),         Seq.iLaplace2D.T1End = tau1(end)*10;end % T1 grid stop
if isemptyfield(Seq.iLaplace2D, 'T2Start'),       Seq.iLaplace2D.T2Start = tau2(1);   end % T2 grid start
if isemptyfield(Seq.iLaplace2D, 'T2End'),         Seq.iLaplace2D.T2End = tau2(end)*10;end % T2 grid stop
if isemptyfield(Seq.iLaplace2D, 'Recovery'),      Seq.iLaplace2D.Recovery = 'Inversion'; end % 'Saturation' or 'Inversion'
if isemptyfield(Seq.iLaplace2D, 'IgnoreFirstEcho'), Seq.iLaplace2D.IgnoreFirstEcho = 0; end % ignore first n Echoes
if isemptyfield(Seq.iLaplace2D, 'LastEchoTrainCorrection')
  % subtract last Echo train from all others
  Seq.iLaplace2D.LastEchoTrainCorrection = ~strcmp(Seq.iLaplace2D.Recovery, 'Decay');
end
if isemptyfield(Seq.iLaplace2D, 'RingFilter'),    Seq.iLaplace2D.RingFilter = 1;      end % mean value of two Echoes

if Seq.iLaplace2D.IgnoreFirstEcho > 0
  amplitude = amplitude(Seq.iLaplace2D.IgnoreFirstEcho+1:end,:);
  tau2 = tau2(Seq.iLaplace2D.IgnoreFirstEcho+1:end);
end

if isemptyfield(Seq.iLaplace2D, 'FullScaleAmplitude'), Seq.iLaplace2D.FullScaleAmplitude = max(abs(amplitude(:))); end % amplitude for scaling results

% Initialization
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
  Seq.iLaplace1D.FitTime = logspace(log10(tau2(1)), log10(tau2(end)), Seq.iLaplace2D.n_tau2Smooth).';
  amplitudeSmooth = zeros(numel(Seq.iLaplace1D.FitTime), numel(tau1));
  FS = ones(1, numel(tau1));
  if ~isfield(Seq, 'iLaplace1D'),             Seq.iLaplace1D = [];     end
  if isemptyfield(Seq.iLaplace1D, 'Problem'), Seq.iLaplace1D.Problem = 'Decay'; end % Type of problem
  if isemptyfield(Seq.iLaplace1D, 'Plot'),    Seq.iLaplace1D.Plot = 0; end  % Plot
  if isemptyfield(Seq.iLaplace1D, 'nSpectrum'), Seq.iLaplace1D.nSpectrum = Seq.iLaplace2D.n_Tcomp; end % number of output grid
  % if isemptyfield(Seq.iLaplace1D, 'SpectrumTimeStart'), Seq.iLaplace1D.SpectrumTimeStart = Seq.iLaplace2D.T2Start; end      % T2 grid start
  % if isemptyfield(Seq.iLaplace1D, 'SpectrumTimeEnd'),   Seq.iLaplace1D.SpectrumTimeEnd = Seq.iLaplace2D.T2End;     end      % T2 grid stop
  % if isemptyfield(Seq.iLaplace1D, 'SpectrumTimeStartCut'), Seq.iLaplace1D.SpectrumTimeStartCut = Seq.iLaplace1D.SpectrumTimeStart*1; end % T2 start cut
  % if isemptyfield(Seq.iLaplace1D, 'SpectrumTimeEndCut'),   Seq.iLaplace1D.SpectrumTimeEndCut = Seq.iLaplace1D.SpectrumTimeEnd; end  % T2 grid stop
  for t = numel(tau1):-1:1
    Seq.iLaplace1D.FullScaleAmplitude = 1;
    [~, ~, Seq] = get_iLaplace1D(tau2, amplitude(:,t), Seq);
    amplitudeSmooth(:,t) = Seq.iLaplace1D.FitAmplitude * Seq.iLaplace1D.FullScaleAmplitude;
    FS(1,t) = Seq.iLaplace1D.FullScaleAmplitude;

    if t==numel(tau1)
      % Subtract last Echo train (~>2*T1) to ensure decay (lsqnonneg)
      Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude = Seq.iLaplace1D.FitAmplitudeAtDataTime*Seq.iLaplace1D.FullScaleAmplitude;
      Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth = Seq.iLaplace1D.FitAmplitude*Seq.iLaplace1D.FullScaleAmplitude;
      amplitude = -bsxfun(@minus, amplitude, Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude);
      Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude = Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth;
      amplitudeSmooth(:,t) = 0;
    end
  end
  amplitude = -amplitudeSmooth ./ FS(end);
  Seq.iLaplace2D.FullScaleAmplitude = Seq.iLaplace2D.FullScaleAmplitude .* FS(end);
  tau2 = Seq.iLaplace1D.FitTime.';
elseif Seq.iLaplace2D.LastEchoTrainCorrection
  Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude = amplitude(:,end);
  amplitude = bsxfun(@minus, amplitude, Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude);
end
tau2 = tau2(:).';

% if Seq.iLaplace2D.n_tau2Smooth
%   tau2Smooth=logspace(log10(tau2(1)),log10(tau2(end)),Seq.iLaplace2D.n_tau2Smooth);
%   IndexPolifitEnd=find(tau2>=tau2(end)/Seq.iLaplace2D.n_tau2Smooth,1, 'first');
%   IndexPolifitSmoothEnd=find(tau2Smooth<tau2(end)/Seq.iLaplace2D.n_tau2Smooth,1, 'last');
%
%   if Seq.iLaplace2D.LastEchoTrainCorrection
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth=zeros(1,numel(tau2Smooth));
%     [p,~,mu] = polyfit(tau2(1:IndexPolifitEnd), Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude(1,1:IndexPolifitEnd), round(max(1,min(6,IndexPolifitEnd/2))));
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth(1,(1:IndexPolifitSmoothEnd)) = polyval(p,tau2Smooth(1:IndexPolifitSmoothEnd),[],mu);
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude(1,(1:IndexPolifitEnd)) = polyval(p,tau2(1:IndexPolifitEnd),[],mu);
%     Seq.iLaplace2D.AmplitudeScalingFactorError=polyval(p,0,[],mu);
%
%     [p,~,mu] = polyfit(tau2(IndexPolifitEnd+1:end), Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude(1,IndexPolifitEnd+1:end), round(max(1,min(16,(numel(tau2)-IndexPolifitEnd)/2))));
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth(1,IndexPolifitSmoothEnd+1:end) = polyval(p,tau2Smooth(IndexPolifitSmoothEnd+1:end),[],mu);
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude(IndexPolifitEnd+1:end)=polyval(p,tau2(IndexPolifitEnd+1:end),[],mu);
%
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth=Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth./Seq.iLaplace2D.AmplitudeScalingFactorError;
%   end
% end
% Seq.iLaplace2D.FullScaleAmplitude=Seq.iLaplace2D.FullScaleAmplitude*Seq.iLaplace2D.AmplitudeScalingFactorError;
% amplitude=amplitude./Seq.iLaplace2D.AmplitudeScalingFactorError;
% if Seq.iLaplace2D.LastEchoTrainCorrection
%   Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude=Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude./Seq.iLaplace2D.AmplitudeScalingFactorError;
%   amplitude=bsxfun(@minus,amplitude,Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude);
%   if Seq.iLaplace2D.n_tau2Smooth
%     Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude=Seq.iLaplace2D.LastEchoTrainCorrectionAmplitudeSmooth;
%   end
% end

% if Seq.iLaplace2D.n_tau2Smooth

%   amplitudeSmooth=zeros(numel(tau1),numel(tau2Smooth));
%   % amplitudeSmoothError=zeros(numel(tau1),numel(tau2));
%   for t=1:numel(tau1)
%     [p,~,mu] = polyfit(tau2(1:IndexPolifitEnd), amplitude(t,1:IndexPolifitEnd), round(max(1,min(6,IndexPolifitEnd/2))));
%     amplitudeSmooth(t,(1:IndexPolifitSmoothEnd)) = polyval(p,tau2Smooth(1:IndexPolifitSmoothEnd),[],mu);
%     % amplitudeSmoothError(t,1:IndexPolifitEnd) = polyval(p,tau2(1:IndexPolifitEnd),[],mu)-amplitude(t,1:IndexPolifitEnd);
%     [p,~,mu] = polyfit(tau2(IndexPolifitEnd+1:end), amplitude(t,IndexPolifitEnd+1:end), round(max(1,min(16,(numel(tau2)-IndexPolifitEnd)/2))));
%     amplitudeSmooth(t,IndexPolifitSmoothEnd+1:end) = polyval(p,tau2Smooth(IndexPolifitSmoothEnd+1:end),[],mu);
%     % amplitudeSmoothError(t,IndexPolifitEnd+1:end) = polyval(p,tau2(IndexPolifitEnd+1:end),[],mu)-amplitude(t,IndexPolifitEnd+1:end);
%     % figure(1)
%     % subplot(2,1,1)
%     % plot(tau2, amplitude(t,:),'-',...
%     %      tau2Smooth,amplitudeSmooth(t,:),'-x')
%     % subplot(2,1,2)
%     % plot(  tau2, amplitudeSmoothError(t,:),'-d')
%     % pause
%   end
%
%   amplitude=amplitudeSmooth;
%   tau2=tau2Smooth;
% end

T1 = logspace(log10(Seq.iLaplace2D.T1Start), log10(Seq.iLaplace2D.T1End), Seq.iLaplace2D.n_Lcomp); % span of T1 values (logarithmic)
T2 = logspace(log10(Seq.iLaplace2D.T2Start), log10(Seq.iLaplace2D.T2End), Seq.iLaplace2D.n_Tcomp); % span of T2 values (logarithmic)

% T1 -Measurement ---------------------------------------------------------
if strcmp(Seq.iLaplace2D.Recovery, 'Inversion')
  % Inversion Recovery
  if Seq.iLaplace2D.LastEchoTrainCorrection
    % Effectively, a decay problem
    L = (-2)*exp(-tau1.' * (1./T1));
    % LR = 1-2*exp(-tau1.' * (1./T1));
  else
    L = 1-2*exp(-tau1.' * (1./T1));
    % LR = 1-2*exp(-tau1.' * (1./T1));
  end
elseif strcmp(Seq.iLaplace2D.Recovery, 'Saturation')
  % Saturation Recovery
  L = 1-1*exp(-tau1.' * (1./T1));
  % LR = 1-1*exp(-tau1.' * (1./T1));
elseif strcmp(Seq.iLaplace2D.Recovery, 'Decay')
  % Decay
  L = exp(-tau1.' * (1./T1));
else
  error('PD:get_iLaplace2D:UnknownRecovery', ...
    'Seq.iLaplace2D.Recovery must be set to ''Inversion'', ''Saturation'', or ''Decay''.');
end

T = exp(-tau2.' * (1./T2));  % T2 (decay) measurement
T = permute(T, [1,3,2]);
TL = zeros(size(T,1), size(L,1), size(T,3), size(L,2));
for t = 1:size(L, 2)
  TL(:,:,:,t) = bsxfun(@times, L(:,t).', T);
  % TLR(:,:,:,t) = bsxfun(@times, LR(:,t).', T);
end
A = reshape(TL, size(T,1)*size(L,1), size(T,3)*size(L,2));
% AR = reshape(TLR, size(T,1)*size(L,1), size(T,3)*size(L,2));
data.iLaplace2D.DataAmplitude = reshape(real(amplitude), [], 1);

if Seq.iLaplace2D.QualityFactor<1e6
  G = [A; repmat((1/Seq.iLaplace2D.QualityFactor^4*sum(abs(A),1)).^0.5,size(A,2),1).*eye(size(A,2))];  % include Regularization in design matrix
  d = [data.iLaplace2D.DataAmplitude; zeros(size(A,2),1)];  % include zeros in signal-vector
else
  G = A;  % no regularization term in design matrix
  d = data.iLaplace2D.DataAmplitude;  % no zeros in signal-vector
end

clear('options');
options.TolX = 1*eps*norm(G,1)*length(G);
[data.iLaplace2D.SpectrumAmplitude, data.iLaplace2D.resnorm, ...
  data.iLaplace2D.FitResidual, data.iLaplace2D.exitflag, ...
  data.iLaplace2D.output, data.iLaplace2D.lambda] = lsqnonneg(G, d,options);  % do regularized non-negative least squares

% data.iLaplace2D.FitAmplitude = AR*data.iLaplace2D.SpectrumAmplitude;
data.iLaplace2D.FitAmplitude = G(1:length(data.iLaplace2D.DataAmplitude),:)*data.iLaplace2D.SpectrumAmplitude;
% data.iLaplace2D.FitAmplitude = data.iLaplace2D.FitAmplitude(1:length(data.iLaplace2D.DataAmplitude));

data.iLaplace2D.FitAmplitude=reshape(data.iLaplace2D.FitAmplitude,size(T,1),size(L,1));
data.iLaplace2D.DataAmplitude=reshape(data.iLaplace2D.DataAmplitude,size(T,1),size(L,1));
data.iLaplace2D.SpectrumAmplitude=reshape(data.iLaplace2D.SpectrumAmplitude,size(T,3),size(L,2));

if Seq.iLaplace2D.LastEchoTrainCorrection
  data.iLaplace2D.FitAmplitude=bsxfun(@plus,data.iLaplace2D.FitAmplitude,Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude);
  data.iLaplace2D.DataAmplitude=bsxfun(@plus,data.iLaplace2D.DataAmplitude,Seq.iLaplace2D.LastEchoTrainCorrectionAmplitude);
end


% data.iLaplace2D.SpectrumAmplitude(or(SpectrumTime<Seq.iLaplace1D.SpectrumTimeStartCut,Seq.iLaplace1D.SpectrumTimeEndCut<SpectrumTime))=0;
% data.iLaplace2D.SpectrumAmplitude(SpectrumTime<diff(DataTime(1:2))/2)=0;
% figure; imagesc(data.iLaplace2D.SpectrumAmplitude)
% figure; imagesc(data.iLaplace2D.DataAmplitude)
% figure; imagesc(data.iLaplace2D.DataAmplitude-data.iLaplace2D.FitAmplitude)
% figure; imagesc(data.iLaplace2D.FitAmplitude)
FS = sum(data.iLaplace2D.SpectrumAmplitude(:));
% data.iLaplace2D.SpectrumAmplitude=data.iLaplace2D.SpectrumAmplitude./FS;
% data.iLaplace2D.DataAmplitude=data.iLaplace2D.DataAmplitude./FS;
% Seq.iLaplace1D.FullScaleAmplitude=Seq.iLaplace1D.FullScaleAmplitude.*FS;



%% plot "measured" signals
%
hf = figure(83);
clf(hf);
hax = axes(hf);
surf(hax, tau2*1000, tau1*1000, data.iLaplace2D.DataAmplitude.', 'LineStyle', 'none');
colormap(hf, 'jet');
title(hax, 'Measured Signal Amplitude');
ylabel(hax, 'Inversion Time [ms]');
xlabel(hax, 'Echo Time [ms]');


%% plot final result
%
hf = figure(85);
clf(hf);
hax = axes(hf);
surf(hax, tau2, tau1, data.iLaplace2D.FitAmplitude.', 'LineStyle', 'none');
colormap(hf, 'jet');
title(hax, 'Fitted Signal Amplitude');
ylabel(hax, 'Inversion Time [s]');
xlabel(hax, 'Echo Time [s]');

hf = figure(86);
clf(hf);
hax = axes(hf);
surf(hax, tau2, tau1, data.iLaplace2D.DataAmplitude.'-data.iLaplace2D.FitAmplitude.', ...
  'LineStyle', 'none');
colormap(hf, 'jet');
title(hax, {'Residuals', 'Measured Signal Amplitude - Fitted Signal Amplitude'});
ylabel(hax, 'Inversion Time [s]');
xlabel(hax, 'Echo Time [s]');

% T1-T2-map
hf = figure(84);
set(hf, 'Name', 'T1-T2-map');
ax(1) = subplot(10,10,[3:10,13:20,23:30,33:40,43:50,53:60,63:70,73:80], 'Parent', hf);
s = data.iLaplace2D.SpectrumAmplitude.';
s(~conv2(double(s~=0), ones(3,3), 'same')) = NaN; % nans are transparent
pcolor(ax(1), T2, T1, s * 100/FS);
set(ax(1), 'XScale', 'log', 'YScale', 'log', 'Tag', 'main');

shading(ax(1), 'interp');
cmap = colormap(ax(1), 'jet');
cmap = [1,1,1; cmap];
colormap(ax(1), cmap);
axis(ax(1), [min(T2) max(T2) min(T1) max(T1)]);

title(ax(1), ['Fitted Spectrum [% FS] (FS = ' num2str(Seq.iLaplace2D.FullScaleAmplitude*1e9) ' nT)']);
set(ax(1), 'XTickLabel', []);
set(ax(1), 'YTickLabel', []);
grid(ax(1), 'on');

ax(2) = subplot(10,10,2:10:72, 'Parent', hf);
semilogy(ax(2), sum(data.iLaplace2D.SpectrumAmplitude.'*100/FS,2), T1, 'LineWidth', 2);
grid(ax(2), 'on');
set(ax(2), 'YTickLabel', [], 'XDir', 'reverse', 'XAxisLocation', 'top', 'XTickLabelRotation', 90, 'Tag', 'tau1_spec');

ax(3) = subplot(10,10,83:90, 'Parent', hf, 'Tag', 'tau2_spec');
semilogx(ax(3), T2, sum(data.iLaplace2D.SpectrumAmplitude.'*100/FS,1), 'LineWidth', 2);
grid(ax(3), 'on');
set(ax(3), 'XTickLabel', [], 'YAxisLocation', 'right');

ax(4) = subplot(10,10,1:10:71, 'Parent', hf);
semilogy(ax(4), cumsum(sum(data.iLaplace2D.SpectrumAmplitude.'*100/FS,2),1), T1, 'LineWidth', 2);
ylabel(ax(4), 'T1 [s]');
xlim(ax(4), [0, 100]);
grid(ax(4), 'on');
set(ax(4), 'XDir', 'reverse', 'XAxisLocation', 'top', 'XTickLabelRotation', 90, 'Tag', 'tau1_accum');

ax(5) = subplot(10,10,93:100, 'Parent', hf, 'YAxisLocation', 'right');
semilogx(ax(5), T2, cumsum(sum(data.iLaplace2D.SpectrumAmplitude.'*100/FS,1),2), 'LineWidth', 2);
xlabel(ax(5), 'T2 [s]');
ylim(ax(5), [0, 100]);
grid(ax(5), 'on');
set(ax(5), 'YAxisLocation', 'right', 'Tag', 'tau2_accum');

linkaxes(ax([1,2,4]), 'y');
linkaxeskeep(ax([1,3,5]), 'x');
set(ax(1), 'XLim', [min(T2(:)), max(T2(:))], 'YLim', [min(T1(:)), max(T1(:))]);

% residual
% Seq.iLaplace2D.res = norm(d - G*s(:),2)^2;

if nargin == 2
  C = data;
  C.iLaplace2D.T1 = T1;
  C.iLaplace2D.T2 = T2;
  T1 = Seq;
else
  C = data.iLaplace2D.SpectrumAmplitude;
end
%
