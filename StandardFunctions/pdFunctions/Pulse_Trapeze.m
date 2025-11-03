function pulseData = Pulse_Trapeze(HW, Center, Pulse, varargin)
%% Pulse shape function for a rectangular RF pulse with a matched bandwidth
%
%   pulseData = Pulse_Rect_Slice(HW, Center, Pulse)
% or:
%   pulseData = Pulse_Rect_Slice(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments,  MaxLength, Frequency, Phase)
% additionally:
%   excitationAngleFactor = Pulse_Rect_Slice(HW, 'Amp')
%   bandwidthFactor = Pulse_Rect_Slice(HW, 'Time')
%
% This function can be used to create rectangular pulses. The bandwidth of these
% pulses is adapted such that the slice thickness corresponds to the FWHM of the
% intensity sinc^2 profile.
%
% INPUT:
%   HW        HW structure
%   Center    The center of the pulse in seconds (tRep).
%   Pulse     A structure with the following fields (if omitted or empty,
%             default values are used):
%     FlipAngle
%             The total (effective) flip angle of the pulse in radians (or the
%             units defined by FlipAngleFullTurn). It is used to set an
%             appropriate pulse amplitude (default: pi/2).
%     FlipAngleFullTurn
%             Value that defines a full turn in FlipAngle units (e.g. 360 for
%             degrees, or 2*pi for radians, default: 2*pi).
%     MaxNumberOfSegments
%             Maximum number of segments of the pulse. Ignored for rectangular
%             pulse which always consists of one segment (default: Inf).
%     MaxLength
%             Maximum length of the pulse in seconds (default: Inf).
%     Frequency
%             Frequency of the rf pulse in Hz (default: HW.fLarmor).
%     Phase   "Local" RF phase of the pulse with respect to the overall sequence
%             in degrees (default: 0).
%     Bandwidth
%             Bandwidth of the pulse in Hz (default:
%             max(1/Pulse.MaxLength, 2e3) where "FactorTime" is the
%             multiplier compared to a simple rectangular pulse with the same
%             flip angle.
%
% OUTPUT:
%   pulseData A structure with the following fields:
%     Start   Column vector with the start times of each component/block in
%             seconds (tRep).
%     Amplitude
%             Column vector with the amplitudes of each component/block in
%             Tesla.
%     Duration
%             Column vector with the durations of each component/block in
%             seconds.
%     Frequency
%             Column vector with the frequencies of each component/block in Hz.
%     Phase
%             Column vector with the phases in degrees of each component/block
%             with respect to the overall sequence.
%
% The additional syntax is used to return the amplitude and bandwidth factors.
% The "excitationAngleFactor" is the factor that must be applied to the
% amplitude of this pulse to have the same excitation angle as a rect pulse of
% the same length. The "bandwidthFactor" is the factor that must be applied to
% the duration of the pulse to have the same bandwidth (FWHM) as a rect pulse.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% Convert from syntax (2) to syntax (1)
if nargin < 3, Pulse = []; end
if ~isstruct(Pulse), Pulse = struct('Bandwidth', Pulse); end
if nargin > 3, Pulse.FlipAngle = varargin{1}; end
if nargin > 4, Pulse.MaxNumberOfSegments = varargin{2}; end
if nargin > 5, Pulse.MaxLength = varargin{3}; end
if nargin > 6, Pulse.Frequency = varargin{4}; end
if nargin > 7, Pulse.Phase = varargin{5}; end

%% default values
bwFactor = 0.442946; % solve sinc(pi*x) = 1/sqrt(2) for x -> x=0.442946
Pulse = set_EmptyField(Pulse, 'FlipAngle', pi/2);
Pulse = set_EmptyField(Pulse, 'FlipAngleFullTurn', 2*pi);
Pulse = set_EmptyField(Pulse, 'MaxNumberOfSegments', 51);
Pulse = set_EmptyField(Pulse, 'MaxLength', Inf);
Pulse = set_EmptyField(Pulse, 'Frequency', HW.fLarmor);
Pulse = set_EmptyField(Pulse, 'Phase', 0);
Pulse = set_EmptyField(Pulse, 'Bandwidth', max(1/Pulse.MaxLength/bwFactor, 2e3)); % FIXME: Is this a sensible default?
Pulse = set_EmptyField(Pulse, 'iDevice', 1);
Pulse = set_EmptyField(Pulse, 'AmpRatio', HW.TX(Pulse.iDevice).AmpRatio);
Pulse = set_EmptyField(Pulse, 'LinearRamp', 0); 

% Q=HW.TX(Pulse.iDevice).CoilQ;                 % quality factor of coil e.g: Q = 55
% D = 1/(2*Q);          % damping 
tau = 1 ./ (2*pi*HW.fLarmor./(2*HW.TX(Pulse.iDevice).CoilQ));  % decay time: time untill amplitude = (1-1/exp(1)) * (final amp) = 63% * (final amp)
% Q=HW.TX(Pulse.iDevice).CoilQDamped;                 % quality factor of coil e.g: Q = 55
if isempty(HW.TX(Pulse.iDevice).CoilQDamped)
  HW.TX(Pulse.iDevice).CoilQDamped = HW.TX(Pulse.iDevice).CoilQ;
end
tauDamped = 1 ./ (2*pi*HW.fLarmor./(2*HW.TX(Pulse.iDevice).CoilQDamped));  % decay time: time untill amplitude = (1-1/exp(1)) * (final amp) = 63% * (final amp)

%% short path
if nargin == 2
  if ischar(Center)
    if strcmp(Center, 'Amp')
      pulseData = 1;
    elseif strcmp(Center, 'Time')
      pulseData = bwFactor*1;
    else
      pulseData = NaN;
    end
  else
    pulseData = NaN;
  end
  return;
end

%% "real" path

if Pulse.MaxNumberOfSegments < 1
  error('Pulse.MaxNumberOfSegments must be >= 1');
end

DurationBW = bwFactor/Pulse.Bandwidth;

AmpDef = (Pulse.FlipAngle/Pulse.FlipAngleFullTurn) ./ (HW.GammaDef/2/pi) ./ DurationBW;
% AmpDef=1;                         % AmpDef=1;
AmpMax=AmpDef*Pulse.AmpRatio;     % AmpMax over drive ratio


t=(0:max([1,tau*50*HW.TX(Pulse.iDevice).fSample]))/HW.TX(Pulse.iDevice).fSample;
nMax=max(1,(Pulse.MaxNumberOfSegments-1)/2); % nMax max number of steps for the Ramp 

tRiseExp=-log(1-AmpDef/AmpMax)*tau;
SlopeRiseMax=AmpMax*(1+exp(-tRiseExp/tau)/tau);
tRiseEndMin=AmpDef/SlopeRiseMax;
tSegment=max(2,ceil(tRiseEndMin/nMax*HW.TX(Pulse.iDevice).fSample))/HW.TX(Pulse.iDevice).fSample;
n=ceil(tRiseEndMin/tSegment);

if AmpDef==AmpMax
  Arise=AmpMax*(1-exp(-t/tau));
  Adamp=-(AmpMax-AmpDef)+(AmpMax)*exp(-t/tauDamped*1);
  if isnan(Arise(1))
    Arise=0;
  else
    Arise=Arise(t<tau*10);
  end
  if isnan(Adamp(1))
    Adamp=0;
  else
    Adamp=Adamp(t<tauDamped*10);
  end
  if sum(t<tau*10)>10000; error('sum(t<tau*5)>10000');end
  [~,IAmin]=min(abs(bsxfun(@minus,Arise,Adamp.')),[],1);
  intPuls=DurationBW*AmpDef;
  intRise=cumsum(Arise.*[0,diff(t(1:numel(Arise)))]);
  intDamp=cumsum(Adamp.*[diff(t(1:numel(Adamp))),0],'reverse');
  intRiseDamp=intRise+intDamp(IAmin);
  [~,IintMin]=min(abs(intPuls-intRiseDamp));
  if IintMin==numel(Arise);
    nAmpDef=ceil((intPuls-intRiseDamp(end))/AmpDef*HW.TX(Pulse.iDevice).fSample/2)*2;
  else
    nAmpDef=0;
  end
  
  AmpAll=[Arise(1:IintMin),AmpDef*ones(1,nAmpDef), Adamp(IAmin(IintMin):end)];
  tPulse=(0:numel(AmpAll))/HW.TX(Pulse.iDevice).fSample;
  intAPulse=cumsum(AmpAll.*diff(tPulse(1:numel(AmpAll)+1)));
  tCenterI=find(intAPulse>=intAPulse(end)/2,1,'first');
  tStart=tPulse(tCenterI);
  
  intARise=sum(Arise(1:IintMin))*diff(t(1:2));
  intAdamp=sum(Adamp(IAmin(IintMin):end))*diff(t(1:2));
  intAmpDef=AmpDef*nAmpDef*diff(t(1:2));
  intExact=DurationBW*AmpDef;
  AmpCorrection=intExact/(intARise+intAmpDef+intAdamp);
  AmpDefCorr=AmpDef*AmpCorrection;
  
%   durationRise=t(IintMin);
%   durationDamp=t(end)-t(IAmin(IintMin));
  durationAmpDef=t(IintMin)+nAmpDef/HW.TX(Pulse.iDevice).fSample;
  durationRise=[];
  durationDamp=[];
  AmpRise=[];
  AmpDamp=[];
  n=0;
  
%   figure(1)
%   subplot(2,1,1)
%   plot(t(1:numel(Arise)),Arise,t(1:numel(Adamp)),Adamp)
%   subplot(2,1,2)
%   t=(0:tau*50*HW.TX(Pulse.iDevice).fSample)/HW.TX(Pulse.iDevice).fSample;
%   plot(t(1:IintMin+numel(Arise)-IAmin(IintMin)+1)-tStart,[Arise(1:IintMin),Adamp(IAmin(IintMin):end)])
%   grid on
    
elseif n==1
  AmpRise=AmpMax;
  Arise=AmpMax*(1-exp(-t/tau));
  AmpDamp=AmpDef-AmpMax;
  Adamp=-(AmpMax-AmpDef)+(AmpMax)*exp(-t/tau*1);

  AiRise=find(Arise<AmpDef,1,'last');
  durationRise=t(AiRise);
  intARise=sum(Arise(1:AiRise))*diff(t(1:2));
  triesAdef=round(intARise/AmpDef*HW.TX(Pulse.iDevice).fSample)/HW.TX(Pulse.iDevice).fSample;

  AiDamp=find(Adamp<0,1,'first');
  if isempty(AiDamp)||AiDamp<=2;AiDamp=3;end
  durationDamp=t(AiDamp);
  intAdamp=sum(Adamp(1:AiDamp))*diff(t(1:2));
  tdampAdef=round(intAdamp/AmpDef*HW.TX(Pulse.iDevice).fSample)/HW.TX(Pulse.iDevice).fSample;
  
  durationAmpDefExact=DurationBW-(triesAdef+tdampAdef);
  durationAmpDef=ceil((durationAmpDefExact-1e-12)*HW.TX(Pulse.iDevice).fSample/2)/HW.TX(Pulse.iDevice).fSample*2;
  AmpCorrection=(durationAmpDefExact+(triesAdef+tdampAdef))/(durationAmpDef+(triesAdef+tdampAdef));
  AmpRise=AmpRise*AmpCorrection;
  AmpDamp=AmpDamp*AmpCorrection;
  AmpDefCorr=AmpDef*AmpCorrection;
  Aip90=max(0,find(durationAmpDef<t+1e-12,1,'first')-1);
  AmpAll=[Arise(1:AiRise), AmpDef*ones(1,Aip90), Adamp(1:AiDamp)];
  tPulse=(0:numel(AmpAll))/HW.TX(Pulse.iDevice).fSample;
  intAPulse=cumsum(AmpAll.*diff(tPulse(1:numel(AmpAll)+1)));
  tCenterI=find(intAPulse>=intAPulse(end)/2,1,'first');
  tStart=tPulse(tCenterI);
  
elseif Pulse.LinearRamp

  AriseExp=AmpMax*(1-exp(-t/tau));
  AiRise=find(AriseExp>AmpDef,1,'first');
  SlopeRise=diff(AriseExp(AiRise-1:AiRise))/diff(t(AiRise-1:AiRise));
  triseEnd=AmpDef/SlopeRise;
  tiend=find(triseEnd>t,1,'last');
  n=floor(tiend/max(2,tiend/(n+1)))-1;
  tRise=round(linspace(0,t(tiend),n+1)*HW.TX(Pulse.iDevice).fSample)/HW.TX(Pulse.iDevice).fSample;
  tdamp=tRise;
  durationRise=round(diff(tRise)*HW.TX(Pulse.iDevice).fSample)/HW.TX(Pulse.iDevice).fSample;
  AmpRise=mean([(tRise(1:end-1)./tRise(end)*AmpDef+(AmpMax-AmpDef));(tRise(2:end)./tRise(end)*AmpDef+(AmpMax-AmpDef))])* 1;
%       AmpRise=AmpRise(1:end-1);
  tRise=tRise(1:end-1);

  durationDamp=durationRise;
  AmpDamp=mean([-tdamp(1:end-1)./tdamp(end)*AmpDef+(AmpDef-(AmpMax-AmpDef));(-tdamp(2:end)./tdamp(end)*AmpDef+(AmpDef-(AmpMax-AmpDef)))])* 1;
%       AmpDamp=AmpDamp(1:end-1);
  tdamp=tdamp(1:end-1);

  triesAdef=sum(durationRise)/2;
  tdampAdef=sum(durationDamp)/2;
  durationAmpDefExact=DurationBW-(triesAdef+tdampAdef);
  durationAmpDef=ceil((durationAmpDefExact-1e-12)*HW.TX(Pulse.iDevice).fSample/2)/HW.TX(Pulse.iDevice).fSample*2;
  AmpCorrection=(durationAmpDefExact+(triesAdef+tdampAdef))/(durationAmpDef+(triesAdef+tdampAdef));
  AmpRise=AmpRise*AmpCorrection;
  AmpDamp=AmpDamp*AmpCorrection;
  AmpDefCorr=AmpDef*AmpCorrection;
  tStart=sum([durationRise,durationAmpDef,durationDamp])/2;
  
else % Ramp exponetial
  AmpRise=AmpMax;
  Arise=AmpMax*(1-exp(-t/tau));
  AiRise=find(Arise<AmpDef,1,'last');
  intRise=cumsum(Arise.*diff(t(1:2)));
  intExact=DurationBW*AmpDef;
  AiRise=min([AiRise,find(intExact/2<intRise,1,'first')]);
  intARise=sum(Arise(1:AiRise))*diff(t(1:2));
  durationRise=t(AiRise);

%   Adamp=AmpMax-(AmpMax+AmpMax)*exp(-t/tau*1);
  AiDamp=AiRise;
  intADamp=intARise;
  n=Pulse.MaxNumberOfSegments-2;
  n=floor(AiDamp/max(2,AiDamp/(n+1)))-1;
  tDamp=round(linspace(t(AiDamp),0,n+1)*HW.TX(Pulse.iDevice).fSample)/HW.TX(Pulse.iDevice).fSample;
  durationDamp=round(abs(diff(tDamp))*HW.TX(Pulse.iDevice).fSample)/HW.TX(Pulse.iDevice).fSample;
  AmpDamp=mean([AmpMax-(AmpMax+AmpMax)*exp(-tDamp(1:end-1)/tau*1) ; AmpMax-(AmpMax+AmpMax)*exp(-tDamp(2:end)/tau*1)]);

  nAmpDef=ceil((intExact-intARise-intADamp)/AmpDef*HW.TX(Pulse.iDevice).fSample/2)*2;
  durationAmpDef=nAmpDef/HW.TX(Pulse.iDevice).fSample;
  AmpAll=[Arise(1:AiRise),AmpDef*ones(1,nAmpDef), Arise(AiRise:-1:1)];
  tPulse=(0:numel(AmpAll))/HW.TX(Pulse.iDevice).fSample;
  intAPulse=cumsum(AmpAll.*diff(tPulse(1:numel(AmpAll)+1)));
  tCenterI=find(intAPulse>=intAPulse(end)/2,1,'first');
  tStart=tPulse(tCenterI);

  intAmpDef=AmpDef*nAmpDef*diff(t(1:2));
  AmpCorrection=intExact/(intARise+intAmpDef+intADamp);
  AmpDefCorr=AmpDef*AmpCorrection;
  AmpRise=AmpRise*AmpCorrection;
  AmpDamp=AmpDamp*AmpCorrection;
  if durationAmpDef==0;durationAmpDef=[];AmpDefCorr=[];end

end


pulseData.Duration     = [durationRise.';durationAmpDef;durationDamp.'];          % duration of rf-pulse
pulseData.Start        = cumsum([0;pulseData.Duration(1:end-1)])-tStart+Center;          % start time of rf-pulse % nan -> no rf-pulse  
pulseData.Frequency    = [zeros(numel(AmpRise),1);zeros(numel(AmpDefCorr),1);zeros(numel(AmpDamp),1)]+Pulse.Frequency;                    % frequency of rf-pulse
pulseData.Phase        = [zeros(numel(AmpRise),1);zeros(numel(AmpDefCorr),1);sign(AmpDamp.')*90-90]+Pulse.Phase;          % phase of rf-pulse
pulseData.Amplitude    = [AmpRise.';AmpDefCorr;abs(AmpDamp.')*1];                 % amplitude of B1p in T


if Pulse.MaxLength < max([abs(pulseData.Start(1)-Center),abs(pulseData.Start(end)+pulseData.Duration(end)-Center)]*2);
  error('maxLength of HF Pulse to short')
end

% pulseData.Start = -BlockLength/2+Center;
% pulseData.Amplitude = gain;
% pulseData.Duration = BlockLength;
% pulseData.Frequency = Pulse.Frequency;
% pulseData.Phase = Pulse.Phase;

end
