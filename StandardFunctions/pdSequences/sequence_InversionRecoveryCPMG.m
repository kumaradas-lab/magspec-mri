function [data, SeqOut] = sequence_InversionRecoveryCPMG(HW, Seq)
%% Inversion RecoveryCPMG method to measure T2/T1 map
%
% This function is deprecated and should no longer be used.
% Use sequence_RecoveryCPMG instead.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2021 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------


% default parameters
if isemptyfield(Seq, 'PlotTR'),  Seq.PlotTR = 0; end
if isemptyfield(Seq, 'Plot'),  Seq.Plot = 0; end

if isemptyfield(Seq, 'FlipPulse'),  Seq.FlipPulse = @Pulse_Rect; end
if isemptyfield(Seq, 'InvertPulse'),  Seq.InvertPulse = @Pulse_Rect; end
if isemptyfield(Seq, 'InvertPulse2'),  Seq.InvertPulse2 = @Pulse_Rect; end

if isemptyfield(Seq, 'tRelax'),  Seq.tRelax = 0.3; end
if isemptyfield(Seq, 'tEcho'),  Seq.tEcho = 5e-3; end
if isemptyfield(Seq, 'nEcho'),  Seq.nEcho = 20; end
if isemptyfield(Seq, 'Tau1'),  Seq.Tau1 = 10e-3; end
if isemptyfield(Seq, 'nTau1'),  Seq.nTau1 = 10; end
if isemptyfield(Seq, 'Tau1Start'),  Seq.Tau1Start = Seq.Tau1; end
if isemptyfield(Seq, 'Tau1End'),  Seq.Tau1End = Seq.Tau1*Seq.nTau1; end
if isemptyfield(Seq, 'Tau1Log'),  Seq.Tau1Log = 0; end
if isemptyfield(Seq, 'tAQEcho'),  Seq.tAQEcho = min(1e-3, Seq.tEcho/8); end
if isemptyfield(Seq, 'fSample'),  Seq.fSample = max(100e3, 7/Seq.tAQEcho); end
if isemptyfield(Seq, 'CLTime'),  Seq.CLTime = 1e-6; end
if isemptyfield(Seq, 'PlotT1'),  Seq.PlotT1 = 1; end
if isemptyfield(Seq, 'PlotT2'),  Seq.PlotT2 = 1; end
if isemptyfield(Seq, 'FitT1'),  Seq.FitT1 = 1; end
if isemptyfield(Seq, 'FitT2'),  Seq.FitT2 = 0; end

if Seq.Tau1Log
  Seq.Tau1 = logspace(log10(Seq.Tau1Start), log10(Seq.Tau1End), Seq.nTau1);
else
  Seq.Tau1 = linspace(Seq.Tau1Start, Seq.Tau1End, Seq.nTau1);
end
Seq.iLaplace2D.Recovery =  'Inversion';
Seq.CalculateFFTOfData = 0;

iDevice = 1;  % FIXME: Support multiple MMRT devices


% AQ
AQ.fSample = 1/(round(HW.RX(iDevice).fSample/Seq.fSample)/HW.RX(iDevice).fSample);
AQ.nSamples = round(Seq.tAQEcho*AQ.fSample);
if AQ.nSamples==0, AQ.nSamples=1; end
AQ.Frequency = HW.fLarmor;


%TX
Seq.t90 = max(1e-6, HW.tFlip90Def * Seq.FlipPulse(HW,'Amp'));
Seq.t90BW = 1/HW.tFlip90Def * Seq.FlipPulse(HW,'Time');

Seq.tInvert = HW.tFlip180Def * Seq.InvertPulse(HW,'Amp');
Seq.tInvert2 = HW.tFlip180Def * Seq.InvertPulse2(HW,'Amp');

pulseInvert =   Seq.InvertPulse(HW, 0, 1/Seq.tInvert*Seq.InvertPulse(HW,'Time'),    pi*HW.tFlip180Def/(HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef),    51, Seq.tInvert,    AQ.Frequency(1), 0);
pulse90 =       Seq.FlipPulse(  HW, 0, Seq.t90BW,                                   pi*HW.tFlip90Def/(HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef),      51, Seq.t90,        AQ.Frequency(1), 0);
pulseInvert2 =  Seq.InvertPulse2(HW, 0, 1/Seq.tInvert2*Seq.InvertPulse2(HW,'Time'),  pi*HW.tFlip180Def/(HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef),    51, Seq.tInvert2,   AQ.Frequency(1), 90);

Seq.tRep=repmat([0,Seq.tEcho/2,repmat(Seq.tEcho,1,Seq.nEcho),Seq.tRelax],1,Seq.nTau1)+reshape([Seq.Tau1(:).';zeros(1+Seq.nEcho+1,length(Seq.Tau1))],1,(Seq.nEcho+3)*length(Seq.Tau1));

TX.Start=nan(max([length(pulseInvert.Start),length(pulse90.Start),length(pulseInvert2.Start)]),Seq.nEcho+3);
TX.Duration=TX.Start;
TX.Amplitude=TX.Start;
TX.Frequency=TX.Start;
TX.Phase=TX.Start;

TX.Start(1:size(pulseInvert.Start,1),1)=pulseInvert.Start;
TX.Duration(1:size(pulseInvert.Start,1),1)=pulseInvert.Duration;
TX.Amplitude(1:size(pulseInvert.Start,1),1)=pulseInvert.Amplitude;
TX.Frequency(1:size(pulseInvert.Start,1),1)=pulseInvert.Frequency;
TX.Phase(1:size(pulseInvert.Start,1),1)=pulseInvert.Phase;

TX.Start(1:size(pulse90.Start,1),2)=pulse90.Start;
TX.Duration(1:size(pulse90.Start,1),2)=pulse90.Duration;
TX.Amplitude(1:size(pulse90.Start,1),2)=pulse90.Amplitude;
TX.Frequency(1:size(pulse90.Start,1),2)=pulse90.Frequency;
TX.Phase(1:size(pulse90.Start,1),2)=pulse90.Phase;

TX.Start(1:size(pulseInvert2.Start,1),2+(1:Seq.nEcho))=repmat(pulseInvert2.Start,1,Seq.nEcho);
TX.Duration(1:size(pulseInvert2.Start,1),2+(1:Seq.nEcho))=repmat(pulseInvert2.Duration,1,Seq.nEcho);
TX.Amplitude(1:size(pulseInvert2.Start,1),2+(1:Seq.nEcho))=repmat(pulseInvert2.Amplitude,1,Seq.nEcho);
TX.Frequency(1:size(pulseInvert2.Start,1),2+(1:Seq.nEcho))=repmat(pulseInvert2.Frequency,1,Seq.nEcho);
TX.Phase(1:size(pulseInvert2.Start,1),2+(1:Seq.nEcho))=repmat(pulseInvert2.Phase,1,Seq.nEcho);

TX.Start=repmat(TX.Start,1,Seq.nTau1);
TX.Duration=repmat(TX.Duration,1,Seq.nTau1);
TX.Amplitude=repmat(TX.Amplitude,1,Seq.nTau1);
TX.Frequency=repmat(TX.Frequency,1,Seq.nTau1);
TX.Phase=repmat(TX.Phase,1,Seq.nTau1);


% AQ
AQ.Start    =   repmat(  ([nan,nan,repmat(Seq.tEcho/2-AQ.nSamples/AQ.fSample/2,1,Seq.nEcho),nan])  ,1,Seq.nTau1);
AQ.ResetPhases=[1,zeros(1,length(Seq.tRep)-1)];


% Gradients
for t = 1:HW.Grad(iDevice).n
  Grad(t).Time = NaN;
  Grad(t).Amp = 0;
  Grad(t).Repeat = [0, ones(1,length(Seq.tRep)-1)];
end


% run measurement
if Seq.Plot || Seq.PlotTR
    [~, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);
else
    [~, SeqOut, data] = set_sequence(HW, Seq, AQ, TX, Grad);
end


% Plot
if SeqOut.Plot
  plot_data_1D(HW, data_1D);
end

if SeqOut.PlotTR
  plot_data_1D_TR(HW, data_1D);
end


Channel = 1;  % FIXME: Add support for multiple acquisition channels?
iAQ = find([SeqOut.AQ(:).Channel] == Channel & [SeqOut.AQ(:).Device] == iDevice, 1, 'first');

% discard data from all but the used channels
data = data(iAQ);

data.SampleEchoTau1=reshape(data.data(~isnan(data.data)), AQ.nSamples, SeqOut.nEcho, SeqOut.nTau1);
data.SampleEchoTau1Time=reshape(data.time_all(~isnan(data.data)), AQ.nSamples, SeqOut.nEcho, SeqOut.nTau1);
data.MeanEchoTau1=squeeze(mean(data.SampleEchoTau1,1));
% data.MeanEchoTau1PhaseOffset=mean(angle(data.MeanEchoTau1(1:min(5,end),1)));              % Phase at first Tau1
data.MeanEchoTau1PhaseOffset=mean(angle(data.MeanEchoTau1(1:round(min(5,end/8)),end)));  % Phase at long Tau1
data.MeanEchoTau1PhaseCorrected=data.MeanEchoTau1*exp(-1i*data.MeanEchoTau1PhaseOffset);
data.EchoTime=cumsum(repmat(SeqOut.tEcho,SeqOut.nEcho,SeqOut.nTau1),1);
data.Tau1Time=repmat(SeqOut.Tau1(:).',SeqOut.nEcho,1);

figure(80);
% imagesc(data.EchoTime(:,1),data.Tau1Time(1,:),real(data.MeanEchoTau1PhaseCorrected))
hp=pcolor(data.EchoTime,data.Tau1Time,real(data.MeanEchoTau1PhaseCorrected));
set(hp, 'ZData',real(data.MeanEchoTau1PhaseCorrected) )
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'linear')
ylabel('Tau1 / sec')
xlabel('Echo time / sec')
shading interp
colorbar

%     figure(6)
%     imagesc(data.EchoTime(:,1),data.Tau1Time(1,:),real(data.MeanEchoTau1PhaseCorrected))
% %     pcolor(data.EchoTime,data.Tau1Time,real(data.MeanEchoTau1PhaseCorrected))
%     xlabel('Tau1 / sec')
%     ylabel('Echo time / sec')
%     shading interp
%     colorbar
%

pause(0.5)
figure(81)
clf
hold off
surface(data.EchoTime,data.Tau1Time,real(data.MeanEchoTau1PhaseCorrected))
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'linear')
ylabel('Tau1 / sec')
xlabel('Echo time / sec')
view(3)
colorbar
grid off
shading interp


%% T1
if SeqOut.FitT1
  t1=nan;
  T1=nan;
  if SeqOut.nTau1>3;
    T1 = fit_exp(real(data.MeanEchoTau1PhaseCorrected(1,:)),data.Tau1Time(1,:),SeqOut.PlotT1,0,0,1);
    t1=T1.tau;
    T1.half=log(2)*T1.tau;
    T1.half1=log(2)*T1.tau1;
    T1.half2=log(2)*T1.tau2;
    data.T1=T1;
  else
    disp('nTau1 too low')
  end
end

%% T2
if SeqOut.FitT2
  fh=figure(82);
  clf(fh)
  ax(1)=subplot(1,1,1, 'Parent',fh);
  if SeqOut.nEcho>3;
    for t=1:SeqOut.nTau1
      % fitting exponential functions to the data
      T2(t) = fit_exp(real(data.MeanEchoTau1PhaseCorrected(:,t)),data.EchoTime(:,t),SeqOut.PlotT2,0,0,0);
      figure(82)
      hold(ax(1), 'all');
      plot(ax(1),T2(t).timeCorrected,T2(t).dataPhaseCorrectedReal,'b');
      % legend(ax(1),'real');
      xlabel(ax(1),'time / s')
      ylabel(ax(1),'amplitude')
      title(ax(1),{[' ', ' ', ' ']})

      plot(ax(1),T2(t).functionTime,T2(t).functionAmpDouble,'b-.');
      legend(ax(1),'real','fit double');

      hold(ax(1), 'off')
    end

    data.T2=T2;
  else
    disp('nEcho too low')
  end
end
SeqOut.iLaplace2D.T1Start=SeqOut.Tau1Start*2;
SeqOut.iLaplace2D.T1End=SeqOut.Tau1End/2;
SeqOut.iLaplace2D.T2Start=SeqOut.tEcho*2;
SeqOut.iLaplace2D.T2End=SeqOut.tEcho*SeqOut.nEcho/2;

end
