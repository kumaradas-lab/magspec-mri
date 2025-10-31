function [Network, SeqOut] = sequence_Network_Gain(HW, Seq, Network)
%% Acquire gain from C2 TX to C1 RxTx channel
%
%   [Network, SeqOut] = sequence_Network_Gain(HW, Seq, Network)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

%% default parameters
iDevice = 1;  % FIXME: Support multiple MMRT devices

AQnSampleLatenzTemp = HW.RX(iDevice).nSampleLatenz;
HW.RX(iDevice).nSampleLatenz = 0;  % latency: number of clock ticks in decimated sampling frequency

if nargin == 1, Seq = struct(); end
if isemptyfield(Seq, 'fSample'), Seq.fSample = HW.RX(iDevice).fSample/1250/1; end
% if isemptyfield(Seq, 'plotSmith'), Seq.plotSmith = true; end
if ~isfield(Seq, 'RxMaxdBm'), Seq.RxMaxdBm = []; end
if ~isfield(Seq, 'GaindB'), Seq.GaindB = []; end
if isemptyfield(Seq, 'Gain')
  if ~isempty(Seq.RxMaxdBm)
    Seq.Gain = HW.RX(iDevice).Amplitude2LnaUin .* HW.RX(iDevice).LNAGain .* HW.RX(iDevice).Uin2VgaUin .* ...
               (1./(sqrt((10^(Seq.RxMaxdBm/10)*1e-3)*HW.RX(iDevice).Rin) .* sqrt(2) .* ...
                HW.RX(iDevice).Uin2VgaUin .* HW.RX(iDevice).VgaUout2MmrtUin .* ...
                HW.RX(iDevice).MmrtUin2AdcUin .* HW.RX(iDevice).AdcUin2Norm)) .* ...
               HW.RX(iDevice).VgaUout2MmrtUin .* HW.RX(iDevice).MmrtUin2AdcUin .* HW.RX(iDevice).AdcUin2Norm;
  elseif ~isempty(Seq.GaindB)
    Seq.Gain = 10^(Seq.GaindB/20) * HW.RX(iDevice).GainMax;
  else
    Seq.Gain = HW.RX(iDevice).GainDef;
  end
end
if isemptyfield(Seq, 'plotGain'), Seq.plotGain = true; end
if ~isfield(Seq, 'plotGainHandle'), Seq.plotGainHandle = []; end
if ~isfield(Seq, 'plotPhaseHandle'), Seq.plotPhaseHandle = []; end
if isemptyfield(Seq, 'plotGainRaw'), Seq.plotGainRaw = false; end
if isemptyfield(Seq, 'tMessung'), Seq.tMessung = 200/Seq.fSample; end
if isemptyfield(Seq, 'nMeasurements'), Seq.nMeasurements = 101; end
if isemptyfield(Seq, 'fSpan'), Seq.fSpan = 1e6; end
if isemptyfield(Seq, 'fCenter'), Seq.fCenter = HW.fLarmor; end
if isemptyfield(Seq, 'fTestStart'), Seq.fTestStart = Seq.fCenter-Seq.fSpan/2; end
if isemptyfield(Seq, 'fTestStop'), Seq.fTestStop = Seq.fCenter+Seq.fSpan/2; end
if isemptyfield(Seq, 'fTest')
  if Seq.nMeasurements >= 2
    Seq.fTest = linspace(Seq.fTestStart,Seq.fTestStop,Seq.nMeasurements);
  else
    Seq.fTest = Seq.fCenter;
  end
end
if ~isfield(Seq, 'TXPowerdBm'), Seq.TXPowerdBm = []; end
if ~isfield(Seq, 'TXPowerW'), Seq.TXPowerW = []; end
if isemptyfield(Seq, 'TXPower')
  if ~isempty(Seq.TXPowerdBm)
    Seq.TXPower = HW.TX(iDevice).dBm2Amp(HW, Seq.TXPowerdBm, 2);
  elseif ~isempty(Seq.TXPowerW)
    Seq.TXPower = (Seq.TXPowerW.*HW.TX(iDevice).Rout)^0.5./HW.TX(iDevice).Amp2Ueff50(2);
  else
    Seq.TXPower = HW.TX(iDevice).AmpDef/1000;
  end
end
if isemptyfield(Seq, 'Cal'), Seq.Cal = 0; end
if isemptyfield(Seq, 'CorrectAQWindowPhase'), Seq.CorrectAQWindowPhase = 1; end
nMessung = 0;

if Seq.nMeasurements <= nMessung
  Seq.CLTime = 100e-3;
  Seq.tRep = (Seq.tMessung+20/Seq.fSample)*Seq.nMeasurements + Seq.CLTime + ...
    10e-3 + HW.TX(iDevice).BlankOffset + HW.TX(iDevice).BlankPostset;
else
  Seq.CLTime = 50e-6;
  Seq.tRep = (Seq.tMessung + get_DeadTimeTX2RX(HW, Seq.fSample) + ...
    10/Seq.fSample + Seq.CLTime + 200e-6 + HW.TX(iDevice).BlankOffset + HW.TX(iDevice).BlankPostset) * ...
    ones(1, Seq.nMeasurements);
end

%% pulse program
if Seq.nMeasurements <= nMessung
  AQ.fSample = Seq.fSample;
  AQ.Start = cumsum(ones(Seq.nMeasurements,1)*(get_DeadTimeTX2RX(HW,AQ.fSample(1))+get_DeadTimeRX2TX(HW,AQ.fSample(1))+Seq.tMessung+1e-6))*ones(1,length(Seq.tRep));
  AQ.nSamples = round(Seq.tMessung*AQ.fSample(1));
  AQ.Frequency = Seq.fTest.';
  AQ.Phase = 0;
  AQ.Gain = Seq.Gain;
  AQ.ResetPhases = 1;
  TX(1).Channel = 2;
  TX(1).BlankOffset = TX.HW.BlankOffset;  % Blank before rf pulse
  TX(1).BlankPostset = TX.HW.BlankPostset;  % Blank after rf pulse
  TX(1).Duration = AQ.nSamples/AQ.fSample(1) + (get_DeadTimeTX2RX(HW,AQ.fSample(1))+get_DeadTimeRX2TX(HW,AQ.fSample(1)));
  TX(1).Start = AQ.Start;  %-get_DeadTimeTX2RX(HW,AQ.fSample(1));
  TX(1).Amplitude = Seq.TXPower;  %./(sin(pi*AQ.Frequency/(HW.TX.fSample))./(pi*AQ.Frequency/(HW.TX.fSample))); %  0 - 1
  TX(1).Frequency = AQ.Frequency;
  TX(1).Phase = 0;

  % TX.setSwitchMode = 1;
  % TX.setNetworkMode = 1;
  % TX.setNetworkForward = [0, 1];

else
  AQ.fSample = Seq.fSample*ones(1,length(Seq.tRep)); % 125e6 /  1 und 4 bis 8192
  % AQ.Start=get_DeadTimeTX2RX(HW,AQ.fSample(1))*ones(1,length(Seq.tRep));
  AQ.Latenz=1./AQ.fSample*HW.RX(iDevice).nSampleLatenz+1./HW.RX(iDevice).fSample*HW.RX(iDevice).nSampleRXLatenz;
  AQ.Start=0*ones(1,length(Seq.tRep))+AQ.Latenz;

  AQ.nSamples=round(Seq.tMessung*AQ.fSample(1))*ones(1,length(Seq.tRep));
  AQ.Frequency=[Seq.fTest];
  AQ.Phase=0;
  AQ.Gain=Seq.Gain;
  AQ.ResetPhases=zeros(1,length(Seq.tRep));
  AQ.ResetPhases=ones(1,length(Seq.tRep));
  AQ.ResetPhases(1)=1;
  % Seq.Gain

  TX(1).Channel=2;
  TX(1).BlankOffset=HW.TX(iDevice).BlankOffset; % Blank vor HF
  TX(1).BlankPostset=HW.TX(iDevice).BlankPostset; % Blank nach HF
  TX(1).Duration=Seq.tMessung+5/Seq.fSample*ones(1,length(Seq.tRep));
  TX(1).Start=0*ones(1,length(Seq.tRep))+HW.TX(iDevice).Latenz;
  TX(1).Amplitude=Seq.TXPower;%./(sin(pi*AQ.Frequency/(HW.TX.fSample))./(pi*AQ.Frequency/(HW.TX.fSample))); %  0 - 1
  TX(1).Frequency=AQ.Frequency;
  TX(1).Phase=0;

  % TX.setSwitchMode=1;
  % TX.setNetworkMode=1;
  % TX.setNetworkForward=[zeros(1,Seq.nMeasurements),ones(1,Seq.nMeasurements)];
end

for t = 1:HW.Grad(iDevice).n
  Grad(t).Time = NaN;
  Grad(t).Amp = 0;
end

%% measurement
% [~, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);
%
% %  plot_data_1D_TR(HW, data_1D);
% %   plot_data_1D(HW, data_1D);
% if Seq.nMeasurements<=nMessung
%   Network.BackRawMean=reshape(mean(data.data(:,:,1),1),[Seq.nMeasurements,1]);
%   Network.Frequency=SeqOut.AQ.Frequency(:,1);
%   Network.ForwardRawMean=reshape(mean(data.data(:,:,2),1),[Seq.nMeasurements,1]);
%
% else
%     Network.BackRawMean=reshape(mean(data.data(:,1,1:Seq.nMeasurements),1),[Seq.nMeasurements,1]);
%     Network.Frequency=SeqOut.AQ.Frequency(1,1:Seq.nMeasurements);
%     Network.ForwardRawMean=reshape(mean(data.data(:,1,Seq.nMeasurements+(1:Seq.nMeasurements)),1),[Seq.nMeasurements,1]);
% end

[dataRaw, SeqOut, dataGain] = set_sequence(HW, Seq, AQ, TX, Grad);
HW.RX(iDevice).nSampleLatenz = AQnSampleLatenzTemp;

%% evaluation
Channel = 1;  % FIXME: Add support for multiple acquisition channels?
iAQ = find([SeqOut.AQ(:).Channel] == Channel & [SeqOut.AQ(:).Device] == iDevice, 1, 'first');

if ~iscell(dataRaw), dataRaw = {dataRaw}; end

Network.dataGain = dataGain;

if Seq.nMeasurements <= nMessung
  Network.AmpRaw = reshape(mean(dataRaw{iAQ}(:,1:Seq.nMeasurements),1), [Seq.nMeasurements,1]);
  Network.FrequencyGain = SeqOut.AQ(iAQ).Frequency(:,1);
else
  Network.AmpRaw = reshape(mean(dataRaw{iAQ}(:,1:Seq.nMeasurements),1), [Seq.nMeasurements,1]);
  Network.FrequencyGain = reshape(SeqOut.AQ(iAQ).Frequency(1,1:Seq.nMeasurements), [Seq.nMeasurements,1]);
end

Network.AmpTXCorr = SeqOut.TX(1).AmplitudeCorr(:);
Network.GainRaw = Network.AmpRaw./SeqOut.TX(1).AmplitudeCorr(:)./SeqOut.AQ(iAQ).Gain;
% Network.Amp2Volt=abs(1./Network.AmpRaw.*(HW.TXPowerW2TXPower(2).*HW.RX(iDevice).Rin).^0.5.*SeqOut.TX(1).AmplitudeCorr(:).*SeqOut.AQ(iAQ).Gain);

if Seq.plotGainRaw
  Network.Gain = Network.GainRaw;
  % Network.Gain = Network.GainRaw./(sin(pi*Network.FrequencyGain/(HW.TX(iDevice).fSample))./(pi*Network.FrequencyGain/(HW.TX(iDevice).fSample)));
  plot_Gain(HW, Network, Seq);
end

CalSet = 0;
switch Seq.Cal
  case 'Cable'
    HW.NetworkCal.Cable = Network;
    CalSet = 1;
end
if CalSet == 1
  NetworkCal = HW.NetworkCal;
  NetworkCal.Seq = Seq;
  NetworkCal.FrequencyGain = Network.FrequencyGain;
  fileID = fopen(HW.NetworkCalPath);
  if fileID ~= -1
    fclose(fileID);
    save(HW.NetworkCalPath, 'NetworkCal');
  else
    if ~strcmp(HW.NetworkCalPath, '')
      save(HW.NetworkCalPath, 'NetworkCal');
    else
      disp('no HW.NetworkCalPath selected');
    end
  end

else
  Network.Gain = get_CalibratedGain(HW, Network);
  [Network.MinGain, Network.iMinGain] = min(Network.Gain);
  Network.fMinGain = Network.FrequencyGain(Network.iMinGain);
  [Network.MaxGain, Network.iMaxGain] = max(Network.Gain);
  Network.fMaxGain = Network.FrequencyGain(Network.iMaxGain);

  if Seq.plotGain
    plot_Gain(HW, Network, Seq);
  end

end

end
