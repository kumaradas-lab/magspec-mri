function [data, SeqOut, mySave, HW] = sequence_Grad_Echo_Slice(HW, mySave, Seq, SliceSelect)
% This function is deprecated and should no longer be used.
% Use sequence_Flash instead.

%----------- Parameter ------------------------------------------------
if ~isfield(Seq, 'WaitSeqenzeTime')
  Seq.WaitSeqenzeTime=[];
end
if isemptyfield(Seq, 'StartSequenceTime')
  if isempty(Seq.WaitSeqenzeTime)
    Seq.WaitSeqenzeTime=0;
  else
    Seq.StartSequenceTime=now*24*3600+Seq.WaitSeqenzeTime;
  end
end
if isemptyfield(Seq, 'plot_k_image')
  Seq.plot_k_image = 0;
end
if isemptyfield(Seq, 'plotPhase')
  Seq.plotPhase = 1;
end

if isemptyfield(Seq, 'RandomTXRXPhase')
  Seq.RandomTXRXPhase = 0;
end
if ~isfield(Seq,'AQPhaseOffset')
  Seq.AQPhaseOffset = [];
end
if ~isfield(Seq, 'TXPhaseOffset')
  Seq.TXPhaseOffset = [];
end
if isemptyfield(Seq, 'TXRXPhaseIncrement')
  Seq.TXRXPhaseIncrement = HW.Constant.GoldenAngle138;
end
if isemptyfield(Seq, 'HzPixMin')
  Seq.HzPixMin = 0;
end

if nargin < 4
  SliceSelect = struct();
end
if isemptyfield(SliceSelect, 'thickness')
  SliceSelect.thickness = 1000000;
end
if isemptyfield(SliceSelect, 'CenterRot')
  SliceSelect.CenterRot = [0.0, 0.0, 0.0];
end
if isemptyfield(SliceSelect, 'alfa')
  SliceSelect.alfa = 0;
end
if isemptyfield(SliceSelect,'phi')
  SliceSelect.phi = 0;
end
if isemptyfield(SliceSelect,'theta')
  SliceSelect.theta = 0;
end


% AQSlice.thickness=SliceSelect.thickness;
if isempty(whos('global','SliceSelect'))

%-----------------------------------------------------------------------

    SliceSelect.Center=SliceSelect.CenterRot([2,3,1]).*[1,-1,1];
    SliceSelect.CenterRot=SliceSelect.Center([3,1,2]).*[1,1,-1];
    [SliceSelect.normV(1),SliceSelect.normV(2),SliceSelect.normV(3)]=sph2cart(SliceSelect.theta,-SliceSelect.phi,1);
    SliceSelect.R=SliceSelect.CenterRot*SliceSelect.normV.';
    SliceSelect.Rauf=SliceSelect.R*SliceSelect.normV;
    SliceSelect.RaufCenter=SliceSelect.CenterRot-SliceSelect.Rauf;
    SliceSelect.CenterRauf=SliceSelect.Rauf-SliceSelect.CenterRot;
    SliceSelect.CenterRaufImage=tpaRotate( SliceSelect.CenterRauf',-SliceSelect.alfa ,-SliceSelect.phi, -SliceSelect.theta)';
    SliceSelect.CenterRaufImage=SliceSelect.CenterRaufImage([3,2,1]).*[1,-1,1];
    SliceSelect.Center2OriginImage=[SliceSelect.CenterRaufImage([1,2]),SliceSelect.R];
else
    clear SliceSelect;
    global SliceSelect;
end

%  [HW,mySave]  = Find_Frequency( HW, mySave);

SliceSelect.Flip=acos(exp(-Seq.tRep/Seq.T1));  %in RAD
% SliceSelect.thickness=AQSlice.thickness; %Achtung
SliceSelect.ReadOS=Seq.ReadOS;
SliceSelect.PhaseOS=Seq.PhaseOS;
AQSlice=SliceSelect;

Seq.tTxSlicemax=max(HW.TX.Amp2FlipPiIn1Sec/(HW.TX.AmpDef*0.9)*AQSlice.Flip/pi*Seq.SlicePulse(HW,'Amp'),  1/(min(HW.Grad.MaxAmp(1:3))*HW.GammaDef/2/pi*AQSlice.thickness)*Seq.SlicePulse(HW,'Time'));

Seq.tRep=Seq.tRep*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
Seq.CLTime=50e-6;
Seq.tInvertmax=HW.TX.Amp2FlipPiIn1Sec/HW.TX.AmpDef;
Seq.TxSliceBW=1/Seq.tTxSlicemax*Seq.SlicePulse(HW,'Time');
Seq.tInvert=0;%Seq.tInvertmax+2*HW.Grad.tEC;
Seq.tSlice=Seq.tTxSlicemax+2*HW.Grad.tRamp+2*HW.Grad.tEC;
Seq.tGrad=Seq.tEcho/2-Seq.tSlice/4;
Seq.tAQmaxt=Seq.tEcho-2*HW.Grad.tRamp-2*HW.Grad.tEC-Seq.tSlice/2;
if Seq.HzPixMin ~= 0 &&  1/Seq.HzPixMin>Seq.tAQmaxt
  warning(['Seq.HzPixMin is too low. ' num2str(ceil(1/Seq.tAQmaxt)) ' Hz is used'])
end

Seq.tAQmax=min(1/Seq.HzPixMin, Seq.tAQmaxt);


AQSlice.nPhase3D=0;
Seq=createSeq(Seq,AQSlice,HW);

% AQPInc=0;
% TXPInc=AQPInc;
if Seq.RandomTXRXPhase
  if isempty(Seq.AQPhaseOffset)
    Seq.AQPhaseOffset=rand(1,AQSlice.nPhase*AQSlice.PhaseOS)*360;
  end
  if isempty(Seq.TXPhaseOffset)
    Seq.TXPhaseOffset=Seq.AQPhaseOffset;
  end
else
  if isempty(Seq.AQPhaseOffset)
    Seq.AQPhaseOffset=zeros(1,AQSlice.nPhase*AQSlice.PhaseOS);
  end
  if isempty(Seq.TXPhaseOffset)
    Seq.TXPhaseOffset=Seq.AQPhaseOffset;
  end
end
% pn=rand(1,AQSlice.nPhase*AQSlice.PhaseOS)*2*pi;
% AQ
AQ.Start=(Seq.tSlice/2+Seq.tEcho-Seq.tAQ/2)*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
AQ.fSample=Seq.fSample*AQSlice.ReadOS;             %125e6 /  1 und 4 bis 8192
AQ.nSamples=AQSlice.nRead*AQSlice.ReadOS;
AQ.Frequency=Seq.fAQCenter;
AQ.Phase=(cumsum(((Seq.PhaseShift/AQSlice.PhaseOS))*ones(1,AQSlice.nPhase*AQSlice.PhaseOS))-mod(Seq.fAQCenter-HW.fLarmor,1./Seq.tEcho).*Seq.tEcho.*360)+mod(cumsum(cumsum(ones(1,AQSlice.nPhase*AQSlice.PhaseOS)*Seq.TXRXPhaseIncrement))-Seq.TXRXPhaseIncrement,360)+Seq.AQPhaseOffset;
AQ.Dur=AQ.nSamples/AQ.fSample;
AQ.ResetPhases=[1,zeros(1,AQSlice.nPhase*AQSlice.PhaseOS-1)];
% AQ.Gain=1; % 0.0032 - 1

% HF TX
TX.BlankOffset=20e-6; % Blank vor HF
pulseData = Seq.SlicePulse(HW, Seq.tSlice/2, Seq.TxSliceBW, AQSlice.Flip, 50,  Seq.tTxSlicemax, Seq.fTxSlice, -90);

TX.Duration=[pulseData.Duration];
TX.Start=[pulseData.Start]; %!!!!!
TX.Amplitude=[pulseData.Amplitude]; %  0 - 1
TX.Frequency=[pulseData.Frequency];
TX.Phase=pulseData.Phase *ones(1,AQSlice.nPhase*AQSlice.PhaseOS)+ones(size(pulseData.Phase))*mod(cumsum(cumsum(ones(1,AQSlice.nPhase*AQSlice.PhaseOS)*(Seq.TXRXPhaseIncrement)))-Seq.TXRXPhaseIncrement,360)+ones(size(pulseData.Phase))*Seq.TXPhaseOffset;

% if size(TX.Start,2)<size(AQ.Start,2)
%     TX.Duration=TX.Duration*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
%     TX.Start=TX.Start*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
%     TX.Amplitude=TX.Amplitude*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
%     TX.Frequency=TX.Frequency*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
%     TX.Phase=TX.Phase*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
% end


GradTime=cumsum([...
    0;...            % slice start
    HW.Grad.tRamp;...
    Seq.tTxSlicemax+2*HW.Grad.tEC;...
    HW.Grad.tRamp;... % slice stop
    HW.Grad.tRamp;... % Grads start
    Seq.tGrad-2*HW.Grad.tRamp;...
    HW.Grad.tRamp;... % Grads stop
%     Seq.tInvert+2*HW.Grad.tEC;...  % invert
%     HW.Grad.tRamp;... % Spoil start
%     Seq.Spoil*Seq.tEcho/2;...
%     HW.Grad.tRamp;... % Spoil stop
    (Seq.tAQmaxt-Seq.tAQmax)/2;
    HW.Grad.tRamp;... % GradRead start
    Seq.tAQmax+2*HW.Grad.tEC;...
    HW.Grad.tRamp;... % GradRead ende
    HW.Grad.tRamp;...
    Seq.tEcho*Seq.Spoil;...
    HW.Grad.tRamp;
    0]);
GradTime(end)=Seq.tRep(1); %!!!!!!!!!!!!!!!

z=zeros(1,AQSlice.nPhase*AQSlice.PhaseOS);
o=ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
for t=1:3;
    Grad(t).Time=GradTime*o;
    Grad(t).Amp=[...
        z;...
%         Seq.AmpSlice(tt)*o;...
%         Seq.AmpSlice(tt)*o;...
%         z;...
%         Seq.AmpPhase(tt,:)+Seq.AmpReadD(tt)+Seq.AmpSpoilD(tt)+Seq.AmpSliceD(tt);...
%         Seq.AmpPhase(tt,:)+Seq.AmpReadD(tt)+Seq.AmpSpoilD(tt)+Seq.AmpSliceD(tt);...

        Seq.AmpSlice(t)*o;...
        Seq.AmpSlice(t)*o;...
        z;...
        -Seq.AmpPhase(t,:)+Seq.AmpReadD(t)+Seq.AmpSliceD(t);...
        -Seq.AmpPhase(t,:)+Seq.AmpReadD(t)+Seq.AmpSliceD(t);...


%         z;...
%         z;...
%         Seq.AmpSpoil(t)*o;...
%         Seq.AmpSpoil(t)*o;...
        z;...
        z;...
        Seq.AmpRead(t)*o;...
        Seq.AmpRead(t)*o;...
        z;...
%         Seq.AmpSpoilP2(t,:);...
%         Seq.AmpSpoilP2(t,:);...
        -Seq.AmpSpoil(t)*o;...
        -Seq.AmpSpoil(t)*o;...

        z;...
        z];
%     tt=tt+1;
end

    Grad(4).Time=nan;
    Grad(4).Amp=0;

%%

[ ~, SeqOut, data, ~] = set_sequence(HW, Seq, AQ, TX, Grad);
SeqOut.AQSlice=AQSlice;
%%
data.k_OS_2D=squeeze(data.data(:,1,find(~isnan(AQ.Start),1)-1+(1:AQSlice.nPhase*AQSlice.PhaseOS)));

data.fft1_OS_2D=squeeze(data.fft1_data(:,1,find(~isnan(AQ.Start),1)-1+(1:AQSlice.nPhase*AQSlice.PhaseOS))); % in Readrichugn FFT Transformierte und Amplitudenkorrigierte Daten des ersten AQ-Fensters pro TR
data.image_OS_2D=fftshift(ifft(fftshift(data.fft1_OS_2D,2),[],2),2); % FFT Transformierte und Amplitudenkorrigierte Daten des ersten AQ-Fensters pro TR
% Ausgeschnitte Bilddaten (PhaseOS und ReadOS)
data.image_2D=data.image_OS_2D(floor(size(data.image_OS_2D,1)/2)-floor(AQSlice.nRead/2)+1:floor(size(data.image_OS_2D,1)/2)+ceil(AQSlice.nRead/2),floor(AQSlice.nPhase*AQSlice.PhaseOS/2)+(-AQSlice.nPhase/2+1:(AQSlice.nPhase/2)));

if Seq.plot_k_image==1;plot_k_image(data.k_OS_2D,data.image_2D,SeqOut.AQSlice,SeqOut.plotPhase); end
SeqOut.data.image_2D=data.image_2D;
