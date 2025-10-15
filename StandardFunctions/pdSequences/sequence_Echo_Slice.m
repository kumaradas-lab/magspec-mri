function [data,SeqOut,mySave, HW]= sequence_Echo_Slice( HW, mySave, Seq, SliceSelect)
% This function is deprecated and should no longer be used.
% Use sequence_SpinEcho instead.

%----------- Parameter ------------------------------------------------
if ~isfield(Seq,'WaitSeqenzeTime');Seq.WaitSeqenzeTime=[];end
if isempty(Seq.WaitSeqenzeTime);Seq.WaitSeqenzeTime=0;end
if Seq.WaitSeqenzeTime ~= 0; Seq.StartSequenceTime=now*24*3600+Seq.WaitSeqenzeTime;    end

if ~isfield(Seq,'InvertPulse');   Seq.InvertPulse=[]; end; if isempty(Seq.InvertPulse);                 Seq.InvertPulse    =   @Pulse_Rect;                    end
if ~isfield(Seq,'SlicePulse');   Seq.SlicePulse=[];  end; if isempty(Seq.SlicePulse);                    Seq.SlicePulse    =   @Pulse_Rect;                    end
if ~isfield(Seq,'tReadDelay');   Seq.tReadDelay=[];  end; if isempty(Seq.tReadDelay);                    Seq.tReadDelay    =   0;                               end
if ~isfield(Seq,'plot_k_image');   Seq.plot_k_image=[];  end; if isempty(Seq.plot_k_image);              Seq.plot_k_image    = 0;                               end
if ~isfield(Seq,'plotPhase');   Seq.plotPhase=[];  end; if isempty(Seq.plotPhase);              Seq.plotPhase    = 1;                               end
if ~isfield(Seq,'HzPixMin');   Seq.HzPixMin=[];  end; if isempty(Seq.HzPixMin);              Seq.HzPixMin    = 0;                               end
if ~isfield(Seq,'Find_Frequency_time');   Seq.Find_Frequency_time=[];  end; if isempty(Seq.Find_Frequency_time);              Seq.Find_Frequency_time    = 1000;                               end
if ~isfield(Seq,'SliceGshift');   Seq.SliceGshift=[];  end; if isempty(Seq.SliceGshift);              Seq.SliceGshift    = 0;                               end

if nargin<4; SliceSelect=[];end
if ~isfield(SliceSelect,'thickness');   SliceSelect.thickness=[];  end; if isempty(SliceSelect.thickness);              SliceSelect.thickness    = 1e12;                               end
if ~isfield(SliceSelect,'CenterRot');   SliceSelect.CenterRot=[];  end; if isempty(SliceSelect.CenterRot);              SliceSelect.CenterRot    =[0.0, 0.0, 0.0];                               end
if ~isfield(SliceSelect,'alfa');   SliceSelect.alfa=[];  end; if isempty(SliceSelect.alfa);              SliceSelect.alfa    =0;                               end
if ~isfield(SliceSelect,'phi');   SliceSelect.phi=[];  end; if isempty(SliceSelect.phi);              SliceSelect.phi    =0;                               end
if ~isfield(SliceSelect,'theta');   SliceSelect.theta=[];  end; if isempty(SliceSelect.theta);              SliceSelect.theta    =0;                               end



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

%  [HW,mySave]  = Find_Frequency( HW, mySave,Seq.Find_Frequency_time);

  if ~isfield(SliceSelect,'Flip');   SliceSelect.Flip=[];  end; if isempty(SliceSelect.Flip);              SliceSelect.Flip    = 90;       end %in Rad
SliceSelect.Flip=SliceSelect.Flip/180*pi;
% SliceSelect.Flip=pi/2;  %in RAD
% SliceSelect.thickness=AQSlice.thickness; %Achtung
SliceSelect.ReadOS=Seq.ReadOS;
SliceSelect.PhaseOS=Seq.PhaseOS;

AQSlice=SliceSelect;

Seq.tTxSlicemax=max(HW.TX.Amp2FlipPiIn1Sec/(HW.TX.AmpDef*0.9)*AQSlice.Flip/pi*Seq.SlicePulse(HW,'Amp'),  1/(min(HW.Grad.MaxAmp(1:3))*HW.GammaDef/2/pi*AQSlice.thickness)*Seq.SlicePulse(HW,'Time'));
useSpoil=double(Seq.Spoil~=0);

Seq.tRep=Seq.tRep*ones(1,SliceSelect.nPhase*SliceSelect.PhaseOS);
Seq.CLTime=50e-6;
Seq.tInvertmax=HW.TX.Amp2FlipPiIn1Sec/HW.TX.AmpDef*Seq.InvertPulse(HW,'Amp');
Seq.TxSliceBW=1/Seq.tTxSlicemax*Seq.SlicePulse(HW,'Time');
Seq.tInvert=Seq.tInvertmax;%+2*HW.Grad.tEC;
Seq.tSlice=Seq.tTxSlicemax+2*HW.Grad.tRamp+2*HW.Grad.tEC;
Seq.tGrad=Seq.tEcho/2-Seq.tSlice/2-Seq.tInvert/2-HW.Grad.tEC;
Seq.tAQmaxt=Seq.tEcho-Seq.tInvert-4*HW.Grad.tRamp*useSpoil-2*HW.Grad.tRamp-Seq.Spoil*Seq.tEcho/2*2*useSpoil-2*2*HW.Grad.tEC;
if Seq.HzPixMin ~= 0;
    if 1/Seq.HzPixMin>Seq.tAQmaxt;
        warning(['Seq.HzPixMin is too low. ' num2str(ceil(1/Seq.tAQmaxt)) ' Hz is used'])
    end
end
Seq.tAQmax=min(1/Seq.HzPixMin, Seq.tAQmaxt);
Seq.tReadGradEmpty=Seq.tAQmaxt-Seq.tAQmax;


Seq=createSeq(Seq,AQSlice,HW);


% AQ
AQ.Start=(Seq.tSlice/2+Seq.tEcho-Seq.tAQ/2-Seq.fftReadoutShift/Seq.fSample/AQSlice.ReadOS/2)*ones(1,AQSlice.nPhase*AQSlice.PhaseOS)+Seq.tReadDelay;
AQ.fSample=Seq.fSample*AQSlice.ReadOS;             %125e6 /  1 und 4 bis 8192
AQ.nSamples=AQSlice.nRead*AQSlice.ReadOS;
AQ.Frequency=Seq.fAQCenter;
AQ.Phase=cumsum(((Seq.PhaseShift/AQSlice.PhaseOS))*ones(1,AQSlice.nPhase*AQSlice.PhaseOS))-mod(Seq.fAQCenter-HW.fLarmor,1./Seq.tEcho).*Seq.tEcho.*360;
AQ.Dur=AQ.nSamples/AQ.fSample;
AQ.ResetPhases=1;
% AQ.Gain=HW.RX.GainDef; % 0.0032 - 1

% HF TX
% TX.BlankOffset=20e-6; % Blank vor HF
pulseData1 = Seq.SlicePulse(HW,     Seq.tSlice/2,               Seq.TxSliceBW,      AQSlice.Flip,   50,  Seq.tTxSlicemax,   Seq.fTxSlice,   -90);
pulseData2 = Seq.InvertPulse(HW,    Seq.tSlice/2+Seq.tEcho/2,   1/Seq.tInvertmax,   pi,             50,  Seq.tInvertmax,    HW.fLarmor,     0);

TX.Duration=[pulseData1.Duration;             pulseData2.Duration]*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
TX.Start=[pulseData1.Start;         pulseData2.Start]*ones(1,AQSlice.nPhase*AQSlice.PhaseOS); %!!!!!
TX.Amplitude=[pulseData1.Amplitude;         pulseData2.Amplitude]*ones(1,AQSlice.nPhase*AQSlice.PhaseOS); %  0 - 1
TX.Frequency=[pulseData1.Frequency; pulseData2.Frequency]*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
TX.Phase=[pulseData1.Phase*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);         pulseData2.Phase*ones(1,AQSlice.nPhase*AQSlice.PhaseOS)+ones(size(pulseData2.Phase))*cumsum(ones(1,AQSlice.nPhase*AQSlice.PhaseOS))*180];

% TX.Duration=[pulseData1.Duration;Seq.tInvertmax];
% TX.Start=[pulseData1.Start;Seq.tSlice/2+Seq.tEcho/2-Seq.tInvertmax/2]; %!!!!!
% TX.Amplitude=[pulseData1.Amplitude; Seq.AmpTxInvert]; %  0 - 1
% TX.Frequency=[pulseData1.Frequency;HW.fLarmor];
% TX.Phase=[pulseData1.Phase;0];

% if size(TX.Start,2)<size(AQ.Start,2)
%     TX.Duration=TX.Duration*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
%     TX.Start=TX.Start*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
%     TX.Amplitude=TX.Amplitude*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
%     TX.Frequency=TX.Frequency*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
%     TX.Phase=TX.Phase*ones(1,AQSlice.nPhase*AQSlice.PhaseOS);
% end

GradTime=cumsum([...
    0-Seq.SliceGshift;...            % slice start
    HW.Grad.tRamp;...
    Seq.tTxSlicemax+2*HW.Grad.tEC;...
    HW.Grad.tRamp+Seq.SliceGshift;... % slice stop
    HW.Grad.tRamp;... % Grads start
    Seq.tGrad-2*HW.Grad.tRamp;...
    HW.Grad.tRamp;... % Grads stop
    Seq.tInvert+2*HW.Grad.tEC+Seq.tReadDelay;...  % invert
    HW.Grad.tRamp*useSpoil;... % Spoil start
    Seq.Spoil*Seq.tEcho/2*useSpoil;...
    HW.Grad.tRamp*useSpoil;... % Spoil stop
    Seq.tReadGradEmpty/2
    HW.Grad.tRamp;... % GradRead start
    Seq.tAQmax+2*HW.Grad.tEC;...
    HW.Grad.tRamp;... % GradRead ende
    Seq.tReadGradEmpty/2
    HW.Grad.tRamp;...
    Seq.tEcho;...
    HW.Grad.tRamp;
    6e-6]);
% GradTime(end)=Seq.tRep(1); %!!!!!!!!!!!!!!!

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
        Seq.AmpPhase(t,:)-Seq.AmpReadD(t)-Seq.AmpSpoilD(t)*useSpoil+Seq.AmpSliceD(t);...
        Seq.AmpPhase(t,:)-Seq.AmpReadD(t)-Seq.AmpSpoilD(t)*useSpoil+Seq.AmpSliceD(t);...


        z;...
        z;...
        Seq.AmpSpoil(t)*o*useSpoil;...
        Seq.AmpSpoil(t)*o*useSpoil;...
        z;...
        z;...
        Seq.AmpRead(t)*o;...
        Seq.AmpRead(t)*o;...
        z;...
        z;...
        Seq.AmpSpoilEnd(t)*o;...
        Seq.AmpSpoilEnd(t)*o;...
        z;...
        z];
%     tt=tt+1;
end

    Grad(4).Time=nan;
    Grad(4).Amp=0;

%%
Seq.AQSlice=AQSlice;
[ ~, SeqOut, data, ~] = set_sequence(HW, Seq, AQ, TX, Grad);
%%
data.k_OS_2D=squeeze(data.data(:,1,find(~isnan(AQ.Start),1)-1+(1:AQSlice.nPhase*AQSlice.PhaseOS)));

data.fft1_OS_2D=squeeze(data.fft1_data(:,1,find(~isnan(AQ.Start),1)-1+(1:AQSlice.nPhase*AQSlice.PhaseOS))); % in Readrichugn FFT Transformierte und Amplitudenkorrigierte Daten des ersten AQ-Fensters pro TR
%     plot_k_image(data.k_OS_2D,data.image_2D/SeqOut.AQSlice.VoxelVolume^2/SeqOut.AQSlice.thickness*HW.RX.AreaCoil,SeqOut.AQSlice,SeqOut.plotPhase);
% data.image_OS_2D=fftshift(ifft(fftshift(data.fft1_OS_2D,2),[],2),2)/SeqOut.AQSlice.VoxelVolume^2/SeqOut.AQSlice.thickness*HW.RX.AreaCoil; % FFT Transformierte und Amplitudenkorrigierte Daten des ersten AQ-Fensters pro TR
data.image_OS_2D=fftshift(ifft(fftshift(data.fft1_OS_2D,2),[],2),2)/HW.B0/SeqOut.AQSlice.VoxelVolume*HW.RX.AreaCoil; % FFT Transformierte und Amplitudenkorrigierte Daten des ersten AQ-Fensters pro TR
% Ausgeschnitte Bilddaten (PhaseOS und ReadOS)
data.image_2D=data.image_OS_2D(floor(size(data.image_OS_2D,1)/2)-floor(AQSlice.nRead/2)+1:floor(size(data.image_OS_2D,1)/2)+ceil(AQSlice.nRead/2),floor(AQSlice.nPhase*AQSlice.PhaseOS/2)+floor(-AQSlice.nPhase/2+1:(AQSlice.nPhase/2)));

if Seq.plot_k_image==1;plot_k_image(data.k_OS_2D,data.image_2D,AQSlice,Seq.plotPhase); end
SeqOut.data.image_2D=data.image_2D;
