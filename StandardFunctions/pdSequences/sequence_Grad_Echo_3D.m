function  [data,SeqOut,mySave, HW]= sequence_Grad_Echo_3D( HW, mySave, Seq, SliceSelect)
% This function is deprecated and should no longer be used.
% Use sequence_Flash instead.

if isemptyfield(Seq, 'get_B0')
  Seq.get_B0 = 0;
end
if isemptyfield(Seq, 'ppmLim')
  Seq.ppmLim = [-50,50];
end

if isemptyfield(Seq, 'RandomTXRXPhase')
  Seq.RandomTXRXPhase = 1;
end
if isemptyfield(Seq, 'AQPhaseOffset')
  Seq.AQPhaseOffset = [];
end
if isemptyfield(Seq, 'TXPhaseOffset')
  Seq.TXPhaseOffset = [];
end
if isemptyfield(Seq, 'HzPixMin')
  Seq.HzPixMin = 0;
end
if isemptyfield(Seq, 'Correct_Phase')
  Seq.Correct_Phase = 0;
end
if isemptyfield(Seq, 'Correct_tAQmax')
  Seq.Correct_tAQmax = 1e-3;
end
if isemptyfield(Seq, 'Correct_nSamples')
  Seq.Correct_nSamples = 20;
end
if isemptyfield(Seq, 'SteadyState_PreShots')
  Seq.SteadyState_PreShots = 10;
end
if isemptyfield(Seq, 'SteadyState_PostShots')
  Seq.SteadyState_PostShots = 10;
end
if isemptyfield(Seq, 'Correct_Plot')
  Seq.Correct_Plot = 0;
end
if isemptyfield(Seq, 'Correct_PlotFrequency')
  Seq.Correct_PlotFrequency = 0;
end
if isemptyfield(Seq, 'WaitCalibrateTime')
  Seq.WaitCalibrateTime = 1;
end
if isemptyfield(Seq, 'CalibrateThanPause')
  Seq.CalibrateThanPause = 0;
end
if isemptyfield(Seq, 'Calibrate')
  Seq.Calibrate = 0;
end
if isemptyfield(Seq, 'averageBreak')
  Seq.averageBreak = Seq.tRep(1);
end
if isemptyfield(Seq, 'average')
  Seq.average = 1;
end
if isemptyfield(Seq, 'Save_loop_kos_3D')
  Seq.Save_loop_kos_3D = 0;
end
if isemptyfield(Seq, 'Save_Global')
  Seq.Save_Global = 0;
end
if isemptyfield(Seq, 'Correct_debug')
  Seq.Correct_debug = 0;
end
if isemptyfield(Seq, 'plotPhase')
  Seq.plotPhase = 0;
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
if isemptyfield(SliceSelect, 'phi')
  SliceSelect.phi = 0;
end
if isemptyfield(SliceSelect, 'theta')
  SliceSelect.theta = 0;
end




AQSlice.nRead=SliceSelect.nRead;
AQSlice.nPhase=SliceSelect.nPhase;
AQSlice.nPhase3D=SliceSelect.nPhase3D;

AQSlice.sizeRead=SliceSelect.sizeRead;
AQSlice.sizePhase=SliceSelect.sizePhase;
AQSlice.sizePhase3D=SliceSelect.sizePhase3D;
AQSlice.thickness=SliceSelect.thickness;

AQSlice.ReadOS=Seq.ReadOS;
AQSlice.PhaseOS=Seq.PhaseOS;
AQSlice.PhaseOS3D=Seq.PhaseOS3D;

AQSlice.alfa=SliceSelect.alfa; %um x Achse
AQSlice.phi=SliceSelect.phi;  %um y Achse
AQSlice.theta=SliceSelect.theta;%um z Achse
AQSlice.CenterRot=SliceSelect.CenterRot;

% if AQSlice.nPhase3D == 1   %%% 2D
%     AQSlice.alfa=0.0*pi; %um x Achse
%     AQSlice.phi=0.5*pi;  %um y Achse
%     AQSlice.theta=0.*pi;%um z Achse
%     AQSlice.CenterRot=[0.000, 0.000, 0.00];
% else                    %%% 3D
%     AQSlice.alfa=0.0*pi; %um x Achse
%     AQSlice.phi=0.0*pi;  %um y Achse
%     AQSlice.theta=0.*pi;%um z Achse
%     AQSlice.CenterRot=[0.000, 0.000, 0.00];
% end


% HW.Grad.TimeDelay=HW.Grad.MMRTTimeOffset+23e-6;

%%%%% Calibration  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Seq.Calibrate>=1
    mySave=[];
    cal=1:1+Seq.Calibrate;
else
    cal=0;
end
% [HW,mySave]  = Find_Frequency_Sweep( HW, mySave);
Seq.StartSequenceTime=[];
HW.ReInit=1;
AQSlice.Flip=acos(exp(-Seq.tRep(1)/Seq.T1));  %in RAD


if AQSlice.nPhase3D>1
    Seq.tTxSlicemax=HW.TX.Amp2FlipPiIn1Sec/HW.TX.AmpDef*AQSlice.Flip/pi;
    Seq.TxSliceBW=1/Seq.tTxSlicemax;




else
%     Seq.TxSliceBWmax=min(HW.Grad.MaxAmp(1:3))*HW.GammaDef/2/pi*AQSlice.thickness;  % smaller BW -> pulse duration (longer than 1/BW*Shaping) limited by maximum gradient amplitude
%     HW.TX.Amp2FlipPiIn1Sec/HW.TX.AmpDef*AQSlice.Flip/pi % pulse duration limited by rf amplifier -> use longer pulse
    Seq.tTxSlicemax=max(HW.TX.Amp2FlipPiIn1Sec/HW.TX.AmpDef*AQSlice.Flip/pi,  1/(min(HW.Grad.MaxAmp(1:3))*HW.GammaDef/2/pi*AQSlice.thickness))*2;
    Seq.TxSliceBW=1/Seq.tTxSlicemax*2;

%     if isempty(whos('global','SliceSelect'))
%         SliceSelect.alfa=0.0*pi; %um x Achse
%         SliceSelect.phi=0.0*pi;  %um y Achse
%         SliceSelect.theta=0.*pi;%um z Achse
%         SliceSelect.CenterRot=[0.00, -0.00, 0.00];
%         SliceSelect.nRead=16*4;
%         SliceSelect.nPhase=16*4;
%         SliceSelect.sizeRead=0.010;
%         SliceSelect.sizePhase=0.010;
%
%     %-----------------------------------------------------------------------
%
%         SliceSelect.Center=SliceSelect.CenterRot([2,3,1]).*[1,-1,1];
%         SliceSelect.CenterRot=SliceSelect.Center([3,1,2]).*[1,1,-1];
%         [SliceSelect.normV(1),SliceSelect.normV(2),SliceSelect.normV(3)]=sph2cart(SliceSelect.theta,-SliceSelect.phi,1);
%         SliceSelect.R=SliceSelect.CenterRot*SliceSelect.normV.';
%         SliceSelect.Rauf=SliceSelect.R*SliceSelect.normV;
%         SliceSelect.RaufCenter=SliceSelect.CenterRot-SliceSelect.Rauf;
%         SliceSelect.CenterRauf=SliceSelect.Rauf-SliceSelect.CenterRot;
%         SliceSelect.CenterRaufImage=tpaRotate( SliceSelect.CenterRauf',-SliceSelect.alfa ,-SliceSelect.phi, -SliceSelect.theta)';
%         SliceSelect.CenterRaufImage=SliceSelect.CenterRaufImage([3,2,1]).*[1,-1,1];
%         SliceSelect.Center2OriginImage=[SliceSelect.CenterRaufImage([1,2]),SliceSelect.R];
%     else
%         global SliceSelect;
%     end
%
%     SliceSelect.thickness=AQSlice.thickness;
%     SliceSelect.ReadOS=AQSlice.ReadOS;
%     SliceSelect.PhaseOS=AQSlice.PhaseOS;
%     SliceSelect.sizePhase3D=AQSlice.sizePhase3D;
%     SliceSelect.PhaseOS3D=AQSlice.PhaseOS3D;
%     SliceSelect.nPhase3D=AQSlice.nPhase3D;
%
%     AQSlice=SliceSelect;

end

    AQSlice.Center=AQSlice.CenterRot([2,3,1]).*[1,-1,1];

    [AQSlice.normV(1),AQSlice.normV(2),AQSlice.normV(3)]=sph2cart(AQSlice.theta,-AQSlice.phi,1);
    AQSlice.R=AQSlice.Center*AQSlice.normV.';
    AQSlice.Rauf=AQSlice.R*AQSlice.normV;
    AQSlice.RaufCenter=AQSlice.Center-AQSlice.Rauf;
    AQSlice.CenterRauf=AQSlice.Rauf-AQSlice.Center;
    AQSlice.CenterRaufImage=tpaRotate( AQSlice.CenterRauf',-AQSlice.alfa ,-AQSlice.phi, -AQSlice.theta)';
    AQSlice.CenterRaufImage=AQSlice.CenterRaufImage([3,2,1]).*[1,-1,1];
    AQSlice.Center2OriginImage=[AQSlice.CenterRaufImage([1,2]),AQSlice.R].*[1,1,-1];


for cal=cal

  switch cal
    case 0
    case 1
      AQSlice.sizePhase=AQSlice.sizePhase+1000;
      AQSlice.sizePhase3D=AQSlice.sizePhase3D+1000;
      AQSlice.nPhase3Dtemp=AQSlice.nPhase3D;
      AQSlice.thickness=AQSlice.thickness+1000;
      AQSlice.nPhase3D=1;
      AQSlice.nPhasetemp=AQSlice.nPhase;
      AQSlice.nPhase=4;
      Seq.averagetemp=Seq.average;
      Seq.average=1;
      mySave.sequence_Echo_3D.GradSystemTimeDelayError=0;
      AmpCorrReadStart=1;
      % RelStepAmpPhaseOffset=0;

    case Seq.Calibrate+1
      AQSlice.sizePhase=AQSlice.sizePhase-1000;
      AQSlice.sizePhase3D=AQSlice.sizePhase3D-1000;
      AQSlice.nPhase3D=AQSlice.nPhase3Dtemp;
      AQSlice.thickness=AQSlice.thickness-1000;
      AQSlice.nPhase=AQSlice.nPhasetemp;
      Seq.Calibrate=0;
      Seq.average=Seq.averagetemp;
      if Seq.CalibrateThanPause == 1
        disp('Exchange sample and hit any key to continue measurement.');
        pause;
      end

  end

  % if ~exist('AmpCorrReadStart','var'), AmpCorrReadStart=0;end
  % if isnan(AmpCorrReadStart), AmpCorrReadStart=0;end
  % AmpCorrReadStart=AmpCorrReadStart(1);
  % if ~exist('RelStepAmpPhaseOffset','var'), RelStepAmpPhaseOffset=0;end
  % RelStepAmpPhaseOffset=RelStepAmpPhaseOffset(1);

  if ~isfield(mySave, 'sequence_Echo_3D')
    mySave.sequence_Echo_3D = [];
  end
  if ~isfield(mySave.sequence_Echo_3D, 'GradSystemTimeDelayError')
    mySave.sequence_Echo_3D.GradSystemTimeDelayError = 0;
  end
  if isnan(mySave.sequence_Echo_3D.GradSystemTimeDelayError)
    mySave.sequence_Echo_3D.GradSystemTimeDelayError = 0;
  end
  mySave.sequence_Echo_3D.GradSystemTimeDelayError=mySave.sequence_Echo_3D.GradSystemTimeDelayError(1);
  HW.Grad.TimeDelay=HW.Grad.SystemTimeDelay+HW.Grad.MMRTTimeOffset+mySave.sequence_Echo_3D.GradSystemTimeDelayError;

%---------------------------------------------------------------




%%%%% Sequenztiming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if AQSlice.nPhase3D>1
    Seq.tSlice=Seq.tTxSlicemax+2*HW.Grad.tEC;
else
    Seq.tSlice=Seq.tTxSlicemax+2*HW.Grad.tRamp+2*HW.Grad.tEC;
end
Seq.onestRep=ones(1,AQSlice.nPhase*AQSlice.PhaseOS*AQSlice.nPhase3D*AQSlice.PhaseOS3D+Seq.SteadyState_PreShots+Seq.SteadyState_PostShots);
Seq.tRep=Seq.tRep(1)*Seq.onestRep;
Seq.TxSliceBW=1/Seq.tTxSlicemax;
Seq.tInvertmax=HW.TX.Amp2FlipPiIn1Sec/HW.TX.AmpDef;


Seq.CLTime=100e-6;                                                                   % Kommando lade Zeit sollte > 0.05 ms sein
% Seq.tInvert=Seq.tInvertmax+1*HW.Grad.tEC;                                            % duration of gradient at inversion pulse
% Seq.tSlice=Seq.tTxSlicemax+2*HW.Grad.tRamp+2*HW.Grad.tEC;                            % Zeit des Slice Gradienten
% Seq.tAktivStart= Seq.tInvert;                                                       % Zeit die vorne vom TR weggeht Seq.tInvert
% Seq.tAktivEnd=  Seq.tEcho-(Seq.CLTime+1e-6+max(HW.BlankOffset,HW.Grad.TimeDelay)+12e-6);% Zeit die hinten vom TR weggeht
% Seq.tEchoInTR= Seq.tEcho/2+Seq.tInvertmax/2;                                        % Echozeitpunkt in TR
% Seq.tSymAktive = min(Seq.tEchoInTR-Seq.tAktivStart,Seq.tAktivEnd-Seq.tEchoInTR);    % Freie Zeit um Echozeitpunkt in TR
% Seq.tSymAktiveStart=Seq.tEchoInTR-Seq.tSymAktive;                                   % start time for gradients in TR
% Seq.tSymAktiveEnd=Seq.tEchoInTR+Seq.tSymAktive;                                     % end time for gradienten in TR
% Seq.tAQmax=(Seq.tSymAktive*2-2*3*HW.Grad.tRamp-2*HW.Grad.tEC)/2;                      % maximum duration of acquisition
% Seq.tGrad=(Seq.tSymAktive*2-2*HW.Grad.tRamp-2*HW.Grad.tEC-Seq.tAQmax)/2;              % duration of phase encoding gradients and read dephase
% Seq.tInvert=Seq.tInvertmax+2*HW.Grad.tEC;
% Seq.tSlice=Seq.tTxSlicemax+2*HW.Grad.tRamp+2*HW.Grad.tEC;
% Seq.tGrad=Seq.tEcho/2-Seq.tSlice/2-Seq.tInvert/2-HW.Grad.tEC;
% Seq.tAQmax=Seq.tEcho-Seq.tInvert-2*3*HW.Grad.tRamp-Seq.Spoil*Seq.tEcho/2*2-2*2*HW.Grad.tEC;
Seq.tInvert=0;%Seq.tInvertmax+2*HW.Grad.tEC;
Seq.tGrad=Seq.tEcho/2-Seq.tSlice/4;
if Seq.Correct_Phase
  if AQSlice.nPhase3D <= 1
    error('AQSlice.nPhase3D and Seq.Correct_Phase');
  end
  Seq.Correct_fSample=HW.RX.fSample/floor((Seq.Correct_tAQmax/(Seq.Correct_nSamples)/(1/HW.RX.fSample))); %Decimierung soll gerade sein, da 2fach oversampling
  Seq.Correct_tAQ=Seq.Correct_nSamples/Seq.Correct_fSample;
  Seq.Correct_tGradShift=(Seq.Correct_tAQ+get_DeadTimeTX2RX(HW,Seq.Correct_fSample));
  Seq.tAQmaxt=Seq.tEcho-2*HW.Grad.tRamp-2*HW.Grad.tEC-Seq.tSlice/2-Seq.Correct_tGradShift*2;
  Seq.Correct_imageAQWindow=2;
else
  Seq.Correct_fSample=[];
  Seq.Correct_tAQ=0;
  Seq.Correct_nSamples=[];
  Seq.tAQmaxt=Seq.tEcho-2*HW.Grad.tRamp-2*HW.Grad.tEC-Seq.tSlice/2;
  Seq.Correct_tGradShift=0;
  Seq.Correct_imageAQWindow=1;
end
if Seq.HzPixMin ~= 0 && 1/Seq.HzPixMin > Seq.tAQmaxt
  warning(['Seq.HzPixMin is too low. ' num2str(ceil(1/Seq.tAQmaxt)) ' Hz is used']);
end

Seq.tAQmax=min(1/Seq.HzPixMin, Seq.tAQmaxt);



[Seq]=createSeq(Seq,AQSlice,HW);

if  Seq.RandomTXRXPhase
  if isempty(Seq.AQPhaseOffset)
    Seq.AQPhaseOffset=rand(1,length(Seq.tRep))*360;
  end
  if isempty(Seq.TXPhaseOffset)
    Seq.TXPhaseOffset=Seq.AQPhaseOffset;
  end
else
  if isempty(Seq.AQPhaseOffset)
    Seq.AQPhaseOffset=zeros(1,length(Seq.tRep));
  end
  if isempty(Seq.TXPhaseOffset)
    Seq.TXPhaseOffset=Seq.AQPhaseOffset;
  end
end


% AQ
if Seq.Correct_Phase
  AQ.Start=[Seq.tSlice+get_DeadTimeTX2RX(HW,Seq.Correct_fSample);...
              (Seq.tSlice/2+Seq.tEcho-Seq.tAQ/2)];
  AQ.fSample=[Seq.Correct_fSample;...
              Seq.fSample*AQSlice.ReadOS];            %125e6 /  1 und 4 bis 8192
  AQ.nSamples=[Seq.Correct_nSamples;...
              AQSlice.nRead*AQSlice.ReadOS];
  AQ.Frequency=Seq.fAQCenter;
  % AQ.Phase=cumsum(((Seq.PhaseShift/AQSlice.PhaseOS))*ones(1,AQSlice.nPhase*AQSlice.PhaseOS))-mod(Seq.fAQCenter-HW.fLarmor,1./Seq.tEcho).*Seq.tEcho.*360;
  AQ.Phase=180+Seq.AQPhaseOffset;
  % AQ.Phase=180;
  AQ.ResetPhases = [1, zeros(1, size(AQ.Phase, 2)-1)];
else
  AQ.Start=(Seq.tSlice/2+Seq.tEcho-Seq.tAQ/2);
  AQ.fSample=Seq.fSample*AQSlice.ReadOS;            %125e6 /  1 und 4 bis 8192
  AQ.nSamples=AQSlice.nRead*AQSlice.ReadOS;
  AQ.Frequency=Seq.fAQCenter;
  % AQ.Phase=cumsum(((Seq.PhaseShift/AQSlice.PhaseOS))*ones(1,AQSlice.nPhase*AQSlice.PhaseOS))-mod(Seq.fAQCenter-HW.fLarmor,1./Seq.tEcho).*Seq.tEcho.*360;
  AQ.Phase=180+Seq.AQPhaseOffset;
  % AQ.Phase=180;
  AQ.ResetPhases = [1, zeros(1, size(AQ.Phase, 2)-1)];
end


% HF TX ------------------------------------------------------------
if AQSlice.nPhase3D>1
    pulseData1 = Pulse_Rect(HW, Seq.tSlice/2, Seq.TxSliceBW, AQSlice.Flip, 20, Seq.tTxSlicemax, Seq.fTxSlice, 90);
else
    pulseData1 = Pulse_RaisedCos(HW, Seq.tSlice/2, Seq.TxSliceBW, AQSlice.Flip, 20, Seq.tTxSlicemax, Seq.fTxSlice, 90);
end
%     pulseData2 = Pulse_Rect(HW, Seq.tSlice/2+Seq.tEcho/2, HW.TX.AmpDef/HW.TX.Amp2FlipPiIn1Sec, pi, 20, Seq.tInvertmax, HW.fLarmor, 0);

TX.Duration=[pulseData1.Duration];
TX.Start=[pulseData1.Start];
TX.Amplitude=[pulseData1.Amplitude]; %  0 - 1
TX.Frequency=[pulseData1.Frequency];
TX.Phase=[pulseData1.Phase]+ones(size(pulseData1.Phase))*Seq.TXPhaseOffset;
% TX.Phase=[pulseData1.Phase];

%%%%%%%%%%%%%%%%%%%%%% timing of gradient for slice selection
if AQSlice.nPhase3D>1

    GradTime=cumsum([...
        0;...            % slice start
        0;...
        Seq.tTxSlicemax+2*HW.Grad.tEC+Seq.Correct_tGradShift;...
        0;... % slice stop
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
else
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
end

%%%%%%%%%%%%%%%%%%%%%% Gradientenzeitpunkte bei den Echos

 GradTime(end)=Seq.tRep(1);             % letzten Zeitpunkt auf das Ende von tRep setzen
% vorletzter Zeitpunkt muss im Aktiven bereich von ersten tRep liegen
if Seq.tRep(1)-(Seq.CLTime+1e-6+max(HW.TX.BlankOffset,HW.Grad.TimeDelay)+12e-6)<GradTime(end-1)
    error(['tRep too short for selected tEcho. tRep min is ' num2str((GradTime(end-1)+(Seq.CLTime+1e-6+max(HW.BlankOffset,HW.Grad.TimeDelay)+12e-6))*1000) ' ms.'])
end
% GradTime(end)=Seq.tRep(2); % letzten Zeitpunkt auf das Ende von tRep setzen

z=zeros(1,AQSlice.nPhase*AQSlice.PhaseOS);  % zeros for padding
o=ones(1,AQSlice.nPhase*AQSlice.PhaseOS);  % ones for padding

%%%%%% Gradienten Kallibrierung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% ende Gradienten Kallibrierung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fft_data3D=zeros(AQSlice.nRead,AQSlice.nPhase,AQSlice.nPhase3D*AQSlice.PhaseOS3D); % Datenspeicher vorbereiten
for t=1:3
    Grad(t).Amp=[];
    Grad(t).Time=[];
end
for tt=1:size(Seq.AmpPhase3D,2)  % loop over 3D phase encoding

%%%%% Zuordung der Gradienten Zeitpunkte und Amplituden
    for t=1:3; % [1, 2, 3]= x y z
        Grad(t).Time=[Grad(t).Time,GradTime*o];
        Grad(t).Amp=[Grad(t).Amp,...
            [...
        z;...
%         Seq.AmpSlice(tt)*o;...
%         Seq.AmpSlice(tt)*o;...
%         z;...
%         Seq.AmpPhase(tt,:)+Seq.AmpReadD(tt)+Seq.AmpSpoilD(tt)+Seq.AmpSliceD(tt);...
%         Seq.AmpPhase(tt,:)+Seq.AmpReadD(tt)+Seq.AmpSpoilD(tt)+Seq.AmpSliceD(tt);...

        Seq.AmpSlice(t)*o;...
        Seq.AmpSlice(t)*o;...
        z;...
%         -Seq.AmpPhase(t,:)+Seq.AmpReadD(t)-Seq.AmpSliceD(t);...
%         -Seq.AmpPhase(t,:)+Seq.AmpReadD(t)-Seq.AmpSliceD(t);...
            -Seq.AmpPhase(t,:)+Seq.AmpReadD(t)+Seq.AmpSliceD(t)+Seq.AmpPhase3D(t,tt);...
            -Seq.AmpPhase(t,:)+Seq.AmpReadD(t)+Seq.AmpSliceD(t)+Seq.AmpPhase3D(t,tt);...


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
        +Seq.AmpSpoil(t)*o;...
        +Seq.AmpSpoil(t)*o;...

        z;...
        z]];


    end

end

for t=1:3
    Grad(t).Amp=[Grad(t).Amp(:,1)*ones(1,Seq.SteadyState_PreShots),Grad(t).Amp,Grad(t).Amp(:,end)*ones(1,Seq.SteadyState_PostShots)];
    Grad(t).Time=[Grad(t).Time(:,1)*ones(1,Seq.SteadyState_PreShots),Grad(t).Time,Grad(t).Time(:,end)*ones(1,Seq.SteadyState_PostShots)];
end

Grad(4).Time=nan;       % B0 nicht verwendet
Grad(4).Amp=0;          % B0 0

%%%%%%%% Startet die Messung %%%%%%%%%%%%
[ ~, SeqOut, data, ~ ] = set_sequence(HW, Seq, AQ, TX, Grad);
SeqOut.AQSlice=AQSlice;

if Seq.Calibrate==0
  % % Start next measurement delayed by Seq.WaitSeqenzeTime
  % Seq.StartSequenceTime=now*24*3600+Seq.WaitSeqenzeTime;
else
  Seq.StartSequenceTime=now*24*3600+Seq.WaitCalibrateTime;
end
if SeqOut.Correct_Phase
    SeqOut.Correct_mywin=cos(linspace(-0.5*pi,0.5*pi,(SeqOut.SteadyState_PreShots+SeqOut.SteadyState_PostShots+1))).'.^2;
    SeqOut.Correct_mywin=SeqOut.Correct_mywin/sum(SeqOut.Correct_mywin(:));
    SeqOut.Correct_foffsetall=conv2(squeeze(mean(diff(unwrap(angle(data.data(1:SeqOut.Correct_nSamples,1,:))))))*SeqOut.AQ.fSample(1)/2/pi,SeqOut.Correct_mywin,'same').';
    SeqOut.Correct_tfoffsetall=SeqOut.AQ.Start(1,:)+SeqOut.AQ.nSamples(1,:)/2*1/SeqOut.AQ.fSample(1,:)+cumsum([0,SeqOut.tRep(1:end-1)]);
    SeqOut.Correct_foffset=conv2(squeeze(mean(diff(unwrap(angle(data.data(1:SeqOut.Correct_nSamples,1,:))))))*SeqOut.AQ.fSample(1)/2/pi,SeqOut.Correct_mywin,'valid').';
    SeqOut.Correct_tfoffset=SeqOut.AQ.Start(1,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)+SeqOut.AQ.nSamples(1,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)/2*1/SeqOut.AQ.fSample(1,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)+cumsum(SeqOut.tRep(SeqOut.SteadyState_PreShots+0:end-SeqOut.SteadyState_PostShots-1))+sum(SeqOut.tRep(1:SeqOut.SteadyState_PreShots-1));
    SeqOut.Correct_tPulse=SeqOut.tSlice/2+cumsum(SeqOut.tRep(SeqOut.SteadyState_PreShots+0:end-SeqOut.SteadyState_PostShots-1))+sum(SeqOut.tRep(1:SeqOut.SteadyState_PreShots-1));
    SeqOut.Correct_foffsetPulse=interp1(SeqOut.Correct_tfoffsetall,SeqOut.Correct_foffsetall,SeqOut.Correct_tPulse);
    SeqOut.Correct_poffsetPulse=(SeqOut.Correct_foffsetPulse+SeqOut.Correct_foffset)./2.*2.*pi.*(SeqOut.Correct_tPulse-SeqOut.Correct_tfoffset); %% +-?
    SeqOut.Correct_tSamples=squeeze(data.time_all(1:SeqOut.AQ.nSamples(2),2,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots));
    SeqOut.Correct_foffsetSamples=zeros(size(SeqOut.Correct_tSamples));
    SeqOut.Correct_foffsetSamples(:)=interp1(SeqOut.Correct_tfoffsetall,SeqOut.Correct_foffsetall,SeqOut.Correct_tSamples(:));
    SeqOut.Correct_poffsetSamples=(ones(size(SeqOut.Correct_foffsetSamples,1),1)*SeqOut.Correct_foffset+SeqOut.Correct_foffsetSamples)./2.*2.*pi.*(SeqOut.Correct_tSamples-ones(size(SeqOut.Correct_foffsetSamples,1),1)*SeqOut.Correct_tfoffset);
    SeqOut.Correct_poffsetBouth=SeqOut.Correct_poffsetSamples-ones(size(SeqOut.Correct_foffsetSamples,1),1)*SeqOut.Correct_poffsetPulse;

    if Seq.Correct_Plot
        figure(29)
        subplot(4,1,1)
        %     plot(squeeze(mean(diff(unwrap(angle(data.data(1:SeqOut.Correct_nSamples,1,SeqOut.SteadyState_PreShots+1:end)))))))
            plot(squeeze((abs(((data.data(1:SeqOut.Correct_nSamples,1,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)))))))
        subplot(4,1,2)
            plot(squeeze(((unwrap(angle(data.data(1:SeqOut.Correct_nSamples,1,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)))))))
        subplot(4,1,3)
            plot(squeeze((diff(unwrap(angle(data.data(1:SeqOut.Correct_nSamples,1,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots))))))*SeqOut.AQ.fSample(1)/2/pi)
        subplot(4,1,4)
            plot(squeeze(mean(diff(unwrap(angle(data.data(1:SeqOut.Correct_nSamples,1,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots))))))*SeqOut.AQ.fSample(1)/2/pi)
            hold on
            plot(conv2(squeeze(mean(diff(unwrap(angle(data.data(1:SeqOut.Correct_nSamples,1,:))))))*SeqOut.AQ.fSample(1)/2/pi,SeqOut.Correct_mywin,'valid'))
            hold off
    else
        if Seq.Correct_PlotFrequency
            figure(29)
            clf('reset')
            plot(squeeze(mean(diff(unwrap(angle(data.data(1:SeqOut.Correct_nSamples,1,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots))))))*SeqOut.AQ.fSample(1)/2/pi)
            hold on
            plot(conv2(squeeze(mean(diff(unwrap(angle(data.data(1:SeqOut.Correct_nSamples,1,:))))))*SeqOut.AQ.fSample(1)/2/pi,SeqOut.Correct_mywin,'valid'),'r')
            hold off
        end
    end



    %%
    if Seq.Correct_Plot

        figure(30)
        subplot(4,1,1)
        %     plot(squeeze(mean(diff(unwrap(angle(data.data(1:SeqOut.Correct_nSamples,1,SeqOut.SteadyState_PreShots+1:end)))))))
            plot(squeeze((abs(((data.data(1:SeqOut.AQ.nSamples(2),2,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)))))))
        subplot(4,1,2)
            plot(squeeze(((unwrap(angle(data.data(1:SeqOut.AQ.nSamples(2),2,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)))))))

        subplot(4,1,3)
            plot(squeeze((diff(unwrap(angle(data.data(SeqOut.AQ.nSamples(2),2,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots))))))*SeqOut.AQ.fSample(2)/2/pi)
        subplot(4,1,4)
            plot(squeeze(mean(diff(unwrap(angle(data.data(1:SeqOut.AQ.nSamples(2),2,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots))))))*SeqOut.AQ.fSample(2)/2/pi)
            hold on
            plot(conv2(squeeze(mean(diff(unwrap(angle(data.data(1:SeqOut.AQ.nSamples(2),2,:))))))*SeqOut.AQ.fSample(2)/2/pi,SeqOut.Correct_mywin,'valid'))

            plot(     squeeze(mean(diff(unwrap(angle(squeeze(data.data(1:SeqOut.AQ.nSamples(2),2,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)).*exp(-1i*SeqOut.Correct_poffsetBouth))))))*SeqOut.AQ.fSample(2)/2/pi,'r')
            plot((SeqOut.SteadyState_PreShots+1):length(data.data(1:SeqOut.AQ.nSamples(2),2,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots))-SeqOut.SteadyState_PostShots,...
                conv2(squeeze(mean(diff(unwrap(angle(squeeze(data.data(1:SeqOut.AQ.nSamples(2),2,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)).*exp(-1i*SeqOut.Correct_poffsetBouth))))))*SeqOut.AQ.fSample(2)/2/pi,SeqOut.Correct_mywin.','valid'),...
                'r')
            plot(conv2(squeeze(mean(diff(unwrap(angle(data.data(1:SeqOut.Correct_nSamples,1,:))))))*SeqOut.AQ.fSample(1)/2/pi,SeqOut.Correct_mywin,'valid'),'g')
            hold off
    end

    data.kos_3D=reshape(squeeze(data.data(1:SeqOut.AQ.nSamples(SeqOut.Correct_imageAQWindow),SeqOut.Correct_imageAQWindow,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)).*exp(-1i*SeqOut.Correct_poffsetBouth),[AQSlice.nRead*AQSlice.ReadOS,AQSlice.nPhase*AQSlice.PhaseOS, AQSlice.nPhase3D*AQSlice.PhaseOS3D]);
    data.fft1os_3D=fftshift(ifft(fftshift(data.kos_3D)))...
                    .*reshape(squeeze(data.cic_corr(1:SeqOut.AQ.nSamples(SeqOut.Correct_imageAQWindow),SeqOut.Correct_imageAQWindow,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)),[AQSlice.nRead*AQSlice.ReadOS,AQSlice.nPhase*AQSlice.PhaseOS, AQSlice.nPhase3D*AQSlice.PhaseOS3D]);




    if Seq.Correct_debug
        data.kos_3D_wo=reshape(squeeze(data.data(:,SeqOut.Correct_imageAQWindow,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)),[AQSlice.nRead*AQSlice.ReadOS,AQSlice.nPhase*AQSlice.PhaseOS, AQSlice.nPhase3D*AQSlice.PhaseOS3D]); % Daten des ersten AQ-Fensters pro TR
        data.fft1os_3D_wo=reshape(squeeze(data.fft1_data(:,SeqOut.Correct_imageAQWindow,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)),[AQSlice.nRead*AQSlice.ReadOS,AQSlice.nPhase*AQSlice.PhaseOS, AQSlice.nPhase3D*AQSlice.PhaseOS3D]); % in Readrichugn FFT Transformierte und Amplitudenkorrigierte Daten des ersten AQ-Fensters pro TR
        data.fft12os_3D_wo=fftshift(ifft(fftshift(data.fft1os_3D_wo,2),[],2),2); % FFT Transformierte und Amplitudenkorrigierte Daten des ersten AQ-Fensters pro TR
        % Ausgeschnitte Bilddaten (PhaseOS und ReadOS)
        data.fft12os3_3D_wo=data.fft12os_3D_wo(   floor(AQSlice.nRead*AQSlice.ReadOS/2)+      (1-floor(AQSlice.nRead/2):ceil(AQSlice.nRead/2)),...
                                            floor(AQSlice.nPhase*AQSlice.PhaseOS/2)+    (1-floor(AQSlice.nPhase/2):ceil(AQSlice.nPhase/2)),...
                                            :);
        data.fftos3_3D_wo=fftshift(ifft(fftshift(data.fft12os3_3D_wo,3),[],3),3);
        data.image_3D_wo=data.fftos3_3D_wo(:,:,floor(AQSlice.nPhase3D*AQSlice.PhaseOS3D/2)+   (1-floor(AQSlice.nPhase3D/2):ceil(AQSlice.nPhase3D/2)));
    else
        SeqOut=rmfield(SeqOut,'Correct_mywin');
        SeqOut=rmfield(SeqOut,'Correct_foffsetall');
        SeqOut=rmfield(SeqOut,'Correct_tfoffsetall');
%         SeqOut=rmfield(SeqOut,'Correct_foffset');
%         SeqOut=rmfield(SeqOut,'Correct_tfoffset');
        SeqOut=rmfield(SeqOut,'Correct_tPulse');
        SeqOut=rmfield(SeqOut,'Correct_foffsetPulse');
        SeqOut=rmfield(SeqOut,'Correct_poffsetPulse');
        SeqOut=rmfield(SeqOut,'Correct_tSamples');
        SeqOut=rmfield(SeqOut,'Correct_foffsetSamples');
        SeqOut=rmfield(SeqOut,'Correct_poffsetSamples');
        SeqOut=rmfield(SeqOut,'Correct_poffsetBouth');
    end
else
    data.kos_3D=reshape(squeeze(data.data(:,SeqOut.Correct_imageAQWindow,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)),[AQSlice.nRead*AQSlice.ReadOS,AQSlice.nPhase*AQSlice.PhaseOS, AQSlice.nPhase3D*AQSlice.PhaseOS3D]); % Daten des ersten AQ-Fensters pro TR
    data.fft1os_3D=reshape(squeeze(data.fft1_data(:,SeqOut.Correct_imageAQWindow,SeqOut.SteadyState_PreShots+1:end-SeqOut.SteadyState_PostShots)),[AQSlice.nRead*AQSlice.ReadOS,AQSlice.nPhase*AQSlice.PhaseOS, AQSlice.nPhase3D*AQSlice.PhaseOS3D]); % in Readrichugn FFT Transformierte und Amplitudenkorrigierte Daten des ersten AQ-Fensters pro TR


end
    data.fft12os_3D=fftshift(ifft(fftshift(data.fft1os_3D,2),[],2),2); % FFT Transformierte und Amplitudenkorrigierte Daten des ersten AQ-Fensters pro TR
    % Ausgeschnitte Bilddaten (PhaseOS und ReadOS)
    data.fft12os3_3D=data.fft12os_3D(   floor(AQSlice.nRead*AQSlice.ReadOS/2)+      (1-floor(AQSlice.nRead/2):ceil(AQSlice.nRead/2)),...
                                        floor(AQSlice.nPhase*AQSlice.PhaseOS/2)+    (1-floor(AQSlice.nPhase/2):ceil(AQSlice.nPhase/2)),...
                                        :);
    data.fftos3_3D=fftshift(ifft(fftshift(data.fft12os3_3D,3),[],3),3);
    data.image_3D=data.fftos3_3D(:,:,floor(AQSlice.nPhase3D*AQSlice.PhaseOS3D/2)+   (1-floor(AQSlice.nPhase3D/2):ceil(AQSlice.nPhase3D/2)));

    if ~Seq.Correct_debug
      data=rmfield(data,'fft1os_3D');
      data=rmfield(data,'fft12os_3D');
      data=rmfield(data,'fft12os3_3D');
      data=rmfield(data,'fftos3_3D');
    end

    %%%%%%%%%%% Kalibrierung der Gradient Offset Zeit und der ersten Read dephase Amplitude
    if SeqOut.Calibrate > 0 && AQSlice.sizePhase>10
      [mySave.sequence_Echo_3D.GradSystemTimeDelayError, AmpCorrReadStart] = ...
        CorrRead(mySave.sequence_Echo_3D.GradSystemTimeDelayError, AmpCorrReadStart, squeeze(data.fft1os_3D), SeqOut, AQ);
    end
end
if AQSlice.nPhase3D>1
    if Seq.Save_Global
        clear global image3D
        global image3D;
        image3D.mySave=mySave;
        image3D.data=data.image_3D;
        image3D.AQSlice=AQSlice;
        image3D.HW=HW;
        %image3D.talker=talker;
        image3D.Seq=SeqOut;
        image3D.xV=image3D.AQSlice.CenterRot(1)+linspace(-image3D.AQSlice.sizeRead/2,image3D.AQSlice.sizeRead/2,image3D.AQSlice.nRead);
        image3D.yV=image3D.AQSlice.CenterRot(2)+linspace(-image3D.AQSlice.sizePhase/2,image3D.AQSlice.sizePhase/2,image3D.AQSlice.nPhase);
        image3D.zV=image3D.AQSlice.CenterRot(3)+linspace(-image3D.AQSlice.sizePhase3D/2,image3D.AQSlice.sizePhase3D/2,image3D.AQSlice.nPhase3D);
        [image3D.x,image3D.y,image3D.z] = ndgrid(image3D.xV,image3D.yV,image3D.zV);
    %     image3D.data=flipdim(image3D.data,3);
        for t=1:size(image3D.data,3)
           image3D.dataYXZ(:,:,t)=image3D.data(:,:,t).';
        end
        if Seq.useSliceomatic

            sliceomatic(abs(image3D.dataYXZ),image3D.xV*1000000,image3D.yV*1000000,image3D.zV*1000000)
            if Seq.plotPhase
                sliceomatic(angle(image3D.dataYXZ),image3D.xV*1000000,image3D.yV*1000000,image3D.zV*1000000)
            end

        end
    end
else
  if ~isfield(Seq, 'plot_k_image')
    Seq.plot_k_image = 1;
  end
  if Seq.plot_k_image == 1 && AQSlice.nPhase3D == 1
    plot_k_image(data.kos_3D, data.image_3D, AQSlice, Seq.plotPhase);
  end
end

end
