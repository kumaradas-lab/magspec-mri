function [HW, mySave, Seq] = Find_ExcitationPulse(HW, mySave, Seq)
%% Find excitation pulse
% It applies an excitation pulse and acquires the FID after some perpendicular pulses.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% input check
if nargin < 1
  error('PD:Find_ExcitationPulse:NoInput', 'At least input argument HW is mandatory.');
end
if nargin < 2
  mySave = [];
end
if nargin < 3
  Seq = struct();
end

iDevice=1;

%% preparations

if isemptyfield(Seq, 'T1') 
  Seq.T1=HW.FindFrequencyPause;
end

if isemptyfield(Seq, 'tExcitation') 
  Seq.tExcitation=HW.tFlip90Def;
end
if isemptyfield(Seq, 'AmpExcitation') 
  Seq.AmpExcitation=HW.TX.Def.AmplitudeCalibrated(HW.TX(iDevice).ChannelDef);
end
if isemptyfield(Seq, 'PaUoutExcitation') 
  Seq.PaUoutExcitation=HW.TX.Def.PaUoutCalibrated(HW.TX(iDevice).ChannelDef);
end
if isemptyfield(Seq, 'FixedExcitationPulseLength') 
  Seq.FixedExcitationPulseLength=false;
end
if isemptyfield(Seq, 'setExcitation') 
  Seq.setExcitation='Amp';
end
if isemptyfield(Seq, 'SearchFlipDeg') 
  Seq.SearchFlipDeg=[90,-90,3*90,5*90];
end
if isemptyfield(Seq, 'SearchFlipDegnStep') 
  Seq.SearchFlipDegnStep=4;
end
if isemptyfield(Seq, 'SearchFlipDegStep') 
  Seq.SearchFlipDegStep=2;
end
if isemptyfield(Seq, 'TrackFrequency') 
  Seq.TrackFrequency=1;
end
if isemptyfield(Seq, 'AQPhaseOffsetDegSaveLevel') 
  Seq.AQPhaseOffsetDegSaveLevel=0.5; % Add a Line to LoadMySystem, if offset is lager. in DEG
end

if isemptyfield(Seq, 'tAQfOffset') 
  Seq.tAQfOffset=5e-3;                      % duration of frequency tracking (pay attention to chemical shift oscillation)
end
if isemptyfield(Seq, 'nPeriods') 
  Seq.nPeriods=1;                           % number of periods
end
if isemptyfield(Seq, 'SearchFlipDegStepsPerPeriod') 
  Seq.SearchFlipDegStepsPerPeriod=4;        % number of AQ windows per period
end
if isemptyfield(Seq, 'tAQFID') 
  Seq.tAQFID=50e-6;        % duration of AQ windows
end

if isemptyfield(Seq, 'Find_Frequency_interval') 
  Seq.Find_Frequency_interval=HW.FindFrequencySweep.maxTime;
end
if ~isnan( Seq.Find_Frequency_interval)
  [HW, mySave] = Find_Frequency_Sweep(HW, mySave, Seq.Find_Frequency_interval, [], 1);  % find magnet frequency
end


%% iteratively search optimum excitation pulse
  switch Seq.setExcitation
    case 'Time'
      tPulse90       = Seq.tExcitation;                % duration of TX pulse in s
      AmpPulse90     = 90/360 /tPulse90/(HW.GammaDef/2/pi);
    case 'Amp'
      AmpPulse90     = Seq.AmpExcitation;
      tPulse90       = 90/360 /AmpPulse90/(HW.GammaDef/2/pi);                % duration of TX pulse in s
    case 'PaUout'
      AmpPulse90     = Seq.PaUoutExcitation * HW.TX(iDevice).PaUout2Amplitude(HW.TX(iDevice).ChannelDef);
      tPulse90       = 90/360 /AmpPulse90/(HW.GammaDef/2/pi);                % duration of TX pulse in s
    otherwise
      tPulse90       = HW.tFlip90Def;                % duration of TX pulse in s
      AmpPulse90     = 90/360 /tPulse90/(HW.GammaDef/2/pi);
  end

if isemptyfield(Seq, 'tPulseRotate') 
  Seq.tPulseRotate=tPulse90*4;        % Amp of rotation pulse
end
if isemptyfield(Seq, 'AmpPulseRotate') 
  Seq.AmpPulseRotate=AmpPulse90;        % Amp of rotation pulse
end

B1=[];
tPulse=[];
FlipDeg=[];
AQPhaseDeg=[];
PaUout=[];
for SearchFlipDeg=Seq.SearchFlipDeg
  testFlipDeg=reshape([-round(Seq.SearchFlipDegnStep/2):-1,1:round(Seq.SearchFlipDegnStep/2)]*Seq.SearchFlipDegStep+SearchFlipDeg,[],1);
  FlipAngleDegPlot=nan(numel(testFlipDeg),1);
  tPulseCorrPlot=nan(numel(testFlipDeg),1);
  fLarmorAll=nan(numel(testFlipDeg),1);
  AQPhaseOffsetDegPlot=nan(numel(testFlipDeg),1);
  PaUoutCorrPlot=nan(numel(testFlipDeg),1);
  AQPhaseOffsetDeg=   0;
  fOffset=0;
  skipPreLopps=0;

  tExcitation = tPulse90*abs(SearchFlipDeg)/90;
  AmpExcitation = AmpPulse90;

  TXVoltage = AmpPulse90 / HW.TX(iDevice).PaUout2Amplitude(HW.TX(iDevice).ChannelDef);
  tPulseTest=tPulse90*abs(SearchFlipDeg)/90;
  tPulseCorr=tPulseTest;
  AmpPulseTest=AmpExcitation;
  AmpPulseCorr=AmpExcitation;
  if Seq.FixedExcitationPulseLength
    fprintf('Searching %.0f%s pulse, with a duration of : %.3f %ss @ about %.6f V\n',SearchFlipDeg, char(176), tPulseTest*1e6, char(181), TXVoltage);
  else
    fprintf('Searching %.0f%s pulse, with a duration of about: %.3f %ss @ %.6f V\n',SearchFlipDeg, char(176), tPulseTest*1e6, char(181), TXVoltage);
  end
  firstLoop=1;
  for tt=[-(1:numel(testFlipDeg)),1:numel(testFlipDeg)]
    t=abs(tt);
    if skipPreLopps && tt<0
      continue
    else
      pause(max(Seq.T1*2,0.5))
    end

    % Parameters used for timing calculations
    Seq.plotSeq   = 0;                           % plot sequence off
    % Seq.average   = 10;                           % plot sequence off

    if Seq.FixedExcitationPulseLength
      AmpPulseTest=AmpExcitation/abs(SearchFlipDeg) * testFlipDeg(t);
    else
      tPulseTest=round(tExcitation/abs(SearchFlipDeg) * testFlipDeg(t)*HW.TX.fSample)/HW.TX.fSample;
      testFlipDeg(t)=tPulseTest/tPulse90*90;
    end
    SearchFlipDegSteps=ceil(Seq.SearchFlipDegStepsPerPeriod*Seq.nPeriods);                             % number of AQ windows of each flip direction

    AQ.fSample    = [repmat(1e6,SearchFlipDegSteps*2+1,1);20e3];                         % sampling rate of AQ window in Hz
    deadTimeTX2RX = get_DeadTimeTX2RX(HW, AQ.fSample(1)); % dead time after TX pulse
    deadTimeRX2TX = get_DeadTimeRX2TX(HW, AQ.fSample(1)); % dead time after TX pulse

    tPulseGrid=(Seq.tPulseRotate/SearchFlipDegSteps*Seq.nPeriods+deadTimeTX2RX+deadTimeRX2TX)+Seq.tAQFID; % Grid of RF-Pules

    % RF transmission parameters
    TX.Duration   = [tPulseTest;repmat(Seq.tPulseRotate/SearchFlipDegSteps*Seq.nPeriods,SearchFlipDegSteps*2,1)];                      % duration of rf-pulse in s
    TX.Start      = -TX.Duration + tPulseGrid*(0:SearchFlipDegSteps*2).' ;                            % start time of rf-pulse in s
    TX.Frequency  = HW.fLarmor;                   % frequency of rf-pulse in Hz
    TX.Phase      = [180*double(AmpPulseTest<0); repmat(90,SearchFlipDegSteps,1);repmat(-90,SearchFlipDegSteps,1)];                            % phase of rf-pulse in degrees
    TX.Amplitude  = [abs(AmpPulseTest);repmat(Seq.AmpPulseRotate,SearchFlipDegSteps*2,1)];
    % TX.Channel      = 2;                            % phase of rf-pulse in degrees
    % Acquisition parameters
    AQ.Start      = [repmat(deadTimeTX2RX,SearchFlipDegSteps*2+1,1) + TX.Start + TX.Duration ; 0];          % acquisition start time in s
    AQ.nSamples   = [repmat(floor((TX.Start(2)-AQ.Start(1)-deadTimeRX2TX)*AQ.fSample(1)),SearchFlipDegSteps*2+1,1);floor(Seq.tAQfOffset*AQ.fSample(end))];                         % number of samples in AQ window
    AQ.Start(end) = AQ.Start(end-1)+AQ.nSamples(end-1)/AQ.fSample(end-1)+deadTimeRX2TX+get_DeadTimeTX2RX(HW, AQ.fSample(end));
    AQ.Frequency  = HW.fLarmor;                   % frequency of AQ window in Hz
    AQ.Phase      = AQPhaseOffsetDeg+180*(mod(SearchFlipDeg/90,4)>2);                            % phase of AQ window in degrees
    % maxInputVoltage = 20e-3;                      % maximum input voltage in V
    % AQ.Gain       = HW.RX(1).Amplitude2Uin / maxInputVoltage;  % amplifier gain
    % AQ.ResetPhases=0;

    Seq.tRep      = ones(1,1)*0.001+AQ.Start(end)+AQ.nSamples(end)/AQ.fSample(end);                       % repetition time in s

    % no Gradient Pulse
    for gt = 1:HW.Grad(iDevice).n
      Grad(gt).Time = NaN;
      Grad(gt).Repeat = [0, ones(1, size(Seq.tRep,2)-1)];
    end    % Sequence parameters

    % % Start measurement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);


    % % Plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_data_1D(HW, data_1D);
    %%

    [fOffset, fOffsetStd] = get_MeanPhaseDiffWeighted(data.data(1:SeqOut.AQ.nSamples(end),end,1),1);
    fOffset =fOffset./(2*pi).*SeqOut.AQ.fSample(end);
    fOffsetStd=fOffsetStd./(2*pi).*SeqOut.AQ.fSample(end);

    if Seq.TrackFrequency
      fLarmorAll(t)=SeqOut.AQ(1).Frequency(1)-fOffset;
      newB0 = fLarmorAll(t) * 2*pi/HW.GammaDef;  % calculate the new magnetic field strength
      if ~isemptyfield(HW.FindFrequencySweep, 'shiftB0') ...
          && HW.FindFrequencySweep.shiftB0
        if isempty(HW.B0Target) || ~isfinite(HW.B0Target)
          warning('PD:Find_Frequency_Sweep:NoB0Target', ...
            ['HW.FindFrequencySweep.shiftB0 is set to true, but HW.B0Target is not set.\n', ...
            'B0 amplitude cannot be shifted without a target. Adjusting HW.fLarmor instead.']);
        else
          newB0shift = HW.Grad(iDevice).AmpOffset(4) + (HW.B0Target - newB0);
          newB0 = HW.B0Target;
          if abs(newB0shift) > HW.Grad(iDevice).HoldShimNormMax(4)*HW.Grad(iDevice).MaxAmp(4)
            desiredB0shift = newB0shift;
            newB0shift = sign(newB0shift) * HW.Grad(iDevice).HoldShimNormMax(4) * HW.Grad(iDevice).MaxAmp(4);
            newB0 = newB0 - (desiredB0shift - newB0shift);
            warning('PD:Find_Frequency_Sweep:B0ShiftExceeded', ...
              ['The required B0 shift (%.1f %cT) is larger than the maximum shift (%.1f %cT).\n', ...
              'Consider setting a different HW.B0Target. Adjusting HW.fLarmor to stay in range.'], ...
              desiredB0shift*1e6, 181, HW.Grad(iDevice).HoldShimNormMax(4)*HW.Grad(iDevice).MaxAmp(4)*1e6, 181);
          end
          HW.Grad(iDevice).AmpOffset(4) = newB0shift;
        end
      end
      HW.B0 = newB0;
    end


    MeanAmp=reshape(mean(data.data(1:SeqOut.AQ.nSamples(1),1:end-1,1),1),[],1);
    errorPhase=reshape(std(unwrap(angle(data.data(1:SeqOut.AQ.nSamples(1),1:end-1,1)),[],1),1)/SeqOut.AQ.nSamples(1).^0.5,[],1);
    % Time=reshape(mean(data.time_of_tRep(1:SeqOut.AQ.nSamples(1),1:end-1,1),1),[],1);
    Time=reshape(SeqOut.TX(1).Start+SeqOut.TX(1).Duration,[],1);


%     figure(5)
%     CMA=mean(MeanAmp(:));
%     plot(Time,real(MeanAmp-CMA)/abs(CMA),Time,imag(MeanAmp-CMA)/abs(CMA),Time,angle(MeanAmp))
%     legend('real', 'imag', 'angle')
%     xlabel( 'time in s')
%     ylabel( 'amp in rel and angle RAD')
%     grid on
%     grid minor

    t1e=Time(SearchFlipDegSteps+1);
    plotShift=t1e*0;
    t1 = Time(1:SearchFlipDegSteps+1)-t1e;
    y = angle(MeanAmp(1:SearchFlipDegSteps+1));
    ypp = max(y)-min(y);                                                   % Range of /y2
    t1Per = (t1(end)-t1(1))/Seq.nPeriods;                                         % Estimate period
    fit = @(b,t1)  b(1).*(sin(2*pi*t1./b(2) + 2*pi*b(3))) + b(4);     % Function to fit
    fcn = @(b) norm(fit(b,t1) - y);                                  % Least-Squares cost function
    options = optimset('TolFun',1e-6,'Tolx',1e-6,'MaxFunEvals',100000);%,'PlotFcns',@optimplotfval);
    [s1,fval1,exitflag1,output1] = fminsearch(fcn, [ypp/2;  t1Per;  0.5/0.05;  mean(y)+ypp/10],options);  % 10 extend step size to 0.5 for phase (5% of 10)  % Minimise Least-Squares
    s1(3)=mod(s1(3),1);
    t1p = linspace(min(t1),max(t1));
    
    fh1=figure(1);
    ax1=subplot(1,1,1,'Parent',fh1);
    errorbar(ax1,t1+plotShift,rad2deg(y),rad2deg(errorPhase(1:SearchFlipDegSteps+1)),'.b'  );
    hold(ax1, 'on');
    plot(ax1,t1p+plotShift,rad2deg(fit(s1,t1p)), 'r')
    xlabel(ax1,'time in s')
    ylabel(ax1,'phase in DEG')
    text(ax1,min(t1+plotShift), rad2deg(mean(y)-ypp/10), sprintf('$y = %.3f\\ sin(2*pi*(%.3f\\cdot t %+.3f)) %+.3f$', [rad2deg(s1(1)); 1./s1(2); s1(3); rad2deg(s1(4))]), 'Interpreter','latex')

    s1Phase0s=mod((-t1e)/s1(2)+s1(3),1);
    s1PhaseEnd=mod(t1(end)/s1(2)+s1(3),1);

    t2s=Time(SearchFlipDegSteps+1);
    t2 = Time(SearchFlipDegSteps+(1:SearchFlipDegSteps+1))-t2s;
    y = angle(MeanAmp(SearchFlipDegSteps+(1:SearchFlipDegSteps+1)));
    ypp = max(y)-min(y);                                                   % Range of /y2
    t2Per = (t2(end)-t2(1))/Seq.nPeriods;                                         % Estimate period
    fit = @(b,t2)  b(1).*(sin(2*pi*t2./b(2) + 2*pi*b(3))) + b(4);     % Function to fit
    fcn = @(b) norm(fit(b,t2) - y);                                  % Least-Squares cost function
    [s2,fval2,exitflag2,output2] = fminsearch(fcn, [ypp/2;  t2Per;  0.5/0.05;  mean(y)+ypp/10],options);                       % Minimise Least-Squares
    s2(3)=mod(s2(3),1);
    t2p = linspace(min(t2),max(t2));

    errorbar(ax1,t2+plotShift,rad2deg(y),rad2deg(errorPhase(SearchFlipDegSteps+(1:SearchFlipDegSteps+1))),'.b');
    plot(ax1,t2p+plotShift,rad2deg(fit(s2,t2p)), 'r')
    grid(ax1, 'on')
    grid(ax1, 'minor')
    legend(ax1,'Data', 'Fitted Curve', 'Location','best')
    text(ax1,min(t2)+plotShift, rad2deg(mean(y)+ypp/10), sprintf('$y = %.3f\\ sin(2*pi*(%.3f\\cdot t %+.3f)) %+.3f$', [rad2deg(s2(1)); 1./s2(2); s2(3); rad2deg(s2(4))]), 'Interpreter','latex')
    hold(ax1, 'off');
    drawnow


    s2PhaseStart=s2(3);
    s2PhaseEnd=mod(t2(end)/s2(2)+s2PhaseStart,1);
    s1PhaseEndMirror=mod(0.5-s1PhaseEnd,1);
    PhaseDiffMirror=mod(s1PhaseEndMirror-s2PhaseStart,1);
    slopdiff=sign(abs(PhaseDiffMirror-0.5)-0.25);
    slopStart=sign(abs(mod(s1Phase0s,1)-0.5)-0.25);

    FlipAngleDeg=rad2deg(slopStart*(s1(1)+slopdiff*s2(1))/2)+SearchFlipDeg;

    if Seq.FixedExcitationPulseLength
      AmpPulseCorr=SeqOut.TX(1).Amplitude(1)/FlipAngleDeg*SearchFlipDeg;
    else
      tPulseCorr=SeqOut.TX(1).Duration(1)/FlipAngleDeg*SearchFlipDeg;
    end

    AQPhaseOffsetDeg=AQPhaseOffsetDeg-s1(4)/pi*180 + 0*(fOffset*(+AQ.Start(1)+AQ.nSamples(1)/AQ.fSample(1)/2)-(TX.Start(1)+TX.Duration(1)/2))*360;

    if tt<0
            PreString=', PreShot';
    else
            PreString='';
    end
    if mean(errorPhase.^2).^0.5*10>min([std(angle(MeanAmp(1:SearchFlipDegSteps+1))),std(angle(MeanAmp(SearchFlipDegSteps+(1:SearchFlipDegSteps+1))))])...
        || abs(fOffset)>10 ...
        || abs(s1(2)-t1Per)/t1Per>0.3 ...
        || abs(s2(2)-t2Per)/t2Per>0.3
      PaUoutCorrPlot(t) = NaN;
      outString = [...
        sprintf('FlipAngle = %.3f%c, tPulseCorr = %.3f %cs, PaUoutCorr = %.3f V, AQPhaseOffset = %.3f%c, (no good conditions)', ...
                FlipAngleDeg, char(176), tPulseCorr*1e6, char(181), PaUoutCorrPlot(t), AQPhaseOffsetDeg, char(176)), ...
        PreString];
      if tt<0
        disp(outString)
        skipPreLopps=0; % if pre loop was with no good conditions
      else
        warning(outString)
      end
      FlipAngleDegPlot(t)=nan;
      tPulseCorrPlot(t)=nan;
      AQPhaseOffsetDegPlot(t)=nan;
    else
      PaUoutCorrPlot(t) = AmpPulseCorr / HW.TX(iDevice).PaUout2Amplitude(HW.TX(iDevice).ChannelDef);
      outString = [...
        sprintf('FlipAngle = %.3f%c, tPulseCorr = %.3f %cs, PaUoutCorr = %.3f V, AQPhaseOffset = %.3f%c', ...
                FlipAngleDeg, char(176), tPulseCorr*1e6, char(181), PaUoutCorrPlot(t), AQPhaseOffsetDeg, char(176)), ...
        PreString];
      disp(outString);
      FlipAngleDegPlot(t)=FlipAngleDeg;
      tPulseCorrPlot(t)=tPulseCorr;

%       AQPhaseOffsetDegPlot(t)=mod(AQPhaseOffsetDeg+SeqOut.TX(1).Phase(1)+180,360)-180;
      AQPhaseOffsetDegPlot(t)=AQPhaseOffsetDeg;
      if tt<0
        skipPreLopps=1; % if pre loop was good
        if abs(SearchFlipDeg)==90
          if   Seq.AmpPulseRotate==AmpPulse90 && Seq.tPulseRotate==tPulse90*4
            Seq.tPulseRotate=tPulseCorr*4;
            Seq.AmpPulseRotate=AmpPulseCorr;
          end
          tPulse90=tPulseCorr; % get better tPulse90
          AmpPulse90=AmpPulseCorr;
        end
        AmpExcitation=AmpPulseCorr;
        tExcitation=tPulseCorr;
      end
    end

    figure(6);clf
    ax(1)=subplot(2,1,1);
    ax(2)=subplot(2,1,2);
    if Seq.FixedExcitationPulseLength
        plot(ax(1),testFlipDeg,testFlipDeg-FlipAngleDegPlot,':x');
        ylabel(ax(1),['FlipAngle offset to ' num2str(SearchFlipDeg,'%.0f') ' DEG']);
        title(ax(1),['Fixed pulse length = ' num2str(mean(tPulseCorrPlot,'omitnan')*1e6,'%.3f') ' ' char(181) 's'])
        plot(ax(2),testFlipDeg,PaUoutCorrPlot*1-mean(PaUoutCorrPlot,'omitnan')*1,':x');
        ylabel(ax(2),['PaUoutCorr to ' num2str(mean(PaUoutCorrPlot,'omitnan')*1,'%.6f') ' V']);
    else
        plot(ax(1),testFlipDeg,testFlipDeg-FlipAngleDegPlot,':x');
        ylabel(ax(1),['FlipAngle offset to ' num2str(SearchFlipDeg,'%.0f') ' DEG']);
        title(ax(1),['Fixed pulse amplitude = ' num2str(mean(PaUoutCorrPlot,'omitnan')*1,'%.3f') ' ' char(181) 'V'])
        plot(ax(2),testFlipDeg,tPulseCorrPlot*1e6-mean(tPulseCorrPlot,'omitnan')*1e6,':x');
        ylabel(ax(2),['tPulseCorr to ' num2str(mean(tPulseCorrPlot,'omitnan')*1e6,'%.3f') ' ' char(181) 's']);
    end
    grid(ax(1), 'on');
    grid(ax(1), 'minor');
    grid(ax(2), 'on');
    grid(ax(2), 'minor');
    linkaxes(ax,'x')

    firstLoop=0;
  end
  if Seq.FixedExcitationPulseLength
    fprintf('tPulse = %.3f %cs, PaUoutCorr mean = %.3f V (STD %.3f %%), AQPhaseOffset = %.3f%c\n', ...
      mean(tPulseCorrPlot(1:tt), 'omitnan') * 1e6, char(181), ...
      mean(PaUoutCorrPlot(1:tt), 'omitnan'), std(PaUoutCorrPlot(1:tt), 'omitnan') / mean(PaUoutCorrPlot(1:tt), 'omitnan') * 100, ...
      mean(AQPhaseOffsetDegPlot, 'omitnan'), char(176));
  else
    fprintf('tPulseCorr mean = %.3f %cs (STD %.3f %%), PaUout = %.3f V, AQPhaseOffset = %.3f%c\n', ...
      mean(tPulseCorrPlot(1:tt), 'omitnan') * 1e6, char(181), std(tPulseCorrPlot(1:tt), 'omitnan') / mean(tPulseCorrPlot(1:tt), 'omitnan') * 100, ...
      mean(PaUoutCorrPlot(1:tt), 'omitnan'), ...
      mean(AQPhaseOffsetDegPlot, 'omitnan'), char(176));
  end
  TXVoltage=mean(PaUoutCorrPlot,'omitnan');

  %% Store calculated values or display them
  Pulse90PO=0;
  if abs(SearchFlipDeg)==90 && ~isnan(mean(tPulseCorrPlot(1:tt),'omitnan'))
    data.tPulse90 = mean(tPulseCorrPlot(1:tt),'omitnan');
    data.savePulseFile = true;

    if isnan(TXVoltage)
      % TXVoltage could not be determined in any of the iteration steps
      data.savePulseFile = false;
    end

    % calculate magnetic flux density of the transmit coil at HW.TX.AmpDef amplitude
    data.B1 = ((pi/2)/data.tPulse90) / HW.GammaDef;
    iDevice = 1;
    PaUout2AmplitudeFieldname = 'PaUout2Amplitude';
    fLarmorPulseDuration = HW.fLarmor;

    newPaUout2Amplitude = HW.TX(iDevice).(PaUout2AmplitudeFieldname);
    newPaUout2Amplitude(HW.TX(iDevice).ChannelDef) = data.B1 / TXVoltage;
    comment = sprintf('%s (%.0f%s pulse duration = %.3f us @ %.3f V @ %.6f MHz) for 0d by %s', ...
      datestr(now, 'yyyy-mm-ddTHH:MM:SS'), SearchFlipDeg, char(176), data.tPulse90*1e6, TXVoltage, fLarmorPulseDuration/1e6, mfilename());

    newCalLine = sprintf('HW.TX(%d).%s = [%.6f, %.6f]*1e-6;', ...
      iDevice, PaUout2AmplitudeFieldname, newPaUout2Amplitude*1e6);

    HW.TX(iDevice).PaUout2Amplitude = newPaUout2Amplitude;

    if abs(mean(AQPhaseOffsetDegPlot,'omitnan'))>Seq.AQPhaseOffsetDegSaveLevel
      newCalLineLMS= sprintf('\nHW = set_AQPhaseOffset(HW, %.3f); %% %s for 0D by %s', mean(AQPhaseOffsetDegPlot,'omitnan'),datestr(now, 'yyyy-mm-ddTHH:MM:SS'), mfilename());
      fid = fopen(fullfile(HW.RootPath, HW.UserPath, 'LoadMySystem.m'), 'a+');
      fwrite(fid, newCalLineLMS);
      fprintf('\nA new line was added to the following file:\n%s%s\n', ...
        fullfile(HW.RootPath, HW.UserPath, 'LoadMySystem.m'), newCalLineLMS);
      Pulse90PO=mean(AQPhaseOffsetDegPlot,'omitnan');
      HW = set_AQPhaseOffset(HW,mean(AQPhaseOffsetDegPlot,'omitnan'));
    end

    if ~isempty(HW.TX(iDevice).CoilName)
      newCalLine = sprintf('if strcmp(HW.TX(%d).CoilName, ''%s''),  %s  end', ...
        iDevice, HW.TX(iDevice).CoilName, newCalLine);
    end

    if iDevice > 1
      newCalLine = sprintf('if numel(HW.TX) >= %d,  %s  end', ...
        iDevice, newCalLine);
    end

    newCalLine = sprintf('%s  %% %s\n', newCalLine, comment);

    if ~isempty(HW.TX(iDevice).PaUout2AmplitudePath) ...
        && (isemptyfield(mySave, 'DummySerial') ...
        || mySave.DummySerial(min(iDevice, numel(mySave.DummySerial))) <= 0)
      addFirstLine = ~exist(HW.TX(iDevice).PaUout2AmplitudePath, 'file');
      if ~exist(fileparts(HW.TX(iDevice).PaUout2AmplitudePath), 'dir')
        mkdir(fileparts(HW.TX(iDevice).PaUout2AmplitudePath));
      end
      fid = fopen(HW.TX(iDevice).PaUout2AmplitudePath, 'a+');
      fid_protect = onCleanup(@() fclose(fid));
      if addFirstLine
        fwrite(fid, ['% factor from voltage amplitude at the coil input to B1+ field strength in T/V', sprintf('\n')]);
      end
      if data.savePulseFile
        fwrite(fid, newCalLine);
        fprintf('\nA new line was added to the following file:\n%s\n%s\n', ...
          HW.TX(iDevice).PaUout2AmplitudePath, newCalLine);
      else
        fwrite(fid, ['% ', newCalLine]);  % add line as comment
      end
      delete(fid_protect);
      [~, name, ~] = fileparts(fopen(fid));
      clear(name);clear('name')

      % save the time of the last RF pulse duration search
      mySave.lastTime_PulseDuration = now*24*3600;
    elseif data.savePulseFile
      fprintf('\n');
      fprintf('\nPlease add the following line to your LoadMySystem.m file:\n%s\n', ...
        newCalLine);
    end
    if ~data.savePulseFile
      fprintf('\n');
      warnStr = ['Determination of pulse length unsuccessful!\n', ...
        'If you want to use the uncertain best guess value anyway, ', ...
        'please manually append or un-comment the following line in your PaUout2AmplitudeCal.m file:\n%s\n'];
      warning('PD:sequence_PulseDuration', warnStr, newCalLine);
    end
  end

  PaUout=[PaUout;mean(PaUoutCorrPlot,'omitnan')]; %#ok<AGROW>
  FlipDeg=[FlipDeg;SearchFlipDeg]; %#ok<AGROW>
  tPulse=[tPulse ; mean(tPulseCorrPlot,'omitnan')]; %#ok<AGROW>
  B1 = [B1; ((SearchFlipDeg/180*pi)/mean(tPulseCorrPlot,'omitnan')) / HW.GammaDef]; %#ok<AGROW>
  AQPhaseDeg=[AQPhaseDeg; (mean(AQPhaseOffsetDegPlot,'omitnan'))-Pulse90PO]; %#ok<AGROW>
  fprintf('%.0f%c pulse duration: %.3f %cs @ %.3f V; B1 = %.3f %cT; B1 fLarmor = %.3f Hz\n', ...
    SearchFlipDeg, char(176), tPulse(end)*1e6, char(181), TXVoltage, B1(end)*1e6, char(181), B1(end)*HW.GammaDef/2/pi);

end


%% display summary
B1end = NaN;
if numel(FlipDeg)>1
  B1end=(((FlipDeg(end)-FlipDeg(end-1))/180*pi) ./ (tPulse(end)-tPulse(end-1)))./ HW.GammaDef;
end
B1rel=B1./B1end.*100;
% PaUoutEnd=(PaUout(end)*tPulse(end)/FlipDeg(end)-PaUout(end-1)*tPulse(end-1)/FlipDeg(end-1));
PaUoutRel=PaUout./PaUout(end).*100;

fprintf('\n');
disp(['FlipDeg   ' num2str(FlipDeg.','%+10.3f') ' DEG']);
disp(['B1        ' num2str(B1.'*1e6,'%+10.3f') ' ' char(181) 'T']);
disp(['B1rel     ' num2str(B1rel.','%+10.3f') ' %']);
disp(['tPulse    ' num2str(tPulse.'*1e6,'%+10.3f') ' ' char(181) 's']);
disp(['Phase     ' num2str(AQPhaseDeg.','%+10.3f') ' DEG']);
disp(['PaUout    ' num2str(PaUout.','%+10.3f') ' V']);
disp(['PaUoutrel ' num2str(PaUoutRel.','%+10.3f') ' %']);
fprintf('\n');

end
