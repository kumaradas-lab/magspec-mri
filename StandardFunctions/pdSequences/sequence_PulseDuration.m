function [SeqOut, mySave] = sequence_PulseDuration(HW, Seq, AQ, TX, Grad, mySave)
%% Find 90 degrees pulse duration and store corresponding calibration value in file
%
%   SeqOut = sequence_PulseDuration(HW, Seq, AQ, TX, Grad, mySave)
%
% This function searches the correct pulse length for a 90 degrees pulse at the
% default amplitude and stores it in HW and mySave. Additionally, it appends the
% corresponding factor for output voltage on the TX channel to B1 field strength
% to the PaUout2AmplitudeCal.m file in the User folder.
% The pulse length is determined by iteratively increasing the pulse length and
% performing a 1D read encoded Spin Echo measurement along the y-axis of the
% sample. The amplitude at the inner half of the acquired 1D images is used for
% calculating the pulse length (i.e. the efficiency of the coil).
% Consistency checks with the apparent 90 degrees and 180 degrees pulse length
% at the center and the inner half of the acquired 1D images are performed to
% evaluate the quality of the results.
%
% After the pulse duration is determined, the execution is paused for
% HW.FindFrequencyPause seconds.
%
%
% INPUT:
%
%   HW
%       HW structure or object
%
%   Seq
%       structure with the following fields (default values are used if omitted
%       or empty):
%
%     doPlot
%         Plot sequence and data (Boolean, default:
%         HW.FindPulseDuration.doPlot = true)
%
%     hParent
%         handle to a figure or uipanel that is used for plotting (default: 320)
%
%     T1Estimated
%         Estimated T1 value of the sample in seconds. A pause of
%         3*Seq.tEstimated is added between each iteration step
%         (default: HW.FindPulseDuration.T1Estimated = 120 ms)
%
%     tPulse90Estimated
%         Estimated 90 degrees pulse length in seconds. This is used to
%         calculate the estimated pulse lengths for a given flip angle.
%         (Default: HW.FindPulseDuration.tPulse90Estimated, depends on the
%         setup).
%
%     excitationFlipAngleStart
%         Smallest flip angle for the scan in degrees (default:
%         HW.FindPulseDuration.excitationFlipAngleStart = 10)
%
%     excitationFlipAngleStop
%         Largest flip angle for the scan in degrees (default:
%         HW.FindPulseDuration.excitationFlipAngleStop = 600)
%
%     excitationFlipAngleSteps
%         Number of the uniformly spaced steps for the scan (default:
%         HW.FindPulseDuration.excitationFlipAngleSteps = 60)
%
%     tEcho
%         Echo time in seconds. If empty, a value sufficiently large for
%         Seq.excitationFlipAngleStop is chosen (default:
%         HW.FindPulseDuration.tEcho = []).
%
%     plotSeq
%         Plot sequence with the given gradients. If empty no pulse program is
%         shown. (Default: HW.FindPulseDuration.plotSeq = []).
%
%     plotSeqTR
%         Plot "stacked" sequence with the given gradients. If empty no
%         "stacked" pulse program is shown. (Default:
%         HW.FindPulseDuration.plotSeqTR = []).
%
%   AQ
%       structure with the settings for the acquisition windows
%
%   TX
%       structure with the settings for the rf pulses
%
%   Grad
%       structure with the settings for the gradients
%
%   mySave
%       structure with data necessary for the timing of frequency searches
%
%
% OUTPUT:
%
%   SeqOut
%       structure as input Seq but with added fields with the actually used
%       values. Additionally, SeqLoop.dataPulseDuration holds the data with the
%       measurement results.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2024 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------


%% Default input
if isemptyfield(Seq, 'PreProcessSequence'), Seq.PreProcessSequence  = 1;  end
if isemptyfield(Seq, 'StartSequence'),      Seq.StartSequence       = 1;  end
if isemptyfield(Seq, 'PollPPGfast'),        Seq.PollPPGfast         = 1;  end
if isemptyfield(Seq, 'GetRawData'),         Seq.GetRawData          = 1;  end
if isemptyfield(Seq, 'PostProcessSequence'),Seq.PostProcessSequence = 1;  end
if isemptyfield(Seq, 'PostProcessLocal'),   Seq.PostProcessLocal    = 1;  end


%% parameters that are also needed when preparation is skipped
iDevice = 1;  % FIXME: Actually support multiple devices
if isemptyfield(Seq, 'useGammaX')
  Seq.useGammaX = false;
end

if isa(HW, 'PD.HWClass')
  if Seq.useGammaX
    % imaging sequences currently only support measurements at the primary nucleus

    % save current settings to restore them on exit
    oldPaUout2Amplitude = HW.TX(iDevice).PaUout2Amplitude;
    oldPaUout2AmplitudeX = HW.TX(iDevice).PaUout2AmplitudeX;
    oldGammaDef = HW.GammaDef;

    guardHW = onCleanup(@() ...
      restoreHWsettings(HW, iDevice, oldPaUout2Amplitude, oldPaUout2AmplitudeX, oldGammaDef));

    HW.TX(iDevice).PaUout2Amplitude = HW.TX(iDevice).PaUout2AmplitudeEstimatedX;
    HW.TX(iDevice).PaUout2AmplitudeX = HW.TX(iDevice).PaUout2AmplitudeEstimatedX;

    HW.GammaDef = HW.GammaX;

  else
    % save current settings to restore them on exit
    oldPaUout2Amplitude = HW.TX(iDevice).PaUout2Amplitude;

    guardHW = onCleanup(@() ...
      restoreHWsettings(HW, iDevice, oldPaUout2Amplitude));

    HW.TX(iDevice).PaUout2Amplitude = HW.TX(iDevice).PaUout2AmplitudeEstimated;
  end

  ampDef = get_TX_Amplitude(HW, ...
    'Uout', HW.TX(iDevice).Def.UoutCalibrated(HW.TX(iDevice).ChannelDef), ...
    'Device', iDevice);

  currentTFlip90Def = 2*pi / HW.GammaDef / ampDef / 4;

else
  if Seq.useGammaX
    warning('PD:sequence_PulseDuration:NoGammaX', ...
      'Pulse duration determination at GammaX is only supported if HW is a PD.HWClass object.');
  end
  currentTFlip90Def = HW.tFlip90Def;
end


%% Preparation
if Seq.PreProcessSequence
  Seq = set_EmptyField(Seq, 'doPlot',                   HW.FindPulseDuration.doPlot);
  Seq = set_EmptyField(Seq, 'T1Estimated',              HW.FindPulseDuration.T1Estimated); % estimated T1 time of the used sample (oil is recommened)
  Seq = set_EmptyField(Seq, 'tPulse90Estimated',        HW.FindPulseDuration.tPulse90Estimated);
  if isempty(Seq.tPulse90Estimated), Seq.tPulse90Estimated = currentTFlip90Def; end
  Seq = set_EmptyField(Seq, 'excitationFlipAngleStart', HW.FindPulseDuration.excitationFlipAngleStart);
  Seq = set_EmptyField(Seq, 'excitationFlipAngleStop',  HW.FindPulseDuration.excitationFlipAngleStop);
  Seq = set_EmptyField(Seq, 'excitationFlipAngleSteps', HW.FindPulseDuration.excitationFlipAngleSteps);
  Seq = set_EmptyField(Seq, 'SteadyState_PreShots90',   HW.FindPulseDuration.SteadyState_PreShots90);
  Seq = set_EmptyField(Seq, 'SteadyState_PostShots90',  0);
  if ~isfield(Seq, 'AQSlice'), Seq.AQSlice = struct(); end
  if isemptyfield(Seq.AQSlice, 'iDevice'), Seq.AQSlice.iDevice = 1; end  % index of device that is used for sending the pulses
  Seq = set_EmptyField(Seq, 'tEcho',                    HW.FindPulseDuration.tEcho); % Echo time in seconds eg. 5e-3
  if isempty(Seq.tEcho)
    % Estimate for dephasing gradient duration: tDephase = tAQ/2 + 2*tRamp + tEC
    % Estimate for tAQ (default below): tAQ = p180 + 5e-3 (????)
    % --> tEcho_min/2 = p90/2 + tDephase + p180/2
    Seq.tEcho = 5e-3 + Seq.tPulse90Estimated*Seq.excitationFlipAngleStop/90*5 + ...
      4*HW.Grad(Seq.AQSlice.iDevice).tRamp + 2*HW.Grad(Seq.AQSlice.iDevice).tEC;
  end
  if ~isfield(Seq, 'plotSeqTR'), Seq.plotSeqTR = HW.FindPulseDuration.plotSeqTR; end  % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
  if ~isfield(Seq, 'plotSeq'), Seq.plotSeq = HW.FindPulseDuration.plotSeq; end  % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
  Seq.AQSlice = set_EmptyField(Seq.AQSlice, 'sizeRead',     HW.FindPulseDuration.AQSlice.sizeRead); % Image size in read in meter (for CSI set to 1e12)
  if isempty(Seq.AQSlice.sizeRead),  Seq.AQSlice.sizeRead = HW.Grad(Seq.AQSlice.iDevice).ImageVol(4) - HW.Grad(Seq.AQSlice.iDevice).ImageVol(3); end
  Seq.AQSlice = set_EmptyField(Seq.AQSlice, 'nRead',        HW.FindPulseDuration.AQSlice.nRead);  % Number of Pixels in read, if nRead>1 nPhase(1)=1 )
  Seq.AQSlice = set_EmptyField(Seq.AQSlice, 'HzPerPixMin',  HW.FindPulseDuration.AQSlice.HzPerPixMin); % Number of Pixels in read, if nRead>1 nPhase(1)=1 )
  if isempty(Seq.AQSlice.HzPerPixMin)
    Seq.AQSlice.HzPerPixMin = 2.5/(Seq.tEcho - Seq.tPulse90Estimated*Seq.excitationFlipAngleStop/90*2 - ...
                                   (4*HW.Grad(Seq.AQSlice.iDevice).tRamp + 2*HW.Grad(Seq.AQSlice.iDevice).tEC));
  end
  Seq.AQSlice = set_EmptyField(Seq.AQSlice, 'alfa',         HW.FindPulseDuration.AQSlice.alfa);   % first rotation around x axis in RAD
  Seq.AQSlice = set_EmptyField(Seq.AQSlice, 'phi',          HW.FindPulseDuration.AQSlice.phi);    % second rotation around y axis in RAD
  Seq.AQSlice = set_EmptyField(Seq.AQSlice, 'theta',        HW.FindPulseDuration.AQSlice.theta);  % third rotation around z axis in RAD


  %%
  Seq.LoopsBreak = Seq.T1Estimated*3;       % Pause between two loop averages in seconds ([]= fast as possible)

  Seq.LoopSaveAllData = 1;                  % Save all data of the Loops
  Seq.LoopPlotAverages = 0;

  Seq.RepetitionTime = Seq.T1Estimated*3;   % Repetition time in seconds; only used if Seq.TurboBreak=[], auto Seq.tRep=[] and Seq.AQSlice(1).TurboBreak= x sec;

  % % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Seq.AQSlice(1).nRead = 64;                % number of pixels in read; if nRead>1, nPhase(1)=1
  Seq.AQSlice(1).nPhase(1) = 1;             % number of pixels in phase(1)
  Seq.AQSlice(1).nPhase(2) = 1;             % number of pixels in phase(2)
  Seq.AQSlice(1).nPhase(3) = 1;             % number of pixels in phase(3) Number of Pixels in for CSI, if nPhase(3)>1 nRead=1 sizeRead=1e12
  % Seq.AQSlice(1).HzPerPixMin = 500;         % bandwith per pixel in Hz (1/HzPerPixMin = duration of AQ)
  % Seq.AQSlice(1).sizeRead = 0.0128;         % image size in read in meter (for CSI set to Inf)
  Seq.AQSlice(1).sizePhase(1) = Inf;        % image size in phase(1) in meter
  Seq.AQSlice(1).sizePhase(2) = Inf;        % image size in phase(2) in meter
  Seq.AQSlice(1).sizePhase(3) = Inf;        % image size in phase(3) in meter
  % Seq.AQSlice(1).thickness = 0.002;         % image thickness in slice direction  used for 2D and 3D! ([] for no Slice) in meter
  % Seq.AQSlice(1).thicknessInversion = Seq.AQSlice(1).thickness+1e12;  % inversion slice thickness in slice direction  used for 2D and 3D! in meter ([] for no Slice)
  % Seq.AQSlice(1).excitationPulse = @Pulse_LoadOwnPulseFromMat;  % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
  Seq.AQSlice = set_EmptyField(Seq.AQSlice, 'excitationPulse', @Pulse_Rect);   % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
  Seq.AQSlice = set_EmptyField(Seq.AQSlice, 'inversionPulse', @Pulse_Rect);   % inversion pulse function

  % % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Seq.AQSlice(1).ReadOS = 16;             % Oversampling read ([] for automatic, recommended >=2) 1...
  Seq.AQSlice(1).PhaseOS(1) = 1;            % Oversampling phase(1)  1...
  Seq.AQSlice(1).PhaseOS(2) = Seq.excitationFlipAngleSteps;  % Oversampling phase(2); used for excitation angle increment
  Seq.AQSlice(1).PhaseOS(3) = 1;            % Oversampling phase(3)  1...; set to 1 if nPhase(3)=1;

  % % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Seq.AQSlice(1).TurboFactor = 1;           % number of echoes per excitation
  % Seq.AQSlice(1).TurboBreak = Seq.RepetitionTime;  % break between last echo and next excitation

  % % Plot        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Seq.AQSlice(1).plotkSpace = 0;            % plot k-space
  Seq.AQSlice(1).plotImage = 0;             % plot image
  Seq.AQSlice(1).plotPhase = 0;             % plot phase of k-space and/or image
  Seq.AQSlice(1).plotFft1_data = 0;         % plot 1D FFT
  Seq.AQSlice(1).plotData = 0;              % plot Echoes
  % Seq.AQSlice(1).raiseFigures = false;      % raise and focus windows when updating plot
  Seq.AQSlice(1).ZeroFillWindowSize = 1.4;  % Zero Fill Window Size (k-Space)
  Seq.AQSlice(1).ZeroFillFactor = 2;        % Zero fill resolution factor
  % Seq.AQSlice(1).plotImageHandle = 111+Loop;  % figure handle of the 2D plot

  % % Some corrections    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Seq.AQSlice(1).sizePhaseSpoil = 10.0001;  % coding length of spoiler
  % HW.Grad.SliceTimeDelay = [0,0,0]*1e-6;    % delay only for the Slice gradient
  % HW.Grad.ReadTimeDelay = [0,0,0]*1e-6;     % delay only for the read gradient
  % HW.Grad.PhaseTimeDelay = [0,0,0]*1e-6;    % delay only for the phase gradient
  Seq.LoopsBreakExactly = Seq.LoopsBreak<2; % Lock Time between two loops
  Seq.LoopPlotAll = Seq.LoopsBreak>=0.4;
  % HW.Function_Run_Drawnow_While_Measurement = 1;

  % Seq.AQSlice(1).ReadOSUsedForImage = 17;   % Number of samples used in CSI
  % Seq.AQSlice(1).AQWindowDirection = [1,-1];  % [1,-1] to reverse even readouts
  % Seq.MaxGradAmpSlice = 0.05;               % limit slice gradient strength
  % Seq.AQSlice(1).sizePhaseSpoil = 0.0001;   % coding length of spoiler
  if ~isfield(Seq.AQSlice(1), 'SpoilFactor')
    Seq.AQSlice(1).SpoilFactor = [1, 1, 1];
  end
  if ~isfield(Seq.AQSlice(1), 'phaseCycling')
    Seq.AQSlice(1).phaseCycling = 0;
  end
  Seq.AQSlice(1).phaseCycling=double(any(Seq.AQSlice(1).phaseCycling)); % zero or one


  Seq.AQSlice(1).excitationFlipAngle = [repmat(Seq.excitationFlipAngleStart, 1, Seq.SteadyState_PreShots90), ... % pre-shots
    reshape(repmat(linspace(Seq.excitationFlipAngleStart, Seq.excitationFlipAngleStop, Seq.excitationFlipAngleSteps),Seq.AQSlice(1).phaseCycling+1,1),1,[]), ...       % actual steps
    repmat(Seq.excitationFlipAngleStop, 1, Seq.SteadyState_PostShots90)] * Seq.tPulse90Estimated/currentTFlip90Def;  % post-shots
  Seq.AQSlice(1).inversionFlipAngle = Seq.AQSlice(1).excitationFlipAngle*2;
  % Seq.tReadoutOffset = -0.2e-3;

  % Seq.CorrectReadRephase = 1;

  oldStartSequence = Seq.StartSequence;
  oldPollPPGfast = Seq.PollPPGfast;
  oldGetRawData = Seq.GetRawData;
  oldPostProcessSequence = Seq.PostProcessSequence;

  Seq.StartSequence = 0;
  Seq.PollPPGfast = 0;
  Seq.GetRawData = 0;
  Seq.PostProcessSequence = 0;
  [Seq, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);
  Seq.StartSequence = oldStartSequence;
  Seq.PollPPGfast = oldPollPPGfast;
  Seq.GetRawData = oldGetRawData;
  Seq.PostProcessSequence = oldPostProcessSequence;
end


%% actual measurement
if Seq.StartSequence || Seq.PollPPGfast || Seq.GetRawData || Seq.PostProcessSequence

  TXVoltage = HW.TX(Seq.AQSlice.iDevice).AmpDef / HW.TX(Seq.AQSlice.iDevice).PaUout2Amplitude(HW.TX(Seq.AQSlice.iDevice).ChannelDef);
  fprintf('Searching pulse duration of 90%c pulse at default amplitude (%.3f V) around %.3f %cs.\n', ...
    176, TXVoltage, Seq.tPulse90Estimated*1e6, 181);

  oldPreProcessSequence = Seq.PreProcessSequence;
  Seq.PreProcessSequence = 0;
  [SeqOut, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);
  SeqOut.PreProcessSequence = oldPreProcessSequence;
else
  SeqOut = Seq;
end


%% evaluate results
iAQ = find([SeqOut.AQ(:).Device] == SeqOut.AQSlice(1).iDevice, 1, 'first');  % FIXME: Support multiple channels?
if SeqOut.PostProcessLocal
  data.Image = mean(SeqOut.data(iAQ).fft1_dataCut(floor(SeqOut.AQSlice(1).nRead*SeqOut.AQSlice(1).ReadOS/2) + (1-floor(SeqOut.AQSlice(1).nRead/2):ceil(SeqOut.AQSlice(1).nRead/2)),:,:,:,:,:),6);
  data.tFlip = currentTFlip90Def/90*SeqOut.AQSlice(1).excitationFlipAngle(SeqOut.SteadyState_PreShots90+1:Seq.AQSlice(1).phaseCycling+1:end-SeqOut.SteadyState_PostShots90);
  data.AmplitudeAbs = mean(abs(data.Image(floor(end/2-end/20):ceil(end/2+end/20),:)), 1);
  [maxA, ~] = max(data.AmplitudeAbs);

  % check range for search of first maximum
  data.savePulseFile = true;
  data.index45 = find(data.AmplitudeAbs>maxA/2, 1, 'first');
  if isempty(data.index45) || data.index45==1
    warning('PD:sequence_PulseDuration', 'Decrease Seq.excitationFlipAngleStart')
    data.savePulseFile = false;
  end

  nRange135 = min(max(2, floor(length(data.tFlip)/10)), 4);
  data.index135 = data.index45 + find(data.AmplitudeAbs(data.index45+1:end)<maxA/2, 1, 'first');
  if isempty(data.index135)
    warning('PD:sequence_PulseDuration', 'Increase Seq.excitationFlipAngleStop')
    data.index135 = min(data.index45+nRange135, length(data.tFlip));
    data.savePulseFile = false;
  elseif data.index135 < data.index45+nRange135
    warning('PD:sequence_PulseDuration', 'Increase Seq.excitationFlipAngleSteps')
    data.index135 = min(data.index45+nRange135, length(data.tFlip));
    data.savePulseFile = false;
  end

  nRange225 = min(max(4, floor(length(data.tFlip)/6)), 8);
  data.index225 = data.index135 + find(data.AmplitudeAbs(data.index135+1:end)>maxA/2, 1, 'first');
  if isempty(data.index225)
    warning('PD:sequence_PulseDuration', 'Increase Seq.excitationFlipAngleStop')
    data.index225 = min(data.index135+nRange225, length(data.tFlip));
    data.savePulseFile = false;
  elseif data.index225 < data.index135+nRange225
    warning('PD:sequence_PulseDuration', 'Increase Seq.excitationFlipAngleSteps')
    data.index225 = min(data.index135+nRange225, length(data.tFlip));
    data.savePulseFile = false;
  end

  [maxA, maxI] = max(data.AmplitudeAbs(data.index45:data.index135));
  maxI = maxI+data.index45-1;

  % FIXME: The following command "corrects" the phase for each pixel in read
  % direction separately. Would it be better to correct all pixels with the
  % average phase from the center? What about phase slope (e.g., due to
  % incorrect gradient delay)?
  data.Image = bsxfun(@times, data.Image, exp(-1i*(0+angle(data.Image(:,maxI)))));

  %   data.Amplitude = mean(data.Image(floor(end/2-end/4):ceil(end/2+end/4),:), 1);
  AmplitudeI = (floor(size(data.Image,1)/2)+1)+(round(-size(data.Image,1)/4+0.5):round(size(data.Image,1)/4-0.5));
  data.Amplitude = mean(data.Image(AmplitudeI,:), 1);
  %   data.AmplitudeCenter = mean(data.Image(floor(end/2-end/20):ceil(end/2+end/20),:), 1);
  AmplitudeCenterI = (floor(size(data.Image,1)/2)+1)+(round(-size(data.Image,1)/32):round(size(data.Image,1)/32));
  data.AmplitudeCenter = mean(data.Image(AmplitudeCenterI,:), 1);

  data.tFlipInterp = interpn(data.tFlip, 8).';
  data.AmplitudeCenterInterp = interpn(data.AmplitudeCenter, 8, 'cubic').';
  data.AmplitudeInterp = interpn(data.Amplitude, 8, 'cubic').';
  data.index45Interp = min([find(data.tFlipInterp>data.tFlip(data.index45), 1, 'first'), length(data.tFlipInterp)]);
  data.index135Interp = min([find(data.tFlipInterp>data.tFlip(data.index135), 1, 'first'), length(data.tFlipInterp)]);
  data.index225Interp = min([find(data.tFlipInterp>data.tFlip(data.index225), 1, 'first'), length(data.tFlipInterp)]);

  [data.maxAmpCenter, data.maxICenter] = max(real(data.AmplitudeCenterInterp(1:data.index225Interp)));
  data.ZCICenter = data.index135Interp-1 + find(real(data.AmplitudeCenterInterp(data.index135Interp:end))<0, 1, 'first');
  data.tFlip90Center =  data.tFlipInterp(data.maxICenter);
  data.tFlip180Center = data.tFlipInterp(data.ZCICenter);
  data.tFlipRatioCenter = data.tFlip180Center/data.tFlip90Center;

  [data.maxA, data.maxI] = max(real(data.AmplitudeInterp(1:data.index225Interp)));
  data.ZCI = data.index135Interp-1 + find(real(data.AmplitudeInterp(data.index135Interp:end))<0, 1, 'first');
  data.tFlip90 = data.tFlipInterp(data.maxI);
  data.tFlip180 = data.tFlipInterp(data.ZCI);
  data.tFlipRatio = data.tFlip180/data.tFlip90;

  if isempty(data.tFlipRatio) || data.tFlipRatio<1.5 || 2.5<data.tFlipRatio
    warning('PD:sequence_PulseDuration', ...
      'Determined 90 degree and 180 degree pulse lengths do not match!')
    data.savePulseFile = false;
  end
  if isempty(data.tFlipRatioCenter) || data.tFlipRatioCenter<1.5 || 2.5<data.tFlipRatioCenter
    warning('PD:sequence_PulseDuration', ...
      'Determined 90 degree and 180 degree pulse lengths at the center do not match!')
    data.savePulseFile = false;
  end

  data.B1 = ((pi/2)/data.tFlip90) / HW.GammaDef;  % calculate magnetic flux density of the transmit coil at HW.TX.AmpDef amplitude

  SeqOut.dataPulseDuration = data;
end

if SeqOut.doPlot
  data = SeqOut.dataPulseDuration;
  if isemptyfield(Seq, 'hParent'), Seq.hParent = 320; end
  if ishghandle(Seq.hParent, 'figure') || isnumeric(Seq.hParent)
    newFigure = ~ishghandle(Seq.hParent, 'figure');
    data.fh = figure(Seq.hParent);
    clf(data.fh, 'reset');
    set(data.fh, 'Name', sprintf('Find 90%c pulse duration', 176), 'NumberTitle', 'off');
    if newFigure
      figPos = get(data.fh, 'Position');
      figPos(1) = figPos(1) - (920-figPos(3))/2;
      figPos(3) = 920;
      set(data.fh, 'Position', figPos);
    end
  elseif ishghandle(Seq.hParent, 'uipanel')
    data.fh = Seq.hParent;
    hKids = get(Seq.hParent, 'Children');
    delete(hKids);
  else
    error('Seq.hParent must be a handle to either a figure or a uipanel.');
  end

  % The following too lines are exact copies from the above section. They are
  % needed here, too, because this function can be called with
  % Seq.PostProcessLocal=false (e.g., by the teach GUI).
  AmplitudeI = (floor(size(data.Image,1)/2)+1)+(round(-size(data.Image,1)/4+0.5):round(size(data.Image,1)/4-0.5));
  AmplitudeCenterI = (floor(size(data.Image,1)/2)+1)+(round(-size(data.Image,1)/32):round(size(data.Image,1)/32));

  maxA = max(max(abs(data.Image(floor(end/2-end/4):ceil(end/2+end/4),:))));
  hax = subplot(2,3,1, 'Parent', data.fh);
  imagesc(data.tFlip*1e6, SeqOut.data(iAQ).Ticks(1).Read*1e3, real(data.Image), 'Parent', hax);
  hold(hax,'on')
  dr=diff(SeqOut.data(iAQ).Ticks(1).Read(1:2))/2*1e3;
  plot(hax, data.tFlip([1,end])*1e6,-dr+SeqOut.data(iAQ).Ticks(1).Read(AmplitudeI([1,1]))*1e3,'k:');
  plot(hax, data.tFlip([1,end])*1e6,+dr+SeqOut.data(iAQ).Ticks(1).Read(AmplitudeI([end,end]))*1e3,'k:');
  plot(hax, data.tFlip([1,end])*1e6,-dr+SeqOut.data(iAQ).Ticks(1).Read(AmplitudeCenterI([1,1]))*1e3,'b:');
  plot(hax, data.tFlip([1,end])*1e6,+dr+SeqOut.data(iAQ).Ticks(1).Read(AmplitudeCenterI([end,end]))*1e3,'b:');
  hold(hax,'off')
  set(hax, 'CLim', [-maxA, maxA]);
  title(hax, 'real');
  ylabel(hax, 'read in mm');
  xlabel(hax, sprintf('pulse length in %cs', 181));
  set(hax, 'YDir', 'normal');

  hax = subplot(2,3,2, 'Parent', data.fh);
  imagesc(data.tFlip*1e6, SeqOut.data(iAQ).Ticks(1).Read*1e3, imag(data.Image), 'Parent', hax);
  hold(hax,'on')
  dr=diff(SeqOut.data(iAQ).Ticks(1).Read(1:2))/2*1e3;
  plot(hax, data.tFlip([1,end])*1e6,-dr+SeqOut.data(iAQ).Ticks(1).Read(AmplitudeI([1,1]))*1e3,'k:');
  plot(hax, data.tFlip([1,end])*1e6,+dr+SeqOut.data(iAQ).Ticks(1).Read(AmplitudeI([end,end]))*1e3,'k:');
  plot(hax, data.tFlip([1,end])*1e6,-dr+SeqOut.data(iAQ).Ticks(1).Read(AmplitudeCenterI([1,1]))*1e3,'b:');
  plot(hax, data.tFlip([1,end])*1e6,+dr+SeqOut.data(iAQ).Ticks(1).Read(AmplitudeCenterI([end,end]))*1e3,'b:');
  hold(hax,'off')
  set(hax, 'CLim', [-maxA, maxA]);
  title(hax, 'imag');
  xlabel(hax, sprintf('pulse length in %cs', 181));
  set(hax, 'YDir', 'normal');

  hax = subplot(2,3,3, 'Parent', data.fh);
  imagesc(abs(data.Image), 'Parent', hax);
  hold(hax,'on')
  dr=1/2;
  plot(hax, [1,numel(data.tFlip)],-dr+AmplitudeI([1,1]),'k:');
  plot(hax, [1,numel(data.tFlip)],+dr+AmplitudeI([end,end]),'k:');
  plot(hax, [1,numel(data.tFlip)],-dr+AmplitudeCenterI([1,1]),'b:');
  plot(hax, [1,numel(data.tFlip)],+dr+AmplitudeCenterI([end,end]),'b:');
  hold(hax,'off')
  set(hax, 'CLim', [-maxA, maxA]);
  title(hax, 'abs');
  ylabel(hax, 'read in samples');
  xlabel(hax, 'excitation flip angle in steps');
  set(hax, 'YDir', 'normal');

  hax = subplot(2,3,4, 'Parent', data.fh);
  plot(hax, ...
    data.tFlip*1e6, real([data.Amplitude;data.AmplitudeCenter])*1e9, ...
    data.tFlipInterp*1e6, real([data.AmplitudeInterp,data.AmplitudeCenterInterp])*1e9, ...
    data.tFlip90Center*1e6, real(data.AmplitudeCenterInterp(data.maxICenter))*1e9, 'x', ...
    data.tFlip90*1e6, real(data.AmplitudeInterp(data.maxI))*1e9, 'x', ...
    data.tFlip180Center*1e6, real(data.AmplitudeCenterInterp(data.ZCICenter))*1e9, 'x', ...
    data.tFlip180*1e6, real(data.AmplitudeInterp(data.ZCI))*1e9, 'x'...
    );
  ylim(hax, [-abs(data.maxAmpCenter),abs(data.maxAmpCenter)]*1.1*1e9);
  grid(hax, 'on');
  title(hax, {sprintf('tFlip90Center = %.1f %cs', data.tFlip90Center*1e6, 181), ...
    sprintf('tFlip90 = %.1f %cs', data.tFlip90*1e6, 181)});
  xlabel(hax, sprintf('pulse length in %cs', 181));
  ylabel(hax, 'amplitude in nT');
  hax = subplot(2,3,5, 'Parent', data.fh);
  plot(hax, data.tFlip*1e6, imag([data.Amplitude;data.AmplitudeCenter])*1e9);
  ylim(hax, [-abs(data.maxAmpCenter),abs(data.maxAmpCenter)]*1.1*1e9);
  grid(hax, 'on');
  title(hax, {sprintf('tFlip180Center = %.1f %cs', data.tFlip180Center*1e6, 181), ...
    sprintf('tFlip180 = %.1f %cs', data.tFlip180*1e6, 181)});
  xlabel(hax, sprintf('pulse length in %cs', 181));
  hax = subplot(2,3,6, 'Parent', data.fh);
  plot(hax, data.tFlip*1e6, abs([data.Amplitude;data.AmplitudeCenter])*1e9, ...
    data.tFlipInterp*1e6, abs([data.AmplitudeInterp,data.AmplitudeCenterInterp])*1e9)
  ylim(hax, [-abs(data.maxAmpCenter),abs(data.maxAmpCenter)]*1.1*1e9);
  grid(hax, 'on');
  title(hax, {['tFlipRatioCenter = ' num2str(data.tFlipRatioCenter)], ['tFlipRatio = ' num2str(data.tFlipRatio)]});
  xlabel(hax, sprintf('pulse length in %cs', 181));
end


end


function restoreHWsettings(HW, iDevice, paUout2Amplitude, paUout2AmplitudeX, gammaDef)
%% function that restores the initial HW settings on exiting the main function

HW.TX(iDevice).PaUout2Amplitude = paUout2Amplitude;

if nargin > 3
  % settings for secondary nucleus
  HW.TX(iDevice).PaUout2AmplitudeX = paUout2AmplitudeX;
  HW.GammaDef = gammaDef;
end

end
