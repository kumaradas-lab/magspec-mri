function [HW, mySave, Seq] = Find_CoilQuality(HW, mySave, Seq)
%% Find parameters that describe the quality of the TX coil
%
%       [HW, mySave, Seq] = Find_CoilQuality(HW, mySave, Seq)
%
% INPUT:
%
%   HW
%       HW object.
%
%   mySave
%       mySave object or structure (optional for storing B0 amplitude
%       information).
%
%   Seq
%       Structure with the following optional fields. All fields are also passed
%       to the function "sequence_RecoveryCPMG".
%
%     getQfactor
%         Boolean value.
%         If true, the Q factor of the rising and falling edges of the TX pulse
%         are determined. For this, a 90 degrees flip angle is determined by
%         varying the pulse duration (keeping the driving voltage constant).
%         Additionally, a small flip angle is determined using the same driving
%         voltage. Finally, the duration and phase of a refocussing pulse is
%         determined at the same driving voltage using a variant of a CPMG
%         experiment.
%         If false, the (apparent) coil efficiency is determined for the 90
%         degrees excitation pulse and the apparent flip angle (and phase) of
%         the 180 degrees refocussing pulse is determined at that (apparent)
%         coil efficiency. For this, a 90 degrees flip angle is determined by
%         varying the driving voltage (keeping the pulse duration constant).
%         Additionally, the duration and phase of a refocussing pulse is
%         determined varying the driving voltage at a set refocusing pulse
%         duration using a variant of a CPMG experiment.
%         (Default: true)
%
% OUTPUT:
%
%   HW
%       HW object.
%
%   mySave
%       mySave object or structure (optional for storing B0 amplitude
%       information).
%
%   Seq
%       Same as the input structure Seq with some additional fields among which:
%
%     Q
%         Structure with results of the coil quality measurement.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2024-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% input check
if nargin < 1
  error('PD:Find_CoilQuality:NoInput', 'At least input argument HW is mandatory.');
end
if nargin < 2
  mySave = [];
end
if nargin < 3
  Seq = struct();
end

%% default parameters
if isemptyfield(Seq, 'getQfactor')
  % true:   measure QFactors of rising and falling ramp
  % false:  measure apparent coil efficiency for 90 degrees pulse and apparent
  %         flip angle and phase
  Seq.getQfactor = true;
end
if isemptyfield(Seq, 'plot')
  Seq.plot = 320;
end
 % pulse duration search settings
if isemptyfield(Seq, 'nStepsPulseAmp')
  % number of steps for the pulse search
  Seq.nStepsPulseAmp = 11;
end
if isemptyfield(Seq, 'minFactorPulseAmp')
  % minimum factor for the currently set pulse length (start loop here)
  Seq.minFactorPulseAmp = 0.9;
end
if isemptyfield(Seq, 'maxFactorPulseAmp')
  % maximum factor for the currently set pulse length (stop loop here)
  Seq.maxFactorPulseAmp = 1.1;
end

% excitation pulse settings
if isemptyfield(Seq, 'excitationPulse')
  % pulse shape function of 90 degree excitation pulse (for step 1)
  Seq.excitationPulse = @Pulse_Rect;
end
if isemptyfield(Seq, 'tExcitation')
  if ~isemptyfield(HW, {'RecoveryCPMG', 'tFlip90Def'})
    Seq.tExcitation = HW.RecoveryCPMG.tFlip90Def * Seq.excitationPulse(HW, 'Amp');
  elseif isempty(HW.RecoveryCPMG.tFlip180Def)
    Seq.tExcitation = HW.tFlip90Def * Seq.excitationPulse(HW, 'Amp');
  else
    Seq.tExcitation = HW.RecoveryCPMG.tFlip180Def/2 * Seq.excitationPulse(HW, 'Amp');
  end
end

% small flip angle settings
% if isemptyfield(Seq, 'smallAngle')
%   % small excitation angle in degrees
%   Seq.smallAngle = 30;
% end
if isemptyfield(Seq, 'tFlipSmall')
  % center duration for small flip angle (step 2)
  Seq.tFlipSmall = Seq.tExcitation / 5;
end


% solid echo pulse settings
if isemptyfield(Seq, 'solidEchoPulse')
  % pulse shape function of solid echo pulse (for step 3)
  Seq.solidEchoPulse = @Pulse_Rect_SolidEcho;
end
if isemptyfield(Seq, 'tauSolidEcho')
  Seq.tauSolidEcho = 10e-6;
end

if isemptyfield(Seq, 'tEcho')
  % echo time (for step 3)
  Seq.tEcho = 0.25e-3*2;
end
if isemptyfield(Seq, 'tAQEcho')
  % acquisition time around echo (for step 3)
  Seq.tAQEcho = 50e-6; %Seq.tRefocus * pi/2;
end
if isemptyfield(Seq, 'tEchoTrain')
  % echo train length (for step 3)
  Seq.tEchoTrain = 50e-3;
end
if isemptyfield(Seq, 'tRefocus')
  % refocusing pulse duration (for step 3)
  if ~isemptyfield(HW, {'RecoveryCPMG', 'tFlip180Def'})
    Seq.tRefocus = HW.RecoveryCPMG.tFlip180Def;
  else
    Seq.tRefocus = 2*Seq.tExcitation;
  end
end
if isemptyfield(Seq, 'refocusingAngle')
  if isemptyfield(HW.RecoveryCPMG, 'refocusingAngle')
    Seq.refocusingAngle = 180;
  else
    Seq.refocusingAngle = HW.RecoveryCPMG.refocusingAngle;
  end
end
if isemptyfield(Seq, 'refocusingPhase')
  if ~isfield(HW.RecoveryCPMG, 'refocusingPhase')
    Seq.refocusingPhase = [];
  else
    Seq.refocusingPhase = HW.RecoveryCPMG.refocusingPhase;
  end
end


% relaxation (T1) properties
if isemptyfield(Seq, 'tRelax')
  % relaxation time after CPMG Echo train in s
  Seq.tRelax = 2;
end

% acceptance criteria
if isemptyfield(Seq, 'relDevMax')
  % tolerance for relative deviation
  Seq.relDevMax = 0.1;
end

% FIXME: Support multiple devices?
iDevice = 1;


%% get current magnet frequency (or shift to target B0 amplitude)
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0, [], 1);  % find magnet frequency


%% step 1: 90 degrees pulse duration using FID amplitude
Seq.nEcho = 0;
Seq.fSample = 1e6;

Seq.nTau1SteadyState = 0;

Seq.nTau1 = 0;
Seq.FitT2 = 0;
Seq.RawData = 0;

Seq.SeqAverage.average = 1;

% slice settings
Seq.thicknessSlice = Inf;  % thickness of slice selected by continuous gradient in m

Seq.plotSeq = 1:3;

if ~Seq.getQfactor
  % Always use same raster of TX amplitudes
  HW.TX(iDevice).PaUout2Amplitude = HW.TX(iDevice).PaUout2AmplitudeEstimated;  % FIXME: Reset when done!
end

Seq.SubtractEmptyMeasurement = 0;  % Run the same measurement without sample to subtract it.

tEcho = Seq.tEcho;  % save and restore the original value for step 2
if Seq.nEcho == 0
  % acquire FID
  Seq.tAQFID = 2e-3;
  Seq.tEcho = 2*Seq.tAQFID;
end

% initialize data for plot and evaluation
meanFID = NaN(1, Seq.nStepsPulseAmp);
factorsFID = linspace(Seq.minFactorPulseAmp, Seq.maxFactorPulseAmp, Seq.nStepsPulseAmp);
fOffsetFIDs = zeros(1, Seq.nStepsPulseAmp);
fOffsetFIDsStd = zeros(1, Seq.nStepsPulseAmp);
newfLarmor = zeros(1, Seq.nStepsPulseAmp);

clear excitationPulse
excitationPulse.FlipAngle = pi/2*factorsFID;
% excitationPulse.MaxLength = HW.RecoveryCPMG.tFlip180Def;

if Seq.getQfactor
  % increment pulse duration and keep constant driving amplitude
  tFlipFID1 = Seq.tExcitation*factorsFID;
  excitationPulse.Bandwidth =  1./tFlipFID1 * Seq.excitationPulse(HW, 'Time');

  exP = Seq.excitationPulse(HW, 0, excitationPulse);
  excPaUout = exP.Amplitude(1)./HW.TX(iDevice).PaUout2Amplitude(HW.TX(iDevice).ChannelDef);

  tFlipFID1_orig = tFlipFID1;
  splitCenterIdx = ceil(numel(tFlipFID1)/2);
  tFlipFID1(1:2:end) = tFlipFID1_orig(1:splitCenterIdx);
  tFlipFID1(2:2:end) = round(tFlipFID1_orig(end:-1:splitCenterIdx+1) * HW.TX.fSample) / HW.TX.fSample;
else
  % increment driving voltage and keep constant pulse duration
  excitationPulse.Bandwidth =  1/(Seq.tExcitation) * Seq.excitationPulse(HW, 'Time');

  exP = Seq.excitationPulse(HW, 0, excitationPulse);
  excPaUout = exP.Amplitude./HW.TX(iDevice).PaUout2Amplitude(HW.TX(iDevice).ChannelDef);

  excPaUout_orig = excPaUout;
  splitCenterIdx = ceil(numel(excPaUout)/2);
  excPaUout(1:2:end) = excPaUout_orig(1:splitCenterIdx);
  excPaUout(2:2:end) = excPaUout_orig(end:-1:splitCenterIdx+1);
  maxPaUout = HW.TX(iDevice).Max.PaUoutCalibrated(HW.TX(iDevice).ChannelDef)*0.999;
  % excPaUout(excPaUout>maxPaUout) = maxPaUout;
  if any(excPaUout>maxPaUout)
    error('PD:Find_CoilQuality:VoltageTooHigh', ...
      'Driving voltage for 90%c pulse (%.3f V) higher than maximum (%.3f V).', ...
      176, max(excPaUout), maxPaUout);
  end
end

if ~isempty(Seq.plot) && ~ishghandle(Seq.plot)
  Seq.plot = figure(Seq.plot);
end

if ishghandle(Seq.plot)
  hKids = get(Seq.plot, 'Children');
  delete(hKids);
  hParent = Seq.plot;

  dataFID = NaN(1, 1);
  timeFID = NaN(1, 1);

  axExc(1) = subplot(2,1,1, 'Parent', hParent);
  plot(axExc(1), timeFID*1e3, abs(dataFID)*1e6, ...
    timeFID*1e3, real(dataFID)*1e6, ...
    timeFID*1e3, imag(dataFID)*1e6);
  hold(axExc(1), 'on');
  title(axExc(1), 'Acquired signal');
  ylabel(axExc(1), ['RX Amplitude in ' char(181) 'V']);
  xlabel(axExc(1), 'Time in ms');
  ylim(axExc(1), [0, Inf]);
  grid(axExc(1), 'on');

  axExc(2) = subplot(2,1,2, 'Parent', hParent);
  if Seq.getQfactor
    hl2 = plot(axExc(2), tFlipFID1*1e6, meanFID, 'x');
    title(axExc(2), 'p90');
    ylabel(axExc(2), ['RX Amplitude in ' char(181) 'V']);
    xlabel(axExc(2), sprintf('flip duration in %cs', char(181)));
    % ylim(ax(2), [0, Inf]);
    grid(axExc(2), 'on');
  else
    hl2 = plot(axExc(2), excPaUout, meanFID, 'x');
    title(axExc(2), 'p90');
    ylabel(axExc(2), ['RX Amplitude in ' char(181) 'V']);
    xlabel(axExc(2), 'TX Amplitude in V');
    % ylim(ax(2), [0, Inf]);
    grid(axExc(2), 'on');
  end

  if ishghandle(Seq.plot, 'figure')
    Seq.plot = figure(Seq.plot);
  end
end

% actual measurement FID
Seq.Reinitialize = 1;
for iStep = reshape([ones(1, length(factorsFID))*numel(factorsFID); 1:length(factorsFID)], 1, [])
  if Seq.getQfactor
    Seq.tExcitation = tFlipFID1(iStep);
    Seq.excitationAngle = HW.GammaDef/2/pi*360 * Seq.tExcitation ...
      * exP.Amplitude(1);
  else
    Seq.excitationAngle = HW.GammaDef/2/pi*360 * Seq.tExcitation ...
      * get_TX_Amplitude(HW, 'PaUout', excPaUout(iStep), 'Device', iDevice);
  end

  [data, SeqOutCoil] = sequence_RecoveryCPMG(HW, Seq);
  Seq.LoopCountStart = SeqOutCoil.LoopCountEnd;
  Seq.Reinitialize = 0;
  Seq.StartSequenceTime = [];

  % meanFID(iStep) = data.MeanEchoTau1PhaseCorrected;
  % meanFID(iStep) = mean(max(abs(data.SampleEchoTau1), [], 1), 2);
  if SeqOutCoil.nEcho > 2
    plot(axExc(1), data.MeanEchoTau1PhaseCorrectedTime*1e3, abs(data.MeanEchoTau1PhaseCorrected)*1e6, '-', ...
      data.MeanEchoTau1PhaseCorrectedTime*1e3, real(data.MeanEchoTau1PhaseCorrected)*1e6, ':', ...
      data.MeanEchoTau1PhaseCorrectedTime*1e3, imag(data.MeanEchoTau1PhaseCorrected)*1e6, '-.');
    meanFID(iStep) = mean(mean(abs(data.MeanEchoTau1PhaseCorrected), 1), 2);
  else
    dataEcho = data.data(~isnan(data.data))*data.Amplitude2Uin(1)/HW.RX.LNAGain;
    timeEcho = data.time_all(~isnan(data.data));
    plot(axExc(1), timeEcho*1e3, abs(dataEcho)*1e6, '-', ...
      timeEcho*1e3, real(dataEcho)*1e6, ':', ...
      timeEcho*1e3, imag(dataEcho)*1e6, '-.');
    % meanFID(iStep) = mean(mean(abs(data.SampleEchoTau1*data.Amplitude2Uin(1)/HW.RX.LNAGain), 1), 2);
    meanFID(iStep) = mean(mean(abs(dataEcho), 1), 2);

    % calculate weighted phase slope for selected samples
    [phaseDiffWeightedMean, phaseDiffWeightedStd] = ...
      get_MeanPhaseDiffWeighted(double(dataEcho), 1, 'omitnan');

    % FIXME: Will the indices for AQ.fSample always be correct here?
    fOffsetFIDs(iStep) = phaseDiffWeightedMean * SeqOutCoil.AQ(1).fSample(1,2) / (2*pi);
    fOffsetFIDsStd(iStep) = phaseDiffWeightedStd * SeqOutCoil.AQ(1).fSample(1,2) / (2*pi);

    newfLarmor(iStep) = HW.fLarmor - double(fOffsetFIDs(iStep)); % calculate new Larmor frequency
    if fOffsetFIDsStd(iStep) < HW.FindFrequencySweep.fOffsetFIDsStdMaxValue ...
        && newfLarmor(iStep) > 0
      % save the new frequency to mySave
      HW.B0 = newfLarmor(iStep) * 2*pi/HW.FindFrequencyGamma;  % calculate the new magnetic field strength
      if ~isemptyfield(HW.FindFrequencySweep, 'shiftB0') ...
          && HW.FindFrequencySweep.shiftB0
        if isempty(HW.B0Target) || ~isfinite(HW.B0Target)
          warning('PD:Find_CoilQuality:NoB0Target', ...
            ['HW.FindFrequencySweep.shiftB0 is set to true, but HW.B0Target is not set.\n', ...
            'B0 amplitude cannot be shifted without a target. Adjusting HW.fLarmor instead.']);
        else
          newB0shift = HW.Grad(iDevice).AmpOffset(4) + (HW.B0Target - HW.B0);
          newB0 = HW.B0Target;
          if abs(newB0shift) > HW.Grad(iDevice).HoldShimNormMax(4)*HW.Grad(iDevice).MaxAmp(4)
            desiredB0shift = newB0shift;
            newB0shift = sign(newB0shift) * HW.Grad(iDevice).HoldShimNormMax(4) * HW.Grad(iDevice).MaxAmp(4);
            newB0 = newB0 - (desiredB0shift - newB0shift);
            warning('PD:Find_CoilQuality:B0ShiftExceeded', ...
              ['The required B0 shift (%.1f %cT) is larger than the maximum shift (%.1f %cT).\n', ...
              'Consider setting a different HW.B0Target. Adjusting HW.fLarmor to stay in range.'], ...
              desiredB0shift*1e6, 181, HW.Grad(iDevice).HoldShimNormMax(4)*HW.Grad(iDevice).MaxAmp(4)*1e6, 181);
          end
          HW.Grad(iDevice).AmpOffset(4) = newB0shift;
          HW.B0 = newB0;
        end
      end
    else
      error('PD:Find_CoilQuality:UncertainFrequency', ...
        'The current Larmor frequency could not be determined with sufficient accuracy.');
    end
  end
  set(hl2, 'YData', abs(meanFID)*1e6);

  % % check interrupt
  % if mriDevice.system.settings.layout.guiType > gui.Type.NOGUI && ...
  %     ~isempty(mriDevice.displayInterface.hGUI) && ...
  %     ~get(mriDevice.displayInterface.hGUI.hButtonTeach, 'Value')
  %   error('Measurement interrupted by user');
  % end
end
hold(axExc(1), 'off');


if Seq.getQfactor
  % pulse duration has been incremented
  tFlip1Interp = linspace(min(tFlipFID1), max(tFlipFID1), 4*Seq.nStepsPulseAmp);
  % quadratic best fit
  A = [tFlipFID1.^2; tFlipFID1; ones(size(tFlipFID1))].' \ meanFID.';
  tFlip90 = -A(2)/(2*A(1));
  maxFID1Amp = [tFlip90^2, tFlip90, 1] * A;
  meanEchoInterp = A.' * [tFlip1Interp.^2; tFlip1Interp; ones(size(tFlip1Interp))];

  hold(axExc(2), 'on');
  plot(axExc(2), tFlip1Interp*1e6, abs(meanEchoInterp)*1e6, '-');
  plot(axExc(2), tFlip90*1e6, maxFID1Amp*1e6, '-o');
  hold(axExc(2), 'off');

  title(axExc(2), ['p90 @ ', num2str(tFlip90*1e6, 4), ' ', char(181), 's']);


  % calculate magnetic flux density of the transmit coil for a 90 degrees pulse
  % with the found pulse duration
  B1 = (1/(tFlip90*4))/(HW.GammaDef/2/pi);

  newPaUout2Amplitude = HW.TX(iDevice).PaUout2Amplitude;
  newPaUout2Amplitude(HW.TX(iDevice).ChannelDef) = B1./excPaUout;
  comment = sprintf('%s (p90 = %.3f %ss @ %.3f V @ %.6f MHz) from FID amplitude by %s', ...
    datestr(now, 'yyyy-mm-ddTHH:MM:SS'), tFlip90*1e6, char(181), ...
    excPaUout, HW.fLarmor/1e6, mfilename());
else
  % driving voltage has been incremented
  excPaUoutInterp = linspace(exP.Amplitude(1), exP.Amplitude(end), 4*Seq.nStepsPulseAmp)./HW.TX.PaUout2Amplitude(HW.TX.ChannelDef);
  % quadratic best fit
  A = [excPaUout.^2; excPaUout; ones(size(excPaUout))].' \ meanFID.';
  paUoutFlip90 = -A(2)/(2*A(1));
  maxFID1Amp = [paUoutFlip90^2, paUoutFlip90, 1] * A;
  meanEchoInterp = A.' * [excPaUoutInterp.^2; excPaUoutInterp; ones(size(excPaUoutInterp))];

  hold(axExc(2), 'on');
  plot(axExc(2), excPaUoutInterp, abs(meanEchoInterp)*1e6, '-');
  plot(axExc(2), paUoutFlip90, maxFID1Amp*1e6, '-o');
  hold(axExc(2), 'off');

  title(axExc(2), ['p90 @ ', num2str(paUoutFlip90, 4), ' V']);


  % calculate magnetic flux density of the transmit coil for a 90 degrees pulse
  % with the selected pulse duration
  B1 = (1/(SeqOutCoil.tExcitation*4))/(HW.GammaDef/2/pi);

  % calculate coil efficiency
  newPaUout2Amplitude = HW.TX(iDevice).PaUout2Amplitude;
  newPaUout2Amplitude(HW.TX(iDevice).ChannelDef) = B1./paUoutFlip90;
  comment = sprintf('%s (p90 = %.3f %cs @ %.3f V @ %.6f MHz) from FID amplitude by %s', ...
    datestr(now, 'yyyy-mm-ddTHH:MM:SS'), Seq.tExcitation*1e6, char(181), ...
    paUoutFlip90, HW.fLarmor/1e6, mfilename());
end

if isempty(HW.TX(iDevice).CoilName)
  newCalLine = sprintf(...
    'HW.TX(%d).PaUout2Amplitude = [%.6f, %.6f]*1e-6;  %% %s\n', ...
    iDevice, newPaUout2Amplitude*1e6, comment);
else
  newCalLine = sprintf(...
    'if strcmp(HW.TX(%d).CoilName, ''%s''),  HW.TX(%d).PaUout2Amplitude = [%.6f, %.6f]*1e-6;  end  %% %s\n', ...
    iDevice, HW.TX.CoilName, iDevice, newPaUout2Amplitude*1e6, comment);
end

% some cross checks to ensure the value can be trusted.
savePulseFile = true;

if A(1) > 0
  warning('PD:sequence_PulseDurationCPMG:NoMax', 'Couldn''t find a maximum.');
  savePulseFile = false;
end

if Seq.getQfactor
  if tFlip90 <= tFlip1Interp(1) || tFlip90 >= tFlip1Interp(end)
    warning('PD:sequence_PulseDurationCPMG:OutOfRange', 'p90 not in range of sweep');
    savePulseFile = false;
  end
else
  if paUoutFlip90 <= excPaUoutInterp(1) || paUoutFlip90 >= excPaUoutInterp(end)
    warning('PD:sequence_PulseDurationCPMG:OutOfRange', 'p90 not in range of sweep');
    savePulseFile = false;
  end
end
SeqOutCoil.relDev = newPaUout2Amplitude(HW.TX.ChannelDef) / HW.TX.PaUout2AmplitudeEstimated(HW.TX.ChannelDef);
% FIXME: The relative deviation tolerance is multiplied by a factor of 2 to
%        work around an issue caused by the expected apparent coil efficiency
%        determined by CPMG echo trains or FIDs respectively.
tolFactor = 2;
if SeqOutCoil.relDev > (1+tolFactor*Seq.relDevMax) ...
    || SeqOutCoil.relDev < (1-tolFactor*Seq.relDevMax)
  warning('PD:sequence_PulseDurationCPMG:OutOfRange', ...
    'Relative deviation to calibrated coil efficiency too high (%.1f %% > %.1f %%).', ...
    (SeqOutCoil.relDev-1)*100, tolFactor*Seq.relDevMax*100);
  savePulseFile = false;
end

if Seq.getQfactor
  paUoutFlip90 = excPaUout;
else
  tFlip90 = SeqOutCoil.tExcitation;
end

fprintf('90%s pulse duration: %.3f %ss @ %.3f V\n', char(176), ...
  SeqOutCoil.tExcitation*1e6, char(181), paUoutFlip90);

% FIXME:  Do we want to write the coil efficiency already to the file when
%         Seq.getQfactor is true?

if ~isempty(HW.TX(iDevice).PaUout2AmplitudePath) ...
    && (isemptyfield(mySave, 'DummySerial') ...
        || mySave.DummySerial(min(iDevice, numel(mySave.DummySerial))) <= 0)
  addFirstLine = ~exist(HW.TX(iDevice).PaUout2AmplitudePath, 'file');
  fid = fopen(HW.TX(iDevice).PaUout2AmplitudePath, 'a+');
  fid_protect = onCleanup(@() fclose(fid));
  if addFirstLine
    fprintf(fid, '%s\n', '% factor from voltage amplitude at the coil input to B1+ field strength in T/V');
  end
  if savePulseFile
    fwrite(fid, newCalLine);
    fprintf('\nA new line was added to the following file:\n%s\n%s\n', ...
      HW.TX(iDevice).PaUout2AmplitudePath, newCalLine);
  else
    fwrite(fid, ['% ', newCalLine]);  % add line as comment
    fprintf('\n');
    warnStr = ['Determination of pulse amplitude unsuccessful!\n', ...
      'If you''d like to use the uncertain best guess value anyway, ', ...
      'please manually append or un-comment the following line in your PaUout2AmplitudeCal.m file:\n%s\n'];
    warning('PD:sequence_PulseDurationCPMG:FID', warnStr, newCalLine);
  end
  delete(fid_protect);
elseif savePulseFile
  fprintf('\n');
  fprintf('\nPlease add the following line to your LoadMySystem.m file:\n%s\n', ...
    newCalLine);
end

if ~savePulseFile
  % if haveGUI
  %   this.SaveDisplay(hParent);
  % end
  %
  % this.result.failure = 'coil efficiency';
  % nsq.dailyCheck.result = this.result.failure;
  % this.result.lastWarnStr = 'detection of 90 degree pulse amplitude unsuccessful';
  % this.result.lastWarnId = 'DailyCheck_solidEcho:Failed90DegAmp';
  % % disp(getReport(ME));
  % if haveGUI
  %   herr = errordlg('Failed to detect coil efficiency.', 'Error', 'modal');
  %   w = warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
  %   warning('off', 'MATLAB:ui:javaframe:PropertyToBeRemoved');
  %   jFrame = get(herr, 'JavaFrame');
  %   warning(w);
  %   jFrame.setFigureIcon(jIcon);
  %   uiwait(herr);
  % end
  return;
elseif ~Seq.getQfactor
  % temporarily set coil efficiency for correct 90° pulses
  HW.TX(iDevice).PaUout2Amplitude = newPaUout2Amplitude;
end

%% get current magnet frequency (or shift to target B0 amplitude)
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0, [], 1);  % find magnet frequency


%% step 2 (optional): small angle pulse duration using FID amplitude
% Search pulse duration for which the rf pulse doesn't reach the end amplitude
if Seq.getQfactor
  % initialize data for plot and evaluation
  ampFID2 = NaN(1, Seq.nStepsPulseAmp);

  % increment pulse duration and keep constant driving amplitude
  clear excitationPulse2
  tFlipFID2 = round(Seq.tFlipSmall * factorsFID * HW.TX.fSample) / HW.TX.fSample;
  excitationPulse2.Bandwidth =  1./tFlipFID2 * Seq.excitationPulse(HW, 'Time');
  excitationPulse2.FlipAngle = tFlipFID2 / Seq.tExcitation * pi/2;

  exP2 = Seq.excitationPulse(HW, 0, excitationPulse2);

  tFlipFID2_orig = tFlipFID2;
  splitCenterIdx = ceil(numel(tFlipFID2)/2);
  tFlipFID2(1:2:end) = tFlipFID2_orig(1:splitCenterIdx);
  tFlipFID2(2:2:end) = tFlipFID2_orig(end:-1:splitCenterIdx+1);

  if ~isempty(Seq.plot)
    if ishghandle(Seq.plot)
      if isgraphics(Seq.plot, 'figure')
        Seq.plot = get(Seq.plot, 'Number') + 1;
      end
    else
      Seq.plot = Seq.plot + 1;
    end
    if ~isgraphics(Seq.plot) || ishghandle(Seq.plot, 'figure')
      Seq.plot = figure(Seq.plot);
    end
  end

  if ishghandle(Seq.plot)
    hKids = get(Seq.plot, 'Children');
    delete(hKids);
    hParent = Seq.plot;

    dataFID = NaN(1, 1);
    timeFID = NaN(1, 1);

    axExc2(1) = subplot(2,1,1, 'Parent', hParent);
    plot(axExc2(1), timeFID*1e3, abs(dataFID)*1e6, ...
      timeFID*1e3, real(dataFID)*1e6, ...
      timeFID*1e3, imag(dataFID)*1e6);
    hold(axExc2(1), 'on');
    title(axExc2(1), 'Acquired signal');
    ylabel(axExc2(1), ['RX Amplitude in ' char(181) 'V']);
    xlabel(axExc2(1), 'Time in ms');
    ylim(axExc2(1), [0, Inf]);
    grid(axExc2(1), 'on');

    axExc2(2) = subplot(2,1,2, 'Parent', hParent);
    hl2 = plot(axExc2(2), tFlipFID2*1e6, ampFID2*1e6, 'x');
    title(axExc2(2), 'pSmall');
    ylabel(axExc2(2), ['RX Amplitude in ' char(181) 'V']);
    xlabel(axExc2(2), sprintf('flip duration in %cs', char(181)));
    % ylim(ax(2), [0, Inf]);
    grid(axExc2(2), 'on');
  end

  % actual measurement FID
  Seq.Reinitialize = 1;
  for iStep = 1:length(factorsFID)
    Seq.tExcitation = tFlipFID2(iStep);
    Seq.excitationAngle = HW.GammaDef/2/pi*360 * Seq.tExcitation ...
      * exP2.Amplitude(iStep);

    [data, SeqOutCoil] = sequence_RecoveryCPMG(HW, Seq);
    Seq.LoopCountStart = SeqOutCoil.LoopCountEnd;
    Seq.Reinitialize = 0;
    Seq.StartSequenceTime = [];

    % meanFID(iStep) = data.MeanEchoTau1PhaseCorrected;
    % meanFID(iStep) = mean(max(abs(data.SampleEchoTau1), [], 1), 2);
    if SeqOutCoil.nEcho > 2
      plot(axExc2(1), data.MeanEchoTau1PhaseCorrectedTime*1e3, abs(data.MeanEchoTau1PhaseCorrected)*1e6, '-', ...
        data.MeanEchoTau1PhaseCorrectedTime*1e3, real(data.MeanEchoTau1PhaseCorrected)*1e6, ':', ...
        data.MeanEchoTau1PhaseCorrectedTime*1e3, imag(data.MeanEchoTau1PhaseCorrected)*1e6, '-.');
      ampFID2(iStep) = mean(mean(abs(data.MeanEchoTau1PhaseCorrected), 1), 2);
    else
      dataEcho = (data.data(:))*data.Amplitude2Uin(1)/HW.RX.LNAGain;
      timeEcho = data.time_of_tRep(:);
      plot(axExc2(1), timeEcho*1e3, abs(dataEcho)*1e6, '-', ...
        timeEcho*1e3, real(dataEcho)*1e6, ':', ...
        timeEcho*1e3, imag(dataEcho)*1e6, '-.');
      ampFID2(iStep) = mean(mean(abs(data.SampleEchoTau1*data.Amplitude2Uin(1)/HW.RX.LNAGain), 1), 2);
    end
    set(hl2, 'YData', abs(ampFID2)*1e6);

    % % check interrupt
    % if mriDevice.system.settings.layout.guiType > gui.Type.NOGUI && ...
    %     ~isempty(mriDevice.displayInterface.hGUI) && ...
    %     ~get(mriDevice.displayInterface.hGUI.hButtonTeach, 'Value')
    %   error('Measurement interrupted by user');
    % end

  end
  hold(axExc2(1), 'off');

  smallAngles = asind(ampFID2/maxFID1Amp);  % flip angle in degrees for each step

  % fit expected pulse amplitude
  tFlip2Interp = linspace(min(tFlipFID2), max(tFlipFID2), 4*Seq.nStepsPulseAmp);

%   % quadratic best fit
%   A2 = [tFlipFID2.^2; tFlipFID2; ones(size(tFlipFID2))].' \ meanFID2.';
%   tFlipSmall = (-A2(2) + sqrt(A2(2)^2 -4*A2(1)*(A2(3)-maxFID1Amp*sind(Seq.smallAngle))))/(2*A2(1));
%   maxFID2Amp = [tFlipSmall^2, tFlipSmall, 1] * A2;
%   meanFID2Interp = [tFlip2Interp.^2; tFlip2Interp; ones(size(tFlip2Interp))].' * A2;

  % linear regression
  A2 = [tFlipFID2; ones(size(tFlipFID2))].' \ ampFID2.';
  tFlipSmall = Seq.tFlipSmall;
  meanFID2Amp = [tFlipSmall, 1] * A2;
  Seq.smallAngle = asind(meanFID2Amp/maxFID1Amp);
  meanFID2Interp = [tFlip2Interp; ones(size(tFlip2Interp))].' * A2;

  hold(axExc2(2), 'on');
  plot(axExc2(2), tFlip2Interp*1e6, abs(meanFID2Interp)*1e6, '-');
  plot(axExc2(2), tFlipSmall*1e6, meanFID2Amp*1e6, '-o');
  hold(axExc2(2), 'off');

  title(axExc2(2), sprintf('p%.2f @ %.2f %cs', Seq.smallAngle, tFlipSmall*1e6, char(181)));
end


%% get current magnet frequency (or shift to target B0 amplitude)
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0, [], 1);  % find magnet frequency


%% step 3: CPMG echo train for 180 degrees pulse amplitude and phase

% use same settings as above as far as possible
Seq.excitationPulse = Seq.solidEchoPulse;  % pulse shape function for excitation pulse
Seq.nEcho = []; 1;
Seq.SeqAverage.average = 1;
Seq.tAQFID = -1;

Seq.tEcho = tEcho;  % restore value before step 1
if Seq.getQfactor
  % FIXME: calculate excitation angle
  Seq.tExcitation = round(tFlip90 * HW.TX.fSample) / HW.TX.fSample;
  Seq.excitationAngle = HW.GammaDef/2/pi*360 * Seq.tExcitation ...
    * exP.Amplitude(1);
  Seq.refocusingAngle = HW.RecoveryCPMG.refocusingAngle;
  Seq.tRefocus = round(Seq.tExcitation/90 * Seq.refocusingAngle * HW.TX.fSample) / HW.TX.fSample;
  Seq.refocusingAngle = HW.GammaDef/2/pi*360 * Seq.tRefocus * exP.Amplitude(1);
else
  Seq.excitationAngle = 90;
end

% Seq.thicknessSlice = 0.05;

Seq.refocusingPhaseOffset = 'alternateIncrement';  % to accumulate errors from inversion angle and phase

Seq.PostProcessSequenceLocal = false;  % skip post-processing with fit_exp

if ~isempty(Seq.plot)
  if ishghandle(Seq.plot)
    if isgraphics(Seq.plot, 'figure')
      Seq.plot = get(Seq.plot, 'Number') + 1;
    end
  else
    Seq.plot = Seq.plot + 1;
  end
  if ~isgraphics(Seq.plot) || ishghandle(Seq.plot, 'figure')
    Seq.plot = figure(Seq.plot);
  end
end


% actual measurement (inversion pulse)
Seq.Reinitialize = 1;
evalPhasePrevious = [];
for iPhase = 1:4
  [data, SeqOutInv] = sequence_RecoveryCPMG(HW, Seq);
  Seq.LoopCountStart = SeqOutInv.LoopCountEnd;
  Seq.Reinitialize = 0;
  Seq.StartSequenceTime = [];

  % [SeqOutInv, data_1D] = get_data_1D(SeqOutInv, data);
  % plot_data_1D(HW, data_1D);

  %% evaluate phases
  szData = size(data.data);
  dataAvg = reshape(data.data, szData(1), szData(2), szData(3)/Seq.SeqAverage.average, Seq.SeqAverage.average);
  dataAvg = mean(dataAvg, 4);
  evalPhase.phaseEchoes = mean(unwrap(angle(dataAvg), 1), 1);

  % FIXME: averages and phase cycle!!!
  evalPhase.phaseEchoesOdd = unwrap(reshape(evalPhase.phaseEchoes(:,:,2:2:end-1), [], 1));
  evalPhase.phaseEchoesEven = unwrap(reshape(evalPhase.phaseEchoes(:,:,3:2:end-1), [], 1));

  % linear fit through phase of echoes
  AOdd = [(1:2:(size(evalPhase.phaseEchoes, 3)-2)).', ones(size(evalPhase.phaseEchoesOdd))];
  evalPhase.fitOdd = AOdd\evalPhase.phaseEchoesOdd;

  AEven = [(2:2:(size(evalPhase.phaseEchoes, 3)-2)).', ones(size(evalPhase.phaseEchoesEven))];
  evalPhase.fitEven = AEven\evalPhase.phaseEchoesEven;

  fitEchoes = [[0; size(evalPhase.phaseEchoes, 3)-2], [1; 1]];

  ampRefocus = deg2rad(SeqOutInv.refocusingAngle) / SeqOutInv.tRefocus / HW.GammaDef;
  paUoutRefocus = ampRefocus / HW.TX(iDevice).PaUout2Amplitude(HW.TX(iDevice).ChannelDef);
  if isgraphics(Seq.plot)
    if iPhase == 1
      hKids = get(Seq.plot, 'Children');
      delete(hKids);
      hParent = Seq.plot;

      haxes = axes('Parent', hParent);
      grid(haxes, 'on');
      xlabel(haxes, 'number of echo');
      ylabel(haxes, 'phase of echo in radians');
    else
      hLines = findobj(haxes, '-depth', 1, 'Type', 'line');
      delete(hLines);
    end

    % hf = figure(40+iPhase); clf(hf)
    % haxes = axes(hf);
    hold(haxes, 'on');
    plot(haxes, AOdd(:,1), evalPhase.phaseEchoesOdd, 'xb');
    plot(haxes, AEven(:,1), evalPhase.phaseEchoesEven, 'xr');
    plot(haxes, fitEchoes(:,1), fitEchoes*evalPhase.fitOdd, '--b', 'LineWidth', 2);
    plot(haxes, fitEchoes(:,1), fitEchoes*evalPhase.fitEven, '--r', 'LineWidth', 2);
    hold(haxes, 'off');
    title(haxes, sprintf('Echo Phase @ p90=%.3f V, p180=%.3f V (= %.2f%c), phase=%.2f%c', ...
      paUoutFlip90, paUoutRefocus, SeqOutInv.refocusingAngle, 176, ...
      SeqOutInv.refocusingPhase-90, 176))
    if iPhase == 1
      % FIXME: hgsave and hgload don't restore legends correctly
      % legend(haxes, {'odd echoes', 'even echoes'});
    end
  end

  evalPhase.phaseErrorPerEcho = rad2deg(evalPhase.fitOdd(1)-evalPhase.fitEven(1));  % additional error in phase of echo per echo

  % % check interrupt
  % if mriDevice.system.settings.layout.guiType > gui.Type.NOGUI && ...
  %     ~isempty(mriDevice.displayInterface.hGUI) && ...
  %     ~get(mriDevice.displayInterface.hGUI.hButtonTeach, 'Value')
  %   error('Measurement interrupted by user');
  % end

  %% incrementally adjust amplitude or duration and phase of inversion pulses
  refocusingAngleCurrent = SeqOutInv.refocusingAngle;
  if iPhase == 1
    % repeat measurement with lower amplitude of inversion pulses
    Seq.refocusingAngle = SeqOutInv.refocusingAngle - 10;
  elseif iPhase == 2
    % calculate phase sensitivity and correct for next step
    AmpDegPerPhase = (SeqPrevious.refocusingAngle - SeqOutInv.refocusingAngle) ...
      / (evalPhase.phaseErrorPerEcho - evalPhasePrevious.phaseErrorPerEcho);
    Seq.refocusingAngle = Seq.refocusingAngle + ...
      AmpDegPerPhase * evalPhase.phaseErrorPerEcho;
  else
    % correct for next step
    Seq.refocusingAngle = Seq.refocusingAngle + ...
      AmpDegPerPhase * evalPhase.phaseErrorPerEcho;
  end
  if Seq.getQfactor
    Seq.tRefocus = round(SeqOutInv.tRefocus * Seq.refocusingAngle / refocusingAngleCurrent * HW.TX.fSample) / HW.TX.fSample;
    Seq.refocusingAngle = HW.GammaDef/2/pi*360 * Seq.tRefocus * exP.Amplitude(1);
  end
  Seq.refocusingPhase = mod(SeqOutInv.refocusingPhase + rad2deg(evalPhase.fitOdd(2)-evalPhase.fitEven(2))/2, 180);


  if ~isfinite(Seq.refocusingAngle) || ~isfinite(Seq.refocusingPhase)
    error('PD:Find_CoilQuality:nonFiniteResult', ...
      'Error calculating apparent refocusing angle or phase for next step.');
  end

  SeqPrevious = SeqOutInv;
  evalPhasePrevious = evalPhase;

end
% phase offset of acquired signal to phase of excitation pulse
Seq.Q.RXphaseOffset = mean([evalPhase.fitOdd(2), evalPhase.fitEven(2)]);
% "apparent" pulse properties (for driving signal)
Seq.Q.apparentRefocusingAngle = Seq.refocusingAngle;
Seq.Q.apparentRefocusingPhase = Seq.refocusingPhase;


%% post-processing of measurements
if Seq.getQfactor
  % calculate time constants from measured durations

  % Assumptions:
  % * "Time profile" while driving coil is strictly exponential.
  % * Flip angle due to the coil ringing is linearly proportional to the
  %   B1-amplitude at the end of the rf pulse.

  % flip angles of the three pulses in degrees
  p1 = Seq.smallAngle;
  p2 = 90;
  p3 = 180;

  % duration of coil driving in seconds
  t1 = tFlipSmall;  % pulse duration of shortest pulse
  t2 = tFlip90;  % pulse duration of medium pulse
  t3 = Seq.tRefocus;  % pulse duration of longest pulse

%   M = p2/p1;  % ratio of flip angles for which one hasn't reached the final amplitude
  M = p3/p1;  % ratio of flip angles for which one hasn't reached the final amplitude
  % ratio of flip angles that are used as "reference" for the minimization
  % (both should(?) have reached the final amplitude)
  N = p3/p2;

%   A = (N*t2-t3)/(N-1);

  clear tau_rise
%   f_min = @(tau_rise) ((A - ...
%     (M*(t1-tau_rise*(1-exp(-t1/tau_rise))) - (t2-tau_rise*(1-exp(-t2/tau_rise)))) ...
%     / (1-M*(1-exp(-t1/tau_rise)))) ...
%     - tau_rise)^2;
%   f_min = @(tau_rise) (...
%     ((t3 - tau_rise*(1-exp(-t3/tau_rise))) - N*(t2 - tau_rise*(1-exp(-t2/tau_rise)))) ...
%     / (N*(1-exp(-t2/tau_rise)) - 1) ...
%     - ...
%     (M*(t1-tau_rise*(1-exp(-t1/tau_rise))) - (t2-tau_rise*(1-exp(-t2/tau_rise)))) ...
%     / (1 - M*(1-exp(-t1/tau_rise)))...
%     )^2;
%   f_min = @(tau_rise) (...
%     ((t3 - tau_rise*(1-exp(-t3/tau_rise))) - N*(t2 - tau_rise*(1-exp(-t2/tau_rise)))) ...
%     / (N*(1-exp(-t2/tau_rise)) - (1-exp(-t3/tau_rise))) ...
%     - ...
%     (M*(t1-tau_rise*(1-exp(-t1/tau_rise))) - (t2-tau_rise*(1-exp(-t2/tau_rise)))) ...
%     / ((1-exp(-t2/tau_rise)) - M*(1-exp(-t1/tau_rise)))...
%     )^2;

  tau_fall_fcn = @(tau_rise, ti1, ti2, flip_ratio) (...
    ((ti1 - tau_rise*(1-exp(-ti1/tau_rise))) - flip_ratio.*(ti2 - tau_rise*(1-exp(-ti2/tau_rise)))) ...
    ./ (flip_ratio.*(1-exp(-ti2/tau_rise)) - (1-exp(-ti1/tau_rise))) ...
    );

  Seq.useAverageSmallAngle = false;  % FIXME: Chose setting and remove alternative?
  if Seq.useAverageSmallAngle
    f_min = @(tau_rise) (...
      ((t3 - tau_rise*(1-exp(-t3/tau_rise))) - N*(t2 - tau_rise*(1-exp(-t2/tau_rise)))) ...
      / (N*(1-exp(-t2/tau_rise)) - (1-exp(-t3/tau_rise))) ...
      - ...
      (M*(t1-tau_rise*(1-exp(-t1/tau_rise))) - (t3-tau_rise*(1-exp(-t3/tau_rise)))) ...
      / ((1-exp(-t3/tau_rise)) - M*(1-exp(-t1/tau_rise)))...
      )^2;
  else
    f_min = @(tau_rise) (...
      numel(smallAngles) * (tau_fall_fcn(tau_rise, Seq.tRefocus, tFlip90, N))^2 ...
      - 2 * tau_fall_fcn(tau_rise, Seq.tRefocus, tFlip90, N) ...
        * sum(tau_fall_fcn(tau_rise, Seq.tRefocus, tFlipFID2, p3./smallAngles)) ...
      + sum((tau_fall_fcn(tau_rise, Seq.tRefocus, tFlipFID2, p3./smallAngles)).^2));
  end
  tau_rise_estimated = tFlipSmall;
  options.TolFun = 1e-9^2;
  options.TolX = 1e-9;
  [tau_rise, Seq.Q.resid, Seq.Q.exitflag, Seq.Q.output] ...
    = fminsearch(f_min, tau_rise_estimated, options);
  Seq.Q.tau_rise = tau_rise;

  if Seq.useAverageSmallAngle
    Seq.Q.tau_fall = ...
      mean(tau_fall_fcn(tau_rise, Seq.tRefocus, [tFlip90, tFlipFID2], p3./[p2, smallAngles]));
  else
    % tau_fall = + (M*(t1-tau_rise*(1-exp(-t1/tau_rise))) - (t2-tau_rise*(1-exp(-t2/tau_rise)))) / (1-M*(1-exp(-t1/tau_rise)))
    Seq.Q.tau_fall = (M*(t1-tau_rise*(1-exp(-t1/tau_rise))) - (t3-tau_rise*(1-exp(-t3/tau_rise)))) ...
      / ((1-exp(-t3/tau_rise)) - M*(1-exp(-t1/tau_rise)));
  end

  B1 = (p3-p2)/360 / (HW.GammaDef/(2*pi) * (t3-t2));

  Seq.Q.PaUout2AmplitudeReal = B1/excPaUout;

  % tau = Q/pi/f0
  Seq.Q.Q_rise = Seq.Q.tau_rise*pi*HW.fLarmor;
  Seq.Q.Q_fall = Seq.Q.tau_fall*pi*HW.fLarmor;

else
  % apparent 180 degrees flip pulse data

  % save calibration results to file
  HW.RecoveryCPMG.refocusingAngle = Seq.refocusingAngle;
  HW.RecoveryCPMG.refocusingPhase = Seq.refocusingPhase;

  newCalLines = {...
    sprintf('HW.RecoveryCPMG.refocusingAngle = %.3f;  %% %s (p180 = %.3f %cs @ %.3f V @%.6f MHz) from CPMG Echo train by %s', ...
    Seq.refocusingAngle, datestr(now, 'yyyy-mm-ddTHH:MM:SS'), SeqOutInv.tRefocus*1e6, 181, ...
    paUoutRefocus, HW.fLarmor/1e6, mfilename()), ...
    sprintf('HW.RecoveryCPMG.refocusingPhase = %.3f;  %% %s phase of refocusing pulse in degrees from CPMG Echo train by %s', ...
    Seq.refocusingPhase, datestr(now, 'yyyy-mm-ddTHH:MM:SS'), mfilename())};

  % some consistency checks
  savePulseFile = true;

  if HW.RecoveryCPMG.refocusingAngle < 145 || HW.RecoveryCPMG.refocusingAngle > 215
    warning('PD:sequence_PulseDurationCPMG:EfficiencyOutOfRange', ...
      'Deviation between efficiency for 90 degree and 180 degree pulses is large.');
    savePulseFile = false;
  end

  if HW.RecoveryCPMG.refocusingPhase < 75 || HW.RecoveryCPMG.refocusingPhase > 105
    warning('PD:sequence_PulseDurationCPMG:PhaseOutOfRange', ...
      'Deviation between effective phase of 90 degree and 180 degree pulses is large.');
    savePulseFile = false;
  end

  if ~isempty(HW.TX(iDevice).PaUout2AmplitudePath) ...
      && (isemptyfield(mySave, 'DummySerial') ...
          || mySave.DummySerial(min(iDevice, numel(mySave.DummySerial))) <= 0)
    fid = fopen(HW.TX(iDevice).PaUout2AmplitudePath, 'a+');
    fid_protect = onCleanup(@() fclose(fid));
    if savePulseFile
      fprintf(fid, '%s\n', newCalLines{:});
    else
      fprintf(fid, '%% %s\n', newCalLines{:});
    end
    delete(fid_protect);
    fprintf('\nNew lines were added to the following file:\n%s\n%s\n', ...
      HW.TX(iDevice).PaUout2AmplitudePath, sprintf('%s\n', newCalLines{:}));
  elseif savePulseFile
    fprintf('\n');
    fprintf('\nPlease add the following lines to your LoadMySystem.m file:\n%s\n', ...
      sprintf('%s\n', newCalLines{:}));
  end
  if ~savePulseFile
    fprintf('\n');
    warnStr = ['Determination of inversion pulse parameters unsuccessful!\n', ...
      'If you''d like to use the uncertain best guess values anyway, ', ...
      'please manually append or un-comment the following lines in your PaUout2AmplitudeCal.m file:\n%s\n'];
    warning('PD:sequence_PulseDurationCPMG', warnStr, sprintf('%s\n', newCalLines{:}));

    % this.result.failure = 'coil efficiency';
    % nsq.dailyCheck.result = this.result.failure;
    % this.result.lastWarnStr = 'detection of 180 degree pulse amplitude or phase unsuccessful';
    % this.result.lastWarnId = 'DailyCheck_solidEcho:Failed180DegAmpPhase';
    % % disp(getReport(ME));
    % if haveGUI
    %   herr = errordlg('Failed to detect 180 degree pulse amplitude or phase.', 'Error', 'modal');
    %   w = warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    %   warning('off', 'MATLAB:ui:javaframe:PropertyToBeRemoved');
    %   jFrame = get(herr, 'JavaFrame');
    %   warning(w);
    %   jFrame.setFigureIcon(jIcon);
    %   uiwait(herr);
    % end
    return;
  end
end

end
