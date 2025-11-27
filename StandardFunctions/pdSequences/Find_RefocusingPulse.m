function [HW, mySave, Seq] = Find_RefocusingPulse(HW, mySave, Seq)
%% CPMG echo train for 180 degrees pulse amplitude and phase
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

if isemptyfield(Seq, 'T1')
  Seq.T1=HW.FindFrequencyPause;
end

if isemptyfield(Seq, 'plot')
  Seq.plot=1;
end

if isemptyfield(Seq, 'Find_Frequency_interval')
  Seq.Find_Frequency_interval=HW.FindFrequencySweep.maxTime;
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
  Seq.tRelax = max(Seq.T1*2,1);
end

% acceptance criteria
if isemptyfield(Seq, 'relDevMax')
  % tolerance for relative deviation
  Seq.relDevMax = 0.1;
end

% FIXME: Support multiple devices?
iDevice = 1;

Seq.nEcho = 0;
Seq.fSample = 1e6;

Seq.nTau1SteadyState = 0;

Seq.nTau1 = 0;
Seq.FitT2 = 0;
Seq.RawData = 0;

Seq.SeqAverage.average = 1;

% slice settings
Seq.thicknessSlice = Inf;  % thickness of slice selected by continuous gradient in m

Seq.excitationPulse = Seq.solidEchoPulse;  % pulse shape function for excitation pulse
Seq.nEcho = [];
Seq.SeqAverage.average = 1;
Seq.tAQFID = -1;


if ~isnan( Seq.Find_Frequency_interval)
  [HW, mySave] = Find_Frequency_Sweep(HW, mySave, Seq.Find_Frequency_interval, [], 1);  % find magnet frequency
end

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

Seq.plotSeq = 1:3;
if Seq.Find_Frequency_interval~=0
  Seq.Find_Frequency_interval=inf;
end
% actual measurement (inversion pulse)
Seq.Reinitialize = 1;
evalPhasePrevious = [];
for iPhase = 1:4
  [data, SeqOutInv] = sequence_RecoveryCPMG(HW, Seq);
  Seq.LoopCountStart = SeqOutInv.LoopCountEnd;
  Seq.Reinitialize = ~Seq.Find_Frequency_interval;
  Seq.StartSequenceTime = [];
  Seq.plotSeq = [];

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

  ampExcitation = deg2rad(SeqOutInv.excitationAngle) / SeqOutInv.tExcitation / HW.GammaDef;
  paUoutExcitation = ampExcitation / HW.TX(iDevice).PaUout2Amplitude(HW.TX(iDevice).ChannelDef);
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
      paUoutExcitation, paUoutRefocus, SeqOutInv.refocusingAngle, 176, ...
      SeqOutInv.refocusingPhase-90, 176))
    if iPhase == 1
      % FIXME: hgsave and hgload don't restore legends correctly
      % legend(haxes, {'odd echoes', 'even echoes'});
    end
  end

  evalPhase.phaseErrorPerEcho = rad2deg(evalPhase.fitOdd(1)-evalPhase.fitEven(1));  % additional error in phase of echo per echo

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
%   if Seq.getQfactor
%     Seq.tRefocus = round(SeqOutInv.tRefocus * Seq.refocusingAngle / refocusingAngleCurrent * HW.TX.fSample) / HW.TX.fSample;
%     Seq.refocusingAngle = HW.GammaDef/2/pi*360 * Seq.tRefocus * exP.Amplitude(1);
%   end
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
end

end
