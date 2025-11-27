function [t1, T1, data, SeqOut] = sequence_Recovery(HW, Seq)
%% Sequence using a inversion or saturation recovery method to measure the T1 time
%
%   [t1, T1, data, SeqOut] = sequence_Recovery(HW, Seq)
%
% The sequence consists of a preparation pulse (180 degrees for inversion or 90
% degrees for saturation). For the snapshot (single shot recovery FLASH), this
% preparation pulse is followed by a series of excitation pulses with low
% excitation angles (typically 3 degrees). Otherwise, the preparation pulse is
% followed by a single 90 degrees excitation pulse with varying flip delays
% (tFlip).
%
% INPUT:
%   HW      HW structure or object.
%   Seq     A structure with the following optional fields. If the fields are
%           omitted or empty, default values are used:
%     recovery
%             String that selects the recovery "mode": 'saturation' or
%             'inversion' (default: 'inversion').
%     T1Estimated
%             Estimated T1 in seconds (default: 140e-3).
%     T1EstimatedMin
%             Minimum estimated T1 (default: Seq.T1Estimated/3).
%     T1EstimatedMax
%             Maximum estimated T1 (default: Seq.T1Estimated*3).
%     fSample Sample frequency in Hz of the acquisition windows after the
%             excitation pulses (default: 100e3).
%     tAQFID  Length of these acquisition windows in seconds (default: 0,
%             corresponds to 1 sample).
%     tRelax  Relaxation time between excitation pulse and next inversion pulse.
%             If greater than 0, a sequence of 180 degrees inversion pulses and
%             90 degrees excitation pulses is used to determine the T1 time. If
%             it is equal to 0, a 180 degrees inversion pulse is followed by a
%             series of flip pulses of Seq.excitationFlipAngle degrees (see
%             below). (Default: 0)
%     preparationPulse
%             Handle to a pulse shape function for the preparation pulse. Enter
%             "@Pulse_" followed by tabulator to get a list of available pulse
%             shape functions (default: @Pulse_Rect_Composite180 for inversion
%             and @Pulse_Rect for saturation).
%     inversionPulse
%             Handle to a pulse shape function for the inversion pulse (only
%             used if preparationPulse is not defined). Enter "@Pulse_" followed
%             by tabulator to get a list of available pulse shape functions
%             (default: @Pulse_Rect_Composite180).
%     saturationPulse
%             Handle to a pulse shape function for the saturation pulse (only
%             used if preparationPulse is not defined). Enter "@Pulse_" followed
%             by tabulator to get a list of available pulse shape functions
%             (default: @Pulse_Rect).
%     preparationFlipAngle
%             Flip angle in degrees of the preparation pulse. (default: 180 for
%             inversion and 90 for saturation).
%     preparationPulsePhase
%             Phase of the preparation pulse(s) in degrees (default: 0).
%     inversionPulsePhase
%             Phase of the inversion pulse(s) in degrees (only used if
%             preparationPulsePhase is not defined). (default: 0).
%     saturationPulsePhase
%             Phase of the saturation pulse(s) in degrees (only used if
%             preparationPulsePhase is not defined). (default: 0).
%     excitationPulse
%             Handle to a pulse shape function for the excitation pulse. Enter
%             "@Pulse_" followed by tabulator to get a list of available pulse
%             shape functions (default: @Pulse_Rect).
%     excitationFlipAngle
%             Flip angle in degrees of the excitation pulse. It is always 90 if
%             Seq.tRelax is greater than 0. The default value for Seq.tRelax ==
%             0 is 5 degrees.
%     excitationPulsePhase
%             Phase of the excitation pulses in degrees (default: 0).
%     tFlip   A vector with the times between the centers of the inversion and
%             excitation pulses in seconds. If it is a scalar, a vector with
%             reasonable flip times is automatically created using the following
%             switches (default: Seq.T1EstimatedMin/3).
%     tFlipLog
%             Boolean value. If true, a logarithmic spacing is used for the
%             times between the inversion and excitation pulses. Otherwise, a
%             linear spacing is used (default: 1).
%     tFlipStart
%             Shortest time between inversion and excitation pulses in seconds
%             (default: Seq.tFlip).
%     tFlipEnd
%             Longest time between inversion and excitation pulses in seconds
%             (default: Seq.T1EstimatedMax*3)
%     tFlipStartFactor
%             Approximate factor that is used for the vector of flip times. In
%             case a logarithmic scaling is used, the flip times increase
%             according to:
%               tFlip(n+1) = tFlipStartFactor * tFlip(n)  (default: exp(1))
%             For linear scaling, the flip times increase according to:
%               tFlip(n) = tFlipStart + n * tFlipStartFactor * tFlipStart
%                                                         (default: 1)
%             The actual factor may deviate such that tFlipEnd is exactly
%             reached.
%     nFids   An integer that defines the number of excitations (instead of the
%             one derived from tFlipStartFactor). Default: Derived from
%             tFlipStartFactor.
%     fitT1   Boolean value. If true, "fit_exp" is used to calculate the T1
%             time (default: 1).
%     SteadyState_PreShot180
%             Number of inversion pulses to reach a steady-state before the
%             actual experiment starts (default: 1 if tRelax > 0; 0 otherwise).
%     AQPhase Phase of the acquisition window in degrees (default: 0).
%     Spoil   Structure defining a spoiler after the preparation pulses with the
%             following fields. If the structure is omitted or empty, no spoiler
%             is used.
%       UseCoordinate
%               Number(s) with coordinate direction of the spoiler where
%               1=x, 2=y, 3=z direction (default: 0 corresponding to no
%               spoiler).
%       tStart  Start of spoiler in seconds after center of inversion pulse
%               (default: 1e-3).
%       Duration
%               Duration of spoiler in seconds (default: 1e-3 but shorter than
%               the time between the centers of the inversion and first
%               excitation pulses minus 1e-3).
%       tRamp   Ramp time of the gradient in seconds (default: 1e-3).
%       Amplitude
%               Amplitude of the spoiler in Tesla/meter (default: 100e-3).
%       PreparationOnly
%               Boolean value. If true, only the preparation pulse is spoiled
%               (default: 0).
%       InversionOnly
%               Boolean value. If true, only the inversion pulse is spoiled
%               (only used if PreparationOnly is undefined). (default: 0).
%       SaturationOnly
%               Boolean value. If true, only the saturation pulse is spoiled
%               (only used if PreparationOnly is undefined). (default: 0).
%     ConsoleOut
%             Boolean value. If true, the result of the T1 measurement is
%             displayed in the "Command Window" (default: 0).
%     plot    Boolean value. If true, the measurement data is plotted on an
%             unwrapped timeline (default: true if Seq.tRelax == 0).
%     plotTR  Boolean value. If true, the measurement data is plotted on a
%             timeline that is wrapped at each tRep (default: true if Seq.tRelax
%             > 0).
%     plotSeq Boolean value. If true, the sequence is plotted on an unwrapped
%             timeline (default: true if Seq.tRelax == 0).
%     plotSeqTR
%             Boolean value. If true, the sequence is plotted on a timeline that
%             is wrapped at each tRep (default: true if Seq.tRelax > 0).
%     fitExp  Structure with options for the T1 fit. For supported fields, see
%             documentation of function "fit_exp".
%
% OUTPUT:
%   t1      Single exponential T1 time (if Seq.FitT1 is true. NaN otherwise.)
%   T1      Structure with results of the exponential fit. See: fit_exp.
%   data    Structure with measured data.
%   SeqOut  Structure with actually used measurement parameters.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% default parameters
if nargin < 2, Seq.Plot = []; end
Seq = set_EmptyField(Seq, 'PreProcessSequence', 1);
Seq = set_EmptyField(Seq, 'StartSequence', 1);
Seq = set_EmptyField(Seq, 'PollPPGfast', 1);
Seq = set_EmptyField(Seq, 'GetRawData', 1);
Seq = set_EmptyField(Seq, 'PostProcessSequence', 1);

Seq = set_EmptyField(Seq, 'recovery', 'Inversion');

if ~isemptyfield(Seq, 'T1Estimated')
  Seq = set_EmptyField(Seq, 'T1EstimatedMin', Seq.T1Estimated/3);
  Seq = set_EmptyField(Seq, 'T1EstimatedMax', Seq.T1Estimated*3);
  Seq = set_EmptyField(Seq, 'tFlip', Seq.T1EstimatedMin/3);
elseif isemptyfield(Seq, 'tFlip') && ~isemptyfield(Seq, 'tFlipStart')
  Seq.tFlip = Seq.tFlipStart;
end
if isemptyfield(Seq, 'T1Estimated')
  if isemptyfield(Seq, 'tFlip')
    Seq.T1Estimated = 140e-3;
  else
    Seq.T1Estimated = Seq.tFlip*3;
  end
end
Seq = set_EmptyField(Seq, 'tFlip', Seq.T1Estimated/3);
Seq = set_EmptyField(Seq, 'tFlipStart', Seq.tFlip);
Seq = set_EmptyField(Seq, 'tRelax', 0);
if Seq.tRelax~=0, Seq.excitationFlipAngle = 90; end
Seq = set_EmptyField(Seq, 'excitationFlipAngle', 5);
Seq = set_EmptyField(Seq, 'excitationPulse', @Pulse_Rect);
Seq = set_EmptyField(Seq, 'excitationPulsePhase', 0);
if strcmpi(Seq.recovery, 'Inversion')
  Seq = set_EmptyField(Seq, 'inversionPulse', @Pulse_Rect_Composite180);
  Seq = set_EmptyField(Seq, 'preparationPulse', Seq.inversionPulse);
  Seq.inversionPulse = Seq.preparationPulse;
  Seq = set_EmptyField(Seq, 'preparationFlipAngle', 180);
  Seq = set_EmptyField(Seq, 'inversionPulsePhase', 0);
  Seq = set_EmptyField(Seq, 'preparationPulsePhase', Seq.inversionPulsePhase);
  Seq.inversionPulsePhase = Seq.preparationPulsePhase;
elseif strcmpi(Seq.recovery, 'Saturation')
  Seq = set_EmptyField(Seq, 'saturationPulse', @Pulse_Rect);
  Seq = set_EmptyField(Seq, 'preparationPulse', Seq.saturationPulse);
  Seq.saturationPulse = Seq.preparationPulse;
  Seq = set_EmptyField(Seq, 'preparationFlipAngle', 90);
  Seq = set_EmptyField(Seq, 'saturationPulsePhase', 0);
  Seq = set_EmptyField(Seq, 'preparationPulsePhase', Seq.saturationPulsePhase);
  Seq.saturationPulsePhase = Seq.preparationPulsePhase;
else
  error('PD:sequence_Recovery:UnknownRecovery', ...
    'Seq.recovery = ''%s'' is not allowed (use ''inversion'' or ''saturation''.', ...
    Seq.recovery);
end
Seq = set_EmptyField(Seq, 'ConsoleOut', 0);
Seq = set_EmptyField(Seq, 'plot', 1:double(Seq.tRelax==0));
Seq = set_EmptyField(Seq, 'plotTR', 1:double(Seq.tRelax>0));
if ~isfield(Seq, 'plotData1D'),  Seq.plotData1D = [];  end
if ~isfield(Seq, 'plotSeq'),  Seq.plotSeq = 1:double(Seq.tRelax==0);  end
if ~isfield(Seq, 'plotSeqTR'),  Seq.plotSeqTR = 1:double(Seq.tRelax>0);  end
Seq = set_EmptyField(Seq, 'fSample', 20e3);
Seq = set_EmptyField(Seq, 'tAQFID', 0);
Seq = set_EmptyField(Seq, 'tFlipLog', 1);
if isemptyfield(Seq, 'tFlipStartFactor')
  if Seq.tFlipLog
    Seq.tFlipStartFactor = exp(1);
  else
    Seq.tFlipStartFactor = 1;
  end
end
if isemptyfield(Seq, 'tFlipEnd')
  if ~isemptyfield(Seq, 'T1EstimatedMax')
    Seq.tFlipEnd = Seq.T1EstimatedMax*3;
  else
    if isemptyfield(Seq, 'nFids')
      Seq.nFids = 3;
    end
    if Seq.tFlipLog
      Seq.tFlipEnd = Seq.tFlipStart * sum(Seq.tFlipStartFactor.^(1:Seq.nFids));
    else
      Seq.tFlipEnd = Seq.tFlipStart * Seq.tFlipStartFactor * Seq.nFids;
    end
  end
end
if isemptyfield(Seq, 'nFids')
  if Seq.tFlipLog
    Seq.nFids = max(2, round((log10(Seq.tFlipEnd)-log10(Seq.tFlipStart))/((log10(Seq.tFlipStart*Seq.tFlipStartFactor)-log10(Seq.tFlipStart))))+1);
  else
    Seq.nFids = max(2, round(Seq.tFlipEnd/(Seq.tFlipStartFactor*Seq.tFlipStart)));
  end
end
Seq = set_EmptyField(Seq, 'fitT1', 1);
Seq = set_EmptyField(Seq, 'SteadyState_PreShots180', double(Seq.tRelax~=0));
Seq = set_EmptyField(Seq, 'AQPhase', 0);

if ~isfield(Seq, 'fitExp'), Seq.fitExp = struct(); end
Seq.fitExp = set_EmptyField(Seq.fitExp, 'hParent', 1);
% Seq.fitExp = set_EmptyField(Seq.fitExp, 'CorrectFrequencyOffset', 1); % Let "fit_exp" chose a sensible default value.
Seq.fitExp = set_EmptyField(Seq.fitExp, 'CorrectFrequencyDrift', 1);
Seq.fitExp = set_EmptyField(Seq.fitExp, 'CorrectPhaseOffset', 1);
Seq.fitExp = set_EmptyField(Seq.fitExp, 'EndOffset', 1); % FIXME: Would it be better to always overwrite this?
Seq.fitExp = set_EmptyField(Seq.fitExp, 'RingFilter', 0);

t1 = NaN;
T1 = NaN;
if Seq.PreProcessSequence
  AQ.fSample = Seq.fSample;
  AQ.Frequency = HW.fLarmor;
  AQ.Phase = Seq.AQPhase;
  iDevice = 1;  % FIXME: Support multiple MMRT devices

  %TX
  Seq.tExcitationPulse = max(50/HW.TX(iDevice).fSample, ...
    HW.tFlip1Def * Seq.excitationFlipAngle * Seq.excitationPulse(HW,'Amp'));

  Seq.tPreparationPulse = HW.tFlip1Def * Seq.preparationFlipAngle * Seq.preparationPulse(HW, 'Amp');

  PreparationPulseData = Seq.preparationPulse(HW, 0, 1/Seq.tPreparationPulse*Seq.preparationPulse(HW,'Time'), ...
    pi*HW.tFlip1Def*Seq.preparationFlipAngle/(HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef), 51, Seq.tPreparationPulse, AQ.Frequency(1), Seq.preparationPulsePhase);
  FlipPulseData = Seq.excitationPulse(HW, 0, 1/(HW.tFlip1Def*Seq.excitationFlipAngle) * Seq.excitationPulse(HW,'Time'), ...
    pi*HW.tFlip1Def*Seq.excitationFlipAngle/(HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef), 51, Seq.tExcitationPulse, AQ.Frequency(1), Seq.excitationPulsePhase);

  if Seq.tFlipStart==0, Seq.tFlipStart = Seq.tPreparationPulse + Seq.tExcitationPulse; end

  if isscalar(Seq.tFlip)
    if Seq.tFlipLog
      Seq.tFlip = logspace(log10(Seq.tFlipStart), log10(Seq.tFlipEnd), Seq.nFids);
    else
      Seq.tFlip = linspace(Seq.tFlipStart, Seq.tFlipEnd, Seq.nFids);
    end
  else
    Seq.nFids = numel(Seq.tFlip);
    Seq.tFlipStart = min(Seq.tFlip);
    Seq.tFlipEnd = max(Seq.tFlip);
  end

  if ~isfield(Seq, 'Spoil'), Seq.Spoil = []; end
  if ~isempty(Seq.Spoil)
    if isemptyfield(Seq.Spoil, 'tStart')
      Seq.Spoil.tStart = max(1e-3, Seq.tPreparationPulse/2 + HW.Grad(iDevice).tEC);
    end
    if isemptyfield(Seq.Spoil, 'Duration')
      Seq.Spoil.Duration =  min(1e-3, ...
        Seq.tFlipStart - Seq.Spoil.tStart - max(1e-3, HW.Grad(iDevice).tEC));
    end
    if isemptyfield(Seq.Spoil, 'tRamp')
      Seq.Spoil.tRamp = HW.Grad(iDevice).tRamp;
    end
    Seq.Spoil = set_EmptyField(Seq.Spoil, 'Amplitude', min(100e-3, min(HW.Grad(iDevice).MaxAmp)));
    Seq.Spoil = set_EmptyField(Seq.Spoil, 'UseCoordinate', 0);
    if strcmpi(Seq.recovery, 'Inversion')
      Seq.Spoil = set_EmptyField(Seq.Spoil, 'InversionOnly', 0);
      Seq.Spoil = set_EmptyField(Seq.Spoil, 'PreparationOnly', Seq.Spoil.InversionOnly);
      Seq.Spoil.InversionOnly = Seq.Spoil.PreparationOnly;
    elseif strcmpi(Seq.recovery, 'Saturation')
      Seq.Spoil = set_EmptyField(Seq.Spoil, 'SaturationOnly', 0);
      Seq.Spoil = set_EmptyField(Seq.Spoil, 'PreparationOnly', 0);
      Seq.Spoil.SaturationOnly = Seq.Spoil.PreparationOnly;
    end

    if ~isempty(Seq.plotSeqTR) && ~(islogical(Seq.plotSeqTR) && isequal(Seq.plotSeqTR, false))
      Seq.plotSeqTR = Seq.Spoil.UseCoordinate;
    end
    if ~isempty(Seq.plotSeq) && ~(islogical(Seq.plotSeq) && isequal(Seq.plotSeq, false))
      Seq.plotSeq = Seq.Spoil.UseCoordinate;
    end
  end


  if Seq.tRelax ~= 0
    Seq.tRep = [repmat(Seq.tFlip(1),1,Seq.SteadyState_PreShots180), Seq.tFlip]+Seq.tRelax;

    n = Seq.nFids+Seq.SteadyState_PreShots180;

    TX.Start = NaN(size(PreparationPulseData.Start,1)+size(FlipPulseData.Start,1), Seq.nFids+Seq.SteadyState_PreShots180);
    TX.Duration = TX.Start;
    TX.Amplitude = TX.Start;
    TX.Frequency = TX.Start;
    TX.Phase = TX.Start;

    TX.Duration(1:size(PreparationPulseData.Start,1),:)=  PreparationPulseData.Duration*  ones(1,n);
    TX.Start(1:size(PreparationPulseData.Start,1),:)=     PreparationPulseData.Start*     ones(1,n);
    TX.Amplitude(1:size(PreparationPulseData.Start,1),:)= PreparationPulseData.Amplitude* ones(1,n);
    TX.Frequency(1:size(PreparationPulseData.Start,1),:)= PreparationPulseData.Frequency* ones(1,n);
    TX.Phase(1:size(PreparationPulseData.Start,1),:)=     PreparationPulseData.Phase*     ones(1,n);

    if Seq.nFids>=1
      TX.Duration(size(PreparationPulseData.Start,1)+(1:size(FlipPulseData.Start,1)),1:end) = FlipPulseData.Duration*   ones(1,n);
      TX.Start(size(PreparationPulseData.Start,1)+(1:size(FlipPulseData.Start,1)),1:end) = FlipPulseData.Start*         ones(1,n) + ...
        repmat([repmat(Seq.tFlip(1),1,Seq.SteadyState_PreShots180),Seq.tFlip],size(FlipPulseData.Start,1),1);
      TX.Amplitude(size(PreparationPulseData.Start,1)+(1:size(FlipPulseData.Start,1)),1:end) = FlipPulseData.Amplitude* ones(1,n);
      TX.Frequency(size(PreparationPulseData.Start,1)+(1:size(FlipPulseData.Start,1)),1:end) = FlipPulseData.Frequency* ones(1,n);
      TX.Phase(size(PreparationPulseData.Start,1)+(1:size(FlipPulseData.Start,1)),1:end) = FlipPulseData.Phase*         ones(1,n);
    end

    % AQ
    AQ.Start    =   [nan(1,Seq.SteadyState_PreShots180),TX.Start(end,Seq.SteadyState_PreShots180+1:end)+TX.Duration(end,Seq.SteadyState_PreShots180+1:end)+get_DeadTimeTX2RX(HW,AQ.fSample(1))];
    AQ.nSamples =   [nan(1,Seq.SteadyState_PreShots180),(floor(Seq.tAQFID*AQ.fSample(1)))*ones(1,Seq.nFids)];
    AQ.ResetPhases= [ones(1,Seq.SteadyState_PreShots180),zeros(1,Seq.nFids)];
  else
    Seq.SteadyState_PreShots180 = 0;
    n = Seq.nFids+1;

    Seq.tRep = [Seq.tFlip(1),diff(Seq.tFlip),Seq.tFlip(1)];

    TX.Start = NaN(max(size(PreparationPulseData.Start,1), size(FlipPulseData.Start,1)), 1+Seq.nFids);
    TX.Duration = TX.Start;
    TX.Amplitude = TX.Start;
    TX.Frequency = TX.Start;
    TX.Phase = TX.Start;

    TX.Duration(1:size(PreparationPulseData.Start,1),1) = PreparationPulseData.Duration;
    TX.Start(1:size(PreparationPulseData.Start,1),1) = PreparationPulseData.Start;
    TX.Amplitude(1:size(PreparationPulseData.Start,1),1) = PreparationPulseData.Amplitude;
    TX.Frequency(1:size(PreparationPulseData.Start,1),1) = PreparationPulseData.Frequency;
    TX.Phase(1:size(PreparationPulseData.Start,1),1) = PreparationPulseData.Phase;

    if Seq.nFids >= 1
      TX.Duration(1:size(FlipPulseData.Start,1),2:end) = FlipPulseData.Duration * ones(1, Seq.nFids);
      TX.Start(1:size(FlipPulseData.Start,1),2:end) = FlipPulseData.Start * ones(1, Seq.nFids);
      TX.Amplitude(1:size(FlipPulseData.Start,1),2:end) = FlipPulseData.Amplitude * ones(1, Seq.nFids);
      TX.Frequency(1:size(FlipPulseData.Start,1),2:end) = FlipPulseData.Frequency * ones(1, Seq.nFids);
      TX.Phase(1:size(FlipPulseData.Start,1),2:end) = FlipPulseData.Phase * ones(1, Seq.nFids);
    end

    % AQ
    AQ.Start    =   ([nan,(TX.Start(length(FlipPulseData.Start),2)+TX.Duration(length(FlipPulseData.Start),2)+get_DeadTimeTX2RX(HW,AQ.fSample(1)))*ones(1,Seq.nFids)]);
    AQ.nSamples =   [nan,(floor(Seq.tAQFID*AQ.fSample(1)))*ones(1,Seq.nFids)];
    AQ.ResetPhases= [1,zeros(1,Seq.nFids)];

  end

  AQ.nSamples(AQ.nSamples<1) = 1;

  % Gradients
  for t = HW.Grad(iDevice).n:-1:1
    if ~isempty(Seq.Spoil) && any(t==Seq.Spoil.UseCoordinate)
      Grad(t).Time = cumsum([Seq.Spoil.tStart; Seq.Spoil.tRamp; Seq.Spoil.Duration-2*Seq.Spoil.tRamp; Seq.Spoil.tRamp]) * ones(1, n);
      Grad(t).Amp = [0; Seq.Spoil.Amplitude; Seq.Spoil.Amplitude; 0] * ones(1, n);
      if Seq.tRelax == 0  %
        if any(Grad(t).Time(end,:) + TX.Start(1,:) + 100e-6 > Seq.tRep) || ...
            Seq.Spoil.PreparationOnly
          Grad(t).Time(:,2:end) = NaN;
          Grad(t).Amp(:,2:end) = NaN;
          Grad(t).Repeat = [0, 0, ones(1, Seq.nFids-1)];
        else
          Grad(t).Repeat = [0, ones(1, Seq.nFids)];
        end
      else
        Grad(t).Repeat = [0, ones(1, n-1)];
      end
    else
      Grad(t).Time = NaN(1, n);
      Grad(t).Amp = zeros(1, n);
      Grad(t).Repeat = [0, ones(1, n)];
    end
  end

  % run measurement
  [~, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);

  if ~(Seq.StartSequence || Seq.PollPPGfast || Seq.GetRawData || Seq.PostProcessSequence)
    return;
  end
end


if ~Seq.PreProcessSequence && (Seq.StartSequence || Seq.PollPPGfast || Seq.GetRawData || Seq.PostProcessSequence)
  % If set_sequence wasn't run with these settings in the above block, run it
  % now.
  [~, SeqOut, data, data_1D] = set_sequence(HW, Seq, Seq.AQ, Seq.TX, Seq.Grad);
end

% Plot
if Seq.plot
  if nargin(@plot_data_1D) < 5
    plot_data_1D(HW, data_1D);
  else
    plot_data_1D(HW, data_1D, [], true, Seq.plotData1D);
  end
end

if Seq.plotTR
  plot_data_1D_TR(HW, data_1D);
end

if Seq.fitT1
  if Seq.nFids >= 3
    try
      if isscalar(Seq.fitExp.CorrectFrequencyDrift) && Seq.fitExp.CorrectFrequencyDrift
        Seq.fitExp.CorrectFrequencyDrift = data(1).time_all(:,1,end-Seq.nFids+1:end);
      end
      if Seq.tRelax ~= 0
        T1 = fit_exp(data(1).data(:,1,end-Seq.nFids+1:end), data(1).time_of_tRep(:,1,end-Seq.nFids+1:end), Seq.fitExp);
      else
        T1 = fit_exp(data(1).data(:,1,end-Seq.nFids+1:end), data(1).time_all(:,1,end-Seq.nFids+1:end), Seq.fitExp);
      end
      Seq.fitExp = T1.fitExpSettings;
      t1 = T1.tau;

      % Console output
      if Seq.ConsoleOut
        fprintf('\nT1   = %10.1f ms\n', T1.tau*1000);
        if Seq.fitExp.DoubleExp
          fprintf('T1_1 = %10.1f ms; weighting factor = %5.1f%%\n', T1.xminDouble(3)*1000, T1.tau1w*100);
          fprintf('T1_2 = %10.1f ms; weighting factor = %5.1f%%\n', T1.xminDouble(5)*1000, T1.tau2w*100);
        end
      end

      data(1).DataAmplitude = T1.dataPhaseCorrected;
      data(1).DataTime = T1.timeCorrected;
      SeqOut.iLaplace1D.Problem = Seq.recovery;
    catch ME
      T1 = NaN;
      t1 = NaN;
      data(1).DataAmplitude = NaN;
      data(1).DataTime = NaN;
      warning(getReport(ME));
    end
  else
    warning('nFids lower than number of free parameters of fit. Must be > 2.');
  end
end

SeqOut.data = data;

end
