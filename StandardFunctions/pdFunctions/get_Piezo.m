function Seq = get_Piezo(Seq, HW)
%% Get (periodic) signal on gradient channel with piezo actuator
%
%     Seq = get_Piezo(Seq, HW)
%
% Seq.Piezo.Frequency                       % Frequency in Hz of Piezo motion
% Seq.Piezo.Start                           % Start of the first Period after PreEmphasis
% Seq.Piezo.Periods                         % Duration of Piezo motion in Periods
% Seq.Piezo.Duration                        % Duration of Piezo motion
% Seq.Piezo.Phase                           % Phase of motion
% Seq.Piezo.Amplitude                       % Amplitude of Piezo current
% Seq.Piezo.Channel                         % Index of Gradient Output
% Seq.Piezo.OffsetCurrent                   % offset current to generate Offset voltage at the Piezo (100 Ohm in parallel)
% Seq.Piezo.OffsetOverdrivePercent
% Seq.Piezo.OffsetRampTime
% Seq.Piezo.OffsetSetTime
% Seq.Piezo.PreEmphasisPeriods
% Seq.Piezo.PreEmphasisOverdrivePercent
% Seq.Piezo.PreEmphasisRelativeAmplitude
% Seq.Piezo.PreEmphasisTime
% Seq.Piezo.PostEmphasisPeriods
% Seq.Piezo.PreEmphasisOverdrivePercent
% Seq.Piezo.PreEmphasisRelativeAmplitude
% Seq.Piezo.PreEmphasisTime
% Seq.Piezo.PostEmphasisPeriods
% Seq.Piezo.PostEmphasisOverdrivePercent
% Seq.Piezo.PostEmphasisRelativeAmplitude
% Seq.Piezo.PostEmphasisTime
% Seq.Piezo.Samples                         % The number of Gradient points maximal used
%
% ------------------------------------------------------------------------------
% (C) Copyright 2017-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%%

if isemptyfield(Seq.Piezo,'Frequency')
  Seq.Piezo.Frequency = HW.Piezo.Frequency;
end
if isemptyfield(Seq.Piezo,'Start')
  Seq.Piezo.Start = HW.Piezo.Start;
end
if isemptyfield(Seq.Piezo,'DisplacementPhase')
  Seq.Piezo.DisplacementPhase = HW.Piezo.DisplacementPhase;
end
if isemptyfield(Seq.Piezo,'Periods')
  Seq.Piezo.Periods = [];
end
if isemptyfield(Seq.Piezo, 'Duration')
  Seq.Piezo.Duration = [];
end
if isempty(Seq.Piezo.Duration)
  if isempty(Seq.Piezo.Periods)
    if isempty(HW.Piezo.Duration)
      Seq.Piezo.Duration = HW.Piezo.Periods/Seq.Piezo.Frequency;
    else
      Seq.Piezo.Duration = HW.Piezo.Duration;
    end
  else
    Seq.Piezo.Duration = Seq.Piezo.Periods/Seq.Piezo.Frequency;
  end
end
if isemptyfield(Seq.Piezo, 'Channel')
  Seq.Piezo.Channel = HW.Piezo.Channel;
end
if isemptyfield(Seq.Piezo, 'DisplacementPerVolt')
  Seq.Piezo.DisplacementPerVolt = HW.Piezo.DisplacementPerVolt;
end
if isemptyfield(Seq.Piezo, 'DisplacementPerAmpere')
  Seq.Piezo.DisplacementPerAmpere = HW.Piezo.DisplacementPerAmpere;
end
if isemptyfield(Seq.Piezo, 'VoltageAmplitudeMax')
  Seq.Piezo.VoltageAmplitudeMax = HW.Piezo.VoltageAmplitudeMax;
end
if isemptyfield(Seq.Piezo, 'CurrentAmplitudeMax')
  Seq.Piezo.CurrentAmplitudeMax = HW.Piezo.CurrentAmplitudeMax;
end
if isemptyfield(Seq.Piezo, 'Capacity')
  % Capacity of piezo
  Seq.Piezo.Capacity = HW.Piezo.Capacity;
end
if isemptyfield(Seq.Piezo, 'Inductance')
  % Inductance of piezo
  Seq.Piezo.Inductance = HW.Piezo.Inductance;
end
if isemptyfield(Seq.Piezo, 'Resistance')
  % Parallel resistance to piezo to get 20 V offset voltage, e.g., 100 Ohm
  Seq.Piezo.Resistance = HW.Piezo.Resistance;
end
if isemptyfield(Seq.Piezo, 'OffsetVoltage')
  % offset voltage, e.g., 20 V
  Seq.Piezo.OffsetVoltage = HW.Piezo.OffsetVoltage;
end
if isemptyfield(Seq.Piezo,'OffsetCurrent')
  % offeset current to get 20 V offset voltage e.g. 0.2 A
  Seq.Piezo.OffsetCurrent = HW.Piezo.OffsetCurrent;
end
if isemptyfield(Seq.Piezo, 'OffsetOverdrivePercent')
  Seq.Piezo.OffsetOverdrivePercent = HW.Piezo.OffsetOverdrivePercent;
end
if isemptyfield(Seq.Piezo, 'OffsetRampTime')
  Seq.Piezo.OffsetRampTime = HW.Piezo.OffsetRampTime;
end
if isemptyfield(Seq.Piezo, 'OffsetSetTime')
  Seq.Piezo.OffsetSetTime = HW.Piezo.OffsetSetTime;
end
if isemptyfield(Seq.Piezo, 'OffsetBeforeAfterSeq')
  Seq.Piezo.OffsetBeforeAfterSeq = HW.Piezo.OffsetBeforeAfterSeq;
end
if isemptyfield(Seq.Piezo, 'PreEmphasisPeriods')
  Seq.Piezo.PreEmphasisPeriods = HW.Piezo.PreEmphasisPeriods;
end
if isemptyfield(Seq.Piezo, 'PreEmphasisOverdrivePercent')
  Seq.Piezo.PreEmphasisOverdrivePercent = HW.Piezo.PreEmphasisOverdrivePercent;
end
if isemptyfield(Seq.Piezo, 'PreEmphasisRelativeAmplitude')
  Seq.Piezo.PreEmphasisRelativeAmplitude = HW.Piezo.PreEmphasisRelativeAmplitude;
end
if isemptyfield(Seq.Piezo, 'PreEmphasisTime')
  Seq.Piezo.PreEmphasisTime = HW.Piezo.PreEmphasisRelativeTime ...
    .* Seq.Piezo.PreEmphasisPeriods ./ Seq.Piezo.Frequency;
end
if isemptyfield(Seq.Piezo, 'PostEmphasisPeriods')
  Seq.Piezo.PostEmphasisPeriods = HW.Piezo.PostEmphasisPeriods;
end
if isemptyfield(Seq.Piezo, 'PostEmphasisOverdrivePercent')
  Seq.Piezo.PostEmphasisOverdrivePercent = HW.Piezo.PostEmphasisOverdrivePercent;
end
if isemptyfield(Seq.Piezo, 'PostEmphasisRelativeAmplitude')
  Seq.Piezo.PostEmphasisRelativeAmplitude = HW.Piezo.PostEmphasisRelativeAmplitude;
end
if isemptyfield(Seq.Piezo, 'PostEmphasisTime')
  Seq.Piezo.PostEmphasisTime = HW.Piezo.PostEmphasisRelativeTime ...
    .* Seq.Piezo.PostEmphasisPeriods ./ Seq.Piezo.Frequency;
end
if isemptyfield(Seq.Piezo, 'UseAtRepetitionTime')
  Seq.Piezo.UseAtRepetitionTime = HW.Piezo.UseAtRepetitionTime;
end
if isemptyfield(Seq.Piezo, 'CLTopBottom')
  Seq.Piezo.CLTopBottom = (Seq.MEG.tEchoExtra > 100e-6);
end
if isemptyfield(Seq.Piezo, 'Samples')
  if isempty(HW.Piezo.Samples)
    Seq.Piezo.Samples = 4*(round(Seq.Piezo.Duration*Seq.Piezo.Frequency)+round(Seq.Piezo.PreEmphasisPeriods)+round(Seq.Piezo.PostEmphasisPeriods))+2*5;
  else
    Seq.Piezo.Samples = HW.Piezo.Samples;
  end
end

if Seq.Piezo.DisplacementPerVolt
  if isemptyfield(Seq.Piezo, 'DisplacementAmplitudeMax')
    Seq.Piezo.DisplacementAmplitudeMax = ...
      Seq.Piezo.VoltageAmplitudeMax .* Seq.Piezo.DisplacementPerVolt;
  end
  Seq.Piezo.Reactance = -1/(2*pi*Seq.Piezo.Frequency*Seq.Piezo.Capacity);
  Seq.Piezo.Impedance = 1/(1/Seq.Piezo.Resistance+1/(1i*Seq.Piezo.Reactance));

  if ~isemptyfield(Seq.Piezo, 'DisplacementAmplitude')
    Seq.Piezo.VoltageAmplitude=abs(Seq.Piezo.DisplacementAmplitude./Seq.Piezo.DisplacementPerVolt);
  end
  if ~isemptyfield(Seq.Piezo, 'VoltageAmplitude')
    Seq.Piezo.DisplacementAmplitude = ...
      abs(Seq.Piezo.DisplacementPerVolt .* Seq.Piezo.VoltageAmplitude);
    Seq.Piezo.CurrentAmplitude = ...
      abs(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance);
  end
  if ~isemptyfield(Seq.Piezo, 'CurrentAmplitude')
    Seq.Piezo.VoltageAmplitude = ...
      abs(Seq.Piezo.CurrentAmplitude .* Seq.Piezo.Impedance);
    Seq.Piezo.DisplacementAmplitude = ...
      abs(Seq.Piezo.DisplacementPerVolt .* Seq.Piezo.VoltageAmplitude);
    Seq.Piezo.CurrentPhase = ...
      angle(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance) / pi * 180;
  end

  MaxRelativeOverdrive = max(...
    max(abs(Seq.Piezo.PreEmphasisRelativeAmplitude)), ...
    max(abs(Seq.Piezo.PostEmphasisRelativeAmplitude)));

  if Seq.Piezo.DisplacementAmplitude*MaxRelativeOverdrive > Seq.Piezo.DisplacementAmplitudeMax
    warning('Piezo:DisplacementAmplitude', ...
      'Seq.Piezo.DisplacementAmplitude too high (%.1f%%): setting to maximum (%.3f)', ...
      Seq.Piezo.DisplacementAmplitude*MaxRelativeOverdrive/Seq.Piezo.DisplacementAmplitudeMax*100, ...
      Seq.Piezo.DisplacementAmplitudeMax/MaxRelativeOverdrive);
    % error(['Seq.Piezo.DisplacementAmplitude too high ( ' num2str( Seq.Piezo.DisplacementAmplitude/Seq.Piezo.DisplacementAmplitudeMax*100,'%10.1f') '% )'])
    Seq.Piezo.DisplacementAmplitude = ...
      Seq.Piezo.DisplacementAmplitudeMax / MaxRelativeOverdrive;
    Seq.Piezo.VoltageAmplitude = ...
      abs(Seq.Piezo.DisplacementAmplitude ./ Seq.Piezo.DisplacementPerVolt);
    Seq.Piezo.CurrentAmplitude = ...
      abs(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance);
    Seq.Piezo.CurrentPhase = ...
      angle(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance) / pi * 180;
  end
  if Seq.Piezo.VoltageAmplitude*MaxRelativeOverdrive > Seq.Piezo.VoltageAmplitudeMax
    warning('Piezo:VoltageAmplitude', ...
      'Seq.Piezo.VoltageAmplitude too high (%.1f%%): setting to maximum (%.3f V)', ...
      Seq.Piezo.VoltageAmplitude*MaxRelativeOverdrive/Seq.Piezo.VoltageAmplitudeMax*100, ...
      Seq.Piezo.VoltageAmplitudeMax/MaxRelativeOverdrive);
    % error(['Seq.Piezo.VoltageAmplitude too high ( ' num2str( Seq.Piezo.VoltageAmplitude/Seq.Piezo.VoltageAmplitudeMax*100,'%10.1f') '% )'])
    Seq.Piezo.VoltageAmplitude = ...
      Seq.Piezo.VoltageAmplitudeMax / MaxRelativeOverdrive;
    Seq.Piezo.DisplacementAmplitude = ...
      abs(Seq.Piezo.DisplacementPerVolt .* Seq.Piezo.VoltageAmplitude);
    Seq.Piezo.CurrentAmplitude = ...
      abs(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance);
    Seq.Piezo.CurrentPhase = ...
      angle(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance) / pi * 180;
  end
  if Seq.Piezo.CurrentAmplitude*MaxRelativeOverdrive > Seq.Piezo.CurrentAmplitudeMax
    warning('Piezo:CurrentAmplitude', ...
      'Seq.Piezo.CurrentAmplitude too high: (%.1f%%): setting to maximum (%.3f A)', ...
      Seq.Piezo.CurrentAmplitude*MaxRelativeOverdrive/Seq.Piezo.CurrentAmplitudeMax*100, ...
      Seq.Piezo.CurrentAmplitudeMax/MaxRelativeOverdrive)
    % error(['Seq.Piezo.CurrentAmplitude too high( ' num2str( Seq.Piezo.CurrentAmplitude/Seq.Piezo.CurrentAmplitudeMax*100,'%10.1f') '% )'])
    Seq.Piezo.CurrentAmplitude = ...
      Seq.Piezo.CurrentAmplitudeMax / MaxRelativeOverdrive;
    Seq.Piezo.VoltageAmplitude = ...
      abs(Seq.Piezo.CurrentAmplitude .* Seq.Piezo.Impedance);
    Seq.Piezo.DisplacementAmplitude = ...
      abs(Seq.Piezo.DisplacementPerVolt .* Seq.Piezo.VoltageAmplitude);
    Seq.Piezo.CurrentPhase = ...
      angle(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance) / pi * 180;
  end

  if HW.Grad.PaCurrentControlled(HW.Piezo.Channel)
    % corrected phase
    Seq.Piezo.Phase = Seq.Piezo.DisplacementPhase + Seq.Piezo.CurrentPhase;
    % amplitude of piezo current
    Seq.Piezo.Amplitude = Seq.Piezo.CurrentAmplitude;
    Seq.Piezo.OffsetAmplitude = Seq.Piezo.OffsetCurrent;
    Seq.Piezo.OffsetVoltage = Seq.Piezo.OffsetCurrent * HW.Piezo.Resistance;
  else
    % not corrected phase
    Seq.Piezo.Phase = Seq.Piezo.DisplacementPhase;
    % amplitude of piezo
    Seq.Piezo.Amplitude = Seq.Piezo.VoltageAmplitude * HW.Piezo.DisplacementPerVolt;
    Seq.Piezo.OffsetAmplitude = Seq.Piezo.OffsetVoltage * HW.Piezo.DisplacementPerVolt;
    Seq.Piezo.OffsetCurrent = Seq.Piezo.OffsetVoltage / HW.Piezo.Resistance;
  end


else % Coil
  if isemptyfield(Seq.Piezo,'DisplacementAmplitudeMax')
    Seq.Piezo.DisplacementAmplitudeMax = ...
      Seq.Piezo.CurrentAmplitudeMax .* Seq.Piezo.DisplacementPerAmpere;
  end
  Seq.Piezo.Reactance=2*pi*Seq.Piezo.Frequency*Seq.Piezo.Inductance;
  Seq.Piezo.Impedance=Seq.Piezo.Resistance+1i*Seq.Piezo.Reactance;

  if ~isemptyfield(Seq.Piezo, 'DisplacementAmplitude')
    Seq.Piezo.CurrentAmplitude = ...
      abs(Seq.Piezo.DisplacementAmplitude ./ Seq.Piezo.DisplacementPerAmpere);
  end
  if ~isemptyfield(Seq.Piezo, 'CurrentAmplitude')
    Seq.Piezo.DisplacementAmplitude = ...
      abs(Seq.Piezo.CurrentAmplitude .* Seq.Piezo.DisplacementPerAmpere);
    Seq.Piezo.VoltageAmplitude = ...
      abs(Seq.Piezo.CurrentAmplitude .* Seq.Piezo.Impedance);
  end
  if ~isemptyfield(Seq.Piezo, 'VoltageAmplitude')
    Seq.Piezo.CurrentAmplitude = ...
      abs(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance);
    Seq.Piezo.DisplacementAmplitude = ...
      abs(Seq.Piezo.CurrentAmplitude .* Seq.Piezo.DisplacementPerAmpere);
    Seq.Piezo.VoltagePhase = ...
      angle(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance) / pi * 180;
  end

  MaxRelativeOverdrive = max(...
    max(abs(Seq.Piezo.PreEmphasisRelativeAmplitude)), ...
    max(abs(Seq.Piezo.PostEmphasisRelativeAmplitude)));

  if Seq.Piezo.DisplacementAmplitude*MaxRelativeOverdrive > Seq.Piezo.DisplacementAmplitudeMax
    warning('Piezo:DisplacementAmplitude', ...
      'Seq.Piezo.DisplacementAmplitude too high (%.1f%%): Setting to maximum (%.3f)', ...
      Seq.Piezo.DisplacementAmplitude*MaxRelativeOverdrive/Seq.Piezo.DisplacementAmplitudeMax*100, ...
      Seq.Piezo.DisplacementAmplitudeMax/MaxRelativeOverdrive)
    % error(['Seq.Piezo.DisplacementAmplitude too high ( ' num2str( Seq.Piezo.DisplacementAmplitude/Seq.Piezo.DisplacementAmplitudeMax*100,'%10.1f') '% )'])
    Seq.Piezo.DisplacementAmplitude = ...
      Seq.Piezo.DisplacementAmplitudeMax / MaxRelativeOverdrive;
    Seq.Piezo.CurrentAmplitude = ...
      abs(Seq.Piezo.DisplacementAmplitude ./ Seq.Piezo.DisplacementPerAmpere);
    Seq.Piezo.VoltageAmplitude = ...
      abs(Seq.Piezo.CurrentAmplitude .* Seq.Piezo.Impedance);
    Seq.Piezo.VoltagePhase = ...
      angle(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance) / pi * 180;
  end
  if Seq.Piezo.VoltageAmplitude*MaxRelativeOverdrive > Seq.Piezo.VoltageAmplitudeMax
    warning('Piezo:VoltageAmplitude', ...
      'Seq.Piezo.VoltageAmplitude too high (%.1f%%): Setting to maximum (%.3f V)', ...
      Seq.Piezo.VoltageAmplitude*MaxRelativeOverdrive/Seq.Piezo.VoltageAmplitudeMax*100, ...
      Seq.Piezo.VoltageAmplitudeMax/MaxRelativeOverdrive);
    % error(['Seq.Piezo.VoltageAmplitude too high ( ' num2str( Seq.Piezo.VoltageAmplitude/Seq.Piezo.VoltageAmplitudeMax*100,'%10.1f') '% )'])
    Seq.Piezo.VoltageAmplitude = ...
      Seq.Piezo.VoltageAmplitudeMax / MaxRelativeOverdrive;
    Seq.Piezo.CurrentAmplitude = ...
      abs(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance);
    Seq.Piezo.DisplacementAmplitude = ...
      abs(Seq.Piezo.CurrentAmplitude .* Seq.Piezo.DisplacementPerAmpere);
    Seq.Piezo.VoltagePhase = ...
      angle(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance) / pi * 180;
  end
  if Seq.Piezo.CurrentAmplitude*MaxRelativeOverdrive > Seq.Piezo.CurrentAmplitudeMax
    warning('Piezo:CurrentAmplitude', ...
      'Seq.Piezo.CurrentAmplitude too high (%.1f%%): Setting to maximum (%.3f A)', ...
      Seq.Piezo.CurrentAmplitude*MaxRelativeOverdrive/Seq.Piezo.CurrentAmplitudeMax*100, ...
      Seq.Piezo.CurrentAmplitudeMax/MaxRelativeOverdrive);
    % error(['Seq.Piezo.CurrentAmplitude too high( ' num2str( Seq.Piezo.CurrentAmplitude/Seq.Piezo.CurrentAmplitudeMax*100,'%10.1f') '% )'])
    Seq.Piezo.CurrentAmplitude = ...
      Seq.Piezo.CurrentAmplitudeMax / MaxRelativeOverdrive;
    Seq.Piezo.DisplacementAmplitude = ...
      abs(Seq.Piezo.CurrentAmplitude .* Seq.Piezo.DisplacementPerAmpere);
    Seq.Piezo.VoltageAmplitude = ...
      abs(Seq.Piezo.CurrentAmplitude .* Seq.Piezo.Impedance);
    Seq.Piezo.VoltagePhase = ...
      angle(Seq.Piezo.VoltageAmplitude ./ Seq.Piezo.Impedance) / pi * 180;
  end

  if HW.Grad.PaCurrentControlled(HW.Piezo.Channel)
    % not corrected phase
    Seq.Piezo.Phase = Seq.Piezo.DisplacementPhase;
    % amplitude of coil
    Seq.Piezo.Amplitude = Seq.Piezo.CurrentAmplitude * HW.Piezo.DisplacementPerAmpere;
    Seq.Piezo.OffsetAmplitude = Seq.Piezo.OffsetCurrent * HW.Piezo.DisplacementPerAmpere;
    Seq.Piezo.OffsetVoltage = Seq.Piezo.OffsetCurrent * HW.Piezo.Resistance;
  else
    % corrected phase
    Seq.Piezo.Phase = Seq.Piezo.DisplacementPhase + Seq.Piezo.VoltagePhase;
    % amplitude of coil Voltage
    Seq.Piezo.Amplitude = Seq.Piezo.VoltageAmplitude;
    Seq.Piezo.OffsetAmplitude = Seq.Piezo.OffsetCurrent * HW.Piezo.Resistance;
    Seq.Piezo.OffsetCurrent = Seq.Piezo.OffsetCurrent;
    HW.Grad.LoadRin(4) = abs(Seq.Piezo.Impedance);
  end


end

Seq.Piezo.Periods = round(Seq.Piezo.Duration*Seq.Piezo.Frequency);
Seq.Piezo.Duration = Seq.Piezo.Periods/Seq.Piezo.Frequency;
Seq.Piezo.PreEmphasisPeriods = round(Seq.Piezo.PreEmphasisPeriods);
Seq.Piezo.PostEmphasisPeriods = round(Seq.Piezo.PostEmphasisPeriods);
if Seq.SingletRep
  Seq.Piezo.SamplesPerPeriode = ...
    floor((Seq.Piezo.Samples-2*16) / ...
          (Seq.Piezo.Periods+Seq.Piezo.PreEmphasisPeriods+Seq.Piezo.PostEmphasisPeriods)/2+1e-12)*2;
else
  Seq.Piezo.SamplesPerPeriode = ...
    floor((Seq.Piezo.Samples-2*16) / ...
          (min([max([Seq.tRep(1), Seq.tEcho-Seq.Piezo.Start]), Seq.Piezo.Duration]) * Seq.Piezo.Frequency + ...
           max([Seq.Piezo.PreEmphasisPeriods, Seq.Piezo.PostEmphasisPeriods]))/2 + 1e-12)*2;
end
Seq.Piezo.SamplesPerPeriode(Seq.Piezo.SamplesPerPeriode>22) = 22;

Seq.Piezo.Phase = mod(Seq.Piezo.Phase, 360);


if Seq.Piezo.SamplesPerPeriode < 4
  error(['Seq.Piezo.Samples too small, min= ' num2str(4*(round(Seq.Piezo.Duration*Seq.Piezo.Frequency)+round(Seq.Piezo.PreEmphasisPeriods)+round(Seq.Piezo.PostEmphasisPeriods))+2*5)])
elseif Seq.Piezo.SamplesPerPeriode < 11
%   figure; plot(linspace(0,2*pi,7),1.095*sin(linspace(0,2*pi,7)),'-x',linspace(0,2*pi,1000),sin(linspace(0,2*pi,1000))); shg;
  Seq.Piezo.AmplitudeCorrected = Seq.Piezo.Amplitude * sqrt(0.5) / sqrt(mean(interp1(linspace(0,2*pi,7),sin(linspace(0,2*pi,7)),linspace(0,2*pi,1001),'linear').^2));
  Seq.Piezo.AmplitudeVector = ...
    repmat(sin([1/6;2/6;4/6;5/6]*2*pi).*Seq.Piezo.AmplitudeCorrected, ...
           [1,Seq.Piezo.Periods+Seq.Piezo.PreEmphasisPeriods+Seq.Piezo.PostEmphasisPeriods+2]);
  Seq.Piezo.TimeLine = ...
    repmat([1/Seq.Piezo.Frequency/6; 2/Seq.Piezo.Frequency/6; 4/Seq.Piezo.Frequency/6; 5/Seq.Piezo.Frequency/6], ...
           [1,Seq.Piezo.Periods+Seq.Piezo.PreEmphasisPeriods+Seq.Piezo.PostEmphasisPeriods+2]) + ...
    repmat(-(Seq.Piezo.PreEmphasisPeriods+1)/Seq.Piezo.Frequency:1/Seq.Piezo.Frequency:(Seq.Piezo.Periods+Seq.Piezo.PostEmphasisPeriods)/Seq.Piezo.Frequency, ...
           [4,1]);
  Seq.Piezo.TimeLine = Seq.Piezo.TimeLine-Seq.Piezo.Phase/360/Seq.Piezo.Frequency;

  Seq.Piezo.AmplitudeVector = Seq.Piezo.AmplitudeVector(Seq.Piezo.TimeLine>=(6e-6-Seq.Piezo.PreEmphasisPeriods/Seq.Piezo.Frequency));
  Seq.Piezo.TimeLine = Seq.Piezo.TimeLine(Seq.Piezo.TimeLine>=(6e-6-Seq.Piezo.PreEmphasisPeriods/Seq.Piezo.Frequency));

  Seq.Piezo.AmplitudeVector = Seq.Piezo.AmplitudeVector(Seq.Piezo.TimeLine+6e-6<=((Seq.Piezo.Periods+Seq.Piezo.PostEmphasisPeriods)/Seq.Piezo.Frequency));
  Seq.Piezo.TimeLine = Seq.Piezo.TimeLine(Seq.Piezo.TimeLine+6e-6<=((Seq.Piezo.Periods+Seq.Piezo.PostEmphasisPeriods)/Seq.Piezo.Frequency));

  Seq.Piezo.AmplitudeVector = [0;Seq.Piezo.AmplitudeVector;0];
  Seq.Piezo.TimeLine = [Seq.Piezo.TimeLine(1)-1/Seq.Piezo.Frequency/6;Seq.Piezo.TimeLine;Seq.Piezo.TimeLine(end)+1/Seq.Piezo.Frequency/6];

  Seq.Piezo.PreEmphasisTimeRegrid = Seq.Piezo.TimeLine(Seq.Piezo.TimeLine<0);
  Seq.Piezo.PreEmphasisAmplitudeRegrid = ...
    interp1(Seq.Piezo.PreEmphasisTime-Seq.Piezo.PreEmphasisPeriods./Seq.Piezo.Frequency, ...
            Seq.Piezo.PreEmphasisRelativeAmplitude, ...
            Seq.Piezo.PreEmphasisTimeRegrid, 'pchip', 'extrap') .* ...
    Seq.Piezo.AmplitudeVector(Seq.Piezo.TimeLine<0);

  Seq.Piezo.PostEmphasisTimeRegrid = Seq.Piezo.TimeLine(Seq.Piezo.TimeLine>Seq.Piezo.Periods/Seq.Piezo.Frequency);
  Seq.Piezo.PostEmphasisAmplitudeRegrid = ...
    interp1(Seq.Piezo.PostEmphasisTime+Seq.Piezo.Periods/Seq.Piezo.Frequency, ...
    Seq.Piezo.PostEmphasisRelativeAmplitude, ...
    Seq.Piezo.PostEmphasisTimeRegrid,'pchip','extrap') .* ...
    Seq.Piezo.AmplitudeVector(Seq.Piezo.TimeLine>Seq.Piezo.Periods/Seq.Piezo.Frequency);

  Seq.Piezo.AmplitudeVector = Seq.Piezo.AmplitudeVector(and(~(Seq.Piezo.TimeLine<0),~(Seq.Piezo.TimeLine>+Seq.Piezo.Periods/Seq.Piezo.Frequency)));
else
  Seq.Piezo.TimeLine = ...
    linspace(-Seq.Piezo.PreEmphasisPeriods/Seq.Piezo.Frequency, ...
             (Seq.Piezo.Periods+Seq.Piezo.PostEmphasisPeriods)/Seq.Piezo.Frequency, ...
             (Seq.Piezo.Periods+Seq.Piezo.PreEmphasisPeriods+Seq.Piezo.PostEmphasisPeriods)*Seq.Piezo.SamplesPerPeriode+1).';
  Seq.Piezo.AmplitudeCorrected = ...
    Seq.Piezo.Amplitude * ...
    sqrt(0.5) / sqrt(mean(interp1(linspace(0,2*pi,Seq.Piezo.SamplesPerPeriode+1), ...
                                  sin(linspace(0,2*pi,Seq.Piezo.SamplesPerPeriode+1)+Seq.Piezo.Phase./180.*pi), ...
                                  linspace(0,2*pi,1001), 'linear').^2));

  Seq.Piezo.AmplitudeVector = ...
    sin(linspace(0, Seq.Piezo.Periods/Seq.Piezo.Frequency, ...
                 Seq.Piezo.Periods*Seq.Piezo.SamplesPerPeriode+1).' .* ...
        2.*pi.*Seq.Piezo.Frequency + Seq.Piezo.Phase./180.*pi) .* Seq.Piezo.AmplitudeCorrected;

  Seq.Piezo.PreEmphasisTimeRegrid = ...
    linspace(-Seq.Piezo.PreEmphasisPeriods./Seq.Piezo.Frequency, 0, ...
             Seq.Piezo.SamplesPerPeriode.*Seq.Piezo.PreEmphasisPeriods+1).';
  Seq.Piezo.PreEmphasisAmplitudeRegrid = ...
    interp1(Seq.Piezo.PreEmphasisTime-Seq.Piezo.PreEmphasisPeriods./Seq.Piezo.Frequency, ...
            Seq.Piezo.PreEmphasisRelativeAmplitude, ...
            Seq.Piezo.PreEmphasisTimeRegrid, 'pchip', 'extrap') .* ...
    sin(Seq.Piezo.PreEmphasisTimeRegrid*2*pi*Seq.Piezo.Frequency+Seq.Piezo.Phase./180.*pi)*Seq.Piezo.AmplitudeCorrected;
  Seq.Piezo.PreEmphasisAmplitudeRegrid=Seq.Piezo.PreEmphasisAmplitudeRegrid(1:end-1);

  Seq.Piezo.PostEmphasisTimeRegrid = ...
    linspace(Seq.Piezo.Periods/Seq.Piezo.Frequency, ...
             (Seq.Piezo.Periods+Seq.Piezo.PostEmphasisPeriods)./Seq.Piezo.Frequency, ...
             Seq.Piezo.SamplesPerPeriode.*Seq.Piezo.PostEmphasisPeriods+1).';
  Seq.Piezo.PostEmphasisAmplitudeRegrid = ...
    interp1(Seq.Piezo.PostEmphasisTime+Seq.Piezo.Periods/Seq.Piezo.Frequency, ...
            Seq.Piezo.PostEmphasisRelativeAmplitude, ...
            Seq.Piezo.PostEmphasisTimeRegrid, 'pchip', 'extrap') .* ...
    sin(Seq.Piezo.PostEmphasisTimeRegrid*2*pi*Seq.Piezo.Frequency+Seq.Piezo.Phase./180.*pi)*Seq.Piezo.AmplitudeCorrected;
  Seq.Piezo.PostEmphasisAmplitudeRegrid = Seq.Piezo.PostEmphasisAmplitudeRegrid(2:end);
end

Seq.Piezo.Samples = numel(Seq.Piezo.TimeLine) + 16 + 16;

Seq.Piezo.Time = NaN(Seq.Piezo.Samples, numel(Seq.tRep));
Seq.Piezo.Amp = NaN(Seq.Piezo.Samples, numel(Seq.tRep));

OffsetRampTime = linspace(0,Seq.Piezo.OffsetRampTime,16).';
OffsetRampAmp = interp1([Seq.Piezo.OffsetRampTime*(-1)/2;0;Seq.Piezo.OffsetRampTime*1/2;Seq.Piezo.OffsetRampTime*2/2;Seq.Piezo.OffsetRampTime*3/2],...
        [0;0;Seq.Piezo.OffsetAmplitude.*(Seq.Piezo.OffsetOverdrivePercent+100)./100;Seq.Piezo.OffsetAmplitude;Seq.Piezo.OffsetAmplitude],...
        OffsetRampTime,'PCHIP');
% figure(234); plot(OffsetRampTime,OffsetRampAmp);

if Seq.Piezo.OffsetBeforeAfterSeq
  Seq.Piezo.OffsetSetTimeBefore=max(Seq.Piezo.OffsetSetTime+(Seq.Piezo.Start+Seq.Piezo.TimeLine(1)),Seq.Piezo.OffsetSetTime);
  Seq.Piezo.OffsetSetTimeAfter = ...
    max(Seq.Piezo.OffsetSetTime + ...
        (Seq.tEcho*(size(Seq.tRepTurboBlock,1)-1-Seq.SteadyState_PostShots180) ...
         - Seq.Piezo.PostEmphasisPeriods/Seq.Piezo.Frequency ...
         - Seq.Piezo.Start - Seq.Piezo.Duration), ...
        Seq.Piezo.OffsetSetTime);
else
  Seq.Piezo.OffsetSetTimeBefore=Seq.Piezo.OffsetSetTime;
  Seq.Piezo.OffsetSetTimeAfter=Seq.Piezo.OffsetSetTime;
end

% time of gradient points in seconds
Seq.Piezo.Time(:,Seq.Piezo.UseAtRepetitionTime) = ...
  (cumsum([...
          (Seq.Piezo.Start+Seq.Piezo.TimeLine(1))-Seq.Piezo.OffsetSetTimeBefore-Seq.Piezo.OffsetRampTime;...
          diff(OffsetRampTime);...
          Seq.Piezo.OffsetSetTimeBefore;...
          diff(Seq.Piezo.TimeLine(1:end));...
          Seq.Piezo.OffsetSetTimeAfter;...
          diff(-OffsetRampTime(end:-1:1))])...
          ) * ones(1, numel(Seq.Piezo.UseAtRepetitionTime));
% amplitude of gradient points
Seq.Piezo.Amp(:,Seq.Piezo.UseAtRepetitionTime) = ...
  ([...                                 
    OffsetRampAmp;...
    Seq.Piezo.PreEmphasisAmplitudeRegrid+Seq.Piezo.OffsetAmplitude;...
    Seq.Piezo.AmplitudeVector+Seq.Piezo.OffsetAmplitude;...
    Seq.Piezo.PostEmphasisAmplitudeRegrid+Seq.Piezo.OffsetAmplitude;...
    OffsetRampAmp(end:-1:1)
    ]) * ones(1, numel(Seq.Piezo.UseAtRepetitionTime));

%  figure(2);plot(Seq.Piezo.Time(:,Seq.Piezo.UseAtRepetitionTime(1)),Seq.Piezo.Amp(:,Seq.Piezo.UseAtRepetitionTime(1)))
if ~Seq.SingletRep

  if isemptyfield(Seq, 'CLTime')
    Seq.CLTime = 20e-6 ...
      + 4/HW.MMRT.fSystem*2*2*(Seq.Piezo.SamplesPerPeriode+4)*ceil(Seq.tEcho*Seq.Piezo.Frequency);
  end

  t2=Seq.CLTime(1)+20e-6;
  t1=Seq.Slice(2).CenterOfPulse-Seq.Slice(2).GradLength/2-t2/2;
  % tend=Seq.Slice(2).CenterOfPulse-Seq.Slice(2).GradLength/2+Seq.tRep(1)-Seq.CLTime(1)-20e-6;
  if Seq.Piezo.CLTopBottom
    t1 = (floor((t1 ...
                 + mod(Seq.Piezo.Phase ...
                       + HW.Grad.TimeDelay(Seq.Piezo.Channel)*(180/Seq.Piezo.Frequency) ...
                       + 90*mod(Seq.MEG.durInvGradientPulses+1,2), 180) /360/Seq.Piezo.Frequency) ...
                 / (0.5/Seq.Piezo.Frequency))*0.5/Seq.Piezo.Frequency ...
                - mod(Seq.Piezo.Phase ...
                      + HW.Grad.TimeDelay(Seq.Piezo.Channel)*(180/Seq.Piezo.Frequency) ...
                      + 90*mod(Seq.MEG.durInvGradientPulses+1,2), 180) /360/Seq.Piezo.Frequency) ...
         + t2/2;
  end
%   tend=t1+Seq.tRep(1)-t2;

  % pezoEndi=numel([Seq.Piezo.OffsetSetTimeAfter;diff(-OffsetRampTime(end:-1:1))]);
  % pezoEnd=Seq.Piezo.Time(end-pezoEndi,Seq.Piezo.UseAtRepetitionTime(1));
  n = 0;
  % loop over all echo trains
  for t = Seq.Piezo.UseAtRepetitionTime.'
    n = n+1;

    % piezo time and amplitude for complete echo train
    Time = Seq.Piezo.Time(:,t);
    Amp = Seq.Piezo.Amp(:,t);

    % initialize matrices for time and amplitude in each tRep
    nAmpSteps = size(Seq.Piezo.Amp,1);
    Seq.Piezo.Time(:,Seq.tRepTurboBlock(:,n)) = NaN(nAmpSteps, numel(Seq.tRepTurboBlock(:,n)));
    Seq.Piezo.Amp(:,Seq.tRepTurboBlock(:,n)) = NaN(nAmpSteps, numel(Seq.tRepTurboBlock(:,n)));

    % times of tRep borders
    tRepEnd = cumsum(Seq.tRep(Seq.tRepTurboBlock(:,n))).';
    tRepStart = tRepEnd-Seq.tRep(Seq.tRepTurboBlock(:,n)).';

    % t1 is calculated such that tRep border can be at piezo maximum or minimum
    TimeStart = tRepStart + t1;
    TimeStart(1) = min(TimeStart(1), Time(1));
    % TimeEnd=tRepStart+tend;

    % t2 is command load time + (arbitrary) margin
    TimeEnd = tRepEnd + t1 - t2;
    % TimeEnd(tRepEnd>=pezoEnd)=(TimeEnd(tRepEnd>=pezoEnd))-t1;
    % TimeStart(tRepEnd>=pezoEnd)=(TimeStart(tRepEnd>=pezoEnd))-t1;
    % TimeEnd(end)=min(TimeEnd(end),max(Time));

    % interpolate piezo amplitude at start and end of tReps
    AmpStart = interp1([Time(1)-1e6;Time;Time(end)+1e6],[Amp(1);Amp;Amp(end)],TimeStart);
    AmpEnd = interp1([Time(1)-1e6;Time;Time(end)+1e6],[Amp(1);Amp;Amp(end)],TimeEnd);
    % figure(234); plot(Time,Amp,TimeStart,AmpStart,'+b',TimeEnd,AmpEnd,'xg');

    % loop over all tReps in echo train
    for tt = 1:find(Time(end)<TimeEnd, 1, 'first')
      % find index of first and last action in tRep. (They might be empty)
      iStart = find((Time>TimeStart(tt))&(Time<TimeEnd(tt)), 1, 'first');
      iEnd = find((Time<TimeEnd(tt))&(Time>TimeStart(tt)), 1, 'last');

      % add points in current tRep
      if Amp(iEnd)
        % If the amplitude at the end of the current tRep is non-zero, add
        % two points with the interpolated end amplitude (to avoid a gradient
        % ramp during the CLTime).
        Seq.Piezo.Time(1:size(iStart:iEnd,2)+3,t+tt-1) = ...
          [TimeStart(tt); Time(iStart:iEnd); TimeEnd(tt); TimeEnd(tt)+1e-6] - tRepStart(tt);
        Seq.Piezo.Amp(1:size(iStart:iEnd,2)+3,t+tt-1) = ...
          [AmpStart(tt); Amp(iStart:iEnd); AmpEnd(tt); AmpEnd(tt)];
      else
        Seq.Piezo.Time(1:size(iStart:iEnd,2)+1,t+tt-1) = ...
          [TimeStart(tt); Time(iStart:iEnd)] - tRepStart(tt);
        Seq.Piezo.Amp(1:size(iStart:iEnd,2)+1,t+tt-1) = ...
          [AmpStart(tt); Amp(iStart:iEnd)];
      end

      if ~isempty(iEnd) && ~isempty(iStart) && (iEnd >= iStart) && ...
          nAmpSteps > 1 && (t+tt-1>1) && isnan(Seq.Piezo.Amp(2,t+tt-2))
        % The point at TimeStart is only necessary if there are any points to
        % set in the current tRep and there was more than one point in the
        % previous tRep.  Remove the point at TimeStart otherwise, to avoid that
        % this (dummy) point extends the tOffset unnecessarily.
        Seq.Piezo.Time(1:end-1,t+tt-1) = Seq.Piezo.Time(2:end,t+tt-1);
        Seq.Piezo.Amp(1:end-1,t+tt-1) = Seq.Piezo.Amp(2:end,t+tt-1);
      end
    end

  %   ttt=Seq.Piezo.Time(:,t:t+tt-1);
  %   tta=Seq.Piezo.Amp(:,t:t+tt-1);
  %   figure(235); plot(ttt(~isnan(ttt)),tta(~isnan(ttt)),'-dk',TimeStart,AmpStart,'+b',TimeEnd,AmpEnd,'xg');
  %   figure(236); plot(Seq.Piezo.Time(:,t:t+tt-1),Seq.Piezo.Amp(:,t:t+tt-1),'-dk',TimeStart,AmpStart,'+b',TimeEnd,AmpEnd,'xg');

  end
  % [Seq.Piezo.Time,Timei]=unique(Seq.Piezo.Time);
  % Seq.Piezo.Amp(:)=Seq.Piezo.Amp(Timei);
  timax=max(sum(~isnan(Seq.Piezo.Time),1),[],2);
  Seq.Piezo.Time=Seq.Piezo.Time(1:timax,:);
  Seq.Piezo.Amp=Seq.Piezo.Amp(1:timax,:);
end

% if or(any(isnan(Seq.Piezo.Time)),any(isnan(Seq.Piezo.Amp)))
%     error('NAN in Seq.Piezo.Amp or Seq.Piezo.Time')
% end

end
