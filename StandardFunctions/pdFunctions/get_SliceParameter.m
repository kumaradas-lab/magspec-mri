function Seq = get_SliceParameter(Seq, HW, idxModule)
%% Create "module" with (slice selective) rf pulse
%
%   Seq = get_SliceParameter(Seq, HW, idxModule)
%
% This function create slice selection gradients (including the respective
% dephase and rephase gradients) derived from the values in Seq.Slice. It also
% creates the corresponding rf pulse.
%
%
% INPUT:
%
% "Seq" is a structure with the field "Slice" that is an (array of) structures
% with the following fields:
%
%   UseCoordinate
%       Number of the coordinate for which the phase encoder is used.
%       (Default: 1)
%
%   Thickness
%       Thickness of slice that is excited in meter. Setting this to Inf (or a
%       very large value, means no slice selection.
%       (Default: Inf)
%
%   CenterOfPulse
%       Center of excitation pulse in seconds on the tRep timeline. (Default: 0)
%
%   MaxGradAmp
%       Maximum gradient amplitude in T/m. (Default: min(HW.Grad.MaxAmp(1:3)) )
%
%   GradAmp
%       Optionally, "fixed" gradient amplitude in T/m. If empty, the slice
%       gradient amplitude is determined from the Pulse.Amplitude or MaxAmpGrad.
%       (Default: [])
%
%   GradSign
%       Sign of the slice selection pulse. (Default: 1)
%
%   tRamp
%       Ramp time of the gradients in seconds. (default: HW.Grad.tRamp)
%
%   GradRephaseLength / GradDephaseLength
%       Length of the rephase/dephase pulse in seconds.
%       (Default: approx. half the slice amplitude, but not shorter than
%       2*tRamp)
%
%   CenterOfRephase / CenterOfDephase
%       Center of rephase/dephase pulse in seconds on the tRep timeline.
%       (default: Such that the ramps of the slice selection pulse and the
%       rephase or dephase pulse overlap)
%
%   GradRephaseSign / GradDephaseSign
%       Sign of the rephase/dephase pulse.
%       (Default: GradRephaseSign = -1, GradDephaseSign = 1)
%
%   tEC
%       Additional time in seconds that the gradients need to reach their set
%       amplitude due to Eddy currents in the pole shoe. (Default: HW.Grad.tEC)
%
%   GradTimeDelay
%       1x3 vector with time delays in seconds for each gradient channel.
%       Positive values lead to the gradient pulses being set earlier.
%       (Default: [0 0 0])
%
%   useAQSlice
%       Index that is used for Seq.AQSlice for the following default parameters.
%       (Default: 1)
%
%   alfa / phi / theta
%       Angles in radians that activley rotate the direction of the selected
%       slice (see UseCoordinate). The first rotation "alfa" is around the
%       x-axis, the second rotation "phi" is around the y-axis and the last
%       rotation "theta" is around the z-axis.
%       (Default: 0 or the corresponding values in
%       Seq.AQSlice(Seq.Slice(t).useAQSlice) )
%
%   UseAtRepetitionTime
%       Indices for the repetition times (tReps) at which the gradient pulses
%       are used. (Mandatory! No default!)
%
%   UseAtRepetitionTimeRephase / UseAtRepetitionTimeDephase
%       Indices for the repetition times where the rephase/dephase pulses should
%       be used. (Default: UseAtRepetitionTime)
%
%   distance
%       Slice coordinate of the origin of the gradient system in the image
%       coordinate system in meter. I.e., a positive value moves the center of
%       the image towards positive slice direction. Calling functions might set
%       this from the value of Seq.AQSlice(1).Center2OriginImage(1).
%       (Default: 0)
%
%   Pulse
%       Structure for the TX pulse with the following fields:
%
%     Function
%         Function handle to the Pulse shape function. (Default: @Pulse_Rect)
%     MaxNumberOfSegments
%         Maximum number of segments for the pulse shape. (Default: 51)
%     Amplitude
%         Amplitude of the TX pulse in Tesla. (Default: HW.TX.AmpDef)
%     FlipAngle
%         Flip angle of the TX pulse in degrees. (Default: 90)
%     Phase
%         Phase of the TX pulse in degrees. (Default: 0)
%     FlipAngleX
%         Flip angle of the TX pulse at the secondary frequency in degrees. See
%         Slice.GammaX. (Default: Slice.Pulse.FlipAngle)
%     PhaseX
%         Phase of the TX pulse at the secondary frequency in degrees. See
%         Slice.GammaX. (Default: Slice.Pulse.Phase)
%
%   Gamma
%       Gyromagnetic ratio of the sample in rad/s/T (default: HW.GammaDef).
%
%   GammaX
%       Gyromagnetic ratio of the sample in rad/s/T for dual frequency pulses
%       (default: []). If set, a pulse is emitted with a bandwidth that selects
%       the same slice thickness and with the same center as the primary pulse.
%
%
% "HW" is the HW object or structure.
%
%
% "idxModule" is an optional input argument with the indices of the slice
% modules for which the pulse program matrices should be created.
% (Default: 1:numel(Seq.Slice))
%
%
% OUTPUT:
%
% To each structure Seq.Phase the following fields are added:
%
%   Grad / GradRephase / GradDephase
%       Structure with the fields "Amp" and "Time" containing the amplitudes and
%       times of the slice selection and rephase/dephase pulses corresponding to
%       the input parameters that can be added to the sequence with "add_Grad".
%
%   TX
%       Structure containing the settings for the (excitation) pulse. It can be
%       added to the sequence with "add_TX".
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default parameters
if nargin < 3 || isempty(idxModule)
  idxModule = 1:numel(Seq.Slice);
end


%% loop over (selected) Seq.Slice "modules"
for t = idxModule
  %% get default or derived readout parameters
  if isemptyfield(Seq.Slice(t), 'useAQSlice'), Seq.Slice(t).useAQSlice = 1; end
  if isemptyfield(Seq.Slice(t), 'UseCoordinate'), Seq.Slice(t).UseCoordinate = 1; end
  if isemptyfield(Seq.Slice(t), 'iDevice')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Slice(t).useAQSlice), 'iDevice')
      Seq.Slice(t).iDevice = 1;
    else
      Seq.Slice(t).iDevice = Seq.AQSlice(Seq.Slice(t).useAQSlice).iDevice;
    end
  end
  if ~isfield(Seq.Slice, 'Pulse'), Seq.Slice(t).Pulse = struct(); end
  if isemptyfield(Seq.Slice(t).Pulse, 'MaxNumberOfSegments'), Seq.Slice(t).Pulse.MaxNumberOfSegments = 51; end
  if isemptyfield(Seq.Slice(t), 'Gamma')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Slice(t).useAQSlice), 'Gamma')
      Seq.Slice(t).Gamma = HW.GammaDef;
    else
      Seq.Slice(t).Gamma = Seq.AQSlice(Seq.Slice(t).useAQSlice).Gamma;
    end
  end
  if ~isfield(Seq.Slice(t), 'GammaX')
    % gamma for simultaneous pulse for X nuclei
    % empty for no dual frequency
    Seq.Slice(t).GammaX = [];
  end

  if isemptyfield(Seq.Slice(t).Pulse, 'FlipAngle'), Seq.Slice(t).Pulse.FlipAngle = 90; end
  if ~isempty(Seq.Slice(t).GammaX) && isemptyfield(Seq.Slice(t).Pulse, 'FlipAngleX')
    Seq.Slice(t).Pulse.FlipAngleX = Seq.Slice(t).Pulse.FlipAngle;
  end

  if isemptyfield(Seq.Slice(t).Pulse, 'Amplitude')
    if isempty(Seq.Slice(t).GammaX)
      Seq.Slice(t).Pulse.Amplitude = HW.TX(Seq.Slice(t).iDevice).AmpDef;
    else
      % For dual frequency pulses, try to stay at the limits set in HW.TX.Def
      % for the sum of the two constituent pulses. In general, this is a
      % non-linear problem when different stages of the amplifier chain have
      % calibration files associated to them. Simplify the calculations by
      % assuming that none of the set limits will be exceeded when the effective
      % "Norm" limit is fulfilled.

      % FIXME: Do we need a better way of setting the limits if dual frequency
      % pulses are used?

      % default amplitudes for each single pulse (in the absense of the other)
      defAmp0 = get_TX_Amplitude(HW, ...
        'Norm', HW.TX(Seq.Slice(t).iDevice).Def.NormCalibrated(HW.TX(Seq.Slice(t).iDevice).ChannelDef), ...
        'Frequency', [Seq.Slice(t).Gamma, Seq.Slice(t).GammaX]/(2*pi)*HW.B0);

      % Scale pulse amplitude such that the "Norm" limit is still respected when
      % both pulses are emitted. That assumes that the ratios at the "Norm"
      % stage scale linearly to the ratios at the "Amplitude" stage. That might
      % not be the case when calibration files are used for intermediate stages.
      Seq.Slice(t).Pulse.Amplitude = 1 / ...
        sum(1 ./ (defAmp0 .* [1, Seq.Slice(t).Pulse.FlipAngle/Seq.Slice(t).Pulse.FlipAngleX]));
    end
  end

  if isemptyfield(Seq.Slice(t).Pulse, 'Function'), Seq.Slice(t).Pulse.Function = @Pulse_Rect; end
  if isemptyfield(Seq.Slice(t), 'tRamp'), Seq.Slice(t).tRamp = HW.Grad(Seq.Slice(t).iDevice).tRamp; end
  if isemptyfield(Seq.Slice(t), 'tEC'), Seq.Slice(t).tEC = HW.Grad(Seq.Slice(t).iDevice).tEC; end
  if isemptyfield(Seq.Slice(t), 'MaxGradAmp'), Seq.Slice(t).MaxGradAmp = min(HW.Grad(Seq.Slice(t).iDevice).MaxAmp(1:3)); end
  if isemptyfield(Seq.Slice(t), 'Thickness'), Seq.Slice(t).Thickness = Inf; end
  if isemptyfield(Seq.Slice(t), 'CenterOfPulse'), Seq.Slice(t).CenterOfPulse = 0; end
  if isemptyfield(Seq.Slice(t), 'GradSign'), Seq.Slice(t).GradSign = 1; end
  if isemptyfield(Seq.Slice(t), 'distance'), Seq.Slice(t).distance = 0; end
  if isempty(Seq.Slice(t).GammaX)
    tmpGammaX = Seq.Slice(t).Gamma;
  else
    tmpGammaX = Seq.Slice(t).GammaX;
  end

  if ~isemptyfield(Seq.Slice(t), 'GradAmp') && isemptyfield(Seq.Slice(t).Pulse, 'Bandwidth')
    Seq.Slice(t).GradLength = ...
      max(1./(min(Seq.Slice(t).GradAmp(:)) * Seq.Slice(t).Gamma/2/pi * Seq.Slice(t).Thickness) * Seq.Slice(t).Pulse.Function(HW, 'Time') ...
          * max(1, Seq.Slice(t).Gamma/tmpGammaX), ...  % "replace" Gamma by GammaX if that requires a longer rf pulse
          Seq.Slice(t).Pulse.Function(HW, 'Length')) ...
        + 2*Seq.Slice(t).tRamp + 2*Seq.Slice(t).tEC;

    Seq.Slice(t).Pulse.MaxLength = Seq.Slice(t).GradLength - 2*Seq.Slice(t).tRamp - 2*Seq.Slice(t).tEC;
    Seq.Slice(t).GradLength = max(Seq.Slice(t).GradLength);

    % calculate bandwidth from gradient amplitude and thickness
    % but limit to maximum rf amplitude
    Seq.Slice(t).Pulse.Bandwidth = ...
      min(abs(Seq.Slice(t).GradAmp) .* Seq.Slice(t).Gamma/(2*pi) .* Seq.Slice(t).Thickness, ...
          HW.TX(Seq.Slice(t).iDevice).AmpDef .* Seq.Slice(t).Gamma/(2*pi) ./ (Seq.Slice(t).Pulse.FlipAngle/360) ...
          * Seq.Slice(t).Pulse.Function(HW, 'Time') / Seq.Slice(t).Pulse.Function(HW, 'Amp')) ...
      / max(1, Seq.Slice(t).Gamma/tmpGammaX);  % "revert" Gamma replacement from above, we need the bandwidth of the Gamma pulse
    Seq.Slice(t).Pulse.MaxLength = Inf;
  else
    Seq.Slice(t).GradLength = ...
      max(max(Seq.Slice(t).Pulse.FlipAngle./(Seq.Slice(t).Pulse.Amplitude * HW.TX(Seq.Slice(t).iDevice).Amp2Deg) * Seq.Slice(t).Pulse.Function(HW, 'Amp'), ...
              1./(min(min(HW.Grad(Seq.Slice(t).iDevice).MaxAmp(1:3)), ...
                      min(Seq.Slice(t).MaxGradAmp(:))) * Seq.Slice(t).Gamma/2/pi * Seq.Slice(t).Thickness) * Seq.Slice(t).Pulse.Function(HW, 'Time')) ...
          * max(1, Seq.Slice(t).Gamma/tmpGammaX), ...  % "replace" Gamma by GammaX if that requires a longer rf pulse
          Seq.Slice(t).Pulse.Function(HW, 'Length')) ...
        + 2*Seq.Slice(t).tRamp + 2*Seq.Slice(t).tEC;

    Seq.Slice(t).Pulse.MaxLength = Seq.Slice(t).GradLength - 2*Seq.Slice(t).tRamp - 2*Seq.Slice(t).tEC;
    Seq.Slice(t).GradLength = max(Seq.Slice(t).GradLength);

    % calculate gradient amplitude from thickness and bandwidth
    Seq.Slice(t).Pulse.Bandwidth = ...
      max(Seq.Slice(t).Pulse.Function(HW, 'Time') ./ (Seq.Slice(t).Pulse.MaxLength ...
          / max(1, Seq.Slice(t).Gamma/tmpGammaX)), ...  % "revert" Gamma replacement from above, we need the bandwidth of the Gamma pulse
          Seq.Slice(t).Pulse.Function(HW, 'Bandwidth'));  % might be a "limited bandwidth" pulse
    % Calculate gradient strength during slice pulse
    Seq.Slice(t).GradAmp = Seq.Slice(t).Pulse.Bandwidth .* (2*pi/Seq.Slice(t).Gamma) ./ Seq.Slice(t).Thickness .* Seq.Slice(t).GradSign;
    Seq.Slice(t).Pulse.MaxLength = max(Seq.Slice(t).Pulse.MaxLength);
  end
  % Calculate frequency of the TX slice pulse. It has an offset if the slice
  % does not go through the center.
  Seq.Slice(t).Pulse.Frequency = (Seq.Slice(t).Gamma/2/pi) .* (HW.B0 - Seq.Slice(t).distance .* Seq.Slice(t).GradAmp);
  if ~isempty(Seq.Slice(t).GammaX)
    Seq.Slice(t).Pulse.FrequencyX = (Seq.Slice(t).GammaX/2/pi) .* (HW.B0 - Seq.Slice(t).distance .* Seq.Slice(t).GradAmp);
  end

  if abs(nargin(Seq.Slice(t).Pulse.Function)) > 4
    maxDuration = 0;
    for tt = 1:numel(Seq.Slice(t).GradAmp)
      TXs(tt) = Seq.Slice(t).Pulse.Function(HW, ...
                                            Seq.Slice(t).CenterOfPulse, ...
                                            Seq.Slice(t).Pulse.Bandwidth(tt), ...
                                            Seq.Slice(t).Pulse.FlipAngle(tt)/180*pi, ...
                                            Seq.Slice(t).Pulse.MaxNumberOfSegments, ...
                                            Seq.Slice(t).Pulse.MaxLength, ...
                                            Seq.Slice(t).Pulse.Frequency(tt), 0);
      maxDuration = max(maxDuration, sum(TXs(tt).Duration));
    end
    if ~isempty(Seq.Slice(t).GammaX)
      for tt = 1:numel(Seq.Slice(t).GradAmp)
        % pulse shape functions should automatically select the matching gamma
        % FIXME: Better be explicit! (off-center pulses!)
        % FIXME: Support different pulse shape functions for pulse at primary
        %        and secondary frequencies.
        TXsX(tt) = Seq.Slice(t).Pulse.Function(HW, ...
                                              Seq.Slice(t).CenterOfPulse, ...
                                              Seq.Slice(t).Pulse.Bandwidth(tt)*Seq.Slice(t).GammaX/Seq.Slice(t).Gamma, ...
                                              Seq.Slice(t).Pulse.FlipAngle(tt)/180*pi, ...
                                              Seq.Slice(t).Pulse.MaxNumberOfSegments, ...
                                              Seq.Slice(t).Pulse.MaxLength, ...
                                              Seq.Slice(t).Pulse.FrequencyX(tt), 0);
        maxDuration = max(maxDuration, sum(TXsX(tt).Duration));
      end
    end
  else
    Pulse = Seq.Slice(t).Pulse;
    Pulse.Phase = 0;
    Pulse.FlipAngleFullTurn = 360; % FlipAngle is passed in degrees
    maxDuration = 0;
    if isemptyfield(Pulse, 'iDevice'), Pulse.iDevice = Seq.Slice(t).iDevice; end
    for tt = 1:numel(Seq.Slice(t).GradAmp)
      Pulse.Bandwidth = Seq.Slice(t).Pulse.Bandwidth(tt);
      Pulse.FlipAngle = Seq.Slice(t).Pulse.FlipAngle(min(tt, end));
      Pulse.Frequency = Seq.Slice(t).Pulse.Frequency(tt);
      TXs(tt) = Seq.Slice(t).Pulse.Function(HW, Seq.Slice(t).CenterOfPulse, Pulse);
      maxDuration = max(maxDuration, sum(TXs(tt).Duration));
    end
    if ~isempty(Seq.Slice(t).GammaX)
      for tt = 1:numel(Seq.Slice(t).GradAmp)
        % pulse shape functions should automatically select the matching gamma
        % FIXME: Better be explicit! (off-center pulses!)
        % FIXME: Support different pulse shape functions for pulse at primary
        %        and secondary frequencies.
        Pulse.Bandwidth = Seq.Slice(t).Pulse.Bandwidth(tt)*Seq.Slice(t).GammaX/Seq.Slice(t).Gamma;
        Pulse.FlipAngle = Seq.Slice(t).Pulse.FlipAngleX(min(tt, end));
        Pulse.Frequency = Seq.Slice(t).Pulse.FrequencyX(tt);
        TXsX(tt) = Seq.Slice(t).Pulse.Function(HW, Seq.Slice(t).CenterOfPulse, Pulse);
        maxDuration = max(maxDuration, sum(TXsX(tt).Duration));
      end
    end
  end
  % return actual maximum length
  Seq.Slice(t).Pulse.MaxLength = maxDuration;

  TXs = correct_PulsePhase(TXs, HW, Seq.Slice(t).iDevice, (Seq.Slice(t).Gamma/2/pi) .* HW.B0);
  Seq.Slice(t).Pulse.CenterOffset = TXs(1).CenterOffset;
  if ~isempty(Seq.Slice(t).GammaX)
    TXsX = correct_PulsePhase(TXsX, HW, Seq.Slice(t).iDevice, (Seq.Slice(t).GammaX/2/pi) .* HW.B0);
  end

  % Calculate integral of the slice gradient
  Seq.Slice(t).GradTimeIntegral = abs(Seq.Slice(t).GradAmp * (Seq.Slice(t).GradLength - Seq.Slice(t).tRamp));
  Seq.Slice(t).GradTimeIntegralRephase = abs(Seq.Slice(t).GradAmp * ((Seq.Slice(t).GradLength-Seq.Slice(t).tRamp)/2 - Seq.Slice(t).Pulse.CenterOffset));
  Seq.Slice(t).GradTimeIntegralDephase = abs(Seq.Slice(t).GradAmp * ((Seq.Slice(t).GradLength-Seq.Slice(t).tRamp)/2 + Seq.Slice(t).Pulse.CenterOffset));
  Seq.Slice(t).GradCenter = Seq.Slice(t).CenterOfPulse - Seq.Slice(t).Pulse.CenterOffset;

  if Seq.Slice(t).GradTimeIntegral==0
    tD = 0.5;
    tR = 0.5;
  else
    tD = Seq.Slice(t).GradTimeIntegralDephase/Seq.Slice(t).GradTimeIntegral;
    tR = Seq.Slice(t).GradTimeIntegralRephase/Seq.Slice(t).GradTimeIntegral;
  end

  if isemptyfield(Seq.Slice(t), 'GradTimeDelay'), Seq.Slice(t).GradTimeDelay = [0,0,0]; end
  if isemptyfield(Seq.Slice(t), 'GradRephaseLength')
    % approx. half the slice amplitude, but not shorter than 2*tRamp
    Seq.Slice(t).GradRephaseLength = max(...
      (Seq.Slice(t).GradLength-Seq.Slice(t).tRamp)*tR + Seq.Slice(t).tRamp, ...
      2*Seq.Slice(t).tRamp + 4/HW.MMRT(Seq.Slice(t).iDevice).fSystem);
  end
  if isemptyfield(Seq.Slice(t), 'CenterOfRephase')
    Seq.Slice(t).CenterOfRephase = Seq.Slice(t).GradLength/2-Seq.Slice(t).Pulse.CenterOffset+Seq.Slice(t).GradRephaseLength/2+Seq.Slice(t).CenterOfPulse-Seq.Slice(t).tRamp;
  end
  if isemptyfield(Seq.Slice(t), 'GradRephaseSign'), Seq.Slice(t).GradRephaseSign = -1; end
  if isemptyfield(Seq.Slice(t), 'GradDephaseLength')
    % approx. half the slice amplitude, but not shorter than 2*tRamp
    Seq.Slice(t).GradDephaseLength =  max(...
      (Seq.Slice(t).GradLength-Seq.Slice(t).tRamp)*tD + Seq.Slice(t).tRamp, ...
      2*Seq.Slice(t).tRamp + 4/HW.MMRT(Seq.Slice(t).iDevice).fSystem);
  end
  if isemptyfield(Seq.Slice(t), 'CenterOfDephase')
    Seq.Slice(t).CenterOfDephase = -Seq.Slice(t).GradLength/2-Seq.Slice(t).Pulse.CenterOffset-Seq.Slice(t).GradDephaseLength/2+Seq.Slice(t).CenterOfPulse+Seq.Slice(t).tRamp;
  end
  if isemptyfield(Seq.Slice(t), 'GradDephaseSign'), Seq.Slice(t).GradDephaseSign = -1; end
  % if isemptyfield(Seq.Slice(t), 'Thickness'), Seq.Slice(t).Thickness = Seq.AQSlice(Seq.Slice(t).useAQSlice).thickness; end
  if isemptyfield(Seq.Slice(t), 'alfa')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Slice(t).useAQSlice), 'alfa')
      Seq.Slice(t).alfa = 0;
    else
      Seq.Slice(t).alfa = Seq.AQSlice(Seq.Slice(t).useAQSlice).alfa;
    end
  end
  if isemptyfield(Seq.Slice(t), 'phi')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Slice(t).useAQSlice), 'phi')
      Seq.Slice(t).phi = 0;
    else
      Seq.Slice(t).phi = Seq.AQSlice(Seq.Slice(t).useAQSlice).phi;
    end
  end
  if isemptyfield(Seq.Slice(t), 'theta')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Slice(t).useAQSlice), 'theta')
      Seq.Slice(t).theta = 0;
    else
      Seq.Slice(t).theta = Seq.AQSlice(Seq.Slice(t).useAQSlice).theta;
    end
  end
  if isemptyfield(Seq.Slice(t), 'angle2Turns')  % conversion factor of the input angles to full turns e.g.: 1/(2*pi)
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Slice(t).useAQSlice), 'angle2Turns')
      Seq.Slice(t).angle2Turns = 1/(2*pi);
    else
      Seq.Slice(t).angle2Turns = Seq.AQSlice(Seq.Slice(t).useAQSlice).angle2Turns;
    end
  end
  if isemptyfield(Seq.Slice(t).Pulse, 'PhaseIncrement'), Seq.Slice(t).Pulse.PhaseIncrement = 0; end
  if isemptyfield(Seq.Slice(t).Pulse, 'Phase'), Seq.Slice(t).Pulse.Phase = 0; end
  if ~isempty(Seq.Slice(t).GammaX) && isemptyfield(Seq.Slice(t).Pulse, 'PhaseX')
    Seq.Slice(t).Pulse.PhaseX = Seq.Slice(t).Pulse.Phase;
  end
  if isemptyfield(Seq.Slice(t), 'UseAtRepetitionTime'), Seq.Slice(t).UseAtRepetitionTime = []; end % !!!
  if isemptyfield(Seq.Slice(t), 'UseAtRepetitionTimeDephase'), Seq.Slice(t).UseAtRepetitionTimeDephase = Seq.Slice(t).UseAtRepetitionTime; end
  if isemptyfield(Seq.Slice(t), 'UseAtRepetitionTimeRephase'), Seq.Slice(t).UseAtRepetitionTimeRephase = Seq.Slice(t).UseAtRepetitionTime; end
  if isemptyfield(Seq.Slice(t), 'Overdrive'), Seq.Slice(t).Overdrive = 0; end
  if isemptyfield(Seq.Slice(t), 'GradTimeIntegralRephaseOffset'), Seq.Slice(t).GradTimeIntegralRephaseOffset = 0; end
  if isemptyfield(Seq.Slice(t), 'GradTimeIntegralDephaseOffset'), Seq.Slice(t).GradTimeIntegralDephaseOffset = 0; end

  if Seq.Slice(t).Thickness <= 1000 && ...
      round((Seq.Slice(t).GradDephaseLength-Seq.Slice(t).tRamp*2)*HW.MMRT(Seq.Slice(t).iDevice).fSystem) < 2
    error('PD:get_SliceParameter:GradDephaseLengthTooShort', ...
      'Seq.Slice(%d).GradDephaseLength too short to fit 2*tRamp', t);
  end
  if Seq.Slice(t).Thickness <= 1000 && ...
      round((Seq.Slice(t).GradRephaseLength-Seq.Slice(t).tRamp*2)*HW.MMRT(Seq.Slice(t).iDevice).fSystem) < 2
    error('PD:get_SliceParameter:GradRephaseLengthTooShort', ...
      'Seq.Slice(%d).GradRephaseLength too short to fit 2*tRamp', t);
  end


  % Calculate gradient amplitudes during slice rephase and dephase
  Seq.Slice(t).GradAmpRephase = (Seq.Slice(t).GradTimeIntegralRephase+Seq.Slice(t).GradTimeIntegralRephaseOffset) ./ ...
                                (Seq.Slice(t).GradRephaseLength-Seq.Slice(t).tRamp);
  Seq.Slice(t).GradAmpDephase = (Seq.Slice(t).GradTimeIntegralDephase+Seq.Slice(t).GradTimeIntegralDephaseOffset) ./ ...
                                (Seq.Slice(t).GradDephaseLength-Seq.Slice(t).tRamp);

  % maxsize=max([numel(Seq.Slice(t).UseAtRepetitionTime);numel(Seq.Slice(t).Pulse.Phase);numel(Seq.Slice(t).Pulse.PhaseIncrement)]);

  % if and(Seq.Slice(t).Pulse.PhaseIncrement==0,maxsize==1)
  %   tempsize=1;
  %   tempPhaseInc=0;
  % else
  %   tempsize=ones(1,size(Seq.tRep,2));
  tempsize = ones(1, numel(Seq.Slice(t).UseAtRepetitionTime));

  if numel(Seq.Slice(t).Pulse.PhaseIncrement)>1
    temp = repmat(Seq.Slice(t).Pulse.PhaseIncrement(:), 1, ...
      numel(Seq.Slice(t).UseAtRepetitionTime)/numel(Seq.Slice(t).Pulse.PhaseIncrement));
    tempPhaseInc = cumsum(temp(:)).';
  else
    tempPhaseInc = cumsum(Seq.Slice(t).Pulse.PhaseIncrement * tempsize);
  end
  if numel(Seq.Slice(t).Pulse.Phase) > 1
    temp = repmat(Seq.Slice(t).Pulse.Phase(:), 1, ...
      numel(Seq.Slice(t).UseAtRepetitionTime)/numel(Seq.Slice(t).Pulse.Phase));
    tempPhase = (temp(:)).';
  else
    tempPhase = Seq.Slice(t).Pulse.Phase * tempsize;
  end

  % end


  %% translate to pulse program matrices
  Seq.Slice(t).TX(1).Start = NaN(size(TXs(1).Start,1), size(Seq.tRep,2));
  Seq.Slice(t).TX(1).Duration = NaN(size(TXs(1).Duration,1), size(Seq.tRep,2));
  Seq.Slice(t).TX(1).Amplitude = NaN(size(TXs(1).Amplitude,1), size(Seq.tRep,2));
  Seq.Slice(t).TX(1).Frequency = NaN(size(TXs(1).Frequency,1), size(Seq.tRep,2));
  Seq.Slice(t).TX(1).Phase = NaN(size(TXs(1).Phase,1), size(Seq.tRep,2));

  Seq.Slice(t).TX(1).Start(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXs.Start], tempsize);
  Seq.Slice(t).TX(1).Duration(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXs.Duration], tempsize);
  Seq.Slice(t).TX(1).Amplitude(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXs.Amplitude], tempsize);
  Seq.Slice(t).TX(1).Frequency(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXs.Frequency], tempsize);
  Seq.Slice(t).TX(1).Phase(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@plus, [TXs.Phase], ones(size(TXs(1).Phase,1),1) * tempPhaseInc + ...
    ones(size(TXs(1).Phase,1),1)*tempPhase);
  clear TXs

  if ~isempty(Seq.Slice(t).GammaX)
    % pulse at X nucleus frequency
    if numel(Seq.Slice(t).Pulse.PhaseX) > 1
      temp = repmat(Seq.Slice(t).Pulse.PhaseX(:), ...
        1, numel(Seq.Slice(t).UseAtRepetitionTime)/numel(Seq.Slice(t).Pulse.PhaseX));
      tempPhase = (temp(:)).';
    else
      tempPhase = Seq.Slice(t).Pulse.PhaseX * tempsize;
    end

    Seq.Slice(t).TX(2).Start = NaN(size(TXsX(1).Start,1), size(Seq.tRep,2));
    Seq.Slice(t).TX(2).Duration = NaN(size(TXsX(1).Duration,1), size(Seq.tRep,2));
    Seq.Slice(t).TX(2).Amplitude = NaN(size(TXsX(1).Amplitude,1), size(Seq.tRep,2));
    Seq.Slice(t).TX(2).Frequency = NaN(size(TXsX(1).Frequency,1), size(Seq.tRep,2));
    Seq.Slice(t).TX(2).Phase = NaN(size(TXsX(1).Phase,1), size(Seq.tRep,2));

    Seq.Slice(t).TX(2).Start(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXsX.Start], tempsize);
    Seq.Slice(t).TX(2).Duration(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXsX.Duration], tempsize);
    Seq.Slice(t).TX(2).Amplitude(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXsX.Amplitude], tempsize);
    Seq.Slice(t).TX(2).Frequency(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXsX.Frequency], tempsize);
    Seq.Slice(t).TX(2).Phase(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@plus, [TXsX.Phase], ones(size(TXsX(1).Phase,1),1) * tempPhaseInc + ...
      ones(size(TXsX(1).Phase,1),1)*tempPhase);
    clear TXsX
  end

  Angle2Deg = Seq.Slice(t).angle2Turns*360;
  [Rx, Ry, Rz] = get_aptDegRotationMatrix(Seq.Slice(t).alfa*Angle2Deg, Seq.Slice(t).phi*Angle2Deg, Seq.Slice(t).theta*Angle2Deg);

  temp = zeros(3,1);
  temp(Seq.Slice(t).UseCoordinate) = 1;
  Seq.Slice(t).GradAmpUnitVector = temp;
  temp = zeros(3,1);
  temp(Seq.Slice(t).UseCoordinate) = 1;
  Seq.Slice(t).GradAmpRephaseUnitVector = temp;
  temp = zeros(3,1);
  temp(Seq.Slice(t).UseCoordinate) = 1;
  Seq.Slice(t).GradAmpDephaseUnitVector = temp;

  Seq.Slice(t).GradAmpUnitVector = Rz*(Ry*(Rx*Seq.Slice(t).GradAmpUnitVector));
  Seq.Slice(t).GradAmpRephaseUnitVector = Rz*(Ry*(Rx*Seq.Slice(t).GradAmpRephaseUnitVector));
  Seq.Slice(t).GradAmpDephaseUnitVector = Rz*(Ry*(Rx*Seq.Slice(t).GradAmpDephaseUnitVector));

  for n = 1:3
    % if Seq.Slice(t).Overdrive == 0
    %   Seq.Slice(t).Grad(n).Amp=zeros(4,size(tempsize,2));
    %   Seq.Slice(t).GradRephase(n).Amp=zeros(4,size(tempsize,2));
    %   Seq.Slice(t).GradDephase(n).Amp=zeros(4,size(tempsize,2));
    %   Seq.Slice(t).GradRephase(n).Amp(2:3,:)=Seq.Slice(t).GradAmpRephase(n);
    %   Seq.Slice(t).GradDephase(n).Amp(2:3,:)=Seq.Slice(t).GradAmpDephase(n);
    %   Seq.Slice(t).Grad(n).Amp=Seq.Slice(t).Grad(n).Amp(1:4,1)*tempsize;
    %   Seq.Slice(t).GradRephase(n).Amp=Seq.Slice(t).GradRephase(n).Amp(1:4,1)*tempsize;
    %   Seq.Slice(t).GradDephase(n).Amp=Seq.Slice(t).GradDephase(n).Amp(1:4,1)*tempsize;
    Seq.Slice(t).Grad(n).Time=nan(4, size(Seq.tRep,2));
    Seq.Slice(t).GradRephase(n).Time=Seq.Slice(t).Grad(n).Time;
    Seq.Slice(t).GradDephase(n).Time=Seq.Slice(t).Grad(n).Time;

    Seq.Slice(t).Grad(n).Amp=nan(4, size(Seq.tRep,2));
    Seq.Slice(t).GradRephase(n).Amp=Seq.Slice(t).Grad(n).Amp;
    Seq.Slice(t).GradDephase(n).Amp=Seq.Slice(t).Grad(n).Amp;

    Seq.Slice(t).Grad(n).Amp(1,Seq.Slice(t).UseAtRepetitionTime) = 0;
    Seq.Slice(t).Grad(n).Amp(2,Seq.Slice(t).UseAtRepetitionTime) = Seq.Slice(t).GradAmpUnitVector(n).*Seq.Slice(t).GradAmp;
    Seq.Slice(t).Grad(n).Amp(3,Seq.Slice(t).UseAtRepetitionTime) = Seq.Slice(t).GradAmpUnitVector(n).*Seq.Slice(t).GradAmp;
    Seq.Slice(t).Grad(n).Amp(4,Seq.Slice(t).UseAtRepetitionTime) = 0;

    Seq.Slice(t).GradRephase(n).Amp(1,Seq.Slice(t).UseAtRepetitionTimeRephase) = 0;
    Seq.Slice(t).GradRephase(n).Amp(2,Seq.Slice(t).UseAtRepetitionTimeRephase) = Seq.Slice(t).GradAmpRephaseUnitVector(n).*Seq.Slice(t).GradAmpRephase.* Seq.Slice(t).GradRephaseSign;
    Seq.Slice(t).GradRephase(n).Amp(3,Seq.Slice(t).UseAtRepetitionTimeRephase) = Seq.Slice(t).GradAmpRephaseUnitVector(n).*Seq.Slice(t).GradAmpRephase.* Seq.Slice(t).GradRephaseSign;
    Seq.Slice(t).GradRephase(n).Amp(4,Seq.Slice(t).UseAtRepetitionTimeRephase) = 0;

    Seq.Slice(t).GradDephase(n).Amp(1,Seq.Slice(t).UseAtRepetitionTimeDephase) = 0;
    Seq.Slice(t).GradDephase(n).Amp(2,Seq.Slice(t).UseAtRepetitionTimeDephase) = Seq.Slice(t).GradAmpDephaseUnitVector(n).*Seq.Slice(t).GradAmpDephase.* Seq.Slice(t).GradDephaseSign;
    Seq.Slice(t).GradDephase(n).Amp(3,Seq.Slice(t).UseAtRepetitionTimeDephase) = Seq.Slice(t).GradAmpDephaseUnitVector(n).*Seq.Slice(t).GradAmpDephase.* Seq.Slice(t).GradDephaseSign;
    Seq.Slice(t).GradDephase(n).Amp(4,Seq.Slice(t).UseAtRepetitionTimeDephase) = 0;


    Seq.Slice(t).Grad(n).Time(:,Seq.Slice(t).UseAtRepetitionTime) = ...
      Seq.Slice(t).GradCenter + [ -Seq.Slice(t).GradLength/2; ...
                                  -Seq.Slice(t).GradLength/2+Seq.Slice(t).tRamp; ...
                                  +Seq.Slice(t).GradLength/2-Seq.Slice(t).tRamp; ...
                                  +Seq.Slice(t).GradLength/2]*tempsize - ...
      Seq.Slice(t).GradTimeDelay(n);
    Seq.Slice(t).GradRephase(n).Time(:,Seq.Slice(t).UseAtRepetitionTimeRephase) = ...
      Seq.Slice(t).CenterOfRephase + [ -Seq.Slice(t).GradRephaseLength/2; ...
                                       -Seq.Slice(t).GradRephaseLength/2+Seq.Slice(t).tRamp; ...
                                       +Seq.Slice(t).GradRephaseLength/2-Seq.Slice(t).tRamp; ...
                                       +Seq.Slice(t).GradRephaseLength/2]*ones(1, numel(Seq.Slice(t).UseAtRepetitionTimeRephase)) - ...
      Seq.Slice(t).GradTimeDelay(n);
    Seq.Slice(t).GradDephase(n).Time(:,Seq.Slice(t).UseAtRepetitionTimeDephase) = ...
      Seq.Slice(t).CenterOfDephase + [ -Seq.Slice(t).GradDephaseLength/2; ...
                                       -Seq.Slice(t).GradDephaseLength/2+Seq.Slice(t).tRamp; ...
                                       +Seq.Slice(t).GradDephaseLength/2-Seq.Slice(t).tRamp; ...
                                       +Seq.Slice(t).GradDephaseLength/2]*ones(1, numel(Seq.Slice(t).UseAtRepetitionTimeDephase)) - ...
      Seq.Slice(t).GradTimeDelay(n);

    % else
    %   error('alt')
    %   Seq.Slice(t).Grad(n).Amp=zeros(8,size(tempsize,2));
    %   Seq.Slice(t).GradRephase(n).Amp=zeros(8,size(tempsize,2));
    %   Seq.Slice(t).GradDephase(n).Amp=zeros(8,size(tempsize,2));
    %   Seq.Slice(t).Grad(n).Amp(2:3,:)=Seq.Slice(t).GradAmp(n)*(1+Seq.Slice(t).Overdrive);
    %   Seq.Slice(t).GradRephase(n).Amp(2:3,:)=Seq.Slice(t).GradAmpRephase(n)*(1+Seq.Slice(t).Overdrive);
    %   Seq.Slice(t).GradDephase(n).Amp(2:3,:)=Seq.Slice(t).GradAmpDephase(n)*(1+Seq.Slice(t).Overdrive);
    %   Seq.Slice(t).Grad(n).Amp(4:5,:)=Seq.Slice(t).GradAmp(n);
    %   Seq.Slice(t).GradRephase(n).Amp(4:5,:)=Seq.Slice(t).GradAmpRephase(n);
    %   Seq.Slice(t).GradDephase(n).Amp(4:5,:)=Seq.Slice(t).GradAmpDephase(n);
    %   Seq.Slice(t).Grad(n).Amp(6:7,:)=Seq.Slice(t).GradAmp(n)*(0-Seq.Slice(t).Overdrive);
    %   Seq.Slice(t).GradRephase(n).Amp(6:7,:)=Seq.Slice(t).GradAmpRephase(n)*(0-Seq.Slice(t).Overdrive);
    %   Seq.Slice(t).GradDephase(n).Amp(6:7,:)=Seq.Slice(t).GradAmpDephase(n)*(0-Seq.Slice(t).Overdrive);
    %   Seq.Slice(t).Grad(n).Amp=Seq.Slice(t).Grad(n).Amp(1:8,1)*tempsize;
    %   Seq.Slice(t).GradRephase(n).Amp=Seq.Slice(t).GradRephase(n).Amp(1:8,1)*tempsize;
    %   Seq.Slice(t).GradDephase(n).Amp=Seq.Slice(t).GradDephase(n).Amp(1:8,1)*tempsize;
    %
    %   Seq.Slice(t).Grad(n).Time=Seq.Slice(t).Pulse.Center...
    %                                                       +[  -Seq.Slice(t).GradLength/2;...
    %                                                           -Seq.Slice(t).GradLength/2+Seq.Slice(t).tRamp*1/3;...
    %                                                           -Seq.Slice(t).GradLength/2+Seq.Slice(t).tRamp*2/3;...
    %                                                           -Seq.Slice(t).GradLength/2+Seq.Slice(t).tRamp*3/3;...
    %                                                           +Seq.Slice(t).GradLength/2-Seq.Slice(t).tRamp*3/3;...
    %                                                           +Seq.Slice(t).GradLength/2-Seq.Slice(t).tRamp*2/3;...
    %                                                           +Seq.Slice(t).GradLength/2-Seq.Slice(t).tRamp*1/3;...
    %                                                           +Seq.Slice(t).GradLength/2]*tempsize-Seq.Slice(t).GradTimeDelay(n);
    %   Seq.Slice(t).GradRephase(n).Time=Seq.Slice(t).CenterOfRephase...
    %                                                       +[  -Seq.Slice(t).GradRephaseLength/2;...
    %                                                           -Seq.Slice(t).GradRephaseLength/2+Seq.Slice(t).tRamp*1/3;...
    %                                                           -Seq.Slice(t).GradRephaseLength/2+Seq.Slice(t).tRamp*2/3;...
    %                                                           -Seq.Slice(t).GradRephaseLength/2+Seq.Slice(t).tRamp*3/3;...
    %                                                           +Seq.Slice(t).GradRephaseLength/2-Seq.Slice(t).tRamp*3/3;...
    %                                                           +Seq.Slice(t).GradRephaseLength/2-Seq.Slice(t).tRamp*2/3;...
    %                                                           +Seq.Slice(t).GradRephaseLength/2-Seq.Slice(t).tRamp*1/3;...
    %                                                           +Seq.Slice(t).GradRephaseLength/2]*tempsize-Seq.Slice(t).GradTimeDelay(n);
    %   Seq.Slice(t).GradDephase(n).Time=Seq.Slice(t).CenterOfDephase...
    %                                                       +[  -Seq.Slice(t).GradDephaseLength/2;...
    %                                                           -Seq.Slice(t).GradDephaseLength/2+Seq.Slice(t).tRamp*1/3;...
    %                                                           -Seq.Slice(t).GradDephaseLength/2+Seq.Slice(t).tRamp*2/3;...
    %                                                           -Seq.Slice(t).GradDephaseLength/2+Seq.Slice(t).tRamp*3/3;...
    %                                                           +Seq.Slice(t).GradDephaseLength/2-Seq.Slice(t).tRamp*3/3;...
    %                                                           +Seq.Slice(t).GradDephaseLength/2-Seq.Slice(t).tRamp*2/3;...
    %                                                           +Seq.Slice(t).GradDephaseLength/2-Seq.Slice(t).tRamp*1/3;...
    %                                                           +Seq.Slice(t).GradDephaseLength/2]*tempsize-Seq.Slice(t).GradTimeDelay(n);
    % end
  end
end

end
