function Seq = get_ReadParameter(Seq, HW)
%% Calculate read encoder
%
%    Seq = get_ReadParameter(Seq, HW)
%
% This function calculates read gradients (and respective rephase and dephase
% gradients) derived from the values in Seq.Read. It also creates settings for
% the corresponding acquisition window.
%
%
% INPUT:
%
% "Seq" is a structure with the field "Read" that is an (array of) structures
% with the following fields:
%   UseCoordinate   Number of the coordinate for which the phase encoder is used
%                   (default: 2).
%   GradLength      Length of the read pulse in seconds (default:
%                   1/Seq.Read(t).HzPerPixelMin+Seq.Read(t).tEC*2+Seq.Read(t).tRamp*2; ).
%   GradSign        Sign of the read pulse (default: 1).
%   GradTimeIntegralOffset
%                   Offset for the time integral of the read pulse in
%                   s*T/m (default: 0).
%   tRamp           Ramp time of the gradients in seconds (default:
%                   HW.Grad.tRamp).
%   tEC             Additional in seconds time that the gradients need to reach
%                   their set amplitude due to Eddy currents in the pole shoe
%                   (default: HW.Grad.tEC).
%   HzPerPixelMin   Minimal Hertz per Pixel in Hz (default:
%                   1/(Seq.Read(t).GradLength-HW.Grad.tEC*2-Seq.Read(t).tRamp*2) )
%   GradTimeDelay   1x3 vector with time delays in seconds for each gradient
%                   channel. Positive values lead to the gradient pulses being
%                   set earlier. (Default: [0 0 0]).
%   Gamma           Gyromagnetic ratio of the sample in rad/s/T (default:
%                   HW.GammaDef).
%   StartWithKSpaceCenter
%                   If not equal to 0, this is a factor for a linear transition
%                   between compensating the read dephase at the center of the
%                   readout (limit -> 0) or at the center of the first acquired
%                   sample (=1). If equal to 0, the acquisition is positioned
%                   such that the k-space center matches the condition for a
%                   discrete Fourier transform. (Default: 0)
%   StartWithKSpaceCenterNoSampleShift
%                   If StartWithKSpaceCenter ~= 0, the acquisition window is not
%                   shifted such that the encoding is correct at the center of
%                   the respective sample, but it is shifted such that the
%                   encoding is correct at the start of each sample.
%                   (Default: 0)
%   StartWithKSpaceCenterNonlinear
%                   Start readout at start of read-out ramp up (default: false).
%   StartWithGradientBeforeExcitationPulse
%                   Start gradient before excitation pulse starts (default:
%                   false).
%   StartWithGradientTime
%                   Time in seconds when the gradient switches on (default: 0).
%   useAQSlice      Index that is used for Seq.AQSlice for the following default
%                   parameters (default: 1).
%   sizeRead        Size of the spectral encoding of the gradient pulses in
%                   meters (default:
%                   Seq.AQSlice(Seq.Read(t).useAQSlice).sizeRead(t))).
%   nRead           Number of points in read direction (default:
%                   Seq.AQSlice(Seq.Read(t).useAQSlice).nRead(t)).
%   ReadOS          Over-sampling factor in read direction (default:
%                   Seq.AQSlice(Seq.Read(t).useAQSlice).ReadOS(t)).
%   SamplingFactor  Integer sampling factor for an additional software
%                   downsampling that mimicks the characteristics of a CIC
%                   filter. (default: 1 or the lowest sampling factor to come
%                   closest to the required sampling rate)
%   alfa / phi / theta
%                   Angles in radians that activley rotate the direction of the
%                   read encoding (see UseCoordinate). The first rotation "alfa"
%                   is around the x-axis, the second rotation "phi" is around
%                   the y-axis and the last rotation "theta" is around the
%                   z-axis.
%   distance        Distance of the center of the encoded space to the center of
%                   the gradient system in meter (default: 0).
%   UseAtRepetitionTime
%                   Repetition times (tReps) at which the gradient pulses are
%                   used (no default!).
%   PhaseIncrement  Increment of the AQ phase (default: 0).
%   Phase           AQ phase in degrees (default: 0).
%   ResetPhases     unused
%   GetData         Boolean value. This can be a scalar or a vector of the size
%                   of UseAtRepetitionTime. If true, the measurement data can be
%                   collected at the end of that tRep. Otherwise, the
%                   measurement data can be collected at the end of the last
%                   tRep of the sequence (default: 0).
%   UseCoordinate   Scalar k-space coordinate along which the read-out gradient
%                   is applied. (Default: 2)
%   kLineDir        Direction in k-space along which the read-out gradient is
%                   applied. This can be a 3x1 vector in which case the same
%                   direction applies to all read-out windows or a 3xN matrix in
%                   which case a direction must be supplied for each read-out
%                   window. The amplitudes of these vectors scale the read-out
%                   gradient amplitude. (Default: unit vector in direction
%                   Seq.Read.UseCoordinate)
% For the rephase and dephase pulses:
%   GradDephaseLengthFactor / GradRephaseLengthFactor
%                   Factor to the default length of the dephase or rephase
%                   gradient pulse, respectively. The default length is such
%                   that the amplitude of the gradient pulse is the same as the
%                   readout gradient amplitude (potentially opposite sign). If
%                   set to 0, the minimum gradient duration (given the ramp
%                   time) is used. Note that using 0 might lead to unpredictable
%                   effects (e.g., due to eddy currents). (Default: 1)
%   GradDephaseLength / GradRephaseLength
%                   Length of the dephase or rephase gradient pulse in seconds
%                   (default: see above GradDephaseLengthFactor /
%                   GradRephaseLengthFactor).
%   CenterOfDephase / CenterOfRephase
%                   Center of dephase or rephase gradient pulse in seconds on
%                   the tRep timeline. (Default: Such that the ramps of the
%                   dephase/rephase pulse and the adjacent readout gradient
%                   coincide.)
%   GradDephaseSign / GradRephaseSign
%                   Sign of the rephase/dephase pulse (default: -1).
%
%
% OUTPUT:
%
% To each structure Seq.Read the following fields are added:
%   Grad / GradRephase / GradDephase
%                   Structure with the fields "Amp" and "Time" containing the
%                   amplitudes and times of the rephase/dephase pulses
%                   corresponding to  the input parameters that can be added to
%                   the sequence with "add_Grad".
%   AQ              Structure containing the settings for the acquisition
%                   window. It can be added to the sequence with "add_AQ".
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%%

if ~isfield(Seq.Read, 'UseCoordinate'), Seq.Read(1).UseCoordinate = 2; end
% if Seq.showSliceRead, Seq.Read(t).UseCoordinate = Seq.Slice(t).UseCoordinate; end
for t = 1:numel(Seq.Read)
  if isemptyfield(Seq.Read(t), 'UseAtRepetitionTime')
    % Element isn't used. Skip this loop, but make sure the created pulse
    % program elements are empty.
    Seq.Read(t).AQ = [];
    Seq.Read(t).Grad = [];
    Seq.Read(t).GradRephase = [];
    Seq.Read(t).GradDephase = [];
    continue;
  end
  if isempty(Seq.Read(t).UseCoordinate), Seq.Read(t).UseCoordinate = 2; end
  if isemptyfield(Seq.Read(t), 'useAQSlice'), Seq.Read(t).useAQSlice = 1; end
  if isemptyfield(Seq.Read(t), 'iDevice')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Read(t).useAQSlice), 'iDevice')
      Seq.Read(t).iDevice = 1;
    else
      Seq.Read(t).iDevice = Seq.AQSlice(Seq.Read(t).useAQSlice).iDevice;
    end
  end
  if isemptyfield(Seq.Read(t), 'tRamp'), Seq.Read(t).tRamp = HW.Grad(Seq.Read(t).iDevice).tRamp; end
  if isemptyfield(Seq.Read(t), 'tEC'), Seq.Read(t).tEC = HW.Grad(Seq.Read(t).iDevice).tEC; end

  if isemptyfield(Seq.Read(t), 'GradSign'), Seq.Read(t).GradSign = 1; end

  if isemptyfield(Seq.Read(t), 'GradTimeDelay'), Seq.Read(t).GradTimeDelay = [0, 0, 0]; end
  if isemptyfield(Seq.Read(t), 'StartWithKSpaceCenterNonlinear')
    Seq.Read(t).StartWithKSpaceCenterNonlinear = 0;
  end
  if isemptyfield(Seq.Read(t), 'StartWithKSpaceCenter')
    Seq.Read(t).StartWithKSpaceCenter = Seq.Read(t).StartWithKSpaceCenterNonlinear;
  end
  if isemptyfield(Seq.Read(t), 'StartWithKSpaceCenterNoSampleShift')
    Seq.Read(t).StartWithKSpaceCenterNoSampleShift = 0;
  end
  if isemptyfield(Seq.Read(t), 'StartWithGradientBeforeExcitationPulse')
    if ~isemptyfield(Seq.Read(t), 'StartWithGradientBeforExcitationPulse')
      % Backwards compatibility for typo in field name
      Seq.Read(t).StartWithGradientBeforeExcitationPulse = Seq.Read(t).StartWithGradientBeforExcitationPulse;
    else
      Seq.Read(t).StartWithGradientBeforeExcitationPulse = 0;
    end
  end
  if isemptyfield(Seq.Read(t), 'Gamma')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Read(t).useAQSlice), 'Gamma')
      Seq.Read(t).Gamma = HW.GammaDef;
    else
      Seq.Read(t).Gamma = Seq.AQSlice(Seq.Read(t).useAQSlice).Gamma;
    end
  end
  if ~isfield(Seq.Read(t), 'GammaX')
    % gamma for simultaneous acquisition at X nucleus frequency
    % empty for no dual frequency
    % FIXME: Return resulting sizeReadX?
    Seq.Read(t).GammaX = [];
  end
  if isemptyfield(Seq.Read(t), 'ReadOS')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Read(t).useAQSlice), 'ReadOS')
      Seq.Read(t).ReadOS = 2;
    else
      Seq.Read(t).ReadOS = Seq.AQSlice(Seq.Read(t).useAQSlice).ReadOS;
    end
  end
  if isemptyfield(Seq.Read(t), 'nRead')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Read(t).useAQSlice), 'nRead')
      error('PD:get_ReadParameter:NoNRead', ...
        'Seq.Read.nRead must be set.');
    else
      Seq.Read(t).nRead = Seq.AQSlice(Seq.Read(t).useAQSlice).nRead;
    end
  end
  if isemptyfield(Seq.Read(t), 'sizeRead')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Read(t).useAQSlice), 'sizeRead')
      error('PD:get_ReadParameter:NoSizeRead', ...
        'Seq.Read.sizeRead must be set.');
    else
      Seq.Read(t).sizeRead = Seq.AQSlice(Seq.Read(t).useAQSlice).sizeRead;
    end
  end
  if isemptyfield(Seq.Read(t), 'SamplingFactor')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Read(t).useAQSlice), 'SamplingFactor')
      % Calculate a sampling factor that allows the CIC decimation factor to be
      % in range.
      CIC_Decimation_without_SamplingFactor = ...
        floor((1/Seq.Read(t).HzPerPixelMin) ...
              / (Seq.Read(t).nRead*Seq.Read(t).ReadOS) ...
              * HW.RX(Seq.Read(t).iDevice).fSample);
      max_CIC_Decimation = min(8000, HW.RX(Seq.Read(1).iDevice).CIC_Decimation_Max);
      if CIC_Decimation_without_SamplingFactor > max_CIC_Decimation
        CIC_Decimation_Rates = ...
          HW.RX(Seq.Read(1).iDevice).CIC_Decimation_Min:max_CIC_Decimation;
        devFromTarget = abs(mod(CIC_Decimation_without_SamplingFactor ./ CIC_Decimation_Rates + 0.5, 1) - 0.5);
        % select smallest sampling factor that "comes closest" to the target
        % sampling rate
        Seq.Read(t).SamplingFactor = CIC_Decimation_without_SamplingFactor ...
          / CIC_Decimation_Rates(find(devFromTarget==min(devFromTarget), 1, 'last'));
      else
        Seq.Read(t).SamplingFactor = 1;
      end
    else
      Seq.Read(t).SamplingFactor = Seq.AQSlice(Seq.Read(t).useAQSlice).SamplingFactor;
    end
  end

  % shift readout if number of samples is even
  if all(mod(Seq.Read(t).nRead*Seq.Read(t).ReadOS,2)==1) || (Seq.Read(t).StartWithKSpaceCenter ~= 0)
    % Center of the readout is at the center of the acquisition window.
    % FIXME: Do we need to make sure that the center of some sample is at the
    %        k-space center if Seq.Read(t).StartWithKSpaceCenter ~= 0?
    Seq.Read(t).fftReadoutShift = 0;
  else
    % Center of the readout is the Seq.Read(t).nRead*Seq.Read(t).ReadOS/2+1
    % sample of the acquisition window.
    % If the readout direction is inverted (GradSign is -1), it is the same
    % sample but counted from the back of the acquisition window.
    Seq.Read(t).fftReadoutShift = Seq.Read(t).GradSign;
  end

  if isemptyfield(Seq.Read(t), 'GradLength')
    % extend gradient length by one pixel to account for a shifted AQ window
    Seq.Read(t).GradLength = (1 + any(Seq.Read(t).fftReadoutShift~=0)./(Seq.Read(t).nRead*Seq.Read(t).ReadOS))/Seq.Read(t).HzPerPixelMin + ...
      Seq.Read(t).tEC*2 + Seq.Read(t).tRamp*2;
  end
  if isemptyfield(Seq.Read(t), 'HzPerPixelMin')
    Seq.Read(t).HzPerPixelMin = (1 + any(Seq.Read(t).fftReadoutShift~=0)./(Seq.Read(t).nRead*Seq.Read(t).ReadOS)) / ...
      (Seq.Read(t).GradLength - Seq.Read(t).tEC*2 - Seq.Read(t).tRamp*2);
  end

  if isemptyfield(Seq.Read(t), 'GradDephaseLengthFactor')
    Seq.Read(t).GradDephaseLengthFactor = 1;
  end
  if isemptyfield(Seq.Read(t), 'GradRephaseLengthFactor')
    Seq.Read(t).GradRephaseLengthFactor = 1;
  end

  if isemptyfield(Seq.Read(t), 'kLineDir')
    Seq.Read(t).kLineDir = zeros(3, 1);
    Seq.Read(t).kLineDir(Seq.Read(t).UseCoordinate) = 1;
  end

  if isemptyfield(Seq.Read(t), 'alfa')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Read(t).useAQSlice), 'alfa')
      Seq.Read(t).alfa = 0;
    else
      Seq.Read(t).alfa = Seq.AQSlice(Seq.Read(t).useAQSlice).alfa;
    end
  end
  if isemptyfield(Seq.Read(t), 'phi')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Read(t).useAQSlice), 'phi')
      Seq.Read(t).phi = 0;
    else
      Seq.Read(t).phi = Seq.AQSlice(Seq.Read(t).useAQSlice).phi;
    end
  end
  if isemptyfield(Seq.Read(t), 'theta')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Read(t).useAQSlice), 'theta')
      Seq.Read(t).theta = 0;
    else
      Seq.Read(t).theta = Seq.AQSlice(Seq.Read(t).useAQSlice).theta;
    end
  end
  if isemptyfield(Seq.Read(t), 'angle2Turns')  % conversion factor of the input angles to full turns e.g.: 1/(2*pi)
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Read(t).useAQSlice), 'angle2Turns')
      Seq.Read(t).angle2Turns = 1/(2*pi);
    else
      Seq.Read(t).angle2Turns = Seq.AQSlice(Seq.Read(t).useAQSlice).angle2Turns;
    end
  end
  if isemptyfield(Seq.Read(t), 'distance'), Seq.Read(t).distance = 0; end
  if isemptyfield(Seq.Read(t), 'Resolution'), Seq.Read(t).Resolution = Seq.Read(t).sizeRead./Seq.Read(t).nRead; end


  % Seq.Read(t).CIC_Decimation = ...
  %   floor((1/Seq.Read(t).HzPerPixelMin) / (Seq.Read(t).nRead*Seq.Read(t).ReadOS + Seq.Read(t).fftReadoutShift) / ...
  %         (1/HW.RX(Seq.Read(t).iDevice).fSample));
  Seq.Read(t).CIC_Decimation = ...
    floor((1/Seq.Read(t).HzPerPixelMin) / (Seq.Read(t).nRead*Seq.Read(t).ReadOS*Seq.Read(t).SamplingFactor) / ...
          (1/HW.RX(Seq.Read(t).iDevice).fSample));
  if Seq.Read(t).CIC_Decimation > HW.RX(Seq.Read(t).iDevice).CIC_Decimation_Max
    error('PD:get_ReadParameter:SamplingRateTooLow', ...
      'Acquisition sampling rate is too low.');
  end
  if Seq.Read(t).CIC_Decimation < HW.RX(Seq.Read(t).iDevice).CIC_Decimation_Min
    error('PD:get_ReadParameter:SamplingRateTooHigh', ...
      'Acquisition sampling rate is too high.');
  end

  Seq.Read(t).fSample = HW.RX(Seq.Read(t).iDevice).fSample / ...
    (Seq.Read(t).CIC_Decimation*Seq.Read(t).SamplingFactor);  % Dezimierung
  Seq.Read(t).AcquisitionTime = Seq.Read(t).nRead*Seq.Read(t).ReadOS / Seq.Read(t).fSample;
  Seq.Read(t).HzPerPixel = 1/Seq.Read(t).AcquisitionTime;

  Seq.Read(t).GradTimeIntegralAQ = 2*pi/Seq.Read(t).Gamma/Seq.Read(t).Resolution .* Seq.Read(t).GradSign;
  Seq.Read(t).GradAmp = Seq.Read(t).GradTimeIntegralAQ./Seq.Read(t).AcquisitionTime;
  Seq.Read(t).GradTimeIntegral = (Seq.Read(t).GradAmp).*(Seq.Read(t).GradLength-Seq.Read(t).tRamp);
  % FIXME: Should "fftReadoutShift" be considered in the limit of
  %        StartWithKSpaceCenter -> 0?
  if Seq.Read(t).StartWithKSpaceCenter == 0
    Seq.Read(t).GradTimeIntegralDephase = Seq.Read(t).GradTimeIntegral/2;
    Seq.Read(t).StartOfReadoutGrad = Seq.Read(t).CenterOfReadout(:).' - Seq.Read(t).GradLength/2;
  else
    Seq.Read(t).GradTimeIntegralDephase = Seq.Read(t).GradAmp * ...
      (Seq.Read(t).AcquisitionTime/2 * ...
      (1 + Seq.Read(t).StartWithKSpaceCenter*(1/(Seq.Read(t).nRead*Seq.Read(t).ReadOS)-1)) + ...
      Seq.Read(t).tEC + Seq.Read(t).tRamp/2);
    Seq.Read(t).StartOfReadoutGrad = Seq.Read(t).CenterOfReadout(:).' - ...
      (Seq.Read(t).AcquisitionTime/2 * ...
      (1 + Seq.Read(t).StartWithKSpaceCenter*(1/(Seq.Read(t).nRead*Seq.Read(t).ReadOS)-1) *(~Seq.Read(t).StartWithKSpaceCenterNoSampleShift)) + ...
      Seq.Read(t).tEC + Seq.Read(t).tRamp);
  end
  Seq.Read(t).GradTimeIntegralRephase = Seq.Read(t).GradTimeIntegral - Seq.Read(t).GradTimeIntegralDephase;

  if Seq.Read(t).GradTimeIntegral == 0
    tD = 0.5;
    tR = 0.5;
  else
    tD = Seq.Read(t).GradTimeIntegralDephase / Seq.Read(t).GradTimeIntegral;
    tR = Seq.Read(t).GradTimeIntegralRephase / Seq.Read(t).GradTimeIntegral;
  end
  tD = tD * Seq.Read(t).GradDephaseLengthFactor;
  tR = tR * Seq.Read(t).GradRephaseLengthFactor;

  if isemptyfield(Seq.Read(t), 'GradTimeIntegralOffset'), Seq.Read(t).GradTimeIntegralOffset = 0; end
  if any(Seq.Read(t).GradTimeIntegralOffset)
    Seq.Read(t).GradTimeDelayOffset = ...
      repmat(Seq.Read(t).GradTimeIntegralOffset(:).', 1, ...
             numel(Seq.Read(t).UseAtRepetitionTime)/numel(Seq.Read(t).GradTimeIntegralOffset)) ...
      ./ abs(Seq.Read(t).GradAmp);
  end
  if isemptyfield(Seq.Read(t), 'GradTimeDelayOffset'), Seq.Read(t).GradTimeDelayOffset = 0; end
  if isemptyfield(Seq.Read(t), 'GradDephaseLength')
    if tD == 0
      % use minimum dephase length
      Seq.Read(t).GradDephaseLength = 2*Seq.Read(t).tRamp + 2/HW.MMRT(Seq.Read(t).iDevice).fSystem;
    else
      Seq.Read(t).GradDephaseLength = (Seq.Read(t).GradLength-Seq.Read(t).tRamp)*tD + Seq.Read(t).tRamp;
    end
  end
  if isemptyfield(Seq.Read(t), 'GradRephaseLength')
    if tR == 0
      % use minimum rephase length
      Seq.Read(t).GradRephaseLength = 2*Seq.Read(t).tRamp + 2/HW.MMRT(Seq.Read(t).iDevice).fSystem;
    else
      Seq.Read(t).GradRephaseLength = (Seq.Read(t).GradLength-Seq.Read(t).tRamp)*tR + Seq.Read(t).tRamp;
    end
  end
  if isemptyfield(Seq.Read(t), 'GradRephaseSign'), Seq.Read(t).GradRephaseSign = -1; end
  if isemptyfield(Seq.Read(t), 'GradDephaseSign'), Seq.Read(t).GradDephaseSign = -1; end
  if isemptyfield(Seq.Read(t), 'CenterOfDephase')
    Seq.Read(t).CenterOfDephase = ...
      Seq.Read(t).StartOfReadoutGrad - Seq.Read(t).GradDephaseLength/2 + Seq.Read(t).tRamp;
  end
  if isemptyfield(Seq.Read(t), 'CenterOfRephase')
    Seq.Read(t).CenterOfRephase = ...
      Seq.Read(t).StartOfReadoutGrad + Seq.Read(t).GradLength + Seq.Read(t).GradRephaseLength/2 - Seq.Read(t).tRamp;
  end
  if isemptyfield(Seq.Read(t), 'PhaseIncrement'), Seq.Read(t).PhaseIncrement = 0; end
  if isemptyfield(Seq.Read(t), 'Phase'), Seq.Read(t).Phase = 0; end
  if ~isfield(Seq.Read, 'ResetPhases'), Seq.Read(t).ResetPhases = []; end
  if ~isfield(Seq.Read, 'UseAtRepetitionTime'), Seq.Read(t).UseAtRepetitionTime = []; end %!!
  if isemptyfield(Seq.Read(t), 'UseAtRepetitionTimeRephase')
    Seq.Read(t).UseAtRepetitionTimeRephase = Seq.Read(t).UseAtRepetitionTime;
  end
  if isemptyfield(Seq.Read(t), 'UseAtRepetitionTimeDephase')
    Seq.Read(t).UseAtRepetitionTimeDephase = Seq.Read(t).UseAtRepetitionTime;
  end
  if isemptyfield(Seq.Read(t), 'GetData'), Seq.Read(t).GetData = 0; end
  if isemptyfield(Seq.Read(t), 'Overdrive'), Seq.Read(t).Overdrive = 0; end
  if isemptyfield(Seq.Read(t), 'StartWithGradientTime'), Seq.Read(t).StartWithGradientTime = 0; end

  if round((Seq.Read(t).GradDephaseLength-Seq.Read(t).tRamp*2)*HW.MMRT(Seq.Read(t).iDevice).fSystem) < 2
    error('PD:get_ReadParameter:GradDephaseLengthTooShort', ...
      'Seq.Read(%d).GradDephaseLength too short to fit 2*tRamp', t);
  end
  if round((Seq.Read(t).GradRephaseLength-Seq.Read(t).tRamp*2)*HW.MMRT(Seq.Read(t).iDevice).fSystem) < 2
    error('PD:get_ReadParameter:GradRephaseLengthTooShort', ...
      'Seq.Read(%d).GradRephaseLength too short to fit 2*tRamp', t);
  end

  Seq.Read(t).GradAmpDephase = Seq.Read(t).GradTimeIntegralDephase ./ ...
    (Seq.Read(t).GradDephaseLength - Seq.Read(t).tRamp) * Seq.Read(t).GradDephaseSign;
  Seq.Read(t).GradAmpRephase = Seq.Read(t).GradTimeIntegralRephase ./ ...
    (Seq.Read(t).GradRephaseLength - Seq.Read(t).tRamp) * Seq.Read(t).GradRephaseSign;
  Seq.Read(t).AcquisitionFrequencyOffset = Seq.Read(t).Gamma/2/pi * Seq.Read(t).distance .* Seq.Read(t).GradAmp;
  Seq.Read(t).AcquisitionFrequency = Seq.Read(t).Gamma/2/pi*HW.B0 + Seq.Read(t).AcquisitionFrequencyOffset;
  FrequencyGridRX = HW.RX(Seq.Read(t).iDevice).fSample / (2^HW.RX(Seq.Read(t).iDevice).DdsPicBits);
  Seq.Read(t).AcquisitionFrequency = round(Seq.Read(t).AcquisitionFrequency/FrequencyGridRX) * FrequencyGridRX;
  if ~isempty(Seq.Read(t).GammaX)
    Seq.Read(t).AcquisitionFrequencyOffsetX = Seq.Read(t).GammaX/2/pi * Seq.Read(t).distance .* Seq.Read(t).GradAmp;
    Seq.Read(t).AcquisitionFrequencyX = HW.fLarmorX + Seq.Read(t).AcquisitionFrequencyOffsetX;
    Seq.Read(t).AcquisitionFrequencyX = round(Seq.Read(t).AcquisitionFrequencyX/FrequencyGridRX) * FrequencyGridRX;
  end
  % maxsize = max([numel(Seq.Read(t).UseAtRepetitionTime); numel(Seq.Read(t).Phase); numel(Seq.Read(t).PhaseIncrement)]);

  % if Seq.Read(t).PhaseIncrement==0 && maxsize==1
  %   tempsize = 1;
  %   tempPhaseInc = 0;
  % else
  tempsize = ones(1, numel(Seq.Read(t).UseAtRepetitionTime));
  if numel(Seq.Read(t).PhaseIncrement) > 1
    % tempPhaseInc = cumsum(Seq.Read(t).PhaseIncrement);
    temp = repmat(Seq.Read(t).PhaseIncrement(:), 1, numel(Seq.Read(t).UseAtRepetitionTime)/numel(Seq.Read(t).PhaseIncrement));
    tempPhaseInc = cumsum(temp(:)).';
  else
    tempPhaseInc = cumsum(Seq.Read(t).PhaseIncrement*tempsize);
  end
  if numel(Seq.Read(t).Phase) > 1
    temp = repmat(Seq.Read(t).Phase(:), 1, numel(Seq.Read(t).UseAtRepetitionTime)/numel(Seq.Read(t).Phase));
    tempPhase = (temp(:)).';
  else
    tempPhase = (Seq.Read(t).Phase*tempsize);
  end
  % end

  tempsize(~Seq.Read(t).UseAtRepetitionTime) = NaN;

  if Seq.Read(t).StartWithKSpaceCenterNonlinear
    Seq.Read(t).nRead = Seq.Read(t).nRead+ceil(Seq.Read(t).tRamp/(Seq.Read(t).ReadOS/Seq.Read(t).fSample));
  end
  % gradAmpTemp = [0, 0, Seq.Read(t).GradAmp, Seq.Read(t).GradAmp];
  % gradTimeTemp = [-1e3, Seq.Read(t).StartOfReadoutGrad, Seq.Read(t).StartOfReadoutGrad+Seq.Read(t).tRamp, Seq.Read(t).GradAmp];

  if Seq.Read(t).StartWithKSpaceCenterNonlinear
    Seq.Read(t).GradTimeIntegralOfSample = ...
        - repmat(Seq.Read(t).GradTimeIntegralDephase, Seq.Read(t).nRead*Seq.Read(t).ReadOS, 1) * (Seq.Read(t).StartWithKSpaceCenter==0) ...
        + linspace(0, 1, Seq.Read(t).nRead*Seq.Read(t).ReadOS).' * Seq.Read(t).GradTimeIntegralAQ * (Seq.Read(t).StartWithKSpaceCenterNonlinear==0) ...
        + [( linspace(0, ...
                      (floor(Seq.Read(t).tRamp*Seq.Read(t).fSample)/Seq.Read(t).fSample)/Seq.Read(t).tRamp, ...
                       floor(Seq.Read(t).tRamp*Seq.Read(t).fSample)).^2).' ...
                      * (Seq.Read(t).GradAmp*Seq.Read(t).tRamp/2); ...
             repmat((Seq.Read(t).GradAmp*Seq.Read(t).tRamp/2), Seq.Read(t).nRead*Seq.Read(t).ReadOS-floor(Seq.Read(t).tRamp*Seq.Read(t).fSample), 1) + ...
             linspace(1/Seq.Read(t).fSample-(Seq.Read(t).tRamp-(floor(Seq.Read(t).tRamp*Seq.Read(t).fSample)/Seq.Read(t).fSample)), ...
                      (Seq.Read(t).nRead*Seq.Read(t).ReadOS-floor(Seq.Read(t).tRamp*Seq.Read(t).fSample))/Seq.Read(t).fSample, ...
                      Seq.Read(t).nRead*Seq.Read(t).ReadOS-floor(Seq.Read(t).tRamp*Seq.Read(t).fSample)).' ...
             * Seq.Read(t).GradAmp ...
           ] * (Seq.Read(t).StartWithKSpaceCenterNonlinear==1);
  else
    Seq.Read(t).GradTimeIntegralOfSample = ...
        -repmat(Seq.Read(t).GradTimeIntegralDephase, Seq.Read(t).nRead*Seq.Read(t).ReadOS, 1) * (Seq.Read(t).StartWithKSpaceCenter==0) ...
        + linspace(0, 1, Seq.Read(t).nRead*Seq.Read(t).ReadOS).' * Seq.Read(t).GradTimeIntegralAQ * (Seq.Read(t).StartWithKSpaceCenterNonlinear==0);
  end


  Seq.Read(t).AQ.Start = NaN(1, size(Seq.tRep,2));
  Seq.Read(t).AQ.fSample = Seq.Read(t).AQ.Start;
  Seq.Read(t).AQ.nSamples = Seq.Read(t).AQ.Start;
  Seq.Read(t).AQ.SamplingFactor = Seq.Read(t).AQ.Start;
  Seq.Read(t).AQ.Frequency = Seq.Read(t).AQ.Start;
  Seq.Read(t).AQ.Phase = Seq.Read(t).AQ.Start;
  Seq.Read(t).AQ.GetData = Seq.Read(t).AQ.Start;
  % Seq.Read(t).AQ.ResetPhases = Seq.Read(t).AQ.Start;


  Seq.Read(t).AQ.Start(:,Seq.Read(t).UseAtRepetitionTime) = ...
    (Seq.Read(t).CenterOfReadout(:).' ...
     - Seq.Read(t).AcquisitionTime/2 * (1 + Seq.Read(t).StartWithKSpaceCenter * (1/(Seq.Read(t).nRead*Seq.Read(t).ReadOS)-1) *(~Seq.Read(t).StartWithKSpaceCenterNoSampleShift)) ...
     - (Seq.Read(t).StartWithKSpaceCenterNonlinear~=0) * (Seq.Read(t).tEC+Seq.Read(t).tRamp) ...
     - (Seq.Read(t).StartWithGradientBeforeExcitationPulse~=0) * (Seq.Read(t).tEC+Seq.Read(t).tRamp) ...
     - Seq.Read(t).fftReadoutShift./Seq.Read(t).fSample/2 ...
    ) .* tempsize;
  Seq.Read(t).AQ.fSample(:,Seq.Read(t).UseAtRepetitionTime) = Seq.Read(t).fSample;
  Seq.Read(t).AQ.nSamples(:,Seq.Read(t).UseAtRepetitionTime) = Seq.Read(t).nRead*Seq.Read(t).ReadOS * tempsize;
  Seq.Read(t).AQ.SamplingFactor(:,Seq.Read(t).UseAtRepetitionTime) = Seq.Read(t).SamplingFactor;
  if numel(Seq.Read(t).AcquisitionFrequency) == 1
    Seq.Read(t).AQ.Frequency(:,Seq.Read(t).UseAtRepetitionTime) = Seq.Read(t).AcquisitionFrequency * tempsize;
  else
    Seq.Read(t).AQ.Frequency(:,Seq.Read(t).UseAtRepetitionTime) = Seq.Read(t).AcquisitionFrequency;
  end
  Seq.Read(t).AQ.Phase(:,Seq.Read(t).UseAtRepetitionTime) = (tempPhaseInc+tempPhase) .* tempsize;
  % Seq.Read(t).AQ.ResetPhases(:,Seq.Read(t).UseAtRepetitionTime) = Seq.Read(t).ResetPhases;
  Seq.Read(t).AQ.GetData(1,Seq.Read(t).UseAtRepetitionTime) = Seq.Read(t).GetData;
  if  Seq.Read(t).StartWithGradientBeforeExcitationPulse
    Seq.Read(t).GradTimeIntegralOfSample = Seq.Read(t).GradTimeIntegralOfSample ...
      + repmat(Seq.Read(t).GradAmp.*Seq.Read(t).AQ.Start(:,Seq.Read(t).UseAtRepetitionTime), Seq.Read(t).nRead*Seq.Read(t).ReadOS, 1);
  end
  if ~isempty(Seq.Read(t).GammaX)
    Seq.Read(t).AQ.FrequencyX = NaN(size(Seq.Read(t).AQ.Start));
    if numel(Seq.Read(t).AcquisitionFrequencyX) == 1
      Seq.Read(t).AQ.FrequencyX(:,Seq.Read(t).UseAtRepetitionTime) = Seq.Read(t).AcquisitionFrequencyX * tempsize;
    else
      Seq.Read(t).AQ.FrequencyX(:,Seq.Read(t).UseAtRepetitionTime) = Seq.Read(t).AcquisitionFrequencyX;
    end

    Seq.Read(t).AQ.PhaseX = Seq.Read(t).AQ.Phase;  % FIXME: Do we need different phase increments at the X frequency?
  end
  Angle2Deg = Seq.Read(t).angle2Turns*360;
  [Rx, Ry, Rz] = get_aptDegRotationMatrix(Seq.Read(t).alfa*Angle2Deg, Seq.Read(t).phi*Angle2Deg, Seq.Read(t).theta*Angle2Deg);

  Seq.Read(t).GradAmpNorm = ...
    bsxfun(@times, reshape(Seq.Read(t).kLineDir, 3, []), ...
                   reshape(Seq.Read(t).GradAmp, 1, []));
  Seq.Read(t).GradAmpDephaseNorm = ...
    bsxfun(@times, reshape(Seq.Read(t).kLineDir, 3, []), ...
                   reshape(Seq.Read(t).GradAmpDephase, 1, []));
  Seq.Read(t).GradAmpRephaseNorm = ...
    bsxfun(@times, reshape(Seq.Read(t).kLineDir, 3, []), ...
                   reshape(Seq.Read(t).GradAmpRephase, 1, []));

  Seq.Read(t).GradAmp = zeros(3, numel(tempsize));
  Seq.Read(t).GradAmpRephase = zeros(3, numel(tempsize));
  Seq.Read(t).GradAmpDephase = zeros(3, numel(tempsize));
  if all(1 == [numel(Seq.Read(t).alfa), numel(Seq.Read(t).phi), numel(Seq.Read(t).theta)])
    Seq.Read(t).GradAmp =        repmat((Rz*Ry*Rx)*Seq.Read(t).GradAmpNorm,        1, numel(tempsize)/size(Seq.Read(t).GradAmpNorm, 2));
    Seq.Read(t).GradAmpRephase = repmat((Rz*Ry*Rx)*Seq.Read(t).GradAmpRephaseNorm, 1, numel(tempsize)/size(Seq.Read(t).GradAmpNorm, 2));
    Seq.Read(t).GradAmpDephase = repmat((Rz*Ry*Rx)*Seq.Read(t).GradAmpDephaseNorm, 1, numel(tempsize)/size(Seq.Read(t).GradAmpNorm, 2));
  else
    for tt = 1:numel(tempsize)
      Seq.Read(t).GradAmp(:,tt) = Rz(:,:,tt) * (Ry(:,:,tt) * (Rx(:,:,tt) * Seq.Read(t).GradAmpNorm(:,tt)));
      Seq.Read(t).GradAmpRephase(:,tt) = Rz(:,:,tt) * (Ry(:,:,tt) * (Rx(:,:,tt) * Seq.Read(t).GradAmpRephaseNorm(:,tt)));
      Seq.Read(t).GradAmpDephase(:,tt) = Rz(:,:,tt) * (Ry(:,:,tt) * (Rx(:,:,tt) * Seq.Read(t).GradAmpDephaseNorm(:,tt)));
    end
  end

  for n = 1:3
    % if Seq.Read(t).Overdrive == 0
    %
    %   Seq.Read(t).Grad(n).Amp = zeros(4, size(tempsize,2));
    %   Seq.Read(t).GradRephase(n).Amp = zeros(4, size(tempsize,2));
    %   Seq.Read(t).GradDephase(n).Amp = zeros(4, size(tempsize,2));
    %   Seq.Read(t).Grad(n).Amp(2:3,:) = Seq.Read(t).GradAmp(n);
    %   Seq.Read(t).GradRephase(n).Amp(2:3,:) = Seq.Read(t).GradAmpRephase(n);
    %   Seq.Read(t).GradDephase(n).Amp(2:3,:) = Seq.Read(t).GradAmpDephase(n);
    %   Seq.Read(t).Grad(n).Amp = Seq.Read(t).Grad(n).Amp(1:4,1) * tempsize;
    %   Seq.Read(t).GradRephase(n).Amp = Seq.Read(t).GradRephase(n).Amp(1:4,1) * tempsize;
    %   Seq.Read(t).GradDephase(n).Amp = Seq.Read(t).GradDephase(n).Amp(1:4,1) * tempsize;
    Seq.Read(t).Grad(n).Amp = NaN(4, size(Seq.tRep,2));
    Seq.Read(t).GradRephase(n).Amp = Seq.Read(t).Grad(n).Amp;
    Seq.Read(t).GradDephase(n).Amp = Seq.Read(t).Grad(n).Amp;

    Seq.Read(t).Grad(n).Amp(1,Seq.Read(t).UseAtRepetitionTime) = 0;
    Seq.Read(t).Grad(n).Amp(2:3,Seq.Read(t).UseAtRepetitionTime) = ...
      [Seq.Read(t).GradAmp(n,:); Seq.Read(t).GradAmp(n,:)];
    Seq.Read(t).Grad(n).Amp(4,Seq.Read(t).UseAtRepetitionTime) = 0;

    Seq.Read(t).GradRephase(n).Amp(1,Seq.Read(t).UseAtRepetitionTimeRephase) = 0;
    Seq.Read(t).GradRephase(n).Amp(2:3,Seq.Read(t).UseAtRepetitionTimeRephase) = ...
      [Seq.Read(t).GradAmpRephase(n,:); Seq.Read(t).GradAmpRephase(n,:)];
    Seq.Read(t).GradRephase(n).Amp(4,Seq.Read(t).UseAtRepetitionTimeRephase) = 0;

    Seq.Read(t).GradDephase(n).Amp(1,Seq.Read(t).UseAtRepetitionTimeDephase) = 0;
    Seq.Read(t).GradDephase(n).Amp(2:3,Seq.Read(t).UseAtRepetitionTimeDephase) = ...
      [Seq.Read(t).GradAmpDephase(n,:); Seq.Read(t).GradAmpDephase(n,:)];
    Seq.Read(t).GradDephase(n).Amp(4,Seq.Read(t).UseAtRepetitionTimeDephase) = 0;

    Seq.Read(t).Grad(n).Time = NaN(4, size(Seq.tRep,2));
    Seq.Read(t).GradRephase(n).Time = Seq.Read(t).Grad(n).Time;
    Seq.Read(t).GradDephase(n).Time = Seq.Read(t).Grad(n).Time;


    if Seq.Read(t).StartWithGradientBeforeExcitationPulse
      Seq.Read(t).Grad(n).Time(:,Seq.Read(t).UseAtRepetitionTime) = ...
        [ (0                                        + Seq.Read(t).StartWithGradientTime); ...
          (Seq.Read(t).tRamp                        + Seq.Read(t).StartWithGradientTime); ...
          (Seq.Read(t).GradLength-Seq.Read(t).tRamp + Seq.Read(t).StartOfReadoutGrad); ...
          (Seq.Read(t).GradLength                   + Seq.Read(t).StartOfReadoutGrad)] ...
        * tempsize ...
        - Seq.Read(t).GradTimeDelay(n) ...
        - repmat(Seq.Read(t).GradTimeDelayOffset(:).', 4, ...
                 numel(Seq.Read(t).UseAtRepetitionTime)/numel(Seq.Read(t).GradTimeDelayOffset));
    else
      Seq.Read(t).Grad(n).Time(:,Seq.Read(t).UseAtRepetitionTime) = ...
        Seq.Read(t).StartOfReadoutGrad ...
        + [ 0; ...
            Seq.Read(t).tRamp; ...
            (Seq.Read(t).GradLength - Seq.Read(t).tRamp); ...
            Seq.Read(t).GradLength] * tempsize ...
        - Seq.Read(t).GradTimeDelay(n) ...
        - repmat(Seq.Read(t).GradTimeDelayOffset(:).', 4, ...
                 numel(Seq.Read(t).UseAtRepetitionTime)/numel(Seq.Read(t).GradTimeDelayOffset));
    end
    Seq.Read(t).GradRephase(n).Time(:,Seq.Read(t).UseAtRepetitionTimeRephase) = ...
      bsxfun(@plus, ...
             repmat(Seq.Read(t).CenterOfRephase(:).', 1, ...
                    numel(Seq.Read(t).UseAtRepetitionTime)/numel(Seq.Read(t).CenterOfRephase)), ...
             [-Seq.Read(t).GradRephaseLength/2; ...
              (-Seq.Read(t).GradRephaseLength/2 + Seq.Read(t).tRamp); ...
              (+Seq.Read(t).GradRephaseLength/2 - Seq.Read(t).tRamp); ...
              +Seq.Read(t).GradRephaseLength/2]) ...
      - Seq.Read(t).GradTimeDelay(n) ...
      - repmat(Seq.Read(t).GradTimeDelayOffset(:).', 4, ...
               numel(Seq.Read(t).UseAtRepetitionTime)/numel(Seq.Read(t).GradTimeDelayOffset));
    Seq.Read(t).GradDephase(n).Time(:,Seq.Read(t).UseAtRepetitionTimeDephase) = ...
      bsxfun(@plus, ...
             repmat(Seq.Read(t).CenterOfDephase(:).', 1, ...
                    numel(Seq.Read(t).UseAtRepetitionTime)/numel(Seq.Read(t).CenterOfDephase)), ...
             [-Seq.Read(t).GradDephaseLength/2; ...
              (-Seq.Read(t).GradDephaseLength/2 + Seq.Read(t).tRamp); ...
              (+Seq.Read(t).GradDephaseLength/2 - Seq.Read(t).tRamp); ...
              +Seq.Read(t).GradDephaseLength/2]) ...
      - Seq.Read(t).GradTimeDelay(n) ...
      - repmat(Seq.Read(t).GradTimeDelayOffset(:).', 4, ...
               numel(Seq.Read(t).UseAtRepetitionTime)/numel(Seq.Read(t).GradTimeDelayOffset));

    % remove pulses without tRep
    Seq.Read(t).Grad(n).Time(:,numel(Seq.tRep)+1:end) = [];
    Seq.Read(t).Grad(n).Amp(:,numel(Seq.tRep)+1:end) = [];
    Seq.Read(t).GradDephase(n).Time(:,numel(Seq.tRep)+1:end) = [];
    Seq.Read(t).GradDephase(n).Amp(:,numel(Seq.tRep)+1:end) = [];
    Seq.Read(t).GradRephase(n).Time(:,numel(Seq.tRep)+1:end) = [];
    Seq.Read(t).GradRephase(n).Amp(:,numel(Seq.tRep)+1:end) = [];

    % else
    %   Seq.Read(t).Grad(n).Amp = zeros(8, size(tempsize,2));
    %   Seq.Read(t).GradRephase(n).Amp = zeros(8, size(tempsize,2));
    %   Seq.Read(t).GradDephase(n).Amp = zeros(8, size(tempsize,2));
    %   Seq.Read(t).Grad(n).Amp(2:3,:) = Seq.Read(t).GradAmp(n) * (1+Seq.Read(t).Overdrive);
    %   Seq.Read(t).GradRephase(n).Amp(2:3,:) = Seq.Read(t).GradAmpRephase(n) * (1+Seq.Read(t).Overdrive);
    %   Seq.Read(t).GradDephase(n).Amp(2:3,:) = Seq.Read(t).GradAmpDephase(n) * (1+Seq.Read(t).Overdrive);
    %   Seq.Read(t).Grad(n).Amp(4:5,:) = Seq.Read(t).GradAmp(n);
    %   Seq.Read(t).GradRephase(n).Amp(4:5,:) = Seq.Read(t).GradAmpRephase(n);
    %   Seq.Read(t).GradDephase(n).Amp(4:5,:) = Seq.Read(t).GradAmpDephase(n);
    %   Seq.Read(t).Grad(n).Amp(6:7,:) = Seq.Read(t).GradAmp(n) * (0-Seq.Read(t).Overdrive);
    %   Seq.Read(t).GradRephase(n).Amp(6:7,:) = Seq.Read(t).GradAmpRephase(n) * (0-Seq.Read(t).Overdrive);
    %   Seq.Read(t).GradDephase(n).Amp(6:7,:) = Seq.Read(t).GradAmpDephase(n) * (0-Seq.Read(t).Overdrive);
    %   Seq.Read(t).Grad(n).Amp = Seq.Read(t).Grad(n).Amp(1:8,1) * tempsize;
    %   Seq.Read(t).GradRephase(n).Amp = Seq.Read(t).GradRephase(n).Amp(1:8,1) * tempsize;
    %   Seq.Read(t).GradDephase(n).Amp = Seq.Read(t).GradDephase(n).Amp(1:8,1) * tempsize;
    %
    %   Seq.Read(t).Grad(n).Time = Seq.Read(t).CenterOfReadout ...
    %     + [ -Seq.Read(t).GradLength/2; ...
    %         (-Seq.Read(t).GradLength/2 + Seq.Read(t).tRamp*1/3); ...
    %         (-Seq.Read(t).GradLength/2 + Seq.Read(t).tRamp*2/3); ...
    %         (-Seq.Read(t).GradLength/2 + Seq.Read(t).tRamp*3/3); ...
    %         (+Seq.Read(t).GradLength/2 - Seq.Read(t).tRamp*3/3); ...
    %         (+Seq.Read(t).GradLength/2 - Seq.Read(t).tRamp*2/3); ...
    %         (+Seq.Read(t).GradLength/2 - Seq.Read(t).tRamp*1/3); ...
    %         +Seq.Read(t).GradLength/2] * tempsize ...
    %     - Seq.Read(t).GradTimeDelay(n);
    %   Seq.Read(t).GradRephase(n).Time = Seq.Read(t).CenterOfRephase ...
    %     + [ -Seq.Read(t).GradRephaseLength/2; ...
    %         (-Seq.Read(t).GradRephaseLength/2 + Seq.Read(t).tRamp*1/3); ...
    %         (-Seq.Read(t).GradRephaseLength/2 + Seq.Read(t).tRamp*2/3); ...
    %         (-Seq.Read(t).GradRephaseLength/2 + Seq.Read(t).tRamp*3/3); ...
    %         (+Seq.Read(t).GradRephaseLength/2 - Seq.Read(t).tRamp*3/3); ...
    %         (+Seq.Read(t).GradRephaseLength/2 - Seq.Read(t).tRamp*2/3); ...
    %         (+Seq.Read(t).GradRephaseLength/2 - Seq.Read(t).tRamp*1/3); ...
    %         +Seq.Read(t).GradRephaseLength/2] * tempsize ...
    %     - Seq.Read(t).GradTimeDelay(n);
    %   Seq.Read(t).GradDephase(n).Time = Seq.Read(t).CenterOfDephase ...
    %     + [ -Seq.Read(t).GradDephaseLength/2; ...
    %         (-Seq.Read(t).GradDephaseLength/2 + Seq.Read(t).tRamp*1/3); ...
    %         (-Seq.Read(t).GradDephaseLength/2 + Seq.Read(t).tRamp*2/3); ...
    %         (-Seq.Read(t).GradDephaseLength/2 + Seq.Read(t).tRamp*3/3); ...
    %         (+Seq.Read(t).GradDephaseLength/2 - Seq.Read(t).tRamp*3/3); ...
    %         (+Seq.Read(t).GradDephaseLength/2 - Seq.Read(t).tRamp*2/3); ...
    %         (+Seq.Read(t).GradDephaseLength/2 - Seq.Read(t).tRamp*1/3); ...
    %         (+Seq.Read(t).GradDephaseLength/2] * tempsize ...
    %     - Seq.Read(t).GradTimeDelay(n);
    % end
  end
end

end
