function Seq = get_PhaseParameter(Seq, HW, idxModule)
%% Calculate phase encoder
%
%    Seq = get_PhaseParameter(Seq, HW, idxModule);
%
% This function calculates rephase and dephase gradients derived from the values
% in Seq.Phase.
%
%
% INPUT:
%
% "Seq" is a structure with the field "Phase" that is an (array of) structures
% with the following fields:
%
%   UseCoordinate
%       Number of the coordinate for which the phase encoder is used.
%       (Default: 3)
%
%   GradRephaseLength / GradDephaseLength
%       Length of the rephase/dephase pulse in seconds.
%       (Default: same as the corresponding read gradient
%       Seq.Read(t).GradRephaseLength/GradRephaseLength)
%
%   CenterOfRephase / CenterOfDephase
%       Center of rephase/dephase pulse in seconds on the tRep timeline.
%       (Default: same as the corresponding read gradient
%       Seq.Read(t).CenterOfRephase/CenterOfDephase)
%
%   GradRephaseSign / GradDephaseSign
%       Sign of the rephase/dephase pulse.
%       (Default: GradRephaseSign = -1, GradDephaseSign = 1)
%
%   GradTimeIntegralRephaseOffset / GradTimeIntegralDephaseOffset
%       Offset for the time integral of the rephase/dephase pulse in s*T/m.
%       (Default: 0)
%
%   tRamp
%       Ramp time of the gradients in seconds.
%       (Default: HW.Grad.tRamp)
%
%   GradTimeDelay
%       1x3 vector with time delays in seconds for each gradient channel.
%       Positive values lead to the gradient pulses being set earlier.
%       (Default: [0 0 0])
%
%   Gamma
%       Gyromagnetic ratio of the sample in rad/s/T.
%       (Default: HW.GammaDef)
%
%   useAQSlice
%       Index that is used for Seq.AQSlice for the following default parameters.
%       (Default: 1)
%
%   sizePhase
%       Size of the spectral encoding of the gradient pulses in meters.
%       (Default: Seq.AQSlice(Seq.Phase(t).useAQSlice).sizePhase(t)))
%
%   nPhase
%       Number of points in phase direction.
%       (Default: Seq.AQSlice(Seq.Phase(t).useAQSlice).nPhase(t))
%
%   PhaseOS
%       Over-sampling factor in phase direction.
%       (Default: Seq.AQSlice(Seq.Phase(t).useAQSlice).PhaseOS(t))
%
%   alfa / phi / theta
%       Angles in radians that activley rotate the direction of the phase
%       encoding (see UseCoordinate). The first rotation "alfa" is around the
%       x-axis, the second rotation "phi" is around the y-axis and the last
%       rotation "theta" is around the z-axis.
%       (Default: 0 or the corresponding values in
%       Seq.AQSlice(Seq.Phase(t).useAQSlice) )
%
%   distance
%       Phase coordinate of the origin of the gradient system in the image
%       coordinate system in meter. I.e., positive value moves the center of the
%       image towards positive phase direction. Calling functions might set this
%       from the value of
%       Seq.AQSlice(1).Center2OriginImage(Seq.Phase(t).UseCoordinate).
%       (Default: 0)
%
%   UseAtRepetitionTime
%       Repetition times (tReps) at which the gradient pulse are used.
%       (No default!)
%
%   UseAtRepetitionTimeRephase / UseAtRepetitionTimeDephase
%       Repetition times where the rephase/dephase pulses should be used.
%       (Mandatory! Default: UseAtRepetitionTime)
%
%   StepIncrement
%       Number of steps before the increment is added for this phase encoder.
%       (Default: 1)
%
%   kLineIncrement
%       Number of steps before the next k-line is encoded taking all phase encoders into account.
%       (Default: 1)
%
%   numKLines
%       Number of k-lines for one complete image.
%       (Default: numel(Seq.Phase(t).UseAtRepetitionTime))
%
%   StepOrder
%       Mapping of elements in "UseAtRepetitionTime",
%       "UseAtRepetitionTimeRephase", or "UseAtRepetitionTimeDephase" to the
%       phase encoding steps.
%       (Default: repelem(1:Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS, 1, ...
%                         Seq.Phase(t).StepIncrement(1)) )
%
%   usedkLines
%       Subset and/or order of k-lines to acquire. These are indices into the
%       k-lines of one complete image in default ordering.
%       (Default: 1:Seq.Phase(t).numKLines)
%
%
% "HW" is the HW object or structure.
%
%
% "idxModule" is an optional input argument with the indices of the phase
% modules for which the pulse program matrices should be created.
% (Default: 1:numel(Seq.Phase))
%
%
% OUTPUT:
%
% To each structure Seq.Phase the following fields are added:
%
%   GradRephase / GradDephase
%       Structure with the fields "Amp" and "Time" containing the amplitudes and
%       times of the rephase/dephase pulses corresponding to the input
%       parameters that can be added to the sequence with "add_Grad".
%
%   AQPhaseShift
%       Phase shift of the read out window for off-center measurements.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default parameters
if nargin < 3 || isempty(idxModule)
  idxModule = 1:numel(Seq.Phase);
end


%% loop over (selected) Seq.Phase "modules"
% if ~isfield(Seq.Phase, 'UseCoordinate'), Seq.Phase(1).UseCoordinate = 3; end
% if Seq.showSlicePhase, Seq.Phase(t).UseCoordinate = Seq.Slice.UseCoordinate; end
for t = idxModule
  %% get default or derived phase encoding parameters
  if isemptyfield(Seq.Phase(t), 'UseCoordinate')
    Seq.Phase(t).UseCoordinate = 3;
  end
  if isemptyfield(Seq.Phase(t), 'useAQSlice')
    % index in Seq.AQSlice
    Seq.Phase(t).useAQSlice = 1;
  end
  if isemptyfield(Seq.Phase(t), 'GradTimeDelay')
    Seq.Phase(t).GradTimeDelay = zeros(1, max(3, Seq.Phase(t).UseCoordinate));
  end
  if isemptyfield(Seq.Phase(t), 'CenterOfRephase')
    Seq.Phase(t).CenterOfRephase = Seq.Read(t).CenterOfRephase;
  end
  if isemptyfield(Seq.Phase(t), 'CenterOfDephase')
    Seq.Phase(t).CenterOfDephase = Seq.Read(t).CenterOfDephase;
  end
  if isemptyfield(Seq.Phase(t), 'GradRephaseLength')
    Seq.Phase(t).GradRephaseLength = Seq.Read(t).GradRephaseLength;
  end
  if isemptyfield(Seq.Phase(t), 'GradDephaseLength')
    Seq.Phase(t).GradDephaseLength = Seq.Read(t).GradDephaseLength;
  end
  if isemptyfield(Seq.Phase(t), 'GradRephaseSign')
    Seq.Phase(t).GradRephaseSign = -1;
  end
  if isemptyfield(Seq.Phase(t), 'GradDephaseSign')
    Seq.Phase(t).GradDephaseSign = 1;
  end
  if isemptyfield(Seq.Phase(t), 'GradTimeIntegralRephaseOffset')
    Seq.Phase(t).GradTimeIntegralRephaseOffset = 0;
  end
  if isemptyfield(Seq.Phase(t), 'GradTimeIntegralDephaseOffset')
    Seq.Phase(t).GradTimeIntegralDephaseOffset = 0;
  end
  if isemptyfield(Seq.Phase(t), 'iDevice')
    if isemptyfield(Seq, 'AQSlice') || ...
        isemptyfield(Seq.AQSlice(Seq.Phase(t).useAQSlice), 'iDevice')
      Seq.Phase(t).iDevice = 1;
    else
      Seq.Phase(t).iDevice = Seq.AQSlice(Seq.Phase(t).useAQSlice).iDevice;
    end
  end
  if isemptyfield(Seq.Phase(t), 'tRamp')
    Seq.Phase(t).tRamp = HW.Grad(Seq.Phase(t).iDevice).tRamp;
  end
  if isemptyfield(Seq.Phase(t), 'Gamma')
    if isemptyfield(Seq.Phase(t), 'Gamma')
      Seq.Phase(t).Gamma = HW.GammaDef;
    else
      Seq.Phase(t).Gamma = Seq.AQSlice(Seq.Phase(t).useAQSlice).Gamma;
    end
  end
  if isemptyfield(Seq.Phase(t), 'sizePhase')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Phase(t).useAQSlice), 'sizePhase')
      error('PD:get_PhaseParameter:NoSizePhase', ...
        'Seq.Phase.sizePhase must be set.');
    else
      Seq.Phase(t).sizePhase = Seq.AQSlice(Seq.Phase(t).useAQSlice).sizePhase(t);
    end
  end
  if isemptyfield(Seq.Phase(t), 'nPhase')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Phase(t).useAQSlice), 'nPhase')
      error('PD:get_PhaseParameter:NoNPhase', ...
        'Seq.Phase.nPhase must be set.');
    else
      Seq.Phase(t).nPhase = Seq.AQSlice(Seq.Phase(t).useAQSlice).nPhase(t);
    end
  end
  if isemptyfield(Seq.Phase(t), 'PhaseOS')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Phase(t).useAQSlice), 'PhaseOS')
      Seq.Phase(t).PhaseOS = 1;
    else
      Seq.Phase(t).PhaseOS = Seq.AQSlice(Seq.Phase(t).useAQSlice).PhaseOS(t);
    end
  end
  if isemptyfield(Seq.Phase(t), 'alfa')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Phase(t).useAQSlice), 'alfa')
      Seq.Phase(t).alfa = 0;
    else
      Seq.Phase(t).alfa = Seq.AQSlice(Seq.Phase(t).useAQSlice).alfa;
    end
  end
  if isemptyfield(Seq.Phase(t), 'phi')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Phase(t).useAQSlice), 'phi')
      Seq.Phase(t).phi = 0;
    else
      Seq.Phase(t).phi = Seq.AQSlice(Seq.Phase(t).useAQSlice).phi;
    end
  end
  if isemptyfield(Seq.Phase(t), 'theta')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(Seq.Phase(t).useAQSlice), 'theta')
      Seq.Phase(t).theta = 0;
    else
      Seq.Phase(t).theta = Seq.AQSlice(Seq.Phase(t).useAQSlice).theta;
    end
  end
  if isemptyfield(Seq.Phase(t), 'angle2Turns')
    % conversion factor of the input angles to full turns  e.g.: 1/(2*pi)
    if isemptyfield(Seq, 'AQSlice') || ...
        isemptyfield(Seq.AQSlice(Seq.Phase(t).useAQSlice), 'angle2Turns')
      Seq.Phase(t).angle2Turns = 1/(2*pi);
    else
      Seq.Phase(t).angle2Turns = Seq.AQSlice(Seq.Phase(t).useAQSlice).angle2Turns;
    end
  end
  if isemptyfield(Seq.Phase(t), 'distance')
    Seq.Phase(t).distance = 0;
  end
  if ~isfield(Seq.Phase(t), 'UseAtRepetitionTime')
    Seq.Phase(t).UseAtRepetitionTime = [];
  end
  if isemptyfield(Seq.Phase(t), 'UseAtRepetitionTimeRephase')
    Seq.Phase(t).UseAtRepetitionTimeRephase = Seq.Phase(t).UseAtRepetitionTime;
  end
  if isemptyfield(Seq.Phase(t), 'UseAtRepetitionTimeDephase')
    Seq.Phase(t).UseAtRepetitionTimeDephase = Seq.Phase(t).UseAtRepetitionTime;
  end
  if isemptyfield(Seq.Phase(t), 'StepIncrement')
    Seq.Phase(t).StepIncrement = 1;
  end
  if isemptyfield(Seq.Phase(t), 'kLineIncrement')
    Seq.Phase(t).kLineIncrement = 1;
  end
  if isemptyfield(Seq.Phase(t), 'StepOrder')
    % encoding order (for one cycle) in this encoding direction
    Seq.Phase(t).StepOrder = ...
      reshape(repmat(1:Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS, ...
                     Seq.Phase(t).StepIncrement(1)/Seq.Phase(t).kLineIncrement(1), 1), ...
              1, []);
  end
  if isemptyfield(Seq.Phase(t), 'numKLines')
    % number of k-lines per one (complete) image
    Seq.Phase(t).numKLines = numel(Seq.Phase(t).UseAtRepetitionTime);
  end
  if isemptyfield(Seq.Phase(t), 'usedkLines')
    % subset and/or order of k-lines to acquire
    Seq.Phase(t).usedkLines = 1:Seq.Phase(t).numKLines;
  end
  if isemptyfield(Seq.Phase(t), 'Overdrive')
    Seq.Phase(t).Overdrive = 0;
  end


  Seq.Phase(t).Resolution = Seq.Phase(t).sizePhase/Seq.Phase(t).nPhase;
  Seq.Phase(t).GradTimeIntegral = pi/Seq.Phase(t).Gamma/Seq.Phase(t).Resolution;

  if any(round((Seq.Phase(t).GradDephaseLength-Seq.Phase(t).tRamp*2)*HW.MMRT(Seq.Phase(t).iDevice).fSystem) < 2)
    error('PD:get_PhaseParameter:GradDephaseLengthTooShort', ...
      'Seq.Phase(%d).GradDephaseLength too short', t);
  end
  if any(round((Seq.Phase(t).GradRephaseLength-Seq.Phase(t).tRamp*2)*HW.MMRT(Seq.Phase(t).iDevice).fSystem) < 2)
    error('PD:get_PhaseParameter:GradRephaseLengthTooShort', ...
      'Seq.Phase(%d).GradRephaseLength too short', t);
  end

  Seq.Phase(t).GradAmpDephase = ...
    (Seq.Phase(t).GradTimeIntegral + Seq.Phase(t).GradTimeIntegralDephaseOffset) ...
    ./ (Seq.Phase(t).GradDephaseLength - Seq.Phase(t).tRamp);
  Seq.Phase(t).GradAmpRephase = ...
    (Seq.Phase(t).GradTimeIntegral + Seq.Phase(t).GradTimeIntegralRephaseOffset) ...
    ./ (Seq.Phase(t).GradRephaseLength - Seq.Phase(t).tRamp);

  if any(abs(Seq.Phase(t).GradAmpDephase/(Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS)) > 0 ...
      & abs(Seq.Phase(t).GradAmpDephase/(Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS)) ...
          < HW.Grad(Seq.Phase(t).iDevice).Dac2Amp(HW.Grad(Seq.Phase(t).iDevice).xyzB(Seq.Phase(t).UseCoordinate)))
    gradAmpDephaseStep = abs(Seq.Phase(t).GradAmpDephase/(Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS));
    gradAmpDephaseStep = min(gradAmpDephaseStep(gradAmpDephaseStep > 0));
    % This can be elevated to an error with
    %   warning('error', 'PD:get_PhaseParameter:smallStep');
    warning('PD:get_PhaseParameter:smallStep', ...
      'Gradient amplitude between phase dephasing steps (%.3fe-6) is smaller than one LSB (%.3fe-6) for %s.', ...
      gradAmpDephaseStep*1e6, ...
      HW.Grad(Seq.Phase(t).iDevice).Dac2Amp(HW.Grad(Seq.Phase(t).iDevice).xyzB(Seq.Phase(t).UseCoordinate))*1e6, ...
      HW.Grad(Seq.Phase(t).iDevice).Name{Seq.Phase(t).UseCoordinate});
  end
  if any(abs(Seq.Phase(t).GradAmpRephase/(Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS)) > 0 ...
      & abs(Seq.Phase(t).GradAmpRephase/(Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS)) ...
          < HW.Grad(Seq.Phase(t).iDevice).Dac2Amp(HW.Grad(Seq.Phase(t).iDevice).xyzB(Seq.Phase(t).UseCoordinate)))
    gradAmpRephaseStep = abs(Seq.Phase(t).GradAmpRephase/(Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS));
    gradAmpRephaseStep = min(gradAmpRephaseStep(gradAmpRephaseStep > 0));
    % This can be elevated to an error with
    %   warning('error', 'PD:get_PhaseParameter:smallStep');
    warning('PD:get_PhaseParameter:smallStep', ...
      'Gradient amplitude between phase rephasing steps (%.3fe-6) is smaller than one LSB (%.3fe-6) for %s.', ...
      gradAmpRephaseStep*1e6, ...
      HW.Grad(Seq.Phase(t).iDevice).Dac2Amp(HW.Grad(Seq.Phase(t).iDevice).xyzB(Seq.Phase(t).UseCoordinate))*1e6, ...
      HW.Grad(Seq.Phase(t).iDevice).Name{Seq.Phase(t).UseCoordinate});
  end

  temp = zeros(3, numel(Seq.Phase(t).GradAmpDephase));
  temp(Seq.Phase(t).UseCoordinate,:) = Seq.Phase(t).GradAmpDephase;
  Seq.Phase(t).GradAmpDephase = temp;
  temp = zeros(3, numel(Seq.Phase(t).GradAmpRephase));
  temp(Seq.Phase(t).UseCoordinate,:) = Seq.Phase(t).GradAmpRephase;
  Seq.Phase(t).GradAmpRephase = temp;

  if Seq.Phase(t).UseCoordinate < 4
    Angle2Deg = Seq.Phase(t).angle2Turns * 360;
    [Rx, Ry, Rz] = get_aptDegRotationMatrix(Seq.Phase(t).alfa*Angle2Deg, Seq.Phase(t).phi*Angle2Deg, Seq.Phase(t).theta*Angle2Deg);
    Seq.Phase(t).GradAmpRephase = Rz*(Ry*(Rx*Seq.Phase(t).GradAmpRephase));
    Seq.Phase(t).GradAmpDephase = Rz*(Ry*(Rx*Seq.Phase(t).GradAmpDephase));
  else
    keyboard();
  end

  % encoding steps for all images in default ordering (chronologically)
  usedkLinesImages = reshape(repmat(Seq.Phase(t).usedkLines, Seq.Phase(t).kLineIncrement, 1), 1, []);
  usedkLinesImages = repmat(usedkLinesImages, 1, ceil(numel(Seq.Phase(t).UseAtRepetitionTime)/numel(usedkLinesImages)));
  Seq.Phase(t).StepOrder = repmat(Seq.Phase(t).StepOrder, 1, ceil(Seq.Phase(t).numKLines/numel(Seq.Phase(t).StepOrder)));
  kLineOrder = Seq.Phase(t).StepOrder(usedkLinesImages);
  % kLineOrder = Seq.Phase(t).usedkLines(Seq.Phase(t).StepOrder);  % original

  if mod(Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS, 2)
    % odd number of lines in k-space
    phaseSteps = linspace(-1+1/(Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS), ...
      1-1/(Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS), ...
      Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS);
  else
    % even number of lines in k-space
    phaseSteps = linspace(-1, ...
      1-2/(Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS), ...
      Seq.Phase(t).nPhase*Seq.Phase(t).PhaseOS);
  end
  phaseSteps = phaseSteps(kLineOrder);

  Seq.Phase(t).GradAmpDephase = bsxfun(@times, Seq.Phase(t).GradAmpDephase, phaseSteps);
  Seq.Phase(t).GradAmpRephase = bsxfun(@times, Seq.Phase(t).GradAmpRephase, phaseSteps);
  Seq.Phase(t).AQPhaseShift = (-Seq.Phase(t).distance)/Seq.Phase(t).sizePhase*pi*Seq.Phase(t).nPhase * phaseSteps;
  Seq.Phase(t).AQPhaseShift = Seq.Phase(t).AQPhaseShift(kLineOrder);


  %% translate to pulse program matrices
  tempsize = ones(1, numel(Seq.Phase(t).UseAtRepetitionTime));

  for n = unique([1:3, Seq.Phase(t).UseCoordinate])
  % if Seq.Phase(t).Overdrive == 0
  %   Seq.Phase(t).GradRephase(n).Amp=nan+zeros(4,size(tempsize,2));
  %   Seq.Phase(t).GradDephase(n).Amp=nan+zeros(4,size(tempsize,2));
  %   Seq.Phase(t).GradRephase(n).Amp([1,4],~isnan(tempsize))=0;
  %   Seq.Phase(t).GradRephase(n).Amp(2,~isnan(tempsize))=Seq.Phase(t).GradAmpRephase(n,Seq.Phase(t).StepOrder);
  %   Seq.Phase(t).GradRephase(n).Amp(3,~isnan(tempsize))=Seq.Phase(t).GradAmpRephase(n,Seq.Phase(t).StepOrder);
  %   Seq.Phase(t).GradDephase(n).Amp([1,4],~isnan(tempsize))=0;
  %   Seq.Phase(t).GradDephase(n).Amp(2,~isnan(tempsize))=Seq.Phase(t).GradAmpDephase(n,Seq.Phase(t).StepOrder);
  %   Seq.Phase(t).GradDephase(n).Amp(3,~isnan(tempsize))=Seq.Phase(t).GradAmpDephase(n,Seq.Phase(t).StepOrder);

    Seq.Phase(t).GradRephase(n).Amp = nan(4, size(Seq.tRep,2));
    Seq.Phase(t).GradDephase(n).Amp = Seq.Phase(t).GradRephase(n).Amp;

    Seq.Phase(t).GradRephase(n).Amp(1,Seq.Phase(t).UseAtRepetitionTimeRephase) = 0;
    Seq.Phase(t).GradRephase(n).Amp(2:3,Seq.Phase(t).UseAtRepetitionTimeRephase) = ...
      [Seq.Phase(t).GradAmpRephase(n,:).*Seq.Phase(t).GradRephaseSign; ...
       Seq.Phase(t).GradAmpRephase(n,:).*Seq.Phase(t).GradRephaseSign];
    Seq.Phase(t).GradRephase(n).Amp(4,Seq.Phase(t).UseAtRepetitionTimeRephase) = 0;

    Seq.Phase(t).GradDephase(n).Amp(1,Seq.Phase(t).UseAtRepetitionTimeDephase) = 0;
    Seq.Phase(t).GradDephase(n).Amp(2:3,Seq.Phase(t).UseAtRepetitionTimeDephase) = ...
      [Seq.Phase(t).GradAmpDephase(n,:).*Seq.Phase(t).GradDephaseSign; ...
       Seq.Phase(t).GradAmpDephase(n,:).*Seq.Phase(t).GradDephaseSign];
    Seq.Phase(t).GradDephase(n).Amp(4,Seq.Phase(t).UseAtRepetitionTimeDephase) = 0;

    Seq.Phase(t).GradRephase(n).Time = NaN(4, size(Seq.tRep, 2));
    Seq.Phase(t).GradDephase(n).Time = Seq.Phase(t).GradRephase(n).Time;

    Seq.Phase(t).GradRephase(n).Time(:,Seq.Phase(t).UseAtRepetitionTimeRephase) = ...
      bsxfun(@plus, Seq.Phase(t).CenterOfRephase(:).' .* tempsize, ...
                    [-Seq.Phase(t).GradRephaseLength/2; ...
                     -Seq.Phase(t).GradRephaseLength/2+Seq.Phase(t).tRamp; ...
                     +Seq.Phase(t).GradRephaseLength/2-Seq.Phase(t).tRamp; ...
                     +Seq.Phase(t).GradRephaseLength/2]) - Seq.Phase(t).GradTimeDelay(n);
    Seq.Phase(t).GradDephase(n).Time(:,Seq.Phase(t).UseAtRepetitionTimeDephase) = ...
      bsxfun(@plus, Seq.Phase(t).CenterOfDephase(:).' .* tempsize, ...
                    [-Seq.Phase(t).GradDephaseLength/2; ...
                     -Seq.Phase(t).GradDephaseLength/2+Seq.Phase(t).tRamp; ...
                     +Seq.Phase(t).GradDephaseLength/2-Seq.Phase(t).tRamp; ...
                     +Seq.Phase(t).GradDephaseLength/2]) - Seq.Phase(t).GradTimeDelay(n);

    % replace Inf by NaN
    Seq.Phase(t).GradRephase(n).Time(isinf(Seq.Phase(t).GradRephase(n).Time)) = NaN;
    Seq.Phase(t).GradDephase(n).Time(isinf(Seq.Phase(t).GradDephase(n).Time)) = NaN;
  % else
  %   Seq.Phase(t).GradRephase(n).Amp=nan+zeros(8,size(tempsize,2));
  %   Seq.Phase(t).GradDephase(n).Amp=nan+zeros(8,size(tempsize,2));
  %   Seq.Phase(t).GradRephase(n).Amp([1,8],~isnan(tempsize))=0;
  %   Seq.Phase(t).GradRephase(n).Amp(2,~isnan(tempsize))=Seq.Phase(t).GradAmpRephase(n,Seq.Phase(t).StepOrder)*(1+Seq.Phase(t).Overdrive);
  %   Seq.Phase(t).GradRephase(n).Amp(3,~isnan(tempsize))=Seq.Phase(t).GradAmpRephase(n,Seq.Phase(t).StepOrder)*(1+Seq.Phase(t).Overdrive);
  %   Seq.Phase(t).GradRephase(n).Amp(4,~isnan(tempsize))=Seq.Phase(t).GradAmpRephase(n,Seq.Phase(t).StepOrder);
  %   Seq.Phase(t).GradRephase(n).Amp(5,~isnan(tempsize))=Seq.Phase(t).GradAmpRephase(n,Seq.Phase(t).StepOrder);
  %   Seq.Phase(t).GradRephase(n).Amp(6,~isnan(tempsize))=Seq.Phase(t).GradAmpRephase(n,Seq.Phase(t).StepOrder)*(0-Seq.Phase(t).Overdrive);
  %   Seq.Phase(t).GradRephase(n).Amp(7,~isnan(tempsize))=Seq.Phase(t).GradAmpRephase(n,Seq.Phase(t).StepOrder)*(0-Seq.Phase(t).Overdrive);
  %   Seq.Phase(t).GradDephase(n).Amp([1,8],~isnan(tempsize))=0;
  %   Seq.Phase(t).GradDephase(n).Amp(2,~isnan(tempsize))=Seq.Phase(t).GradAmpDephase(n,Seq.Phase(t).StepOrder)*(1+Seq.Phase(t).Overdrive);
  %   Seq.Phase(t).GradDephase(n).Amp(3,~isnan(tempsize))=Seq.Phase(t).GradAmpDephase(n,Seq.Phase(t).StepOrder)*(1+Seq.Phase(t).Overdrive);
  %   Seq.Phase(t).GradDephase(n).Amp(4,~isnan(tempsize))=Seq.Phase(t).GradAmpDephase(n,Seq.Phase(t).StepOrder);
  %   Seq.Phase(t).GradDephase(n).Amp(5,~isnan(tempsize))=Seq.Phase(t).GradAmpDephase(n,Seq.Phase(t).StepOrder);
  %   Seq.Phase(t).GradDephase(n).Amp(6,~isnan(tempsize))=Seq.Phase(t).GradAmpDephase(n,Seq.Phase(t).StepOrder)*(0-Seq.Phase(t).Overdrive);
  %   Seq.Phase(t).GradDephase(n).Amp(7,~isnan(tempsize))=Seq.Phase(t).GradAmpDephase(n,Seq.Phase(t).StepOrder)*(0-Seq.Phase(t).Overdrive);
  %
  %
  %
  %   Seq.Phase(t).GradRephase(n).Time=Seq.Phase(t).CenterOfRephase...
  %                                                       +[  -Seq.Phase(t).GradRephaseLength/2;...
  %                                                           -Seq.Phase(t).GradRephaseLength/2+Seq.Phase(t).tRamp*1/3;...
  %                                                           -Seq.Phase(t).GradRephaseLength/2+Seq.Phase(t).tRamp*2/3;...
  %                                                           -Seq.Phase(t).GradRephaseLength/2+Seq.Phase(t).tRamp*3/3;...
  %                                                           +Seq.Phase(t).GradRephaseLength/2-Seq.Phase(t).tRamp*3/3;...
  %                                                           +Seq.Phase(t).GradRephaseLength/2-Seq.Phase(t).tRamp*2/3;...
  %                                                           +Seq.Phase(t).GradRephaseLength/2-Seq.Phase(t).tRamp*1/3;...
  %                                                           +Seq.Phase(t).GradRephaseLength/2]*tempsize-Seq.Phase(t).GradTimeDelay(n);
  %   Seq.Phase(t).GradDephase(n).Time=Seq.Phase(t).CenterOfDephase...
  %                                                       +[  -Seq.Phase(t).GradDephaseLength/2;...
  %                                                           -Seq.Phase(t).GradDephaseLength/2+Seq.Phase(t).tRamp*1/3;...
  %                                                           -Seq.Phase(t).GradDephaseLength/2+Seq.Phase(t).tRamp*2/3;...
  %                                                           -Seq.Phase(t).GradDephaseLength/2+Seq.Phase(t).tRamp*3/3;...
  %                                                           +Seq.Phase(t).GradDephaseLength/2-Seq.Phase(t).tRamp*3/3;...
  %                                                           +Seq.Phase(t).GradDephaseLength/2-Seq.Phase(t).tRamp*2/3;...
  %                                                           +Seq.Phase(t).GradDephaseLength/2-Seq.Phase(t).tRamp*1/3;...
  %                                                           +Seq.Phase(t).GradDephaseLength/2]*tempsize-Seq.Phase(t).GradTimeDelay(n);
  % end
  end
end

end
