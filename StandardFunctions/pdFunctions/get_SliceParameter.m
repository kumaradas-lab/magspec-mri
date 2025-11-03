function Seq = get_SliceParameter(Seq, HW)
%% Calculate slice select pulse
%
%   Seq = get_SliceParameter(Seq, HW);
%
% This function calculates slice gradients (and respective rephase and dephase
% spoilers) derived from the values in Seq.Slice. It also creates the
% corresponding excitation pulse.
%
% INPUT:
% "Seq" is a structure with the field "Slice" that is an (array of) structures
% with the following fields:
%   UseCoordinate   Number of the coordinate for which the phase encoder is used
%                   (default: 1).
%   Thickness       Thickness of slice that is excited in meter (default: 1e12).
%   CenterOfPulse   Center of excitation pulse in seconds on the tRep timeline
%                   (default: 0).
%   MaxGradAmp      Maximum gradient amplitude in T/m (default:
%                   min(HW.Grad.MaxAmp(1:3)) ).
%   GradSign        Sign of the slice selection pulse (default: 1).
%   tRamp           Ramp time of the gradients in seconds (default:
%                   HW.Grad.tRamp).
%   tEC             Additional in seconds time that the gradients need to reach
%                   their set amplitude due to Eddy currents in the pole shoe
%                   (default: HW.Grad.tEC).
%   GradTimeDelay   1x3 vector with time delays in seconds for each gradient
%                   channel. Positive values lead to the gradient pulses being
%                   set earlier. (Default: [0 0 0]).
%   useAQSlice      Index that is used for Seq.AQSlice for the following default
%                   parameters (default: 1).
%   alfa / phi / theta
%                   Angles in radians that activley rotate the direction of the
%                   selected slice (see UseCoordinate). The first rotation
%                   "alfa" is around the x-axis, the second rotation "phi" is
%                   around the y-axis and the last rotation "theta" is around
%                   the z-axis.
%   UseAtRepetitionTime
%                   Repetition times (tReps) at which the gradient pulse are
%                   used (mandatory! no default!).
%   UseAtRepetitionTimeRephase / UseAtRepetitionTimeDephase
%                   Repetition times where the rephase/dephase pulses should be
%                   used (mandatory! default: UseAtRepetitionTime).
%   distance        Displacement of the selected slice to the center of the
%                   gradient system in image coordinate system in meter
%                   (default: 0). I.e. a positive value moves the selected slice
%                   towards negative slice direction.
%   Pulse           Structure for the TX pulse with the following fields:
%     Function        Function handle to the Pulse shape function. (default:
%                     @Pulse_Rect).
%     MaxNumberOfSegments
%                     Maximum number of segments for the pulse shape (default:
%                     51)
%     Amplitude       Amplitude of the TX pulse in Tesla (default:
%                     HW.TX.AmpDef).
%     FlipAngle       Flip angle in degrees (default: 90).
% For the spoilers:
%   GradRephaseLength / GradDephaseLength
%                   Length of the rephase/dephase pulse in seconds (default: ).
%   CenterOfRephase / CenterOfDephase
%                   Center of rephase/dephase pulse in seconds on the tRep
%                   timeline (default: ).
%   GradRephaseSign / GradDephaseSign
%                   Sign of the rephase/dephase pulse (default: GradRephaseSign
%                   = -1, GradDephaseSign = 1).
%
% OUTPUT:
% To each structure Seq.Phase the following fields are added:
%   Grad / GradRephase / GradDephase
%                   Structure with the fields "Amp" and "Time" containing the
%                   amplitudes and times of the slice selection and
%                   rephase/dephase pulses corresponding to the input parameters
%                   that can be added to the sequence with "add_Grad".
%   TX              Structure containing the settings for the (excitation)
%                   pulse. It can be added to the sequence with "add_TX".
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%%
for t = 1:numel(Seq.Slice)
  if isemptyfield(Seq.Slice(t), 'useAQSlice'), Seq.Slice(t).useAQSlice = 1; end
  if isemptyfield(Seq.Slice(t), 'UseCoordinate'), Seq.Slice(t).UseCoordinate = 1; end
  if isemptyfield(Seq.Slice(t), 'iDevice')
    if isemptyfield(Seq.AQSlice(Seq.Slice(t).useAQSlice), 'iDevice')
      Seq.Slice(t).iDevice = 1;
    else
      Seq.Slice(t).iDevice = Seq.AQSlice(Seq.Slice(t).useAQSlice).iDevice;
    end
  end
  if ~isfield(Seq.Slice, 'Pulse'), Seq.Slice(t).Pulse = struct(); end
  if isemptyfield(Seq.Slice(t).Pulse, 'MaxNumberOfSegments'), Seq.Slice(t).Pulse.MaxNumberOfSegments = 51; end
  if isemptyfield(Seq.Slice(t).Pulse, 'Amplitude'), Seq.Slice(t).Pulse.Amplitude = HW.TX(Seq.Slice(t).iDevice).AmpDef; end
  if isemptyfield(Seq.Slice(t).Pulse, 'FlipAngle'), Seq.Slice(t).Pulse.FlipAngle = 90; end
  if isemptyfield(Seq.Slice(t).Pulse, 'Function'), Seq.Slice(t).Pulse.Function = @Pulse_Rect; end
  if isemptyfield(Seq.Slice(t), 'tRamp'), Seq.Slice(t).tRamp = HW.Grad(Seq.Slice(t).iDevice).tRamp; end
  if isemptyfield(Seq.Slice(t), 'tEC'), Seq.Slice(t).tEC = HW.Grad(Seq.Slice(t).iDevice).tEC; end
  if isemptyfield(Seq.Slice(t), 'MaxGradAmp'), Seq.Slice(t).MaxGradAmp = min(HW.Grad(Seq.Slice(t).iDevice).MaxAmp(1:3)); end
  if isemptyfield(Seq.Slice(t), 'Thickness'), Seq.Slice(t).Thickness = 1e12; end
  if isemptyfield(Seq.Slice(t), 'CenterOfPulse'), Seq.Slice(t).CenterOfPulse = 0; end
  if isemptyfield(Seq.Slice(t), 'GradSign'), Seq.Slice(t).GradSign = 1; end
  if isemptyfield(Seq.Slice(t), 'distance'), Seq.Slice(t).distance = 0; end
  if isemptyfield(Seq.Slice(t), 'Gamma')
    if isemptyfield(Seq.AQSlice(Seq.Slice(t).useAQSlice), 'Gamma')
      Seq.Slice(t).Gamma = HW.GammaDef;
    else
      Seq.Slice(t).Gamma = Seq.AQSlice(Seq.Slice(t).useAQSlice).Gamma;
    end
  end

  Seq.Slice(t).GradLength = (max(max(Seq.Slice(t).Pulse.FlipAngle./(Seq.Slice(t).Pulse.Amplitude * HW.TX(Seq.Slice(t).iDevice).Amp2Deg) * Seq.Slice(t).Pulse.Function(HW, 'Amp'), ...
                                  1./(min(min(HW.Grad(Seq.Slice(t).iDevice).MaxAmp(1:3), min(Seq.Slice(t).MaxGradAmp(:)))) * Seq.Slice(t).Gamma/2/pi * Seq.Slice(t).Thickness) * Seq.Slice(t).Pulse.Function(HW, 'Time')), ...
                              Seq.Slice(t).Pulse.Function(HW,'Length'))...
                              +Seq.Slice(t).tRamp*2+Seq.Slice(t).tEC*2);

  Seq.Slice(t).Pulse.MaxLength = Seq.Slice(t).GradLength-Seq.Slice(t).tRamp*2-Seq.Slice(t).tEC*2;
  Seq.Slice(t).Pulse.Bandwidth = max(Seq.Slice(t).Pulse.Function(HW,'Time')./Seq.Slice(t).Pulse.MaxLength,...
                                   Seq.Slice(t).Pulse.Function(HW,'Bandwidth'));
  Seq.Slice(t).GradLength = max(Seq.Slice(t).GradLength);
  Seq.Slice(t).Pulse.MaxLength = max(Seq.Slice(t).Pulse.MaxLength);
  % Calculate gradient strength during slice pulse
  Seq.Slice(t).GradAmp = Seq.Slice(t).Pulse.Bandwidth .* (2*pi/Seq.Slice(t).Gamma) ./ Seq.Slice(t).Thickness .* Seq.Slice(t).GradSign;
  % Calculate frequency of the TX slice pulse, it has an offset if the Slice does not touch the Center
  Seq.Slice(t).Pulse.Frequency = HW.fLarmor - (Seq.Slice(t).Gamma/2/pi) .* Seq.Slice(t).distance .* Seq.Slice(t).GradAmp;

  if abs(nargin(Seq.Slice(t).Pulse.Function)) > 4
    for tt = 1:numel(Seq.Slice(t).GradAmp)
      TXs(tt) = Seq.Slice(t).Pulse.Function(HW, ...
                                            Seq.Slice(t).CenterOfPulse, ...
                                            Seq.Slice(t).Pulse.Bandwidth(tt), ...
                                            Seq.Slice(t).Pulse.FlipAngle(tt)/180*pi, ...
                                            Seq.Slice(t).Pulse.MaxNumberOfSegments, ...
                                            Seq.Slice(t).Pulse.MaxLength, ...
                                            Seq.Slice(t).Pulse.Frequency(tt), 0);
    end
  else
    Pulse = Seq.Slice(t).Pulse;
    Pulse.Phase = 0;
    Pulse.FlipAngleFullTurn = 360; % FlipAngle is passed in degrees
    if isemptyfield(Pulse, 'iDevice'), Pulse.iDevice = Seq.Slice(t).iDevice; end
    for tt = 1:numel(Seq.Slice(t).GradAmp)
      Pulse.Bandwidth = Seq.Slice(t).Pulse.Bandwidth(tt);
      Pulse.FlipAngle = Seq.Slice(t).Pulse.FlipAngle(tt);
      Pulse.Frequency = Seq.Slice(t).Pulse.Frequency(tt);
      TXs(tt) = Seq.Slice(t).Pulse.Function(HW, Seq.Slice(t).CenterOfPulse, Pulse);
    end
  end

  TXs = correct_PulsePhase(TXs, HW, Seq.Slice(t).iDevice);
  Seq.Slice(t).Pulse.CenterOffset = TXs(1).CenterOffset;

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
  if isemptyfield(Seq.Slice(t), 'GradRephaseLength'), Seq.Slice(t).GradRephaseLength = (Seq.Slice(t).GradLength-Seq.Slice(t).tRamp)*tR+Seq.Slice(t).tRamp; end
  if isemptyfield(Seq.Slice(t), 'CenterOfRephase')
    Seq.Slice(t).CenterOfRephase = Seq.Slice(t).GradLength/2-Seq.Slice(t).Pulse.CenterOffset+Seq.Slice(t).GradRephaseLength/2+Seq.Slice(t).CenterOfPulse-Seq.Slice(t).tRamp;
  end
  if isemptyfield(Seq.Slice(t),'GradRephaseSign'), Seq.Slice(t).GradRephaseSign = -1; end
  if isemptyfield(Seq.Slice(t),'GradDephaseLength'), Seq.Slice(t).GradDephaseLength = (Seq.Slice(t).GradLength-Seq.Slice(t).tRamp)*tD+Seq.Slice(t).tRamp; end
  if isemptyfield(Seq.Slice(t), 'CenterOfDephase')
    Seq.Slice(t).CenterOfDephase = -Seq.Slice(t).GradLength/2-Seq.Slice(t).Pulse.CenterOffset-Seq.Slice(t).GradDephaseLength/2+Seq.Slice(t).CenterOfPulse+Seq.Slice(t).tRamp;
  end
  if isemptyfield(Seq.Slice(t), 'GradDephaseSign'), Seq.Slice(t).GradDephaseSign = -1; end
  if isemptyfield(Seq.Slice(t), 'Thickness'), Seq.Slice(t).Thickness = Seq.AQSlice(Seq.Slice(t).useAQSlice).thickness; end
  if isemptyfield(Seq.Slice(t), 'alfa'), Seq.Slice(t).alfa = Seq.AQSlice(Seq.Slice(t).useAQSlice).alfa; end
  if isemptyfield(Seq.Slice(t), 'phi'), Seq.Slice(t).phi = Seq.AQSlice(Seq.Slice(t).useAQSlice).phi; end
  if isemptyfield(Seq.Slice(t), 'theta'), Seq.Slice(t).theta = Seq.AQSlice(Seq.Slice(t).useAQSlice).theta; end
  if isemptyfield(Seq.Slice(t), 'angle2Turns')  % conversion factor of the input angles to full turns e.g.: 1/(2*pi)
    if isemptyfield(Seq.AQSlice(Seq.Slice(t).useAQSlice), 'angle2Turns')
      Seq.Slice(t).angle2Turns = 1/(2*pi);
    else
      Seq.Slice(t).angle2Turns = Seq.AQSlice(Seq.Slice(t).useAQSlice).angle2Turns;
    end
  end
  if isemptyfield(Seq.Slice(t).Pulse, 'PhaseIncrement'), Seq.Slice(t).Pulse.PhaseIncrement = 0; end
  if isemptyfield(Seq.Slice(t).Pulse, 'Phase'), Seq.Slice(t).Pulse.Phase = 0; end
  if isemptyfield(Seq.Slice(t), 'UseAtRepetitionTime'), Seq.Slice(t).UseAtRepetitionTime = []; end % !!!
  if isemptyfield(Seq.Slice(t), 'UseAtRepetitionTimeDephase'), Seq.Slice(t).UseAtRepetitionTimeDephase = Seq.Slice(t).UseAtRepetitionTime; end
  if isemptyfield(Seq.Slice(t), 'UseAtRepetitionTimeRephase'), Seq.Slice(t).UseAtRepetitionTimeRephase = Seq.Slice(t).UseAtRepetitionTime; end
  if isemptyfield(Seq.Slice(t), 'Overdrive'), Seq.Slice(t).Overdrive = 0; end
  if isemptyfield(Seq.Slice(t), 'GradTimeIntegralRephaseOffset'), Seq.Slice(t).GradTimeIntegralRephaseOffset = 0; end
  if isemptyfield(Seq.Slice(t), 'GradTimeIntegralDephaseOffset'), Seq.Slice(t).GradTimeIntegralDephaseOffset = 0; end

  if Seq.Slice(t).Thickness <= 1000 && Seq.Slice(t).GradDephaseLength-Seq.Slice(t).tRamp*2 < 2/HW.MMRT(Seq.Slice(t).iDevice).fSystem
    error(['Seq.Slice(' num2str(t) ').GradDephaseLength too short']);
  end
  if Seq.Slice(t).Thickness <= 1000 && Seq.Slice(t).GradRephaseLength-Seq.Slice(t).tRamp*2 < 2/HW.MMRT(Seq.Slice(t).iDevice).fSystem
    error(['Seq.Slice(' num2str(t) ').GradRephaseLength too short']);
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
      temp=repmat(Seq.Slice(t).Pulse.PhaseIncrement(:),1,numel(Seq.Slice(t).UseAtRepetitionTime)/numel(Seq.Slice(t).Pulse.PhaseIncrement));
      tempPhaseInc=cumsum(temp(:)).';
  else
      tempPhaseInc=cumsum(Seq.Slice(t).Pulse.PhaseIncrement*tempsize);
  end
  if numel(Seq.Slice(t).Pulse.Phase)>1
      temp=repmat(Seq.Slice(t).Pulse.Phase(:),1,numel(Seq.Slice(t).UseAtRepetitionTime)/numel(Seq.Slice(t).Pulse.Phase));
      tempPhase=(temp(:)).';
  else
      tempPhase=(Seq.Slice(t).Pulse.Phase*tempsize);
  end

  % end
  % tempsize(~Seq.Slice(t).UseAtRepetitionTime)=nan;

  Seq.Slice(t).TX.Start = NaN(size(TXs(1).Start,1), size(Seq.tRep,2));
  Seq.Slice(t).TX.Duration = NaN(size(TXs(1).Duration,1), size(Seq.tRep,2));
  Seq.Slice(t).TX.Amplitude = NaN(size(TXs(1).Amplitude,1), size(Seq.tRep,2));
  Seq.Slice(t).TX.Frequency = NaN(size(TXs(1).Frequency,1), size(Seq.tRep,2));
  Seq.Slice(t).TX.Phase = NaN(size(TXs(1).Phase,1), size(Seq.tRep,2));

  Seq.Slice(t).TX.Start(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXs.Start], tempsize);
  Seq.Slice(t).TX.Duration(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXs.Duration], tempsize);
  Seq.Slice(t).TX.Amplitude(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXs.Amplitude], tempsize);
  Seq.Slice(t).TX.Frequency(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@times, [TXs.Frequency], tempsize);
  Seq.Slice(t).TX.Phase(:,Seq.Slice(t).UseAtRepetitionTime) = bsxfun(@plus, [TXs.Phase], ones(size(TXs(1).Phase,1),1) * tempPhaseInc + ...
    ones(size(TXs(1).Phase,1),1)*tempPhase);
  clear TXs

  Angle2Deg=Seq.Slice(t).angle2Turns*360;
  [Rx, Ry, Rz] = get_aptDegRotationMatrix(Seq.Slice(t).alfa*Angle2Deg, Seq.Slice(t).phi*Angle2Deg, Seq.Slice(t).theta*Angle2Deg);

  temp=zeros(3,1);
  temp(Seq.Slice(t).UseCoordinate)=1;
  Seq.Slice(t).GradAmpUnitVector=temp;
  temp=zeros(3,1);
  temp(Seq.Slice(t).UseCoordinate)=1;
  Seq.Slice(t).GradAmpRephaseUnitVector=temp;
  temp=zeros(3,1);
  temp(Seq.Slice(t).UseCoordinate)=1;
  Seq.Slice(t).GradAmpDephaseUnitVector=temp;

  Seq.Slice(t).GradAmpUnitVector=Rz*(Ry*(Rx*Seq.Slice(t).GradAmpUnitVector));
  Seq.Slice(t).GradAmpRephaseUnitVector=Rz*(Ry*(Rx*Seq.Slice(t).GradAmpRephaseUnitVector));
  Seq.Slice(t).GradAmpDephaseUnitVector=Rz*(Ry*(Rx*Seq.Slice(t).GradAmpDephaseUnitVector));

  for n = 1:3
    % if Seq.Slice(t).Overdrive==0;
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
    Seq.Slice(t).Grad(n).Amp(2:3,Seq.Slice(t).UseAtRepetitionTime) = ...
      [Seq.Slice(t).GradAmpUnitVector(n); Seq.Slice(t).GradAmpUnitVector(n)] * ...
      (Seq.Slice(t).GradAmp .* tempsize);
    Seq.Slice(t).Grad(n).Amp(4,Seq.Slice(t).UseAtRepetitionTime) = 0;

    Seq.Slice(t).GradRephase(n).Amp(1,Seq.Slice(t).UseAtRepetitionTimeRephase) = 0;
    Seq.Slice(t).GradRephase(n).Amp(2:3,Seq.Slice(t).UseAtRepetitionTimeRephase) = ...
      [Seq.Slice(t).GradAmpRephaseUnitVector(n); Seq.Slice(t).GradAmpRephaseUnitVector(n)] * ...
      (Seq.Slice(t).GradAmpRephase .* tempsize .* Seq.Slice(t).GradRephaseSign);
    Seq.Slice(t).GradRephase(n).Amp(4,Seq.Slice(t).UseAtRepetitionTimeRephase) = 0;

    Seq.Slice(t).GradDephase(n).Amp(1,Seq.Slice(t).UseAtRepetitionTimeDephase) = 0;
    Seq.Slice(t).GradDephase(n).Amp(2:3,Seq.Slice(t).UseAtRepetitionTimeDephase) = ...
      [Seq.Slice(t).GradAmpDephaseUnitVector(n); Seq.Slice(t).GradAmpDephaseUnitVector(n)] * ...
      (Seq.Slice(t).GradAmpDephase .* tempsize .* Seq.Slice(t).GradDephaseSign);
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
                                       +Seq.Slice(t).GradRephaseLength/2]*tempsize - ...
      Seq.Slice(t).GradTimeDelay(n);
    Seq.Slice(t).GradDephase(n).Time(:,Seq.Slice(t).UseAtRepetitionTimeDephase) = ...
      Seq.Slice(t).CenterOfDephase + [ -Seq.Slice(t).GradDephaseLength/2; ...
                                       -Seq.Slice(t).GradDephaseLength/2+Seq.Slice(t).tRamp; ...
                                       +Seq.Slice(t).GradDephaseLength/2-Seq.Slice(t).tRamp; ...
                                       +Seq.Slice(t).GradDephaseLength/2]*tempsize - ...
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
