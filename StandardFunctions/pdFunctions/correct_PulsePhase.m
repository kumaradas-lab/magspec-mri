function TXShape = correct_PulsePhase(TXShape, HW, iDevice, fDDS)
%% Correct phase of a slice selective rf pulses with offset frequency
%
%   TXShape = correct_PulsePhase(TXShape, HW, iDevice, fDDS)
%
% The phase of the (shaped) rf pulses in the structure TXShape is corrected such
% that the "effective phase offset" at the center of the pulse is as set in
% TXShape.Phase.
% Additionally, the "dummy pulses" that are added before and after the rf pulse
% are used to shift the phase of the DDS such that the effective phase after the
% pulse is as if no off-frequency pulse was emitted. It is assumed that the
% reference frequency is HW.fLarmor (or HW.fLarmorX for the second frequency of
% dual frequency pulses).
%
%
% INPUT:
%
%   TXShape
%         TX structure which contains only the pulse for which the phase should
%         be corrected.
%
%   HW
%         HW object (or structure)
%
%   iDevice
%         Index of the used device in case multiple MRT devices are connected.
%         (Default: 1)
%
%   fDDS
%         Reference frequency of the DDS. (Default: HW.fLarmor)
%
%
% OUTPUT:
%
%   TXShape
%         TX structure which contains the pulse from the input structure TXShape
%         potentially extended by "dummy pulses" before and after the actual rf
%         pulse that are used to shift the phase of the DDS.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if nargin < 3, iDevice = 1; end
if nargin < 4, fDDS = HW.fLarmor; end

for t = 1:numel(TXShape)
  FrequencyGridTX = HW.TX(iDevice).fSample / 2^HW.TX(iDevice).DdsPicBits;
  TXShape(t).Frequency = round(TXShape(t).Frequency./FrequencyGridTX) .* FrequencyGridTX;
  fDDS = round(fDDS./FrequencyGridTX) .* FrequencyGridTX;
  % apply settings to all time steps
  if size(TXShape(t).Frequency, 1) == 1
    TXShape(t).Frequency = repmat(TXShape(t).Frequency, size(TXShape(t).Start,1), 1);
  end
  if size(TXShape(t).Duration, 1) == 1
    TXShape(t).Duration = repmat(TXShape(t).Duration, size(TXShape(t).Start,1), 1);
  end
  if size(TXShape(t).Phase, 1) == 1
    TXShape(t).Phase = repmat(TXShape(t).Phase, size(TXShape(t).Start,1), 1);
  end

  if all(TXShape(t).Frequency(:) == fDDS)
    % no correction necessary
    if isemptyfield(TXShape(t), 'Start')
      TXShape(t).Center = 0;
      TXShape(t).CenterOffset = 0;
    else
      if isemptyfield(TXShape(t), 'Center')
        TXShape(t).Center = TXShape(t).Start(1,:) + sum(TXShape(t).Duration,1)/2;
      end
      if isemptyfield(TXShape(t), 'CenterOffset')
        TXShape(t).CenterOffset = TXShape(t).Center - (TXShape(t).Start(1,:) + sum(TXShape(t).Duration,1)/2);
      end
    end
    continue;
  end

  tGrid = 1/HW.TX(iDevice).fSample;
  tMinPulse = tGrid*2;  % duration of shortest possible pulse
  % PhaseGridTX = 360 / 2^HW.TX(iDevice).DdsPofBits;

  % txStart = TXShape(t).Start;
  % txDuration = TXShape(t).Duration;
  TXShape(t).Start = round(TXShape(t).Start/tGrid)*tGrid;
  TXShape(t).Duration = round(TXShape(t).Duration/tGrid)*tGrid;
  TXShape(t).Duration(1:end-1,:) = TXShape(t).Start(2:end,:) - TXShape(t).Start(1:end-1,:);

  % if any(abs(txStart-TXShape(t).Start) > 0.99*0.5/HW.TX(iDevice).fSample)
  %   if any(TXShape(t).Frequency(~isnan(TXShape(t).Start))~=round(HW.fLarmor/FrequencyGridTX)*FrequencyGridTX)
  %     error(' TX frequency is not HW.fLarmor and start times are not on 1/HW.TX.fSample Grid')
  %   end
  % end
  % if any(abs(txDuration-TXShape(t).Duration) > 0.001*1/HW.TX(iDevice).fSample)
  %   if any(TXShape(t).Frequency(~isnan(TXShape(t).Start))~=round(HW.fLarmor/FrequencyGridTX)*FrequencyGridTX)
  %     error(' TX frequency is not HW.fLarmor and duration is not on 1/HW.TX.fSample Grid')
  %   end
  % end

  if ~all(all((round(TXShape(t).Start(2:end,:)*HW.TX(iDevice).fSample)-round(TXShape(t).Start(1:end-1,:)*HW.TX(iDevice).fSample)) == ...
              round(TXShape(t).Duration(1:end-1,:) * HW.TX(iDevice).fSample)))
    error('TX pulses are not on 1/HW.TX.fSample grid')
  end

  if isemptyfield(TXShape(t), 'Center')
    TXShape(t).Center = TXShape(t).Start(1,:) + sum(TXShape(t).Duration,1)/2;
  end
  % if isemptyfield(TXShape(t), 'CenterOffset'), TXShape(t).CenterOffset = TXShape(t).Center-TXShape(t).Start(1,:) + sum(TXShape(t).Duration,1)/2; end
  TXShape(t).CenterOffset = TXShape(t).Center - (TXShape(t).Start(1,:)+sum(TXShape(t).Duration,1)/2);

  % apply settings to all tReps
  sz2 = size(TXShape(t).Frequency, 2);
  if size(TXShape(t).Duration, 2) == 1
    TXShape(t).Duration = TXShape(t).Duration * ones(1, sz2);
  end
  if size(TXShape(t).Start, 2) == 1
    TXShape(t).Start = TXShape(t).Start * ones(1, sz2);
  end
  if size(TXShape(t).Center, 2) == 1
    TXShape(t).Center = TXShape(t).Center * ones(1, sz2);
  end

  % Find index with element that is at the center of the pulse
  pci = sum(TXShape(t).Start <= repmat(TXShape(t).Center, size(TXShape(t).Start,1), 1), 1);
  if any(pci <= 0), error('Pulse center before pulse'); end

  tDur = TXShape(t).Duration;
  % set all durations after the pulse center to 0
  tDur(TXShape(t).Start>repmat(TXShape(t).Center,size(TXShape(t).Start,1),1)) = 0;
  % change duration of the element at the center to the duration until the center
  % (linear index into 2d matrix)
  tDur((0:size(TXShape(t).Start,2)-1) * size(TXShape(t).Start,1) + pci) = ...
    TXShape(t).Center - TXShape(t).Start((0:size(TXShape(t).Start,2)-1) * size(TXShape(t).Start,1) + pci);
  % Adjusted duration (to the center of the pulse) must not be longer than the
  % original duration.
  if any(tDur((0:size(TXShape(t).Start,2)-1) * size(TXShape(t).Start,1) + pci) > ...
      TXShape(t).Duration((0:size(TXShape(t).Start,2)-1) * size(TXShape(t).Start,1) + pci) + 1e-12)
    error('Pulse center after pulse');
  end

  % Get the frequency offset to set for one DDS clock to compensate the phase
  % offset in the DDS phase accumulator.
  FrequencyOffsetOfCorrectPhasePulse = (mod(0.5 + ...
                                            sum(...
                                                TXShape(t).Duration ...     % rounded duration of pulse segments
                                                .*(TXShape(t).Frequency ... % rounded frequency of pulse segments
                                                   -  fDDS )...             % get offset frequency
                                                , 1)...
                                            , ones(1, sz2))-0.5)...  % get +/- of a half turn (full turns are not subtracted)
                                       ./(tMinPulse);
  PhaseOffsetAtPulseCenter = (mod(0.5 + ...
                                  sum(...
                                      tDur ...                        % rounded duration of pulse segments
                                      .*( TXShape(t).Frequency ...    % rounded frequncy of pulse segments
                                         -  fDDS  ))...               % get rounded offset frequency
                                  ...             % phase offset from start of pulse to center
                                  , ones(1, sz2))-0.5) ...  % get +/- of a half turn (full turns are not subtracted)
                             * 360;                   % get the phase in deg
  % PhaseOffsetAtPulseCenter = (mod(0.5+...
  %                                 sum(...
  %                                     TXShape(t).Duration  ...                       % rounded duration of pulse segments
  %                                     .*( TXShape(t).Frequency ...    % rounded frequncy of pulse segments
  %                                        -  fLarmor  ))...   % get rounded offset frequency
  %                                     ./2 ...             % phase offset to the center
  %                                ,ones(1,s2))-0.5)...             % get +/- of a half turn (full turns are not subtracted)
  %                             * 360;                   % get the phase in deg

  % Prepend "infinitely short" pulse with HW.fLarmor
  % Append "infinitely short" pulse that compensates the phase offset due to the
  % offset frequency followed by another "infinitely short" pulse that resets to
  % HW.fLarmor
  % All of these "infinitely short" pulses have zero amplitude and are only for
  % the DDS phase accumulator.
  TXShape(t).Start = [ TXShape(t).Start(1,:)-tMinPulse; ...
                       TXShape(t).Start; ...
                       TXShape(t).Start(end,:)+TXShape(t).Duration(end,:); ...
                       TXShape(t).Start(end,:)+TXShape(t).Duration(end,:)+tMinPulse]+1.111e-12;  % round up a 0.5*tGrid

  TXShape(t).Frequency = [ fDDS*ones(1, sz2); ...
                           TXShape(t).Frequency; ...
                           fDDS-FrequencyOffsetOfCorrectPhasePulse; ...
                           fDDS*ones(1, sz2)];

  TXShape(t).Duration = [ tMinPulse*ones(1, sz2); ...
                          TXShape(t).Duration; ...
                          tMinPulse*ones(1, sz2); ...
                          tMinPulse*ones(1, sz2)];

  if size(TXShape(t).Phase, 2) == 1
    TXShape(t).Phase = TXShape(t).Phase * ones(1, size(PhaseOffsetAtPulseCenter,2));
  end

  TXShape(t).Phase = [ TXShape(t).Phase(1,:); ...
                       TXShape(t).Phase; ...
                       TXShape(t).Phase(end,:); ...
                       TXShape(t).Phase(end,:)] - ...
                     ones(size(TXShape(t).Phase,1)+3, 1) * PhaseOffsetAtPulseCenter;

  TXShape(t).Amplitude = [ zeros(1, sz2);...
                           TXShape(t).Amplitude;...
                           zeros(1, sz2);...
                           zeros(1, sz2)];

end

end
