function [HW, Seq, AQ, TX, Grad] = prepare_DigitalIO_Trigger(HW, Seq, AQ, TX, Grad)
%% Prepare function that switches the Digital output for a trigger
%
%   [HW, Seq, AQ, TX, Grad] = prepare_DigitalIO_Trigger(HW, Seq, AQ, TX, Grad)
%
% This prepare function adds digital output signals to a Spin Echo or FLASH
% measurement that can be used to trigger an external device synchronized to the
% MR pulse sequence.
%
%
% INPUT:
%
%   The input parameter Seq might contain a field "triggerSettings" which is a
%   structure with the following (optional) fields. All of these fields may be
%   vectors to define multiple triggers:
%
%     digiOutChannel
%         Digital output channel for the trigger. (Default: 1)
%
%     onTime
%         Time in seconds from the center of the excitation pulse to when the
%         trigger is switched on. (Default: -3)
%
%     keepOn
%         Boolean value to select whether the trigger should be kept on during
%         the entire pulse sequence (true) or whether the trigger should be
%         switched off before the actual pulse sequence starts (false).
%         (Default: true)
%
%     alwaysOn
%         Only applies if keepOn is true. If set to false, switch trigger off
%         between echo train. If true, keep trigger on during and in between all
%         echo trains. (Default: false)
%
%     UseAtRepetitionTime
%         Vector with the tReps in which the trigger signal will be added.
%         (Default: Seq.P90tReps or, if that field doesn't exist, 1)
%
%     referenceTime
%         Reference time for the end of the trigger signal. The trigger pulses
%         are added to the tReps listed in UseAtRepetitionTime.
%         (Default: Seq.tEcho/2 corresponding to the center of the excitation
%         pulse in sequence_Spin_Echo)
%
%     offTime
%         Time in seconds before "referenceTime" (i.e., the center of the
%         excitation pulse) if keepOn is false, or the center of the last rf
%         pulse if keepOn is true when the trigger is switched off.
%
%     Device
%         If multiple MMRT devices are connected, select the device that should
%         emit the trigger signal. (Default: 1)
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2022-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input parameters
if isemptyfield(Seq, 'triggerSettings')
  Seq.triggerSettings = struct();
end

if isemptyfield(Seq.triggerSettings, 'digiOutChannel')
  Seq.triggerSettings.digiOutChannel = 1;
end

if isemptyfield(Seq.triggerSettings, 'onTime')
  Seq.triggerSettings.onTime = -3;
end
if isemptyfield(Seq.triggerSettings, 'keepOn')
  Seq.triggerSettings.keepOn = true;
end
if isemptyfield(Seq.triggerSettings, 'alwaysOn')
  Seq.triggerSettings.alwaysOn = false;
end
if isemptyfield(Seq.triggerSettings, 'offTime')
  if Seq.triggerSettings.keepOn
    Seq.triggerSettings.offTime = 1;
  else
    Seq.triggerSettings.offTime = -1;
  end
end

if isemptyfield(Seq.triggerSettings, 'UseAtRepetitionTime')
  if isemptyfield(Seq, 'P90tReps')
    Seq.triggerSettings.UseAtRepetitionTime = 1;
  else
    % use the tReps with excitation pulses in sequence_Spin_Echo by default
    Seq.triggerSettings.UseAtRepetitionTime = Seq.P90tReps;
  end
end
if isemptyfield(Seq.triggerSettings, 'referenceTime')
  if isemptyfield(Seq, 'tEcho')
    Seq.triggerSettings.referenceTime = 0;
  else
    % use center of excitation pulse in sequence_Spin_Echo by default
    Seq.triggerSettings.referenceTime = Seq.tEcho/2;
  end
end

if isemptyfield(Seq.triggerSettings, 'Device')
  Seq.triggerSettings.Device = 1;
end


%% expand for multiple triggers
triggerFieldnames = fieldnames(Seq.triggerSettings);
numelTriggerFields = cellfun(@(x) numel(Seq.triggerSettings.(x)), triggerFieldnames);
numTriggers = max(numelTriggerFields);

if any((numelTriggerFields ~= numTriggers) & (numelTriggerFields ~= 1))
  error('PD:prepare_DigitalIO_Trigger', 'Number of settings for trigger are inconsistent');
end

for iField = 1:numel(triggerFieldnames)
  if numelTriggerFields(iField) ~= numTriggers
    Seq.triggerSettings.(triggerFieldnames{iField}) = ...
      repmat(Seq.triggerSettings.(triggerFieldnames{iField}), 1, numTriggers);
  end
end


%% create signal to digital output for each trigger
[DigiOutTrigger(1:numTriggers).SetTime] = deal(NaN(1+any(~Seq.triggerSettings.keepOn), size(Seq.tRep,2)));
[DigiOutTrigger(1:numTriggers).SetValue] = deal(NaN(size(DigiOutTrigger(1).SetTime)));
for iTrigger = 1:numTriggers
  if Seq.triggerSettings.keepOn(iTrigger)
    if ~isemptyfield(Seq, {'AQSlice', 'EPIFactor'})  ... % Check for a field that doesn't exist for Spin Echo measurements.
        || Seq.triggerSettings.alwaysOn(iTrigger)
      % For FLASH or alwaysOn
      % keep trigger on during entire pulse program
      DigiOutTrigger(iTrigger).SetTime(1,1) = Seq.triggerSettings.onTime(iTrigger) + Seq.triggerSettings.referenceTime(iTrigger);
      DigiOutTrigger(iTrigger).SetValue(1,1) = 2^(Seq.triggerSettings.digiOutChannel(iTrigger)-1);
      if isscalar(Seq.triggerSettings.UseAtRepetitionTime) && ~Seq.triggerSettings.keepOn(iTrigger)
        DigiOutTrigger(iTrigger).SetTime(2,1) = Seq.triggerSettings.offTime(iTrigger);
        DigiOutTrigger(iTrigger).SetValue(2,1) = 0;
      else
        DigiOutTrigger(iTrigger).SetTime(1,end) = Seq.triggerSettings.offTime(iTrigger);
        DigiOutTrigger(iTrigger).SetValue(1,end) = 0;
      end
    else
      % For Spin Echo without alwaysOn
      % switch on before echo trains
      DigiOutTrigger(iTrigger).SetTime(1,Seq.triggerSettings.UseAtRepetitionTime) = Seq.triggerSettings.onTime(iTrigger) + Seq.triggerSettings.referenceTime(iTrigger);
      DigiOutTrigger(iTrigger).SetValue(1,Seq.triggerSettings.UseAtRepetitionTime) = 2^(Seq.triggerSettings.digiOutChannel(iTrigger)-1);
      if isscalar(Seq.triggerSettings.UseAtRepetitionTime) && ~Seq.triggerSettings.keepOn(iTrigger)
        % switch off after echo trains
        DigiOutTrigger(iTrigger).SetTime(2,Seq.triggerSettings.UseAtRepetitionTime) = Seq.triggerSettings.offTime(iTrigger);
        DigiOutTrigger(iTrigger).SetValue(2,Seq.triggerSettings.UseAtRepetitionTime) = 0;
      else
        % switch off after echo trains
        DigiOutTrigger(iTrigger).SetTime(1,[Seq.triggerSettings.UseAtRepetitionTime(2:end)-1,end]) = Seq.triggerSettings.offTime(iTrigger);
        DigiOutTrigger(iTrigger).SetValue(1,[Seq.triggerSettings.UseAtRepetitionTime(2:end)-1,end]) = 0;
      end
    end
    % extend last tRep to make room for switch off signal
    Seq.tRep(end) = max(Seq.tRep(end), DigiOutTrigger(iTrigger).SetTime(1,end)+0.1e-3);
  else
    % switch trigger on before (initial) excitation
    if isemptyfield(Seq, {'AQSlice', 'EPIFactor'})  % Check for a field that doesn't exist for Spin Echo measurements.
      % Spin Echo
      % Add trigger before each turbo block.
      DigiOutTrigger(iTrigger).SetTime(1,Seq.triggerSettings.UseAtRepetitionTime) = Seq.triggerSettings.onTime(iTrigger) + Seq.triggerSettings.referenceTime(iTrigger);
      DigiOutTrigger(iTrigger).SetTime(2,Seq.triggerSettings.UseAtRepetitionTime) = Seq.triggerSettings.offTime(iTrigger) + Seq.triggerSettings.referenceTime(iTrigger);
      DigiOutTrigger(iTrigger).SetValue(1,Seq.triggerSettings.UseAtRepetitionTime) = 2^(Seq.triggerSettings.digiOutChannel(iTrigger)-1);
      DigiOutTrigger(iTrigger).SetValue(2,Seq.triggerSettings.UseAtRepetitionTime) = 0;
    else
      % FLASH
      % FIXME: Is this sensible?
      DigiOutTrigger(iTrigger).SetTime(1,1) = Seq.triggerSettings.onTime(iTrigger) + Seq.triggerSettings.referenceTime(iTrigger);
      DigiOutTrigger(iTrigger).SetTime(2,1) = Seq.triggerSettings.offTime(iTrigger) + Seq.triggerSettings.referenceTime(iTrigger);
      DigiOutTrigger(iTrigger).SetValue(1,1) = 2^(Seq.triggerSettings.digiOutChannel(iTrigger)-1);
      DigiOutTrigger(iTrigger).SetValue(2,1) = 0;
    end
  end
end


%% combine triggers into one output
DigiOutTriggerAll = DigiOutTrigger(1);
for iTrigger = 2:numTriggers
  DigiOutTriggerAll = add_DigitalIO(DigiOutTriggerAll, DigiOutTrigger(iTrigger));
end


%% add to digital output
if isemptyfield(Seq, 'DigitalIO')
  Seq.DigitalIO(Seq.triggerSettings.Device) = DigiOutTriggerAll;
else
  Seq.DigitalIO(Seq.triggerSettings.Device) = ...
    add_DigitalIO(Seq.DigitalIO(Seq.triggerSettings.Device), DigiOutTriggerAll);
end

end
