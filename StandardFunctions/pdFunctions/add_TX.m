function TX = add_TX(TX1, TX2)
%% Merge rf pulses from two TX structures into one structure
%
%   TX = add_TX(TX1, TX2)
%
% The rf pulses of both input structures are added such that they lie within the
% corresponding tRep of the output structure. Pulses in TX1 and TX2 must not
% overlap.
% Both input TX structures can be vectors (= rf channels and/or devices). In any
% case, the input structures must be sorted equivalently such that the channel
% and device of the respective elements match.
%
% INPUT:
%
% The data in the following fields of TX1 and TX2 are added together (i.e.
% concatenated and sorted by Start):
%
%   Start
%       Matrix (time steps x tReps) with the start times of the rf (block)
%       pulses in seconds.
%
%   Duration
%       Matrix (time steps x tReps) with the duration of the rf (block) pulses
%       in seconds.
%
%   Amplitude
%       Matrix (time steps x tReps) with the amplitude of the rf (block) pulses
%       in units as given by HW.TX.PaUout2Amplitude.
%
%   Frequency
%       Matrix (time steps x tReps) with the frequency of the rf (block) pulses
%       in Hz.
%
%   Phase
%       Matrix (time steps x tReps) with the phase of the rf (block) pulses in
%       degrees.
%
% Optional fields are:
%
%   Channel
%       Scalar with the channel number. If the channels in TX1 and TX2 don't
%       match, a warning is issued and the channel from TX2 takes precedence.
%
%   Device
%       Scalar with the device number. If the device number in TX1 and TX2
%       doesn't match, a warning is issued and the device number from TX2 takes
%       precedence.
%
% All other input fields are ignored and discarded.
%
%
% OUTPUT:
%
% A TX structure where the above fields are "merged".
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2013-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% merge structure elements that are present in both inputs
for t = 1:min(numel(TX1), numel(TX2))
  TX(t).Start = NaN(size(TX1(t).Start,1) + size(TX2(t).Start,1), ...
    max(size(TX1(t).Start,2), size(TX2(t).Start,2)));
  TX(t).Duration = TX(t).Start;
  TX(t).Amplitude = TX(t).Start;
  TX(t).Frequency = TX(t).Start;
  TX(t).Phase = TX(t).Start;

  % FIXME: Support "mixed" output channels.
  if isfield(TX1(t), 'Channel')  TX(t).Channel = TX1(t).Channel;  end
  if isfield(TX2(t), 'Channel')
    if isfield(TX(t), 'Channel') && TX(t).Channel ~= TX2(t).Channel
      warning('PD:add_TX:ChannelDisagree', ...
        'The channel of merged TX structures must agree.');
    end
    TX(t).Channel = TX2(t).Channel;
  end
  if isfield(TX1(t), 'Device'),  TX(t).Device = TX1(t).Device;  end
  if isfield(TX2(t), 'Device')
    if isfield(TX(t), 'Device') && TX(t).Device ~= TX2(t).Device
      warning('PD:add_TX:DeviceDisagree', ...
        'The device of merged TX structures must agree.');
    end
    TX(t).Device = TX2(t).Device;
  end

  % FIXME: What about the fields "BlankOffset", "BlankPostset", "Repeat"?

  % expand to all tReps
  if size(TX1(t).Start, 2) == 1
    TX1(t).Start = TX1(t).Start * ones(1, max(size(TX2(t).Start,2),1));
    TX1(t).Duration = TX1(t).Duration * ones(1, max(size(TX2(t).Start,2),1));
    TX1(t).Amplitude = TX1(t).Amplitude * ones(1, max(size(TX2(t).Start,2),1));
    TX1(t).Frequency = TX1(t).Frequency * ones(1, max(size(TX2(t).Start,2),1));
    TX1(t).Phase = TX1(t).Phase * ones(1, max(size(TX2(t).Start,2),1));
  end
  if size(TX2(t).Start, 2) == 1
    TX2(t).Start = TX2(t).Start * ones(1, max(size(TX1(t).Start,2),1));
    TX2(t).Duration = TX2(t).Duration * ones(1, max(size(TX1(t).Start,2),1));
    TX2(t).Amplitude = TX2(t).Amplitude * ones(1, max(size(TX1(t).Start,2),1));
    TX2(t).Frequency = TX2(t).Frequency * ones(1, max(size(TX1(t).Start,2),1));
    TX2(t).Phase = TX2(t).Phase * ones(1, max(size(TX1(t).Start,2),1));
  end

  % expand all fields to the number of rows of TX.Start
  if size(TX1(t).Duration, 1) == 1
    TX1(t).Duration = ones(size(TX1(t).Start,1),1) * TX1(t).Duration;
  end
  if size(TX1(t).Amplitude, 1) == 1
    TX1(t).Amplitude = ones(size(TX1(t).Start,1),1) * TX1(t).Amplitude;
  end
  if size(TX1(t).Frequency, 1) == 1
    TX1(t).Frequency = ones(size(TX1(t).Start,1),1) * TX1(t).Frequency;
  end
  if size(TX1(t).Phase, 1) == 1
    TX1(t).Phase = ones(size(TX1(t).Start,1),1) * TX1(t).Phase;
  end

  if size(TX2(t).Duration, 1) == 1
    TX2(t).Duration = ones(size(TX2(t).Start,1),1) * TX2(t).Duration;
  end
  if size(TX2(t).Amplitude, 1) == 1
    TX2(t).Amplitude = ones(size(TX2(t).Start,1),1) * TX2(t).Amplitude;
  end
  if size(TX2(t).Frequency, 1) == 1
    TX2(t).Frequency = ones(size(TX2(t).Start,1),1) * TX2(t).Frequency;
  end
  if size(TX2(t).Phase, 1) == 1
    TX2(t).Phase = ones(size(TX2(t).Start,1),1) * TX2(t).Phase;
  end

  % concatenate fields from both input structures
  TX(t).Start = [TX1(t).Start; TX2(t).Start];
  TX(t).Duration = [TX1(t).Duration; TX2(t).Duration];
  TX(t).Amplitude = [TX1(t).Amplitude; TX2(t).Amplitude];
  TX(t).Frequency = [TX1(t).Frequency; TX2(t).Frequency];
  TX(t).Phase = [TX1(t).Phase; TX2(t).Phase];

  % sort all fields by start time
  [TX(t).Start, IX] = sort(TX(t).Start, 1);
  IX = IX + ones(size(IX,1),1) * (cumsum(size(IX,1)*ones(1,size(IX,2)),2) - size(IX,1));
  TX(t).Duration = TX(t).Duration(IX);
  TX(t).Amplitude = TX(t).Amplitude(IX);
  TX(t).Frequency = TX(t).Frequency(IX);
  TX(t).Phase = TX(t).Phase(IX);

  % remove trailing rows with NaN
  row = max(sum(~isnan(TX(t).Start),1));
  TX(t).Start = TX(t).Start(1:row,:);
  TX(t).Duration = TX(t).Duration(1:row,:);
  TX(t).Amplitude = TX(t).Amplitude(1:row,:);
  TX(t).Frequency = TX(t).Frequency(1:row,:);
  TX(t).Phase = TX(t).Phase(1:row,:);

end


%% add structure elements that only occur in one of the two input arguments
for t = (numel(TX2)+1):numel(TX1)
  if isfield(TX1(t), 'Channel'),  TX(t).Channel = TX1(t).Channel;  end
  if isfield(TX1(t), 'Device'),  TX(t).Device = TX1(t).Device;  end
  if isfield(TX1(t), 'Start'),  TX(t).Start = TX1(t).Start;  end
  if isfield(TX1(t), 'Frequency'),  TX(t).Frequency = TX1(t).Frequency;  end
  if isfield(TX1(t), 'Duration'),  TX(t).Duration = TX1(t).Duration;  end
  if isfield(TX1(t), 'Amplitude'),  TX(t).Amplitude = TX1(t).Amplitude;  end
  if isfield(TX1(t), 'Phase'),  TX(t).Phase = TX1(t).Phase;  end
  if isfield(TX1(t), 'BlankOffset'),  TX(t).BlankOffset = TX1(t).BlankOffset;  end
  if isfield(TX1(t), 'BlankPostset'),  TX(t).BlankPostset = TX1(t).BlankPostset;  end
  if isfield(TX1(t), 'Repeat'),  TX(t).Repeat = TX1(t).Repeat;  end
end

for t = (numel(TX1)+1):numel(TX2)
  if isfield(TX2(t), 'Channel'),  TX(t).Channel = TX2(t).Channel;  end
  if isfield(TX2(t), 'Device'),  TX(t).Device = TX2(t).Device;  end
  if isfield(TX2(t), 'Start'),  TX(t).Start = TX2(t).Start;  end
  if isfield(TX2(t), 'Frequency'),  TX(t).Frequency = TX2(t).Frequency;  end
  if isfield(TX2(t), 'Duration'),  TX(t).Duration = TX2(t).Duration;  end
  if isfield(TX2(t), 'Amplitude'),  TX(t).Amplitude = TX2(t).Amplitude;  end
  if isfield(TX2(t), 'Phase'),  TX(t).Phase = TX2(t).Phase;  end
  if isfield(TX2(t), 'BlankOffset'),  TX(t).BlankOffset = TX2(t).BlankOffset;  end
  if isfield(TX2(t), 'BlankPostset'),  TX(t).BlankPostset = TX2(t).BlankPostset;  end
  if isfield(TX2(t), 'Repeat'),  TX(t).Repeat = TX2(t).Repeat;  end
end


end
