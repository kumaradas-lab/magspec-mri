function AQ = add_AQ(AQ1, AQ2)
%% Add two AQ structures
%
%   AQ = add_AQ(AQ1, AQ2)
%
% The acquisition windows of both input structures are added such that they lie
% within the corresponding tRep of the output structure. Both input AQs can be
% vectors (= acquisition channels) in which case fields of the corresponding
% channels are added.
%
%
% INPUT:
%
% The data in the following fields of AQ1 and AQ2 are added together (i.e.
% concatenated and sorted by Start):
%
%   Start
%       Matrix (time steps x tReps) with the start times of the acquisition
%       windows in seconds.
%
%   fSample
%       Matrix (time steps x tReps) with the sample frequencies in of the
%       acquisition windows in Hz.
%
%   Frequency
%       Matrix (time steps x tReps) with the acquisition frequencies (for the
%       mixer) in Hz.
%
%   Phase
%       Matrix (time steps x tReps) with the phases of the acquisition windows
%       in degrees.
%
%   nSamples
%       Matrix (time steps x tReps) with the number of samples in the
%       acquisition windows.
%
% Optional fields are:
%
%   ResetPhases
%       Vector (1 x tRep) with Boolean values that indicate whether the DDS
%       should be reset at the beginning of the respective tRep. That
%       essentially re-references all phase relationships.
%
%   GetData
%       Vector (1 x tReps) with Boolean values that indicate whether a packet
%       end signal is required at the end of the tRep. If any of the input AQs
%       sets this to 1, the resulting AQ also sets this to 1 for the
%       corresponding tRep (i.e. an "or" operation).
%
%   GetDataAtEnd
%       Scalar Boolean value that indicates whether a packet end signal is
%       requested at the end of the sequence. If any of the input AQs set this
%       to 0, this value is set to 0 in the resulting AQ as well (i.e. an "and"
%       operation).
%
%
% OUTPUT:
%
% An AQ structure where the above fields are "merged".
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% check input for optional fields
hasResetPhases = isfield(AQ1, 'ResetPhases') || isfield(AQ2, 'ResetPhases');
hasGetData = isfield(AQ1, 'GetData') || isfield(AQ2, 'GetData');
hasGetDataAtEnd = isfield(AQ1, 'GetDataAtEnd') || isfield(AQ2, 'GetDataAtEnd');
if hasResetPhases
  if ~isempty(AQ1) && ~isfield(AQ1, 'ResetPhases'), [AQ1.ResetPhases] = deal(0); end
  if ~isempty(AQ2) && ~isfield(AQ2, 'ResetPhases'), [AQ2.ResetPhases] = deal(0); end
end
if hasGetData
  if ~isempty(AQ1) && ~isfield(AQ1, 'GetData'), [AQ1.GetData] = deal(0); end
  if ~isempty(AQ2) && ~isfield(AQ2, 'GetData'), [AQ2.GetData] = deal(0); end
end
if hasGetDataAtEnd
  if ~isempty(AQ1) && ~isfield(AQ1, 'GetDataAtEnd'), [AQ1.GetDataAtEnd] = deal(1); end
  if ~isempty(AQ2) && ~isfield(AQ2, 'GetDataAtEnd'), [AQ2.GetDataAtEnd] = deal(1); end
end

if ~isempty(AQ1) && ~isfield(AQ1, 'SamplingFactor'),  [AQ1.SamplingFactor] = deal(1);  end
if ~isempty(AQ2) && ~isfield(AQ2, 'SamplingFactor'),  [AQ2.SamplingFactor] = deal(1);  end


%% loop over acquisition channels
for t = 1:min(numel(AQ1),numel(AQ2))
  sizeAQ1 = size(AQ1(t).Start);
  sizeAQ2 = size(AQ2(t).Start);

  AQ(t).Start = nan(sizeAQ1(1)+sizeAQ2(1), max(sizeAQ1(2), sizeAQ2(2)));
  AQ(t).fSample = AQ(t).Start;
  AQ(t).Frequency = AQ(t).Start;
  AQ(t).Phase = AQ(t).Start;
  AQ(t).nSamples = AQ(t).Start;
  % AQ(t).Gain=AQ.Start;
  if hasResetPhases, AQ(t).ResetPhases = AQ1.ResetPhases; end
  if hasGetData, AQ(t).GetData = AQ1.GetData; end
  if hasGetDataAtEnd, AQ(t).GetDataAtEnd = AQ1.GetDataAtEnd; end

  % Expand structure with one tRep to number of tReps in other structure
  if size(AQ1(t).Start, 2) == 1
    ones2AQ2 = ones(1, max(sizeAQ2(2),1));
    AQ1(t).Start = AQ1(t).Start * ones2AQ2;
    AQ1(t).fSample = AQ1(t).fSample * ones2AQ2;
    AQ1(t).Frequency = AQ1(t).Frequency * ones2AQ2;
    AQ1(t).Phase = AQ1(t).Phase * ones2AQ2;
    AQ1(t).nSamples = AQ1(t).nSamples * ones2AQ2;
    AQ1(t).SamplingFactor = AQ1(t).SamplingFactor * ones2AQ2;
    if hasResetPhases, AQ1(t).ResetPhases = AQ1(t).ResetPhases * ones2AQ2; end
    if hasGetData, AQ1(t).GetData = AQ1(t).GetData * ones2AQ2; end
    % AQ1(t).Gain=AQ1(t).Gain*ones(1,max(size(AQ2(t).Start,2),1));
  end
  if size(AQ2(t).Start, 2) == 1
    ones2AQ1 = ones(1, max(sizeAQ1(2),1));
    AQ2(t).Start = AQ2(t).Start * ones2AQ1;
    AQ2(t).fSample = AQ2(t).fSample * ones2AQ1;
    AQ2(t).Frequency = AQ2(t).Frequency * ones2AQ1;
    AQ2(t).Phase = AQ2(t).Phase * ones2AQ1;
    AQ2(t).nSamples = AQ2(t).nSamples * ones2AQ1;
    AQ2(t).SamplingFactor = AQ2(t).SamplingFactor * ones2AQ1;
    if hasResetPhases, AQ2(t).ResetPhases = AQ2(t).ResetPhases * ones2AQ1; end
    if hasGetData, AQ2(t).GetData = AQ2(t).GetData * ones2AQ1; end
    % AQ2(t).Gain=AQ2(t).Gain*ones2AQ1;
  end

  % Expand values that are set only once per tRep to each step in the tRep
  ones1AQ1 = ones(sizeAQ1(1), 1);
  if size(AQ1(t).fSample,1)==1, AQ1(t).fSample = ones1AQ1 * AQ1(t).fSample; end
  if size(AQ1(t).Frequency,1)==1, AQ1(t).Frequency = ones1AQ1 * AQ1(t).Frequency; end
  if size(AQ1(t).Phase,1)==1, AQ1(t).Phase = ones1AQ1 * AQ1(t).Phase; end
  if size(AQ1(t).nSamples,1)==1, AQ1(t).nSamples = ones1AQ1 * AQ1(t).nSamples; end
  if size(AQ1(t).SamplingFactor,1)==1, AQ1(t).SamplingFactor = ones1AQ1 * AQ1(t).SamplingFactor; end
  % if size(AQ1(t).Gain,1)==1, AQ1(t).Gain=ones1AQ1*AQ1(t).Gain;end         % 1D

  ones1AQ2 = ones(sizeAQ2(1), 1);
  if size(AQ2(t).fSample,1)==1, AQ2(t).fSample = ones1AQ2 * AQ2(t).fSample; end
  if size(AQ2(t).Frequency,1)==1, AQ2(t).Frequency = ones1AQ2 * AQ2(t).Frequency; end
  if size(AQ2(t).Phase,1)==1, AQ2(t).Phase = ones1AQ2 * AQ2(t).Phase; end
  if size(AQ2(t).nSamples,1)==1, AQ2(t).nSamples = ones1AQ2 * AQ2(t).nSamples; end
  if size(AQ2(t).SamplingFactor,1)==1, AQ2(t).SamplingFactor = ones1AQ2 * AQ2(t).SamplingFactor; end
  % if size(AQ2(t).Gain,1)==1, AQ2(t).Gain=ones1AQ2*AQ2(t).Gain;end         % 1D

  % Concatenate fields in both structures
  AQ(t).Start = [AQ1(t).Start; AQ2(t).Start];
  AQ(t).fSample = [AQ1(t).fSample; AQ2(t).fSample];
  AQ(t).Frequency = [AQ1(t).Frequency; AQ2(t).Frequency];
  AQ(t).Phase = [AQ1(t).Phase; AQ2(t).Phase];
  AQ(t).nSamples = [AQ1(t).nSamples; AQ2(t).nSamples];
  AQ(t).SamplingFactor = [AQ1(t).SamplingFactor; AQ2(t).SamplingFactor];
  if hasResetPhases
    % GetData 1 overrules GetData 0 (the default)
    AQ(t).ResetPhases = bsxfun(@plus, AQ1(t).ResetPhases, AQ2(t).ResetPhases);
    AQ(t).ResetPhases(AQ(t).ResetPhases > 0) = 1; % FIXME: @or doesn't work with NaNs
  end
  if hasGetData
    % GetData 1 overrules GetData 0 (the default)
    AQ(t).GetData = bsxfun(@plus, AQ1(t).GetData, AQ2(t).GetData);
    AQ(t).GetData(AQ(t).GetData > 0) = 1; % FIXME: @or doesn't work with NaNs
  end
  if hasGetDataAtEnd
    % GetDataAtEnd 0 overrules GetDataAtEnd 1 (the default)
    AQ(t).GetDataAtEnd = AQ1(t).GetDataAtEnd * AQ2(t).GetDataAtEnd;
  end
  % AQ(t).Gain=[AQ1(t).Gain;AQ2(t).Gain];

  % Sort according to start times
  [AQ(t).Start, IX] = sort(AQ(t).Start,1);
  IX = IX + ones(size(IX,1),1)*(cumsum(size(IX,1)*ones(1,size(IX,2)),2)-size(IX,1));
  AQ(t).fSample = AQ(t).fSample(IX);
  AQ(t).Frequency = AQ(t).Frequency(IX);
  AQ(t).Phase = AQ(t).Phase(IX);
  AQ(t).nSamples = AQ(t).nSamples(IX);
  AQ(t).SamplingFactor = AQ(t).SamplingFactor(IX);
  % AQ(t).Gain = AQ(t).Gain(IX);

  % Remove time steps that are "empty" (NaN) in all tReps
  [row,col] = find(sum(~isnan(AQ(t).Start),2), 1, 'last');
  firstEmpty = max(row)+1;
  AQ(t).Start(firstEmpty:end,:) = [];
  AQ(t).fSample(firstEmpty:end,:) = [];
  AQ(t).Frequency(firstEmpty:end,:) = [];
  AQ(t).Phase(firstEmpty:end,:) = [];
  AQ(t).nSamples(firstEmpty:end,:) = [];
  AQ(t).SamplingFactor(firstEmpty:end,:) = [];
  % AQ(t).Gain(firstEmpty:end,:) = [];
end

% Add elements of AQ with more AQ channels
for t = (numel(AQ2)+1):numel(AQ1), AQ(t) = AQ1(t); end  % AQ1 larger than AQ2
for t = (numel(AQ1)+1):numel(AQ2), AQ(t) = AQ2(t); end  % AQ2 larger than AQ1

end
