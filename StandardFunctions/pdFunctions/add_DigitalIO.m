function IO = add_DigitalIO(IO1, IO2)
%% Add two digital IO structures
%
%   IO = add_DigitalIO(IO1, IO2)
%
% The digital IO (output) channels of both input structures are added such that
% they lie within the corresponding tRep of the output structure. Both input IOs
% can be vectors (= multiple ports) in which case fields of the corresponding
% ports are added.
% It is assumed that all bits are low at the beginning of the tRep and that the
% settings in IO1 and IO2 are independent.
%
% INPUT:
% The data in the following fields of IO1 and IO2 are merged together (i.e.
% concatenated and sorted by Start):
%   SetTime     Matrix (time steps x tReps) with the times when the digital IO
%               is changed in seconds.
%   SetValue    Matrix (time steps x tReps) with the value that is set at time
%               "SetTime".
%
% OUTPUT:
% An IO structure where the above fields are "merged".
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% FIXME: The "Repeat" fields aren't correctly merged. Probably no longer needed.
%        Handled by .NET.

%% loop over IO ports
for t = 1:min(numel(IO1), numel(IO2))
  if isemptyfield(IO1(t), 'SetTime'),  IO1(t).SetTime = NaN;  end
  if isemptyfield(IO2(t), 'SetTime'),  IO2(t).SetTime = NaN;  end
  if isemptyfield(IO1(t), 'SetValue'),  IO1(t).SetValue = NaN;  end
  if isemptyfield(IO2(t), 'SetValue'),  IO2(t).SetValue = NaN;  end

  IO(t) = IO1(t);
  sizeIO1 = size(IO1(t).SetTime);
  sizeIO2 = size(IO2(t).SetTime);

  % clear values that don't have a corresponding time
  IO1(t).SetValue(isnan(IO1(t).SetTime)) = NaN;
  IO2(t).SetValue(isnan(IO2(t).SetTime)) = NaN;

  if all(isnan(IO1(t).SetValue(:)))
    IO(t).SetTime = IO2(t).SetTime;
    IO(t).SetValue = IO2(t).SetValue;
  elseif ~all(isnan(IO2(t).SetValue(:)))
    % remove time steps that don't change the value and sort in ascending order
    noChange = [false(1, sizeIO1(2)); diff(IO1(t).SetValue, 1, 1) == 0];
    IO1(t).SetTime(noChange) = NaN;
    IO1(t).SetValue(noChange) = NaN;
    [IO1(t).SetTime, idx] = sort(IO1(t).SetTime, 1);
    IO1(t).SetValue(repmat(((1:sizeIO1(2))-1)*sizeIO1(1), sizeIO1(1), 1) + idx) = IO1(t).SetValue;
    noChange = [false(1, sizeIO2(2)); diff(IO2(t).SetValue, 1, 1) == 0];
    IO2(t).SetTime(noChange) = NaN;
    IO2(t).SetValue(noChange) = NaN;
    [IO2(t).SetTime, idx] = sort(IO2(t).SetTime, 1);
    IO2(t).SetValue(repmat(((1:sizeIO2(2))-1)*sizeIO2(1), sizeIO2(1), 1) + idx) = IO2(t).SetValue;

    %   IO(t).SetTime = nan(sizeIO1(1)+sizeIO2(1), max(sizeIO1(2), sizeIO2(2))); %#ok<*AGROW>
    %   IO(t).SetValue = IO(t).SetTime;

    % Expand structure with one tRep to number of tReps in other structure
    if sizeIO1(2) == 1
      ones2IO2 = ones(1, max(sizeIO2(2), 1));
      IO1(t).SetTime = IO1(t).SetTime * ones2IO2;
      IO1(t).SetValue = IO1(t).SetValue * ones2IO2;
    end
    if sizeIO2(2) == 1
      ones2IO1 = ones(1, max(sizeIO1(2), 1));
      IO2(t).SetTime = IO2(t).SetTime * ones2IO1;
      IO2(t).SetValue = IO2(t).SetValue * ones2IO1;
    end

    % Expand values that are set only once per tRep to each step in the tRep
    ones1IO1 = ones(sizeIO1(1), 1);
    % if size(IO1(t).SetTime,1)==1, IO1(t).SetTime = ones1IO1 * IO1(t).SetTime; end
    if size(IO1(t).SetValue,1)==1, IO1(t).SetValue = ones1IO1 * IO1(t).SetValue; end

    ones1IO2 = ones(sizeIO2(1), 1);
    % if size(IO2(t).SetTime,1)==1, IO2(t).SetTime = ones1IO2 * IO2(t).SetTime; end
    if size(IO2(t).SetValue,1)==1, IO2(t).SetValue = ones1IO2 * IO2(t).SetValue; end

    % replace NaNs with first finite preceeding value
    % assume that the state at the beginning is the state at the end
    isNaN1 = isnan(IO1(t).SetValue);
    idxPreceeding1 = reshape(cumsum(~isNaN1(:)), size(isNaN1));
    finiteValues = IO1(t).SetValue(~isNaN1);
    idxPreceeding1(idxPreceeding1 < 1) = numel(finiteValues);
    IO1(t).SetValue(isNaN1) = finiteValues(idxPreceeding1(isNaN1));
    isNaN2 = isnan(IO2(t).SetValue);
    idxPreceeding2 = reshape(cumsum(~isNaN2(:)), size(isNaN2));
    finiteValues = IO2(t).SetValue(~isNaN2);
    idxPreceeding2(idxPreceeding2 < 1) = numel(finiteValues);
    IO2(t).SetValue(isNaN2) = finiteValues(idxPreceeding2(isNaN2));

    % Concatenate fields in both structures
    IO(t).SetTime = [IO1(t).SetTime; IO2(t).SetTime];
    IO(t).SetValue = [IO1(t).SetValue; IO2(t).SetValue];
    isSecond = [zeros(size(IO1(t).SetValue)); ones(size(IO2(t).SetValue))];

    % Sort according to SetTimes
    [IO(t).SetTime, IX] = sort(IO(t).SetTime,1);
    IX = IX + ones(size(IX,1),1)*(cumsum(size(IX,1)*ones(1,size(IX,2)),2)-size(IX,1));
    IO(t).SetValue = IO(t).SetValue(IX);
    isSecond = isSecond(IX);

    % >0: from 1 to 2, <0: from 2 to 1
    switchStruct1 = reshape(diff([0; isSecond(:)], 1, 1), size(isSecond));
    switchStruct2 = reshape(diff([1; isSecond(:)], 1, 1), size(isSecond));

    % bit-wise addition (bitor)
    idx1 = reshape(cumsum(~isSecond(:), 1), size(isSecond));
    idx2 = reshape(cumsum(isSecond(:), 1), size(isSecond));
    addToIO1 = (reshape(cumsum(switchStruct1(:)), size(switchStruct1)) > 0) & (idx1 > 0);
    addToIO2 = (reshape(cumsum(-switchStruct2(:)), size(switchStruct2)) > 0) & (idx2 > 0);
    IO2plus1 = IO(t).SetValue;
    IO1plus2 = IO2plus1;
    if ~isempty(IO2plus1(addToIO1))
      IO_a = IO2plus1(addToIO1);
      IO_b = IO1(t).SetValue(idx1(addToIO1));
      IO2plus1(addToIO1) = bitor(IO_a(:), IO_b(:));
    end
    if ~isempty(IO1plus2(addToIO2))
      IO_a = IO1plus2(addToIO2);
      IO_b = IO2(t).SetValue(idx2(addToIO2));
      IO1plus2(addToIO2) = bitor(IO_a(:), IO_b(:));
    end
    IO(t).SetValue = max(IO2plus1, IO1plus2);
    IO(t).SetValue(isnan(IO(t).SetTime)) = NaN;

    % eliminate concurrent values
    concurrent = [abs(diff(IO(t).SetTime, 1, 1))<eps; false(1, size(IO(t).SetTime, 2))];
    falseLine = false(1, size(IO(t).SetTime, 2));
    selectMin = switchStruct1 < 0 & switchStruct2 < 0 & concurrent;
    elimFirstMin = IO(t).SetValue(selectMin) < IO(t).SetValue([falseLine;selectMin(1:end-1,:)]);
    selectMinIdx = find(selectMin);

    selectMax = concurrent & ~selectMin;
    elimFirstMax = IO(t).SetValue(selectMax) > IO(t).SetValue([falseLine;selectMax(1:end-1,:)]);
    selectMaxIdx = find(selectMax);

    IO(t).SetValue([selectMinIdx+elimFirstMin; selectMaxIdx+elimFirstMax]) = NaN;
  end


  % sort eliminated timesteps to the end
  IO(t).SetTime(isnan(IO(t).SetValue)) = NaN;
  [IO(t).SetTime, IX] = sort(IO(t).SetTime,1);
  IX = IX + ones(size(IX,1),1)*(cumsum(size(IX,1)*ones(1,size(IX,2)),2)-size(IX,1));
  IO(t).SetValue = IO(t).SetValue(IX);


  % Remove time steps that are "empty" (NaN) in all tReps
  emptyRow = all(isnan(IO(t).SetTime), 2);
  IO(t).SetTime(emptyRow,:) = [];
  IO(t).SetValue(emptyRow,:) = [];

  % IO fields mustn't be empty
  if isempty(IO(t).SetTime), IO(t).SetTime = NaN(1, max([size(IO1(t).SetTime, 2), size(IO2(t).SetTime, 2)])); end
  if isempty(IO(t).SetValue), IO(t).SetValue = NaN(1, max([size(IO1(t).SetTime, 2), size(IO2(t).SetTime, 2)])); end
end

% Add elements of IO with more IO ports
for t = (numel(IO2)+1):numel(IO1), IO(t)=IO1(t); end  % IO1 larger than IO2
for t = (numel(IO1)+1):numel(IO2), IO(t)=IO2(t); end  % IO2 larger than IO1

end
