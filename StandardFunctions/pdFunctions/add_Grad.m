function Grad = add_Grad(Grad1, Grad2)
%% Add gradient pulses from two Grad structures into a single one
%
%   Grad = add_Grad(Grad1, Grad2)
%
%
% INPUT:
%
%   Grad1, Grad2
%         Two input Grad structures (see manual) that should be combined into
%         one. Both input structures might be vectors of potentially different
%         size. The fields in the corresponding elements are "added" on top of
%         each other. Points which would not change the resulting gradient pulse
%         shape are removed.
%
%
% OUTPUT:
%
%   Grad
%         Output Grad structure with the combination of the two input Grad
%         structures.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2013-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% default input
[Grad1(arrayfun(@(x) isemptyfield(x, 'Shim'), Grad1)).Shim] = deal(0);
[Grad2(arrayfun(@(x) isemptyfield(x, 'Shim'), Grad2)).Shim] = deal(0);

% loop over common gradient elements
for t = min(numel(Grad1), numel(Grad2)):-1:1
  if all(abs(Grad1(t).Amp(~isnan(Grad1(t).Time)))<1e-12)
    % no pulses in Grad1
    Grad(t).Time = Grad2(t).Time;
    Grad(t).Amp = Grad2(t).Amp;
    Grad(t).Shim = sum([Grad1(t).Shim, Grad2(t).Shim]);
    continue;
  end
  if all(abs(Grad2(t).Amp(~isnan(Grad2(t).Time)))<1e-12)
    % no pulses in Grad2
    Grad(t).Time = Grad1(t).Time;
    Grad(t).Amp = Grad1(t).Amp;
    Grad(t).Shim = sum([Grad1(t).Shim, Grad2(t).Shim]);
    continue;
  end

  numtRep = max(size(Grad1(t).Time,2), size(Grad2(t).Time,2));
  % combine pulses in both Grad structures
  Grad(t).Time = NaN(size(Grad1(t).Time,1) + size(Grad2(t).Time,1), numtRep);
  Grad(t).Amp = NaN(size(Grad1(t).Time,1) + size(Grad2(t).Time,1), numtRep);
  if size(Grad1(t).Time,2) == 1
    % apply pulse(s) in Grad1 to all tReps in Grad2
    Grad1(t).Time = Grad1(t).Time * ones(1,max(size(Grad2(t).Time,2),1));
    Grad1(t).Amp = Grad1(t).Amp * ones(1,max(size(Grad2(t).Time,2),1));
  end
  if size(Grad2(t).Time,2) == 1
    % apply pulse(s) in Grad2 to all tReps in Grad1
    Grad2(t).Time = Grad2(t).Time * ones(1,max(size(Grad1(t).Time,2),1));
    Grad2(t).Amp = Grad2(t).Amp * ones(1,max(size(Grad1(t).Time,2),1));
  end
  if ~any(~isnan(Grad1(t).Time(:))) && ~any(~isnan(Grad2(t).Time(:)))
    % all time values in Grad1 and Grad2 are none
    Grad(t).Time = [];
    Grad(t).Amp = [];
    continue;
  end

  % create timeline with all tReps (each tRep separated by 1000 seconds)
  Grad1(t).Time = bsxfun(@plus, Grad1(t).Time, (1:numtRep)*1000);
  Grad2(t).Time = bsxfun(@plus, Grad2(t).Time, (1:numtRep)*1000);

  % indices for tRep of each data point
  temp1_tRep = cumsum(ones(size(Grad1(t).Time)), 2);
  temp2_tRep = cumsum(ones(size(Grad2(t).Time)), 2);
  temp1_tRep = temp1_tRep(~isnan(Grad1(t).Time(:)));
  temp2_tRep = temp2_tRep(~isnan(Grad2(t).Time(:)));

  if numel(temp1_tRep) == 0
    % no pulses in Grad1 --> use data from Grad2
    TimeAll = Grad2(t).Time(~isnan(Grad2(t).Time(:)));
    AmpAll = Grad2(t).Amp(~isnan(Grad2(t).Time(:)));
    temp_tRepAll = temp2_tRep(:);
  elseif numel(temp2_tRep) == 0
    % no pulses in Grad2 --> use data from Grad1
    TimeAll = Grad1(t).Time(~isnan(Grad1(t).Time(:)));
    AmpAll = Grad1(t).Amp(~isnan(Grad1(t).Time(:)));
    temp_tRepAll = temp1_tRep(:);
  else
    % combine data from both Grad1 and Grad2
    % remove data where time is NaN
    Amp1 = Grad1(t).Amp(~isnan(Grad1(t).Time(:)));
    Amp2 = Grad2(t).Amp(~isnan(Grad2(t).Time(:)));
    Time1 = Grad1(t).Time(~isnan(Grad1(t).Time(:)));
    Time2 = Grad2(t).Time(~isnan(Grad2(t).Time(:)));
    % linearly interpolate Grad1 at times from Grad2 and vice versa
    % (copy first and last amplitude value to very large times for extrapolation)
    Amp1at2 = interp1([-1e12;Time1;1e12], [Amp1(1);Amp1;Amp1(end)], Time2, 'linear');
    Amp2at1 = interp1([-1e12;Time2;1e12], [Amp2(1);Amp2;Amp2(end)], Time1, 'linear');
    % calculate sum on both timelines respectively
    Amp22 = Amp1at2 + Amp2;
    Amp11 = Amp2at1 + Amp1;
    % combine data from both inputs
    AmpAll = [Amp11; Amp22];
    TimeAll = [Time1; Time2];
    temp_tRepAll = [temp1_tRep(:); temp2_tRep(:)];
  end

  % remove duplicate times
  [Time12, IT12] = unique(TimeAll);  % FIXME: Do we need to check if amplitude in both inputs matches?
  Amp12 = AmpAll(IT12);
  temp_tRep = temp_tRepAll(IT12);

  % convert back to times of tRep
  Time12 = Time12 - temp_tRep*1000;

  % count number of (valid) data points in each tRep
  count_tRep = histc(temp_tRep, 1:numtRep);
  % count_start = 0;

  % create logical matrix for indexing into rectangular Grad.Time and Grad.Amp
  logic_use = false(size(Grad(t).Time));
  for n = 1:numtRep
    % Grad(t).Time(1:count_tRep(n),n) = Time12(count_start+(1:count_tRep(n)));
    % Grad(t).Amp(1:count_tRep(n),n) = Amp12(count_start+(1:count_tRep(n)));
    logic_use(1:count_tRep(n),n) = true;
    % count_start = count_start+count_tRep(n);
  end
  Grad(t).Time(logic_use) = Time12;
  Grad(t).Amp(logic_use) = Amp12;


  % remove unnecessary steps
  % duplicate amplitude at start of each tRep
  % (necessary for detection of unnecessary steps at beginning)
  Grad(t).Time = [Grad(t).Time(1,:)-1000; Grad(t).Time; NaN(2, numtRep)];
  Grad(t).Amp = [Grad(t).Amp(1,:); Grad(t).Amp; NaN(2, numtRep)];

  dummy_values = false(size(Grad(t).Time));
  dummy_values(1,:) = true;

  % duplicate last valid amplitude in each tRep
  % (necessary for detection of unnecessary steps at end)
  lastValidIdx = sum(~isnan(Grad(t).Amp) & ~isnan(Grad(t).Time), 1);
  tRepIdx = 1:numtRep;
  lastValidIdx(lastValidIdx==0) = 1;
  linIdxLastValid = sub2ind(size(Grad(t).Time), lastValidIdx, tRepIdx);
  Grad(t).Time(linIdxLastValid+1) = Grad(t).Time(linIdxLastValid)+1000;
  Grad(t).Amp(linIdxLastValid+1) = Grad(t).Amp(linIdxLastValid);
  dummy_values(linIdxLastValid+1) = true;

  % If the slope doesn't change (i.e. curvature is 0), we can remove that point.
  slope = diff(Grad(t).Amp, 1, 1) ./ diff(Grad(t).Time, 1, 1);
  curve = diff(slope, 1, 1);  % no need to scale with the time steps because only look for ==0

  Grad(t).Time(dummy_values) = [];
  Grad(t).Amp(dummy_values) = [];
  Grad(t).Time = reshape(Grad(t).Time, [], numtRep);
  Grad(t).Amp = reshape(Grad(t).Amp, [], numtRep);

  % Move subsequent actions in each tRep upwards whereever an action can be
  % removed.
  % Indices with actions in each tRep we want to keep.
  % (Matlab keeps the last one if an index is used multiple times.)
  idxAct = bsxfun(@minus, sum(curve~=0,1)+1, cumsum(curve~=0,1,'reverse'));
  % convert to linear indices
  linIdx = sub2ind(size(Grad(t).Time), idxAct, repmat(1:numtRep, size(Grad(t).Time, 1), 1));
  % move valid actions to front in each tRep
  Grad(t).Time(linIdx) = Grad(t).Time;
  Grad(t).Amp(linIdx) = Grad(t).Amp;
  % delete invalid actions at end of each tRep
  iLastValid = max(idxAct, [], 1);
  Grad(t).Time(bsxfun(@gt, (1:size(Grad(t).Time,1)).', iLastValid)) = NaN;


  % for n=1:max(size(Grad1(t).Time,2),size(Grad2(t).Time,2))
  %   if size(Grad1(t).Time,2)==0
  %     TimeAll=Grad2(t).Time(1:find(~isnan(Grad2(t).Time(:,n)),1,'last'),n);
  %     AmpAll=Grad2(t).Amp(1:find(~isnan(Grad2(t).Time(:,n)),1,'last'),n);
  %   elseif size(Grad2(t).Time,2)==0
  %     TimeAll=Grad1(t).Time(1:find(~isnan(Grad1(t).Time(:,n)),1,'last'),n);
  %     AmpAll=Grad1(t).Amp(1:find(~isnan(Grad1(t).Time(:,n)),1,'last'),n);
  %   else
  %     Amp1=Grad1(t).Amp(1:find(~isnan(Grad1(t).Time(:,n)),1,'last'),n);
  %     Amp2=Grad2(t).Amp(1:find(~isnan(Grad2(t).Time(:,n)),1,'last'),n);
  %     Time1=Grad1(t).Time(1:find(~isnan(Grad1(t).Time(:,n)),1,'last'),n);
  %     Time2=Grad2(t).Time(1:find(~isnan(Grad2(t).Time(:,n)),1,'last'),n);
  %     % Amp1at2=interp1([-1e12;Time1;1e12],[Amp1(1);Amp1;Amp1(end)],Time2,'linear');
  %     % Amp2at1=interp1([-1e12;Time2;1e12],[Amp2(1);Amp2;Amp2(end)],Time1,'linear');
  %     Amp1at2=interp1q([-1e12;Time1;1e12],[Amp1(1);Amp1;Amp1(end)],Time2);
  %     Amp2at1=interp1q([-1e12;Time2;1e12],[Amp2(1);Amp2;Amp2(end)],Time1);
  %
  %     Amp22=Amp1at2+Amp2;
  %     Amp11=Amp2at1+Amp1;
  %     AmpAll=[Amp11;Amp22];
  %     TimeAll=[Time1;Time2];
  %   end
  %   [Time12,IT12] = unique(TimeAll);  % attention: no discontinuities in amplitude
  %   Amp12=AmpAll(IT12);
  %
  %   Grad(t).Time(1:size(Time12,1),n)=Time12;
  %   Grad(t).Amp(1:size(Time12,1),n)=Amp12;
  % end

  % remove rows at the end without any action in all tReps
  [row, ~] = find(sum(~isnan(Grad(t).Time),2), 1, 'last');
  Grad(t).Time = Grad(t).Time(1:max(row),:);
  Grad(t).Amp = Grad(t).Amp(1:max(row),:);

  % use sum of shim
  Grad(t).Shim = sum([Grad1(t).Shim,Grad2(t).Shim]);
end

% copy Grad1 vector elements without matching partner in Grad2
for t = (numel(Grad2)+1):numel(Grad1)
  Grad(t).Time = Grad1(t).Time;
  Grad(t).Amp = Grad1(t).Amp;
  Grad(t).Shim = sum(Grad1(t).Shim);
end
% copy Grad2 vector elements without matching partner in Grad1
for t = (numel(Grad1)+1):numel(Grad2)
  Grad(t).Time = Grad2(t).Time;
  Grad(t).Amp = Grad2(t).Amp;
  Grad(t).Shim = sum(Grad2(t).Shim);
end

end
