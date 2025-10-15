function [Seq, data_1D, data_S_AQs_TRs] = PostProcessing_Filter(nargout_set_sequence, ...
  Seq, raw_data_AQs, num_averages, Channel, iDevice, nucleusX)
%% Sort raw data into structures and do some post-processing (incl. Frequency Filter)
%
%   [Seq, data_1D, data_S_AQs_TRs] = PostProcessing_Filter(nargout_set_sequence, ...
%     Seq, raw_data_AQs, num_averages, Channel, iDevice, nucleusX)
%
% This function uses zeroFill_image as a "band-pass" filter for the measured
% data.
%
% See also: get_data, get_data_1d, Postprocessing
%
% ------------------------------------------------------------------------------
% (C) Copyright 2020-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if nargin < 7, nucleusX = false; end

%% find corresponding AQ structure
iAQ = find([Seq.AQ(:).Channel] == Channel & [Seq.AQ(:).Device] == iDevice, 1, 'first');

%% raw data
if isempty(raw_data_AQs)
  % create arrays with NaN
  if all(Seq.AQ(iAQ).fSample(Seq.AQ(iAQ).nSamples(:)>0) == Seq.HW.RX(iDevice).fSample)
    % un-mixed raw data
    raw_data_AQs = NaN(max(Seq.AQ(iAQ).nSamples(:)), sum(Seq.AQ(iAQ).nSamples(:)>0));
  else
    % raw data from mixer
    raw_data_AQs = NaN(max(Seq.AQ(iAQ).nSamples(:)), sum(Seq.AQ(iAQ).nSamples(:)>0)) + 1i*NaN;
  end
end

%% data (samples x AQ windows x tReps)
if nargout_set_sequence > 2
  [Seq, data_S_AQs_TRs] = get_data(Seq, raw_data_AQs, num_averages, Channel, iDevice, nucleusX);

  % Filter data
  if isemptyfield(Seq, 'PostProcessing'), Seq.PostProcessing = struct(); end
  if isemptyfield(Seq.PostProcessing, 'BW'), Seq.PostProcessing.BW = HW.fLarmor/8; end  % frequency bandwidth in Hz
  if isemptyfield(Seq.PostProcessing, 'offset'), Seq.PostProcessing.offset = 0; end  % frequency offset in Hz
  if isemptyfield(Seq.PostProcessing, 'onlyValid'), Seq.PostProcessing.onlyValid = 0; end  % remove re-sampled samples at border that are probably invalid

  iAQ_dev = find([Seq.AQ(:).Device] == iDevice, 1, 'first');
  for itRep = 1:size(Seq.AQ(iAQ_dev).nSamples, 2)
    for iAQ = 1:size(Seq.AQ(iAQ_dev).nSamples, 1)
      if isnan(Seq.AQ(iAQ_dev).Start(iAQ,itRep))
        continue;
      end
      removeSamples = Seq.PostProcessing.onlyValid*ceil(Seq.AQ(iAQ_dev).fSample(iAQ,itRep)/Seq.PostProcessing.BW);

      DataOfEcho = data_S_AQs_TRs.data(1:Seq.AQ(iAQ_dev).nSamples(iAQ,itRep),iAQ,itRep);
      TimeOfEcho = data_S_AQs_TRs.time_all(1:Seq.AQ(iAQ_dev).nSamples(iAQ,itRep),iAQ,itRep);
      % FIXME: Which is the reference sample with fixed phase?
      % DataOfEchoShift = DataOfEcho .* exp(1i*2*pi*Seq.PostProcessing.offset * (TimeOfEcho-Seq.tEcho));
      DataOfEchoShift = DataOfEcho .* exp(1i*2*pi*Seq.PostProcessing.offset * TimeOfEcho);
      DataOfEchoBW = zeroFill_image(DataOfEchoShift, size(DataOfEchoShift), ...
        [2*Seq.PostProcessing.BW/Seq.AQ(iAQ_dev).fSample(iAQ,itRep), Inf, Inf]);
      if 0
        % replace the invalid samples by complex NaN
        data_S_AQs_TRs.data(1:removeSamples,iAQ,itRep) = NaN+1i*NaN;
        data_S_AQs_TRs.data((removeSamples+1):(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples),iAQ,itRep) = ...
          DataOfEchoBW((removeSamples+1):(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples),iAQ,itRep);
        data_S_AQs_TRs.data(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples:end,iAQ,itRep) = NaN+1i*NaN;
      else
        % remove invalid samples and move valid samples to front
        data_S_AQs_TRs.data(1:(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples),iAQ,itRep) = ...
          DataOfEchoBW((removeSamples+1):(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples));
        data_S_AQs_TRs.data(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples+1:end,iAQ,itRep) = NaN+1i*NaN;
        data_S_AQs_TRs.time_all(1:(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples),iAQ,itRep) = ...
          TimeOfEcho((removeSamples+1):(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples));
        data_S_AQs_TRs.time_all(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples+1:end,iAQ,itRep) = NaN;
        data_S_AQs_TRs.time_of_tRep(1:(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples),iAQ,itRep) = ...
          data_S_AQs_TRs.time_of_tRep((removeSamples+1):(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples),iAQ,itRep);
        data_S_AQs_TRs.time_of_tRep(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples+1:end,iAQ,itRep) = NaN;
      end
    end
  end
%   % reduce size of resulting array
%   allNaN = all(all(isnan(data_S_AQs_TRs.data), 2), 3);
%   firstAllNaN = find(~allNaN, 1, 'last')+1;
%   data_S_AQs_TRs.data(firstAllNaN:end,:,:) = [];
%   data_S_AQs_TRs.time_all(firstAllNaN:end,:,:) = [];
%   data_S_AQs_TRs.time_of_tRep(firstAllNaN:end,:,:) = [];
else
  data_S_AQs_TRs = struct('time_of_tRep', NaN, 'time_all', NaN, ...
    'data', NaN, 'averages', [], 'overflow', [], ...
    'device', NaN, 'channel', NaN, 'nucleusX', NaN, ...
    'f_fft1_data', [], 'fft1_data', [], ...
    'cic_N', NaN, 'cic_M', NaN, 'cic_corr', [], 'timestamp', [], ...
    'Amplitude2Norm', [], 'Amplitude2Uin', [], 'WindowPhaseOffset', []);
end

%% data_1D (all samples divided by NaN at AQ window borders)
if nargout_set_sequence > 3
  [Seq, data_1D] = get_data_1D(Seq, data_S_AQs_TRs, Channel, iDevice, nucleusX);
else
  data_1D = struct('time_of_tRep', NaN, 'time_all', NaN, ...
    'data', NaN, 'AqFrequency', NaN, 'channel', NaN, 'device', NaN, 'nucleusX', NaN);
end

end
