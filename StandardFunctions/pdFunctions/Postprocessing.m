function [Seq, data_1D, data_S_AQs_TRs] = Postprocessing(nargout_set_sequence, ...
  Seq, raw_data_AQs, num_averages, Channel, iDevice, nucleusX)
%% Sort raw data into structures and do some post-processing
%
%   [Seq, data_1D, data_S_AQs_TRs] = Postprocessing(nargout_set_sequence, ...
%     Seq, raw_data_AQs, num_averages, Channel, iDevice, nucleusX)
%
% See also: get_data, get_data_1d
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2020 Pure Devices GmbH, Wuerzburg, Germany
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
else
  data_S_AQs_TRs = struct('time_of_tRep', NaN, 'time_all', NaN, ...
  'data', NaN, 'averages', [], 'channel', NaN, 'device', NaN, 'nucleusX', NaN, ...
  'f_fft1_data', [], 'fft1_data', [], ...
  'cic_N', NaN, 'cic_M', NaN, 'cic_corr', [], ...
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
