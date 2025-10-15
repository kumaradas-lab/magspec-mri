function Seq = load_Seq(fName)
%% Load Seq structure and measurement data from file
%
%     Seq = load_Seq(fName)
%
%
% INPUT:
%
%   fName
%           Optional string with the file name that stores the sequence and
%           measurement data. If omitted, the user is asked to select a file
%           interactively.
%
%
% OUTPUT:
%
%   Seq
%         Structure with sequence data and measurement data (stored in
%         Seq.data).
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
%
% See also: save_Seq

%% input check
if nargin < 1 || isempty(fName)
  [fName, fPath] = uigetfile('*.mat', 'Select File');
  if isnumeric(fName), return; end
  fName = fullfile(fPath, fName);
end

if ~exist(fName, 'file')
  if exist([fName '.mat'], 'file')
    fName = [fName '.mat'];
  else
    error('PD:load_Seq:UnknownFile', 'Cannot find file "%s".', fName);
  end
end

%% load file
averages = [];  % Default value might be overwritten when loading the file.

load(fName);

if ~exist('Seq', 'var') || ~exist('data', 'var') || ~exist('raw_data', 'var')
  error('PD:load_Seq:InvalidFile', 'Could not load data from file "%s".', fName);
end

%% reconstruct data
% FIXME: Loading the raw data of multiple AQ channels isn't handled.

if Seq.PostProcessSequence
  oldCorrectAQWindowPhase = {Seq.AQ(:).CorrectAQWindowPhase};
  [Seq.AQ(:).CorrectAQWindowPhase] = deal(false);
  [Seq, data_S_AQs_TRs] = get_data(Seq, raw_data, averages, 1);
  Seq.data = data_S_AQs_TRs;
  [Seq.AQ(:).CorrectAQWindowPhase] = deal(oldCorrectAQWindowPhase{:});
end

dataFields = fieldnames(data);

for iField = 1:numel(dataFields)
  Seq.data.(dataFields{iField}) = data.(dataFields{iField});
end

end
