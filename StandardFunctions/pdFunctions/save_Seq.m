function save_Seq(Seq, fName, force, channel)
%% Reduce data duplication and save Seq structure to file
%
%     save_Seq(Seq, fName, force, channel)
%
%
% INPUT:
%
%   Seq
%         Structure with sequence data and measurement data (stored in
%         Seq.data).
%
%   fName
%         Optional string with filename. If omitted, a file with auto-generated
%         name in the current directory is used.
%
%   force
%         If true, don't ask user if they want to over-write an existing file.
%         (Default: false)
%
%   channel
%         For systems with multiple AQ channels, the number of the channel that
%         is used. (Default: 1)
%
%
% OUTPUT:
%
%   none
%
% ------------------------------------------------------------------------------
% (C) Copyright 2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
%
% See also: load_Seq

%% input check
if nargin < 1 || isempty(Seq)
  error('PD:save_Seq:NoInput', 'Function "save_Seq" must be called with at least input "Seq".');
end
if isemptyfield(Seq, 'data') || isemptyfield(Seq.data, 'data')
  error('PD:save_Seq:NoData', 'First input argument must be a Seq structure containing measurement data.');
end

%% default input
if nargin < 2 || isempty(fName), fName = sprintf('%s_Seq', datestr(now, 'yyyymmddHHMMSS')); end
if nargin < 3 || isempty(force), force = false; end
if nargin < 4 || isempty(channel), channel = 1; end

% FIXME: Saving the raw data of multiple AQ channels isn't handled completely.

if ~force && (exist(fName, 'file') || exist([fName '.mat'], 'file'))
  res = questdlg(sprintf(['A file with the name "%s" already exists.\n' ...
    'Do you want to over-write the existing file?'], fName), ...
    'Over-write File', 'Yes', 'No', 'Yes');
  if ~strcmp(res, 'Yes')
    return;
  end
end

%% prepare structure and data for saving
% get "original" raw data array (still including AQ window phase correction)
AQs = ~isnan(Seq.AQ(channel).Start(:));  % FIXME: Can we rely on this? Should we better use nSamples > 0?
raw_data = Seq.data.data(:,AQs);

if isfield(Seq.data, 'averages') && isequal(size(Seq.data.data), size(Seq.data.averages))
  averages = Seq.data.averages(:,AQs);
else
  averages = [];
end

% store data fields that cannot be restored otherwise
data = struct();
if isfield(Seq.data, 'WindowPhaseOffset')
  data.WindowPhaseOffset = Seq.data.WindowPhaseOffset;
end
if isfield(Seq.data, 'tImageZ')
  data.tImageZ = Seq.data.tImageZ;
end
if isfield(Seq.data, 'StartSequenceTime')
  data.StartSequenceTime = Seq.data.StartSequenceTime;
end
if isfield(Seq.data, 'fCenter')
  data.fCenter = Seq.data.fCenter;
end

% remove data structure from Seq
Seq.data = [];
Seq = rmfield(Seq, 'data');

%% save data in file
Seq; data; raw_data; averages; %#ok<VUNUS>
save(fName, 'Seq', 'data', 'raw_data', 'averages');

end
