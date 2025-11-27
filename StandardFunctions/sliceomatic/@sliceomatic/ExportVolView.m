function ExportVolView(this, settingsIn)
%% Export image in VolView format
%
%     sliceomatic.ExportVolView(settings)
%
%
% INPUT:
%
%   If the function is called without input arguments or settings.showDialog is
%   true, a dialog is shown where a user can set properties for exporting the
%   image to VolView.
%
%   settings
%       Optional structure with the following optional fields. If this structure
%       is passed as an input argument, no dialog is shown and the video is
%       exported with the settings as supplied in this argument.
%
%     showDialog
%         Structure with the following optional fields:
%       enable
%           Boolean value to indicate whether a user dialog to set video export
%           properties should be shown.
%       showFilename
%           Boolean value to indicate whether the filename of the exported video
%           can be changed in the dialog.
%
%     filename
%         String with the filename of the exported VolView data file. A file of
%         the same name but with the appended extension .vvi (containing the
%         metadata for VolView) is created at the same location as the data
%         file. If this string doesn't include a path, the files are stored in
%         the current folder.
%         (Default: 'sliceomatic')
%
%     scalingFactor
%         Scaling factor from figure coordinates to meters.
%         (Default: 1e-3)
%
%
% OUTPUT:
%
%   none
%
%
% -----------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% -----------------------------------------------------------------------

if nargin < 2
  settings.showDialog.enable = true;
else
  settings = settingsIn;
end

if isemptyfield(settings, {'showDialog', 'enable'})
  settings.showDialog.enable = false;
end
if isemptyfield(settings, {'showDialog', 'showFilename'})
  settings.showDialog.showFilename = true;
end

if settings.showDialog.enable
  settings.abort = true;
  uiwait(getSettings(this));
  if settings.abort
    return;
  end
end

if isemptyfield(settings, 'filename')
  settings.filename = 'sliceomatic';
end
if isemptyfield(settings, 'scalingFactor') || ~isfinite(settings.scalingFactor)
  settings.scalingFactor = 1e-3;
end

% create Seq structure with data from sliceomatic figure
d = getappdata(this.hParent, 'sliceomatic');

Seq.data.Image = permute(d.data, [2, 1, 3]);
Seq.data.ImageZ = Seq.data.Image;
Seq.AQSlice(1).ReadOS = 1;
Seq.AQSlice(1).PhaseOS = [1, 1, 1];

Seq.AQSlice(1).sizeRead = ...
  (diff(this.xmesh([1,end])) + diff(this.xmesh(1:2))) ...
  * settings.scalingFactor;
Seq.AQSlice(1).nRead = size(Seq.data.Image, 1);
Seq.AQSlice(1).sizePhase(1) = ...
  (diff(this.ymesh([1,end])) + diff(this.ymesh(1:2))) ...
  * settings.scalingFactor;
Seq.AQSlice(1).sizePhase(2) = ...
  (diff(this.zmesh([1,end])) + diff(this.zmesh(1:2))) ...
  * settings.scalingFactor;
Seq.AQSlice(1).sizePhase(3) = 1;
Seq.AQSlice(1).nPhase = [size(Seq.data.Image,2), size(Seq.data.Image,3), 1];

Seq.AQSlice(1).Center2OriginImage(1) = mean(this.xmesh([1,end])) ...
  + diff(this.xmesh(1:2))/2 * mod(numel(this.xmesh) + 1, 2);
Seq.AQSlice(1).Center2OriginImage(2) = mean(this.ymesh([1,end])) ...
  + diff(this.ymesh(1:2))/2 * mod(numel(this.ymesh) + 1, 2);
Seq.AQSlice(1).Center2OriginImage(3) = mean(this.zmesh([1,end])) ...
  + diff(this.zmesh(1:2))/2 * mod(numel(this.zmesh) + 1, 2);
Seq.AQSlice(1).Center2OriginImage = Seq.AQSlice(1).Center2OriginImage ...
  * settings.scalingFactor;

% export to VolView format
Image3D2VolView(Seq, settings.filename);

  function hfSettings = getSettings(this)
    %% open figure with settings for video

    posSli = get(this.hFigure, 'Position');

    szSettings = [370, 85+55*settings.showDialog.showFilename];
    centerSli = posSli(1:2) + posSli(3:4)/2;
    posSett = [centerSli(1:2)-szSettings/2, szSettings];
    hfSettings = figure(...
      'Name', 'VolView Export Settings', ...
      'Position', posSett, ...
      'Menubar', 'none', ...
      'Toolbar', 'none', ...
      'NumberTitle', 'off');


    baseLine = szSettings(2);
    if settings.showDialog.showFilename
      % file name
      baseLine = baseLine - 30;
      baseLine = baseLine - 25;
      uicontrol(...
        'Style', 'text', ...
        'String', 'Filename:', ...
        'Position', [20, baseLine, 150, 18], ...
        'HorizontalAlignment', 'left');
      uiHandles.filename = uicontrol(...
        'Style', 'edit', ...
        'String', 'sliceomatic', ...
        'Position', [180, baseLine, 170, 22], ...
        'HorizontalAlignment', 'left');

      uicontrol(...
        'Style', 'pushbutton', ...
        'String', 'Choose filename...', ...
        'Position', [180, baseLine+25, 170, 22], ...
        'Callback', @(h,e) FilenameCallback(uiHandles.filename));
    end

    % scaling factor
    baseLine =  baseLine - 25;
    uicontrol(...
      'Style', 'text', ...
      'String', 'Scaling factor:', ...
      'Position', [20, baseLine, 150, 18], ...
      'HorizontalAlignment', 'left');
    uiHandles.scalingFactor = uicontrol(...
      'Style', 'edit', ...
      'String', '1e-3', ...
      'Position', [180, baseLine, 170, 22], ...
      'HorizontalAlignment', 'left');


    %%
    baseLine =  baseLine - 45;
    uicontrol(...
      'Style', 'pushbutton', ...
      'String', 'Export VolView', ...
      'Position', [20, baseLine, 150, 30], ...
      'HorizontalAlignment', 'left', ...
      'Callback', @(h,e) ExportButtonCallback(hfSettings, uiHandles));

    uicontrol(...
      'Style', 'pushbutton', ...
      'String', 'Abort', ...
      'Position', [200, baseLine, 150, 30], ...
      'HorizontalAlignment', 'left', ...
      'Callback', @(h,e) close(hfSettings));

  end


  function FilenameCallback(huiFilename)
    %% file selection dialog
    [fileName, pathName] = uiputfile('*', 'Select file name', 'sliceomatic');
    if ~(isnumeric(fileName) && fileName == 0)
      set(huiFilename, 'String', fullfile(pathName, fileName));
    end
  end


  function ExportButtonCallback(hfSettings, uiHandles)
    %% get settings from ui fields and close dialog
    if settings.showDialog.showFilename
      settings.filename = get(uiHandles.filename, 'String');
    end
    settings.scalingFactor = str2double(get(uiHandles.scalingFactor, 'String'));

    settings.abort = false;
    close(hfSettings);
  end

end
