function ExportAvi(this, settingsIn)
%% Export video with moving slices or surfaces
%
%     sliceomatic.ExportAvi(settings)
%
%
% INPUT:
%
%   If the function is called without input arguments or settings.showDialog is
%   true, a dialog is shown where a user can set properties for exporting the
%   video.
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
%         String with the filename of the exported video. If this doesn't
%         include a path, the video is stored in the current folder.
%         (Default: 'sliceomatic.avi')
%
%     axis
%         Integer that selects the direction in which a slice is moved (1: x, 2:
%         y, 3: z) or 4 for moving iso values. Note that these integers refer to
%         the axis direction in the figure. That does not necessarily correspond
%         to the axis direction of the displayed object.
%         (Default: z)
%
%     numFrames
%         Total number of frames in the video. The total number of frames might
%         be doubled when also exporting the backwards movement of the slice or
%         iso surface (see below).
%         (Default: 180)
%
%     framesPerSec
%         Number of frames per second in the video. The duration of the video
%         in seconds is numFrames/framesPerSec.
%         (Default: 60)
%
%     startVal
%         Coordinate or value of the slice or iso surface in the first frame of
%         the video. If the value exceeds the limits of the respective axis or
%         the data, it is clipped to respective limit.
%         (Default: lowest coordinate or value of the respective axis or data)
%
%     stopVal
%         Coordinate or value of the slice or iso surface in the last frame of
%         the video (not taking a potential backwards movement into account).
%         stopVal can be larger or lower than startVal. If the value exceeds the
%         limits of the respective axis or the data, it is clipped to respective
%         limit.
%         (Default: highest coordinate or value of the respective axis or data)
%
%     loopSlice
%         Boolean value to indicate whether the video should include additional
%         frames in which the slice or the iso surface moves back from the
%         stopVal to the startVal. If this field is true, the duration of the
%         video doubles.
%         (Default: true)
%
%     autoPlay
%         Boolean value to indicate whether the exported video is opened in the
%         default application when the export is complete. This is only
%         supported on Windows. On other platforms, the value of this field is
%         ignored.
%         (Default: true on Windows, false otherwise)
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
  settings.filename = 'sliceomatic.avi';
end
if numel(settings.filename) < 4 || ~strcmp(settings.filename(end-3:end), '.avi')
  settings.filename = [settings.filename, '.avi'];
end
if isemptyfield(settings, 'axis')
  settings.axis = 3;
end
if isemptyfield(settings, 'numFrames') || ~isfinite(settings.numFrames)
  settings.numFrames = 180;
end
if isemptyfield(settings, 'framesPerSec') || ~isfinite(settings.framesPerSec)
  settings.framesPerSec = 60;
end
if settings.axis < 4
  limits = get(this.hAxes, [char(87+settings.axis), 'Lim']);
else
  limits = get(this.hAxes, 'CLim');
end
if isemptyfield(settings, 'startVal') || isnan(settings.startVal)
  settings.startVal = limits(1);
end
if isemptyfield(settings, 'stopVal') || isnan(settings.stopVal)
  settings.stopVal = limits(2);
end
if settings.startVal < settings.stopVal
  settings.startVal = max(settings.startVal, limits(1));
  settings.stopVal = min(settings.stopVal, limits(2));
else
  settings.startVal = min(settings.startVal, limits(2));
  settings.stopVal = max(settings.stopVal, limits(1));
end

if isemptyfield(settings, 'loopSlice')
  settings.loopSlice = true;
end
if isemptyfield(settings, 'autoPlay')
  settings.autoPlay = ispc();
end

% create vector with the positions of the slice in the video
% y_pos = linspace(min(SeqLoop.data.Ticks(2).PhaseZ), max(SeqLoop.data.Ticks(2).PhaseZ), num_frames);
pos = linspace(settings.startVal, settings.stopVal, settings.numFrames);
if settings.loopSlice
  % move slice forth and back in the video
  pos = [pos, fliplr(pos)];
end

if ~ishghandle(this.hParent, 'figure')
  % copy sliceomatic to (hidden) figure
  hf = figure('Visible', 'off');
  protectFigure = onCleanup(@() close(hf));
  hsl = this.copyobj(hf);
else
  hsl = this;
end

% add slice at initial position
switch settings.axis
  case 1
    hSlice = hsl.AddSliceX(pos(1));
  case 2
    hSlice = hsl.AddSliceY(pos(1));
  case 3
    hSlice = hsl.AddSliceZ(pos(1));
  case 4
    hSlice = hsl.AddIsoSurface(pos(1));
  otherwise
    error('PD:sliceomatic:aviExport:UnsupportedAxes', ...
      'The value of the axes setting must be 1 (x), 2 (y), 3 (z), or 4 (val).');
end

% initialize video writer
vidWriter = VideoWriter(settings.filename);
vidWriter.Quality = 80;
vidWriter.FrameRate = settings.framesPerSec;

% show waitbar
hwb = waitbar(0, 'Exporting video...');
hwb_protect = onCleanup(@() close(hwb));

open(vidWriter);
for iFrame = 1:numel(pos)
  % move slice for current frame
  waitbar(iFrame/numel(pos), hwb);
  if settings.axis == 4
    hsl.MoveIsoSurface(hSlice, pos(iFrame));
  else
    hsl.MoveSlice(hSlice, pos(iFrame));
  end
  drawnow();

  % get frame and store in video
  currFrame = getframe(hsl.hFigure);
  vidWriter.writeVideo(currFrame);
end
close(vidWriter);

if settings.axis == 4
  hsl.DeleteIso(hSlice)
else
  hsl.DeleteSlice(hSlice)
end
delete(hwb_protect);

if ispc() && settings.autoPlay
  % open created video in default application
  winopen(settings.filename);
end

  function hfSettings = getSettings(this)
    %% open figure with settings for video

    posSli = get(this.hFigure, 'Position');

    szSettings = [370, 215+55*settings.showDialog.showFilename];
    centerSli = posSli(1:2) + posSli(3:4)/2;
    posSett = [centerSli(1:2)-szSettings/2, szSettings];
    hfSettings = figure(...
      'Name', 'Video Export Settings', ...
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
        'String', 'sliceomatic.avi', ...
        'Position', [180, baseLine, 170, 22], ...
        'HorizontalAlignment', 'left');

      uicontrol(...
        'Style', 'pushbutton', ...
        'String', 'Choose filename...', ...
        'Position', [180, baseLine+25, 170, 22], ...
        'Callback', @(h,e) FilenameCallback(uiHandles.filename));
    end

    % slice direction
    baseLine = baseLine - 30;
    axesLabels = cell(1,3);
    for iAxis = 1:3
      axesLabels{iAxis} = get(get(this.hAxes, [char(87+iAxis), 'Label']), 'String');
    end
    uicontrol(...
      'Style', 'text', ...
      'String', 'Move surface or slice along:', ...
      'Position', [20, baseLine, 150, 22], ...
      'HorizontalAlignment', 'left');
    uiHandles.axis = uicontrol(...
      'Style', 'popupmenu', ...
      'String', [axesLabels, {'ISO surface'}], ...
      'Value', 3, ...
      'Position', [180, baseLine-4, 170, 28]);

    % number of frames
    baseLine =  baseLine - 25;
    uicontrol(...
      'Style', 'text', ...
      'String', 'Number of frames:', ...
      'Position', [20, baseLine, 150, 18], ...
      'HorizontalAlignment', 'left');
    uiHandles.numFrames = uicontrol(...
      'Style', 'edit', ...
      'String', '180', ...
      'Position', [180, baseLine, 170, 22], ...
      'HorizontalAlignment', 'left');

    % frames per second
    baseLine =  baseLine - 25;
    uicontrol(...
      'Style', 'text', ...
      'String', 'Frames per second:', ...
      'Position', [20, baseLine, 150, 18], ...
      'HorizontalAlignment', 'left');
    uiHandles.framesPerSec = uicontrol(...
      'Style', 'edit', ...
      'String', '60', ...
      'Position', [180, baseLine, 170, 22], ...
      'HorizontalAlignment', 'left');

    % limits = get(this.hAxes, [char(87+iAxis), 'Lim']);
    % start value
    baseLine =  baseLine - 25;
    uicontrol(...
      'Style', 'text', ...
      'String', 'Start value:', ...
      'Position', [20, baseLine, 150, 18], ...
      'HorizontalAlignment', 'left');
    uiHandles.startVal = uicontrol(...
      'Style', 'edit', ...
      'String', '-Inf', ... % num2str(limits(1)), ...
      'Position', [180, baseLine, 170, 22], ...
      'HorizontalAlignment', 'left');

    % stop value
    baseLine =  baseLine - 25;
    uicontrol(...
      'Style', 'text', ...
      'String', 'Stop value:', ...
      'Position', [20, baseLine, 150, 18], ...
      'HorizontalAlignment', 'left');
    uiHandles.stopVal = uicontrol(...
      'Style', 'edit', ...
      'String', 'Inf', ... % num2str(limits(2)), ...
      'Position', [180, baseLine, 170, 22], ...
      'HorizontalAlignment', 'left');

    % move slice forth and back in the video
    baseLine =  baseLine - 25;
    uiHandles.loopSlice = uicontrol(...
      'Style', 'checkbox', ...
      'String', 'Also move back', ...
      'Value', true, ...
      'Position', [20, baseLine, 170, 22], ...
      'HorizontalAlignment', 'left');

    % autoplay
    uiHandles.autoPlay = uicontrol(...
      'Style', 'checkbox', ...
      'String', 'Auto play', ...
      'Value', true, ...
      'Position', [180, baseLine, 170, 22], ...
      'HorizontalAlignment', 'left');
    if ~ispc()
      set(uiHandles.autoPlay, ...
        'Value', false, ...
        'Enable', 'off');
    end

    %%
    baseLine =  baseLine - 45;
    uicontrol(...
      'Style', 'pushbutton', ...
      'String', 'Export Video', ...
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
    [fileName, pathName] = uiputfile('*.avi', 'Select file name', 'sliceomatic.avi');
    if ~(isnumeric(fileName) && fileName == 0)
      set(huiFilename, 'String', fullfile(pathName, fileName));
    end
  end


  function ExportButtonCallback(hfSettings, uiHandles)
    %% get settings from ui fields and close dialog
    if settings.showDialog.showFilename
      settings.filename = get(uiHandles.filename, 'String');
    end
    settings.axis = get(uiHandles.axis, 'Value');
    settings.numFrames = str2double(get(uiHandles.numFrames, 'String'));
    settings.framesPerSec = str2double(get(uiHandles.framesPerSec, 'String'));
    settings.startVal = str2double(get(uiHandles.startVal, 'String'));
    settings.stopVal = str2double(get(uiHandles.stopVal, 'String'));
    settings.loopSlice = get(uiHandles.loopSlice, 'Value');
    settings.autoPlay = get(uiHandles.autoPlay, 'Value');

    settings.abort = false;
    close(hfSettings);
  end

end
