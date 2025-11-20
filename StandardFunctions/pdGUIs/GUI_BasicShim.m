function varargout = GUI_BasicShim(varargin)
%% GUI_BASICSHIM MATLAB code for GUI_BasicShim.fig
%      GUI_BASICSHIM, by itself, creates a new GUI_BASICSHIM or raises the existing
%      singleton*.
%
%      H = GUI_BASICSHIM returns the handle to a new GUI_BASICSHIM or the handle to
%      the existing singleton*.
%
%      GUI_BASICSHIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_BASICSHIM.M with the given input arguments.
%
%      GUI_BASICSHIM('Property','Value',...) creates a new GUI_BASICSHIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_BasicShim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_BasicShim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_BasicShim

% Last Modified by GUIDE v2.5 29-Mar-2017 10:12:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @GUI_BasicShim_OpeningFcn, ...
  'gui_OutputFcn',  @GUI_BasicShim_OutputFcn, ...
  'gui_LayoutFcn',  [], ...
  'gui_Callback',   []);
if nargin && ischar(varargin{1})
  gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
  [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
  gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

end

% --- Executes just before GUI_BasicShim is made visible.
function GUI_BasicShim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_BasicShim (see VARARGIN)

% Choose default command line output for GUI_BasicShim
handles.output = hObject;

%% initialize handles.data structure
[handles.data.plot(1:3).hold] = deal(false);
[handles.data.plot(1:3).data] = deal([]);
handles.data.sliderFrequencyOldValue = 0.5;
handles.data.sliderT90OldValue = 0.5;
handles.data.sliderShimXOldValue = 0.5;
handles.data.sliderShimYOldValue = 0.5;
handles.data.sliderShimZOldValue = 0.5;
handles.data.sliderChanged = false;
handles.data.sliderFirstChange = true;

%% Update handles structure
guidata(hObject, handles);

%% add handler to uitable
% uitablepeer = findjobj(hObject, '-nomenu', 'class', 'uitablepeer');
% set(uitablepeer, 'MouseClickedCallback', {@uitable_hold_MouseClickHandler, handles});

% UIWAIT makes GUI_BasicShim wait for user response (see UIRESUME)
% uiwait(handles.GUI_BasicShim);

end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_BasicShim_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end

function GUI_BasicShim_ResizeFcn(hObject, eventdata, handles)
%% Executes when GUI_BasicShim is resized
% hObject    handle to setupDevice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get window size and position
windowPosition = get(hObject, 'Position');

%% check for minimum window size
% window height
if windowPosition(4) < 655
  windowPosition(4) = 655;
  set(hObject, 'Position', windowPosition);
  GUI_BasicShim_ResizeFcn(hObject, eventdata, handles);
  return
end
% window width
if windowPosition(3) < 740
  windowPosition(3) = 740;
  set(hObject, 'Position', windowPosition);
  GUI_BasicShim_ResizeFcn(hObject, eventdata, handles);
  return
end

%% adjust for height changes
windowHeight = windowPosition(4);
% left controls
pushbutton_loadSystem_Position = get(handles.pushbutton_loadSystem, 'Position');
pushbutton_loadSystem_Position(2) = windowHeight - 50;
set(handles.pushbutton_loadSystem, 'Position', pushbutton_loadSystem_Position);
pushbutton_loadHW_Position = get(handles.pushbutton_loadHW, 'Position');
pushbutton_loadHW_Position(2) = windowHeight - 50;
set(handles.pushbutton_loadHW, 'Position', pushbutton_loadHW_Position);
uipanel_parameter_Position = get(handles.uipanel_parameter, 'Position');
uipanel_parameter_Position(2) = windowHeight - 300;
set(handles.uipanel_parameter, 'Position', uipanel_parameter_Position);
uipanel_FID_Position = get(handles.uipanel_FID, 'Position');
uipanel_FID_Position(2) = windowHeight - 470;
set(handles.uipanel_FID, 'Position', uipanel_FID_Position);
uipanel_autoConfig_Position = get(handles.uipanel_autoConfig, 'Position');
uipanel_autoConfig_Position(2) = windowHeight - 600;
set(handles.uipanel_autoConfig, 'Position', uipanel_autoConfig_Position);

% right graphics
% keep history the same size, change height of FID
axes_FID_Position = get(handles.axes_FID, 'Position');
axes_FID_Position(4) = windowHeight - 317;
% set(handles.axes_FID, 'Position', axes_FID_Position)

%% adjust for width changes
windowWidth = windowPosition(3);
% keep left controls at fixed width

% right graphics
axes_FID_Position(3) = windowWidth - 555;
set(handles.axes_FID, 'Position', axes_FID_Position)
set(handles.axes_FID2, 'Position', axes_FID_Position)
axes_fLarmor_Position = get(handles.axes_fLarmor, 'Position');
axes_fLarmor_Position(3) = windowWidth - 555;
set(handles.axes_fLarmor, 'Position', axes_fLarmor_Position)

%% move window on screen
monitorPosition = get(0, 'MonitorPositions');
windowPositionOld = windowPosition;
% assume that the area of all monitors is rectangular
if windowPosition(1) < min(monitorPosition(:,1))
  windowPosition(1) = min(monitorPosition(:,1));
end
if windowPosition(2) < min(monitorPosition(:,2))
  windowPosition(2) = min(monitorPosition(:,2));
end
% move on screen
if windowPosition(1) + windowPosition(3) - 1 > max(monitorPosition(:,3))
  windowPosition(1) = max(monitorPosition(:,3)) - windowPosition(3);
end
% reduce size
if windowPosition(3) > max(monitorPosition(:,3)) - min(monitorPosition(:,1))
  windowPosition(3) = max(monitorPosition(:,3)) - min(monitorPosition(:,1));
end
% move on screen
if windowPosition(2) + windowPosition(4) - 1 > ...
    max(monitorPosition(:,4))
  windowPosition(2) = max(monitorPosition(:,4)) - windowPosition(4);
end
% reduce size
if windowPosition(4) > max(monitorPosition(:,4)) - min(monitorPosition(:,2))
  windowPosition(4) = max(monitorPosition(:,4)) - min(monitorPosition(:,2));
end

%{
screenSize = get(0, 'ScreenSize');
windowPositionOld = windowPosition;
if windowPosition(1) < 0
    windowPosition(1) = 0;
end
if windowPosition(2) < 0
    windowPosition(2) = 0;
end
if windowPosition(1) + windowPosition(3) > screenSize(3)
    windowPosition(1) = screenSize(3) - windowPosition(3);
end
if windowPosition(2) + windowPosition(4) > screenSize(4)
    windowPosition(2) = screenSize(4) - windowPosition(4);
end
%}
if any(windowPosition - windowPositionOld)
  set(hObject, 'Position', windowPosition);
  GUI_BasicShim_ResizeFcn(hObject, eventdata, handles);
  return
end


end


function pushbutton_loadSystem_Callback(hObject, eventdata, handles)
%% Executes on button press in pushbutton_loadSystem.
% hObject    handle to pushbutton_loadSystem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% load system data and set it up
evalin('base', 'LoadSystem');

%% update values in GUI
set(handles.edit_frequency, 'String', num2str(evalin('base', 'HW.fLarmor'), '%.0f'));
set(handles.edit_t90,       'String', num2str(evalin('base', 'HW.tFlip90Def*1e6'), '%.0f'));
set(handles.edit_shimX,     'String', num2str(evalin('base', 'HW.MagnetShim(1)*1e3'), '%.4f'));
set(handles.edit_shimY,     'String', num2str(evalin('base', 'HW.MagnetShim(2)*1e3'), '%.4f'));
set(handles.edit_shimZ,     'String', num2str(evalin('base', 'HW.MagnetShim(3)*1e3'), '%.4f'));

%% enable buttons
set(handles.pushbutton_measureFID,      'Enable', 'on');
set(handles.checkbox_autoRunFID,        'Enable', 'on');
set(handles.pushbutton_frequencySweep,  'Enable', 'on');
set(handles.pushbutton_P90,             'Enable', 'on');
set(handles.pushbutton_shim,            'Enable', 'on');
set(handles.pushbutton_frequencyFID,    'Enable', 'on');

end


function pushbutton_loadHW_Callback(hObject, eventdata, handles)
%% Executes on button press in pushbutton_loadHW.
% hObject    handle to pushbutton_loadHW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% update values in GUI
set(handles.edit_frequency, 'String', num2str(evalin('base', 'HW.fLarmor'), '%.0f'));
set(handles.edit_t90,       'String', num2str(evalin('base', 'HW.tFlip90Def*1e6'), '%.0f'));
set(handles.edit_shimX,     'String', num2str(evalin('base', 'HW.MagnetShim(1)*1e3'), '%.4f'));
set(handles.edit_shimY,     'String', num2str(evalin('base', 'HW.MagnetShim(2)*1e3'), '%.4f'));
set(handles.edit_shimZ,     'String', num2str(evalin('base', 'HW.MagnetShim(3)*1e3'), '%.4f'));

%% enable buttons
set(handles.pushbutton_measureFID,      'Enable', 'on');
set(handles.checkbox_autoRunFID,        'Enable', 'on');
set(handles.pushbutton_frequencySweep,  'Enable', 'on');
set(handles.pushbutton_P90,             'Enable', 'on');
set(handles.pushbutton_shim,            'Enable', 'on');
set(handles.pushbutton_frequencyFID,    'Enable', 'on');


end


function slider_frequency_CreateFcn(hObject, eventdata, handles)
%% Executes during object creation, after setting all properties.
% hObject    handle to slider_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
  set(hObject, 'BackgroundColor', [.9 .9 .9]);
end

%% add listener to slider
addlistener(hObject, 'ContinuousValueChange', @slider_frequency_OnSlide);

end

function slider_frequency_Callback(hObject, eventdata, handles)
%% Executes after movement of frequency slider
% hObject    handle to slider_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get config parameters (better use config file?)
settings.sliderRange = 5e6;
settings.sliderOldValueName = 'sliderFrequencyOldValue';
settings.editCallback = @edit_frequency_Callback;
textbox.format = '%.0f';
textbox.unitFactor = 1;

%% update GUI and system
hEdit = handles.edit_frequency;

sliderCallbackHelper(hObject, hEdit, eventdata, handles, settings, textbox);

end

function slider_frequency_OnSlide(hObject, eventdata)
%% Executes on dragging of slider for frequency
% hObject    handle to slider_shimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

%% get handles from figure
handles = guidata(hObject);

%%
% This callback is executed once when the slider is clicked after it was
% dragged.
% Workaround: Only execute the second time.
if handles.data.sliderFirstChange
  handles.data.sliderFirstChange = false;
  guidata(hObject, handles);
  return
end

%% get config parameters (better use config file?)
settings.sliderRange = 5e6;
settings.sliderOldValueName = 'sliderFrequencyOldValue';
settings.editCallback = @edit_frequency_Callback;
textbox.format = '%.0f';
textbox.unitFactor = 1;

%% update GUI and system
hEdit = handles.edit_frequency;
hSlider = handles.slider_frequency;

sliderOnSlideHelper(hSlider, hEdit, eventdata, handles, settings, textbox);

end


function edit_frequency_CreateFcn(hObject, eventdata, handles)
%% Executes during object creation, after setting all properties.
% hObject    handle to edit_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
  set(hObject, 'BackgroundColor', 'white');
end

end

function edit_frequency_Callback(hObject, eventdata, handles)
%% Executes when edit box looses focus
% hObject    handle to edit_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get config parameters (better use config file)
threshold.name = 'frequency';
threshold.min = 10e6;
threshold.max = 25e6;
threshold.unit = 'MHz'; % for error message
threshold.unitFactor = 1e-6; % for error message

%% read value from GUI
newValue = str2double(get(hObject, 'String'));

%% check thresholds
updatedValue = checkThreshold(newValue, threshold);

if updatedValue ~= newValue
  set(hObject, 'String', num2str(updatedValue, '%.0f'));
end

%% set variables in base workspace
HW = evalin('base', 'HW'); % get current variable from base
HW.fLarmor = newValue;
assignin('base', 'HW', HW); % write updated variable to base

end


function slider_t90_CreateFcn(hObject, eventdata, handles)
%% Executes during object creation, after setting all properties.
% hObject    handle to slider_t90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
  set(hObject, 'BackgroundColor', [.9 .9 .9]);
end

%% add listener to slider
addlistener(hObject, 'ContinuousValueChange', @slider_t90_OnSlide);

end

function slider_t90_Callback(hObject, eventdata, handles)
%% Executes after movement of slider for 90 degrees pulse length
% hObject    handle to slider_t90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get config parameters (better use config file?)
settings.sliderRange = 100;
settings.sliderOldValueName = 'sliderShimXOldValue';
settings.editCallback = @edit_t90_Callback;
textbox.format = '%.2f';
textbox.unitFactor = 1e-6;

%% update GUI and system
hEdit = handles.edit_t90;

sliderCallbackHelper(hObject, hEdit, eventdata, handles, settings, textbox);

end

function slider_t90_OnSlide(hObject, eventdata)
%% Executes on dragging of slider for 90 degrees pulse length
% hObject    handle to slider_shimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

%% get handles from figure
handles = guidata(hObject);

%%
% This callback is executed once when the slider is clicked after it was
% dragged.
% Workaround: Only execute the second time.
if handles.data.sliderFirstChange
  handles.data.sliderFirstChange = false;
  guidata(hObject, handles);
  return
end

%% get config parameters (better use config file?)
settings.sliderRange = 100;
settings.sliderOldValueName = 'sliderShimXOldValue';
settings.editCallback = @edit_t90_Callback;
textbox.format = '%.2f';
textbox.unitFactor = 1e-6;

%% update GUI and system
hEdit = handles.edit_t90;
hSlider = handles.slider_t90;

sliderOnSlideHelper(hSlider, hEdit, eventdata, handles, settings, textbox);

end

function edit_t90_CreateFcn(hObject, eventdata, handles)
%% Executes during object creation, after setting all properties.
% hObject    handle to edit_t90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
  set(hObject, 'BackgroundColor', 'white');
end

end

function edit_t90_Callback(hObject, eventdata, handles)
%% Executes when edit box looses focus
% hObject    handle to edit_t90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get config parameters (better use config file)
threshold.name = sprintf('time for 90%c puls', 176);
threshold.min = 1;
threshold.max = 1000;
threshold.unit = [char(181) 's']; % for error message
threshold.unitFactor = 1; % for error message

%% read value from GUI
newValue = str2double(get(hObject, 'String'));

%% check thresholds
updatedValue = checkThreshold(newValue, threshold);

if updatedValue ~= newValue
  set(hObject, 'String', num2str(updatedValue, '%.0f'))
end

%% set variables in base workspace
% HW = evalin('base', 'HW'); % get current variable from base
% HW.fLarmor = newValue;
% assignin('base', 'HW', HW); % write updated variable to base

end


function slider_shimX_CreateFcn(hObject, eventdata, handles)
%% Executes during object creation, after setting all properties.
% hObject    handle to slider_shimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
  set(hObject, 'BackgroundColor', [.9 .9 .9]);
end

%% add listener to slider
addlistener(hObject, 'ContinuousValueChange', @slider_shimX_OnSlide);

end

function slider_shimX_Callback(hObject, eventdata, handles)
%% Executes after movement of slider for shim in X
% hObject    handle to slider_shimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get config parameters (better use config file?)
settings.sliderRange = 1e-1;
settings.sliderOldValueName = 'sliderShimXOldValue';
settings.editCallback = @edit_shimX_Callback;
textbox.format = '%.4f';
textbox.unitFactor = 1e-3;

%% update GUI and system
hEdit = handles.edit_shimX;

sliderCallbackHelper(hObject, hEdit, eventdata, handles, settings, textbox);

end

function slider_shimX_OnSlide(hObject, eventdata)
%% Executes on dragging of slider for shim in X
% hObject    handle to slider_shimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

%% get handles from figure
handles = guidata(hObject);

%%
% This callback is executed once when the slider is clicked after it was
% dragged.
% Workaround: Only execute the second time.
if handles.data.sliderFirstChange
  handles.data.sliderFirstChange = false;
  guidata(hObject, handles);
  return
end

%% get config parameters (better use config file?)
settings.sliderRange = 1e-1;
settings.sliderOldValueName = 'sliderShimXOldValue';
settings.editCallback = @edit_shimX_Callback;
textbox.format = '%.4f';
textbox.unitFactor = 1e-3;

%% update GUI and system
hEdit = handles.edit_shimX;
hSlider = handles.slider_shimX;

sliderOnSlideHelper(hSlider, hEdit, eventdata, handles, settings, textbox);

end


function edit_shimX_CreateFcn(hObject, eventdata, handles)
%% Executes during object creation, after setting all properties.
% hObject    handle to edit_shimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
  set(hObject, 'BackgroundColor', 'white');
end

end

function edit_shimX_Callback(hObject, eventdata, handles)
%% Executes when edit box looses focus
% hObject    handle to edit_shimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get config parameters (better use config file)
threshold.name = 'Shim X';
threshold.min = -10;
threshold.max = 10;
threshold.unit = 'mT/m'; % for error message
threshold.unitFactor = 1; % for error message

%% read value from GUI
newValue = str2double(get(hObject, 'String'));

%% check thresholds
updatedValue = checkThreshold(newValue, threshold);

if updatedValue ~= newValue
  set(hObject, 'String', num2str(updatedValue, '%.4f'));
end

%% set variables in base workspace
HW = evalin('base', 'HW'); % get current variable from base
HW.MagnetShim(1) = newValue/1e3;
assignin('base', 'HW', HW); % write updated variable to base

end


function slider_shimY_CreateFcn(hObject, eventdata, handles)
%% Executes during object creation, after setting all properties.
% hObject    handle to slider_shimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
  set(hObject, 'BackgroundColor', [.9 .9 .9]);
end

%% add listener to slider
addlistener(hObject, 'ContinuousValueChange', @slider_shimY_OnSlide);

end

function slider_shimY_Callback(hObject, eventdata, handles)
%% Executes after movement of slider for shim in Y
% hObject    handle to slider_shimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get config parameters (better use config file?)
settings.sliderRange = 1e-1;
settings.sliderOldValueName = 'sliderShimYOldValue';
settings.editCallback = @edit_shimY_Callback;
textbox.format = '%.4f';
textbox.unitFactor = 1e-3;

%% update GUI and system
hEdit = handles.edit_shimY;

sliderCallbackHelper(hObject, hEdit, eventdata, handles, settings, textbox);

end

function slider_shimY_OnSlide(hObject, eventdata)
%% Executes on dragging of slider for shim in Y
% hObject    handle to slider_shimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

%% get handles from figure
handles = guidata(hObject);

%%
% This callback is executed once when the slider is clicked after it was
% dragged.
% Workaround: Only execute the second time.
if handles.data.sliderFirstChange
  handles.data.sliderFirstChange = false;
  guidata(hObject, handles);
  return
end

%% get config parameters (better use config file?)
settings.sliderRange = 1e-1;
settings.sliderOldValueName = 'sliderShimYOldValue';
settings.editCallback = @edit_shimY_Callback;
textbox.format = '%.4f';
textbox.unitFactor = 1e-3;

%% update GUI and system
hEdit = handles.edit_shimY;
hSlider = handles.slider_shimY;

sliderOnSlideHelper(hSlider, hEdit, eventdata, handles, settings, textbox);

end

function edit_shimY_CreateFcn(hObject, eventdata, handles)
%% Executes during object creation, after setting all properties.
% hObject    handle to edit_shimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
  set(hObject, 'BackgroundColor', 'white');
end

end

function edit_shimY_Callback(hObject, eventdata, handles)
%% Executes when edit box looses focus
% hObject    handle to edit_shimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get config parameters (better use config file)
threshold.name = 'Shim Y';
threshold.min = -10;
threshold.max = 10;
threshold.unit = 'mT/m'; % for error message
threshold.unitFactor = 1; % for error message

%% read value from GUI
newValue = str2double(get(hObject, 'String'));

%% check thresholds
updatedValue = checkThreshold(newValue, threshold);

if updatedValue ~= newValue
  set(hObject, 'String', num2str(updatedValue, '%.4f'));
end

%% set variables in base workspace
HW = evalin('base', 'HW'); % get current variable from base
HW.MagnetShim(2) = newValue/1e3;
assignin('base', 'HW', HW); % write updated variable to base

end

function slider_shimZ_CreateFcn(hObject, eventdata, handles)
%% Executes during object creation, after setting all properties.
% hObject    handle to slider_shimZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
  set(hObject, 'BackgroundColor', [.9 .9 .9]);
end

%% add listener to slider
addlistener(hObject, 'ContinuousValueChange', @slider_shimZ_OnSlide);

end

function slider_shimZ_Callback(hObject, eventdata, handles)
%% Executes after movement of slider for shim in Z
% hObject    handle to slider_shimZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get config parameters (better use config file?)
settings.sliderRange = 1e-1;
settings.sliderOldValueName = 'sliderShimZOldValue';
settings.editCallback = @edit_shimZ_Callback;
textbox.format = '%.4f';
textbox.unitFactor = 1e-3;

%% update GUI and system
hEdit = handles.edit_shimZ;

sliderCallbackHelper(hObject, hEdit, eventdata, handles, settings, textbox);

end

function slider_shimZ_OnSlide(hObject, eventdata)
%% Executes on dragging of slider for shim in Z
% hObject    handle to slider_shimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

%% get handles from figure
handles = guidata(hObject);

%%
% This callback is executed once when the slider is clicked after it was
% dragged.
% Workaround: Only execute the second time.
if handles.data.sliderFirstChange
  handles.data.sliderFirstChange = false;
  guidata(hObject, handles);
  return
end

%% get config parameters (better use config file?)
settings.sliderRange = 1e-1;
settings.sliderOldValueName = 'sliderShimZOldValue';
settings.editCallback = @edit_shimZ_Callback;
textbox.format = '%.4f';
textbox.unitFactor = 1e-3;

%% update GUI and system
hEdit = handles.edit_shimZ;
hSlider = handles.slider_shimZ;

sliderOnSlideHelper(hSlider, hEdit, eventdata, handles, settings, textbox);

end

function edit_shimZ_CreateFcn(hObject, eventdata, handles)
%% Executes during object creation, after setting all properties.
% hObject    handle to edit_shimZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
  set(hObject, 'BackgroundColor', 'white');
end

end

function edit_shimZ_Callback(hObject, eventdata, handles)
%% Executes when edit box looses focus
% hObject    handle to edit_shimZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get config parameters (better use config file)
threshold.name = 'Shim Z';
threshold.min = -10;
threshold.max = 10;
threshold.unit = 'mT/m'; % for error message
threshold.unitFactor = 1; % for error message

%% read value from GUI
newValue = str2double(get(hObject, 'String'));

%% check thresholds
updatedValue = checkThreshold(newValue, threshold);

if updatedValue ~= newValue
  set(hObject, 'String', num2str(updatedValue, '%.4f'));
end

%% set variables in base workspace
HW = evalin('base', 'HW'); % get current variable from base
HW.MagnetShim(3) = newValue/1e3;
assignin('base', 'HW', HW); % write updated variable to base

end


function checkbox_realImag_Callback(hObject, eventdata, handles)
%% Executes on button press in checkbox_realImag.
% hObject    handle to checkbox_realImag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_realImag

end


function checkbox_plotPhase_Callback(hObject, eventdata, handles)
%% Executes on button press in checkbox_plotPhase.
% hObject    handle to checkbox_plotPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plotPhase
end


function pushbutton_measureFID_Callback(hObject, eventdata, handles)
%% Executes on button press in pushbutton_measureFID.
% hObject    handle to pushbutton_measureFID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% start measurement
showFID(handles)

end


function checkbox_autoRunFID_Callback(hObject, eventdata, handles)
%% Executes on button press in checkbox_autoRunFID.
% hObject    handle to checkbox_autoRunFID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% start measurement
if get(hObject, 'Value')
  % fprintf('starting showFID at %s\n', datestr(now));
  showFID(handles);
  % dbstack
  % fprintf('stopped showFID at %s\n', datestr(now));
end

end


function edit_repRateFID_CreateFcn(hObject, eventdata, handles)
%% Executes during object creation, after setting all properties.
% hObject    handle to edit_repRateFID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
  set(hObject, 'BackgroundColor', 'white');
end

end

function edit_repRateFID_Callback(hObject, eventdata, handles)
%% Executes when edit box looses focus
% hObject    handle to edit_repRateFID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_repRateFID as text
%        str2double(get(hObject,'String')) returns contents of edit_repRateFID as a double

end


function uitable_hold_CellEditCallback(hObject, eventdata, handles)
%% Executes when entered data in editable cell(s) in uitable_hold.
% hObject    handle to uitable_hold (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%  Indices: row and column indices of the cell(s) edited
%  PreviousData: previous data for the cell(s) edited
%  EditData: string(s) entered by the user
%  NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%  Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%% check whether last column was selected
data = get(hObject, 'Data');
lastColIndex = size(data, 2);
if isempty(eventdata.Indices) || eventdata.Indices(2) ~= lastColIndex
  return
end
isSelected = [data{:,lastColIndex}];
selectedRow = find(isSelected, 1, 'first');
if isempty(selectedRow)
  return
end

%% reset selector in last column
data{isSelected,lastColIndex} = false;
set(handles.uitable_hold, 'Data', data);

%% check whether there is data in this row
if any( cellfun(@isempty, data(isSelected,:)) )
  return
end

%% set values in GUI and call callback manually
set(handles.edit_frequency, 'String', num2str(data{selectedRow,2}, '%.0f'));
edit_frequency_Callback(handles.edit_frequency, [], handles);
set(handles.edit_t90, 'String', num2str(data{selectedRow,3}, '%.0f'));
edit_t90_Callback(handles.edit_t90, [], handles);
set(handles.edit_shimX, 'String', num2str(data{selectedRow,4}, '%.4f'));
edit_shimX_Callback(handles.edit_shimX, [], handles);
set(handles.edit_shimY, 'String', num2str(data{selectedRow,5}, '%.4f'));
edit_shimY_Callback(handles.edit_shimY, [], handles);
set(handles.edit_shimX, 'String', num2str(data{selectedRow,6}, '%.4f'));
edit_shimZ_Callback(handles.edit_shimZ, [], handles);

end


function pushbutton_frequencySweep_Callback(hObject, eventdata, handles)
%% Executes on button press in pushbutton_frequencySweep
% hObject    handle to pushbutton_frequencySweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% parameters
maxtime =   1;              % Maximal time in seconds since the last Find_Frequency. i.e. if set to 10, find freq will only be run, if 10 seconds have passed since the last search
span    =   4e6;            % searching span of frequency around fLarmor (HW.B0*HW.FindFrequencyGamma/2/pi)
doplot  =   1;              % Plot sequence and data
tPulse90=   str2double(get(handles.edit_t90, 'String'))/1e6;  % duration of the 90 degrees pulse

%% check auto-run and disable if necessary
isAutoRun = get(handles.checkbox_autoRunFID, 'Value');
if isAutoRun
  set(handles.checkbox_autoRunFID, 'Value', false);
  pause(1);
end

%% get variables from base workspace
HW = evalin('base', 'HW');
mySave = evalin('base', 'mySave');

%% actual sweep
[HW, mySave]  = Find_Frequency_Sweep(HW, mySave, maxtime, span, doplot, tPulse90, 101);

disp(['New frequency: ' num2str(HW.fLarmor) ' Hz.']);

%% update variables in base workspace
assignin('base', 'HW', HW);
assignin('base', 'mySave', mySave);

%% update value in GUI
set(handles.edit_frequency, 'String', num2str(HW.fLarmor, '%.0f'));

msgbox(sprintf('Frequency updated to %.0f Hz.', HW.fLarmor));

%% re-set auto-run
if isAutoRun
  pause(1);
  set(handles.checkbox_autoRunFID, 'Value', true);
  checkbox_autoRunFID_Callback(handles.checkbox_autoRunFID, [], handles);
end

end


function pushbutton_frequencyFID_Callback(hObject, eventdata, handles)
%% Executes on button press in pushbutton_frequencyFID
% hObject    handle to pushbutton_frequencyFID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% parameter
maxtime         =   00;           % Maximal time in seconds since the last Find_Frequency.
iterations      =   2;            % number of iterations used for freqeuency determination
tAQ             =   0.5e-3;       % acquisition time

%% check auto-run and disable if necessary
isAutoRun = get(handles.checkbox_autoRunFID, 'Value');
if isAutoRun
  set(handles.checkbox_autoRunFID, 'Value', false);
  pause(0.5);
end

%% get variables from base workspace
HW = evalin('base', 'HW');
mySave = evalin('base', 'mySave');

%% actual measurements
[HW, mySave] = Find_Frequency_FID(HW, mySave, maxtime, iterations, tAQ);

disp(['New frequency: ' num2str(HW.fLarmor) ' Hz.']);

%% update variables in base workspace
assignin('base', 'HW', HW);
assignin('base', 'mySave', mySave);

%% update value in GUI
set(handles.edit_frequency, 'String', num2str(HW.fLarmor, '%.0f'));

msgbox(sprintf('Frequency updated to %.0f Hz.', HW.fLarmor));

%% re-set auto-run
if isAutoRun
  pause(2);
  set(handles.checkbox_autoRunFID, 'Value', true);
  checkbox_autoRunFID_Callback(handles.checkbox_autoRunFID, [], handles);
end

end


function pushbutton_P90_Callback(hObject, eventdata, handles)
%% Executes on button press in pushbutton_P90
% hObject    handle to pushbutton_P90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% parameters
maxtime         =   0;                  % Maximal time in seconds since the last Find_Pulseduration.
doplot          =   1;                  % Plot sequence and data
iterations      =   2;                  % number of iterations used for pulselenght determination / 0: Set the tPulse90 without searching.
tPulse90        =   str2double(get(handles.edit_t90, 'String'))/1e6;  % duration of the 90 degrees pulse
T1              =   0.5;                % estimated T1 value of the used sample

%% check auto-run and disable if necessary
isAutoRun = get(handles.checkbox_autoRunFID, 'Value');
if isAutoRun
  set(handles.checkbox_autoRunFID, 'Value', false);
  pause(1);
end

%% get variables from base workspace
HW = evalin('base', 'HW');
mySave = evalin('base', 'mySave');

%% actual measurements
[HW, mySave] = Find_PulseDuration(HW, mySave, maxtime, doplot, iterations, tPulse90, T1);

%LoadSystem;                 % load new parameters from file.
disp(['New pulse duration: ' num2str(HW.tFlip90Def) ' s.']);

%% update variables in base workspace
assignin('base', 'HW', HW);
assignin('base', 'mySave', mySave);

%% update value in GUI
set(handles.edit_t90, 'String', num2str(HW.tFlip90Def*1e6, '%.0f'));

msgbox(sprintf('Time for 90%c puls updated to %.0f %cs.', 176, HW.tFlip90Def*1e6, 181));

%% re-set auto-run
if isAutoRun
  set(handles.checkbox_autoRunFID, 'Value', true);
  checkbox_autoRunFID_Callback(handles.checkbox_autoRunFID, [], handles);
end

end


function pushbutton_shim_Callback(hObject, eventdata, handles)
%% Executes on button press in pushbutton_shim
% hObject    handle to pushbutton_shim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% parameters
maxtime             = 10;                 % Minimal time in seconds since the last Find_Frequency_Sweep.
doplot              = 1;                  % Plot sequence and data
Seq.iterations      = 100;                % number of iterations used for shim determination
Seq.tEcho           = 50e-3;              % echo time
Seq.T1              = 0.1;                % estimated T1 value of the used sample (generates a sleep(Seq.T1))
Seq.ShimStart       = [0,0,0,0];          % start with these shim values
Seq.ShimStep        = [0.02,0.02,0.02,0.000001];   % stepwidth
Seq.nEchos          = 1;                  % use FID or Echo for shimming (0 for Fid 1,2,3... for nEcho)

Seq.useSliceSelect  = 0;                  % use slice gradient
SliceSelect.alfa    = 0.5;                % around x Achse
SliceSelect.phi     = 0.5;                % around y Achse
SliceSelect.theta   = 0.0;                % around z Achse
% SliceSelect.thickness= 0.001;             % Slice thickness

%% check auto-run and disable if necessary
isAutoRun = get(handles.checkbox_autoRunFID, 'Value');
if isAutoRun
  set(handles.checkbox_autoRunFID, 'Value', false);
  pause(1);
end

%% get variables from base workspace
HW = evalin('base', 'HW');
mySave = evalin('base', 'mySave');

%% actual measurements
[HW, mySave] = Find_Shim(HW, mySave, maxtime, doplot, Seq );

disp(['New shim values: ' num2str(HW.MagnetShim) ' T/m.']);

%% update variables in base workspace
assignin('base', 'HW', HW);
assignin('base', 'mySave', mySave);

%% update value in GUI
set(handles.edit_shimX, 'String', num2str(HW.MagnetShim(1)*1e3, '%.4f'));
set(handles.edit_shimY, 'String', num2str(HW.MagnetShim(2)*1e3, '%.4f'));
set(handles.edit_shimZ, 'String', num2str(HW.MagnetShim(3)*1e3, '%.4f'));

msgbox([sprintf('Values for shim updated to:\n') ...
  sprintf('\t%.4f mT/m\n', HW.MagnetShim(:)*1e3)]);

%% re-set auto-run
if isAutoRun
  set(handles.checkbox_autoRunFID, 'Value', true);
  checkbox_autoRunFID_Callback(handles.checkbox_autoRunFID, [], handles);
end

end


function pushbutton_resetHistory_Callback(hObject, eventdata, handles)
%% Executes on button press in pushbutton_resetHistory
% hObject    handle to pushbutton_resetHistory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% clear history data from handles
handles.data.fLarmor = [];

%% save updated handles in GUI
guidata(hObject, handles);

%% clear plot
cla(handles.axes_fLarmor);

end



%% inline functions
function newValue = checkThreshold(newValue, threshold)
%% check whether value is within range and display message

if ~isnumeric(newValue) || ~isfinite(newValue)
  errordlg(sprintf('Error setting %s', threshold.name));
  return
elseif newValue < threshold.min
  errordlg(...
    sprintf('Set %s is lower than threshold. Setting to %f %s.', ...
    threshold.name, threshold.min*threshold.unitFactor, threshold.unit));
  newValue = threshold.min;
elseif newValue > threshold.max
  errordlg(...
    sprintf('Set %s is higher than threshold. Setting to %f %s.', ...
    threshold.name, threshold.max*threshold.unitFactor, threshold.unit));
  newValue = threshold.max;
end


end

function showFID(handles)

ax2 = axes('Position', get(handles.axes_FID, 'Position'), ...
  'XAxisLocation','top',...
  'YAxisLocation','right',...
  'Color','none',...
  'XColor','k','YColor','k');

evalin('base', 'HW.tRepInit = 300e-3;')


while true
  tStart = tic;

  %% get variables from base workspace
  HW = evalin('base', 'HW');
  Seq = evalin('base', 'Seq');
  AQ = evalin('base', 'AQ');
  TX = evalin('base', 'TX');
  Grad = evalin('base', 'Grad');

  %% set parameters
  % Parameters used for timing calculations
  Seq.p90         = str2double(get(handles.edit_t90, 'String'))*1e-6;  % duration 1. TX pulse
  Seq.plotSeq     = [];                           % plot sequence off

  % RF transmission parameters
  TX.Start        = 0;                            % start time of rf-pulse
  TX.Duration     = Seq.p90;                      % duration of rf-pulse
  TX.Frequency    = HW.fLarmor;                   % frequency of rf-pulse
  TX.Phase        = 0;                            % phase of rf-pulse
  TX.Amplitude    = HW.TX.AmpDef; % amplitude in Volt

  % Acquisition parameters
  AQ.Start        = 100e-6 + Seq.p90;             % acquisition start time
  AQ.fSample      = 50e3;                         % sampling rate AQ window
  AQ.nSamples     = 1024*2;                       % number of samples in AQ window
  AQ.Frequency    = HW.fLarmor;                   % frequency of AQ window
  AQ.Phase        = 0;                            % phase of AQ window
  AQ.Gain         = HW.RX.GainDef; % HW.RX.Amplitude2Uin / 20e-3;  % maximum input voltage

  % Sequence parameters
  Seq.tRep        = AQ.Start + AQ.nSamples / AQ.fSample + 1e-3;                       % repetition time

  %% update variables in base workspace
  assignin('base', 'Seq', Seq);
  assignin('base', 'TX', TX);
  assignin('base', 'AQ', AQ);

  %% Start measurement
  [ Raw, SeqOut, data, data_1D ] = set_sequence(HW, Seq, AQ, TX, Grad);

  evalin('base', 'HW.tRepInit = 0.01;');

  %% check hold states
  tableData = get(handles.uitable_hold, 'Data');
  isHoldTable = cell2mat(tableData(1:3,1))';
  isHoldData = [handles.data.plot(:).hold];
  dataToSave = find(isHoldTable & ~isHoldData, 1, 'first');
  handles = guidata(handles.uitable_hold);
  if ~isempty(dataToSave)
    handles.data.plot(dataToSave).hold = true;
    handles.data.plot(dataToSave).data = data.data*data.Amplitude2Uin(1);
    handles.data.plot(dataToSave).fLarmor = HW.fLarmor;
    handles.data.plot(dataToSave).p90 = Seq.p90;
    handles.data.plot(dataToSave).shimX = HW.MagnetShim(1);
    handles.data.plot(dataToSave).shimY = HW.MagnetShim(2);
    handles.data.plot(dataToSave).shimZ = HW.MagnetShim(3);
    tableData{dataToSave,2} = handles.data.plot(dataToSave).fLarmor;
    tableData{dataToSave,3} = handles.data.plot(dataToSave).p90*1e6;
    tableData{dataToSave,4} = handles.data.plot(dataToSave).shimX*1e3;
    tableData{dataToSave,5} = handles.data.plot(dataToSave).shimY*1e3;
    tableData{dataToSave,6} = handles.data.plot(dataToSave).shimZ*1e3;
    set(handles.uitable_hold, 'Data', tableData);
  end
  dataToDelete = ~isHoldTable & isHoldData;
  if any(dataToDelete)
    handles.data.plot(dataToDelete).hold = false;
    handles.data.plot(dataToDelete).data = [];
    tableData{dataToDelete,2} = '';
    tableData{dataToDelete,3} = '';
    tableData{dataToDelete,4} = '';
    tableData{dataToDelete,5} = '';
    tableData{dataToDelete,6} = '';
    set(handles.uitable_hold, 'Data', tableData);
  end
  guidata(handles.edit_frequency, handles);

  %% plot saved data
  lineColors = {'r', 'g', 'k'};
  cla(handles.axes_FID);
  cla(handles.axes_FID2);
  linkaxes([handles.axes_FID handles.axes_FID2], 'x')
  % set(handles.axes_FID2, 'XTick', [])
  set(handles.axes_FID2, 'YTickMode', 'auto')
  set(handles.axes_FID2, 'Color', 'none')
  % set(handles.axes_FID, 'Color', 'none')
  for iPlot = 1:length(isHoldData)
    if handles.data.plot(iPlot).hold
      plot(handles.axes_FID, ...
        data.time_of_tRep, abs(handles.data.plot(iPlot).data), ...
        ['-' lineColors{iPlot}]);
      hold(handles.axes_FID, 'on');
    end
  end

  %% Plot results
  plot(handles.axes_FID, data.time_of_tRep, abs(data.data)*data.Amplitude2Uin(1));
  title(handles.axes_FID, 'Acquired signal');
  ylabel(handles.axes_FID, 'Amplitude in V');
  xlabel(handles.axes_FID, 'Time in s');
  maxFIDPlot = max(abs(data.data)*data.Amplitude2Uin(1));
  maxFIDPlot = ceil(maxFIDPlot/10^(round(log10(maxFIDPlot))-1))*10^(round(log10(maxFIDPlot))-1);
  set(handles.axes_FID, 'YLim', [0 maxFIDPlot]);

  %% plot complex components
  if get(handles.checkbox_realImag, 'Value')
    hold(handles.axes_FID, 'all');
    plot(handles.axes_FID, data.time_of_tRep, real(data.data)*data.Amplitude2Uin(1), 'Color', [0 .5 0]);
    plot(handles.axes_FID, data.time_of_tRep, imag(data.data)*data.Amplitude2Uin(1), 'Color', [.6 .3 0]);
    set(handles.axes_FID, 'YLim', [-maxFIDPlot maxFIDPlot]);
  end

  %% plot phase of signal
  if get(handles.checkbox_plotPhase, 'Value')
    % FIDIdx = find(data.time_of_tRep < 2.5e-3);
    FIDIdxLast = find(abs(data.data) < 0.02*abs(data.data(1)), 1, 'first');
    if isempty(FIDIdxLast)
      FIDIdx = 1:SeqOut.AQ.nSamples(1)-1;
    else
      FIDIdx = 1:(FIDIdxLast-1);
    end
    fOffsetFIDs = diff(unwrap(angle(data.data(1:SeqOut.AQ.nSamples(1),:)))) * SeqOut.AQ.fSample(1)/2/pi;  % get frequency offset between each sample
    fOffsetFID = mean(fOffsetFIDs(FIDIdx), 1);                                                                   % get the mean frequeny offset of each AQ window

    data.dataPC = data.data(:,1) .* ...
      exp(-1i * 2 * pi * fOffsetFID * (data.time_of_tRep(:) - data.time_of_tRep(1))); % correct the linear phase error
    PhaseOffsetRad = mean(unwrap(angle(data.dataPC(1:10))));    % get the phase offset
    data.dataPCo = data.dataPC(:) .* exp(-1i * PhaseOffsetRad);    % correct the phase offset
    data.dataPCo = data.dataPC(:);
    % hold(handles.axes_FID(1), 'on');

    line(data.time_of_tRep(FIDIdx), unwrap(angle(data.dataPCo(FIDIdx))), 'Color', 'black', 'Parent', handles.axes_FID2);
    set(handles.axes_FID2, 'YAxisLocation', 'right', 'YLim', [-1 1] + PhaseOffsetRad);
  end
  grid(handles.axes_FID, 'on');
  drawnow

  %% plot Larmor frequency over runs
  % get Larmor frequency from spectral response
  absFFT1_data = abs(data.fft1_data);
  fLarmor = data.f_fft1_data(absFFT1_data==max(absFFT1_data));
  % save Larmor frequency in handles
  handles = guidata(handles.uitable_hold);
  if ~isfield(handles.data, 'fLarmor')
    handles.data.fLarmor = fLarmor(1);
  else
    handles.data.fLarmor(end+1) = fLarmor(1);
  end
  guidata(handles.uitable_hold, handles);
  % actual plot
  plot(handles.axes_fLarmor, 1:length(handles.data.fLarmor), handles.data.fLarmor - handles.data.fLarmor(1));
  title(handles.axes_fLarmor, 'Larmor frequency over runs');
  ylabel(handles.axes_fLarmor, sprintf('Larmor frequency - %.0f Hz', handles.data.fLarmor(1)));
  xlabel(handles.axes_fLarmor, 'Number of run');

  %% check whether loop is on
  if ~get(handles.checkbox_autoRunFID, 'Value')
    evalin('base', 'HW.tRepInit = 0.1;');
    return
  end

  %% add delay
  repRate = str2double(get(handles.edit_repRateFID, 'String'));
  pause(repRate - toc(tStart));
end

end


function handles = sliderCallbackHelper(hSlider, hEdit, eventdata, handles, settings, textbox)
%% helper function for all slider callbacks (executes after movement of slider)
%
%   handles = sliderCallbackHelper(hSlider, hEdit, eventdata, handles, settings, textbox)
%
% Updates the corresponding text fields and the settings in the HW
% structure.
%

%%
if handles.data.sliderChanged
  % Slider was dragged and therefore updated in sliderOnSlideHelper.
  % This Callback executes on Mouse button release. Do not change any
  % values nothing here.
  handles.data.sliderChanged = false;
else
  %% check whether there is something to be changed
  newSliderPos = get(hSlider, 'Value');
  if newSliderPos == 0.5 % continue only if there is something to change
    return
  end

  %% slider range
  if ismember('control', get(gcf, 'currentModifier'))
    % if control key is held, slider is more precise
    settings.sliderRange = settings.sliderRange / 10;
  end

  %% get new frequency
  currentValue = str2double(get(hEdit, 'String'));

  newValue = currentValue + ...
    (newSliderPos - 0.5) * settings.sliderRange;

  %% set value in text box
  set(hEdit, 'String', num2str(newValue, textbox.format));

  %% re-set slider to middle position
  set(hSlider, 'Value', 0.5);
  drawnow
end

%% re-set slider to middle position
set(hSlider, 'Value', 0.5);
handles.data.(settings.sliderOldValueName) = 0.5;
handles.data.sliderFirstChange = true;
guidata(hSlider, handles);
drawnow

%% callback of edit box is not executed automatically
settings.editCallback(hEdit, eventdata, handles);

end


function handles = sliderOnSlideHelper(hSlider, hEdit, eventdata, handles, settings, textbox)
%% helper function for all slider 'ContinuousValueChange' events (executes on dragging of slider)
%
%   handles = sliderOnSlideHelper(hSlider, hEdit, eventdata, handles, settings, textbox)
%
% Updates the corresponding text fields and the settings in the HW
% structure.
%


%% check whether there is something to be changed
newSliderPos = get(hSlider, 'Value');
% continue only if there is something to change
if newSliderPos == handles.data.(settings.sliderOldValueName)
  return
end

%% slider range
if ismember('control', get(handles.setupDevice, 'currentModifier'))
  % if control key is held, slider is more precise
  settings.sliderRange = settings.sliderRange / 10;
end

%% get new frequency
currentValue = str2double(get(hEdit, 'String'));

newValue = currentValue + ...
  (newSliderPos - handles.data.(settings.sliderOldValueName)) * settings.sliderRange;

%% save data to figure
handles.data.(settings.sliderOldValueName) = newSliderPos;
handles.data.sliderChanged = true;
guidata(hEdit, handles)

%% set value in text box
set(hEdit, 'String', num2str(newValue, textbox.format))
drawnow

%% callback of edit box is not executed automatically
settings.editCallback(hEdit, eventdata, handles);

end
