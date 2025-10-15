function varargout = GUI_getSNR(varargin)
% GUI_GETSNR MATLAB code for GUI_getSNR.fig
%      GUI_GETSNR, by itself, creates a new GUI_GETSNR or raises the existing
%      singleton*.
%
%      H = GUI_GETSNR returns the handle to a new GUI_GETSNR or the handle to
%      the existing singleton*.
%
%      GUI_GETSNR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GETSNR.M with the given input arguments.
%
%      GUI_GETSNR('Property','Value',...) creates a new GUI_GETSNR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_getSNR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_getSNR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_getSNR

% Last Modified by GUIDE v2.5 04-Oct-2016 14:28:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_getSNR_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_getSNR_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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

% --- Executes just before GUI_getSNR is made visible.
function GUI_getSNR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_getSNR (see VARARGIN)

% Choose default command line output for GUI_getSNR
handles.output = hObject;

handles.noiseValues = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_getSNR wait for user response (see UIRESUME)
% uiwait(handles.figure1);

end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_getSNR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end

% --- Executes on button press in pushbutton_measureSNR.
function pushbutton_measureSNR_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_measureSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[~,~,Fnoise] = get_SNR(handles.HW);

handles.noiseValues = [handles.noiseValues; Fnoise];
if ~ishghandle(hObject)
    return
end
guidata(hObject, handles)

% set value
set(handles.edit_SNR, 'String', num2str(Fnoise, '%.1f'))

if Fnoise > 30
    set(handles.edit_SNR, 'BackgroundColor', [255 0 0]/255)
elseif Fnoise > 6
    set(handles.edit_SNR, 'BackgroundColor', [255 204 0]/255)
else
    set(handles.edit_SNR, 'BackgroundColor', [0 255 0]/255)
end

% sliding averages:
avgNum = str2double(get(handles.edit_slideAvgNum, 'String'));

if avgNum >= length(handles.noiseValues)
    avgSNR = mean(handles.noiseValues);
else
    avgSNR = mean(handles.noiseValues(end-avgNum:end));
end

% set value
set(handles.edit_avgSNR, 'String', num2str(avgSNR, '%.1f'))

if avgSNR > 30
    set(handles.edit_avgSNR, 'BackgroundColor', [255 0 0]/255)
elseif avgSNR > 6
    set(handles.edit_avgSNR, 'BackgroundColor', [255 204 0]/255)
else
    set(handles.edit_avgSNR, 'BackgroundColor', [0 255 0]/255)
end



end

% --- Executes on button press in checkbox_autoSNR.
function checkbox_autoSNR_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_autoSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_autoSNR

isActive = get(hObject, 'Value');

while isActive
    tSNR_loop = tic;
    handles = guidata(hObject);
    pushbutton_measureSNR_Callback(hObject, eventdata, handles)
    if ~ishghandle(hObject)
        break;
    end
    isActive = get(hObject, 'Value');
    measureInterval = str2double(get(handles.edit_measureInterval, 'String'));
    pause(measureInterval-toc(tSNR_loop));
end

end


% --- Executes on button press in pushbutton_loadHW.
function pushbutton_loadHW_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadHW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~evalin('base', 'exist(''HW'', ''var'')')
    errordlg('HW not found. Run "LoadSystem" first!')
    return
end

handles.HW = evalin('base', 'HW');
guidata(hObject, handles)

% enable buttons and controls
set(handles.pushbutton_measureSNR, 'Enable', 'on')
set(handles.checkbox_autoSNR, 'Enable', 'on')
set(handles.pushbutton_measureAvgSNR, 'Enable', 'on')


end


% --- Executes on button press in pushbutton_measureAvgSNR.
function pushbutton_measureAvgSNR_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_measureAvgSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

avgNum = str2double(get(handles.edit_slideAvgNum, 'String'));
measureInterval = str2double(get(handles.edit_measureInterval, 'String'));

for iLoop = 1:avgNum
    tSNR_loop = tic;
    pushbutton_measureSNR_Callback(hObject, eventdata, handles)
    pause(measureInterval-toc(tSNR_loop));
end


end


function edit_slideAvgNum_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slideAvgNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_slideAvgNum as text
%        str2double(get(hObject,'String')) returns contents of edit_slideAvgNum as a double

slideNum = str2double(get(hObject, 'String'));

if isnan(slideNum) || (mod(slideNum, 1) ~= 0) || slideNum < 1
    errordlg('Value must be a positive integer!')
    set(hObject, 'String', '5')
    return
end

end

function edit_measureInterval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_measureInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_measureInterval as text
%        str2double(get(hObject,'String')) returns contents of edit_measureInterval as a double

measureIntervalStr = get(hObject, 'String');
measureIntervalStr(measureIntervalStr==',') = '.'; % replace "," with "."
measureInterval = str2double(measureIntervalStr);

if isnan(measureInterval) || measureInterval < 1
    errordlg('Value must be a positive scalar!')
    set(hObject, 'String', '1.0')
    return
end
set(hObject, 'String', measureIntervalStr)


end

% --- Executes during object creation, after setting all properties.
function edit_measureInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_measureInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end
