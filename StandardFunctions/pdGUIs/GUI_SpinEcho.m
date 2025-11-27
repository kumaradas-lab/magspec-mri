function varargout = GUI_SpinEcho(varargin)
% GUI_SPINECHO M-file for GUI_SpinEcho.fig
%      GUI_SPINECHO, by itself, creates a new GUI_SPINECHO or raises the existing
%      singleton*.
%
%      H = GUI_SPINECHO returns the handle to a new GUI_SPINECHO or the handle to
%      the existing singleton*.
%
%      GUI_SPINECHO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SPINECHO.M with the given input arguments.
%
%      GUI_SPINECHO('Property','Value',...) creates a new GUI_SPINECHO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_SpinEcho_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_SpinEcho_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% Edit the above text to modify the response to help GUI_SpinEcho

% Last Modified by GUIDE v2.5 23-Apr-2012 14:20:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_SpinEcho_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_SpinEcho_OutputFcn, ...
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


% --- Executes just before GUI_SpinEcho is made visible.
function GUI_SpinEcho_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_SpinEcho (see VARARGIN)

% Choose default command line output for GUI_SpinEcho
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using GUI_SpinEcho.
if strcmp(get(hObject,'Visible'),'off')
%     plot(rand(5));
end

if evalin('base', 'exist(''HW'', ''var'')')
  HW = evalin('base', 'HW');

  set(handles.edit_shimX, 'String', num2str(HW.MagnetShim(1)*1e3, '%.3f'));
  set(handles.edit_shimY, 'String', num2str(HW.MagnetShim(2)*1e3, '%.3f'));
  set(handles.edit_shimZ, 'String', num2str(HW.MagnetShim(3)*1e3, '%.3f'));
end

% UIWAIT makes GUI_SpinEcho wait for user response (see UIRESUME)
% uiwait(handles.figure_SpinEcho);

end


% --- Outputs from this function are returned to the command line.
function varargout = GUI_SpinEcho_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end


function figure_SpinEcho_SizeChangeFcn(hObject, eventdata, handles)
% executes on size change

% keep controls in upper right corner
figPos = get(hObject, 'Position');
controlsPos = get(handles.panel_controls_content, 'Position');
controlsPos(1:2) = figPos(3:4)-controlsPos(3:4);
set(handles.panel_controls_content, 'Position', controlsPos);

% keep popupmenu close to upper border
popupPos = get(handles.popupmenu_axes_large, 'Position');
popupPos(2) = figPos(4) - 32;
set(handles.popupmenu_axes_large, 'Position', popupPos);

% resize large axes to fit available space
axesPos = get(handles.axes_large, 'Position');
axesPos(3:4) = figPos(3:4) - [(1120-595), (732-620)];
set(handles.axes_large, 'Position', axesPos);

end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

end


% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure_SpinEcho)

end


% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure_SpinEcho,'Name') '?'],...
                     ['Close ' get(handles.figure_SpinEcho,'Name') '...'],...
                     'Yes', 'No', 'Yes');
if ~strcmp(selection, 'Yes')
    return;
end

delete(handles.figure_SpinEcho)

end


% --- Executes on selection change in popupmenu_axes_large.
function popupmenu_axes_large_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_axes_large (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_axes_large contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_axes_large
if ~handles.Running
    handles.refreshPlot=1;
    guidata(hObject, handles);
    pushbutton_update_Callback(hObject, eventdata, handles)
end

end


% --- Executes during object creation, after setting all properties.
function popupmenu_axes_large_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_axes_large (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

% set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});

end


function edit_p90_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p90 as text
%        str2double(get(hObject,'String')) returns contents of edit_p90 as a double

end


% --- Executes during object creation, after setting all properties.
function edit_p90_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit_p180_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p180 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p180 as text
%        str2double(get(hObject,'String')) returns contents of edit_p180 as a double

end


% --- Executes during object creation, after setting all properties.
function edit_p180_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p180 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit_tEcho_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tEcho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tEcho as text
%        str2double(get(hObject,'String')) returns contents of edit_tEcho as a double

end


% --- Executes during object creation, after setting all properties.
function edit_tEcho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tEcho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit_nEchoes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nEchoes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nEchoes as text
%        str2double(get(hObject,'String')) returns contents of edit_nEchoes as a double

end


% --- Executes during object creation, after setting all properties.
function edit_nEchoes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nEchoes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on selection change in popupmenu_axes_small.
function popupmenu_axes_small_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_axes_small (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_axes_small contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_axes_small
if ~handles.Running
    handles.refreshPlot=1;
    guidata(hObject, handles);
    pushbutton_update_Callback(hObject, eventdata, handles)
end

end


% --- Executes during object creation, after setting all properties.
function popupmenu_axes_small_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_axes_small (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in pushbutton_update.
function pushbutton_update_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'Update'),  handles.Update = 1;  end
handles.Update = 1;
if ~isfield(handles, 'Running'),  handles.Running = 0;  end
if ~isfield(handles, 'Cal'),     handles.Cal.Cal = 0;  end
if ~isfield(handles.Cal, 'Cal'), handles.Cal.Cal = 0;  end
%if isfield(handles, 'talker'), talker = handles.talker;  end
mySave = [];
if isfield(handles, 'mySave'),  mySave = handles.mySave;  end
if ~isfield(handles, 'refreshPlot'),  handles.refreshPlot = 0;  end
handles.Running = 1;
guidata(hObject, handles);
% FIXME: Support multiple MMRT devices
if isemptyfield(handles, 'iDevice'),  handles.iDevice = 1;  end

if evalin('base', 'exist(''HW'', ''var'')')
  HW = evalin('base', 'HW');
end
if exist('HW', 'var') && isa(HW, 'PD.HWClass')
  ResetStructs;
else
  if isfield(HW, 'MagnetShim')
    oldShim = HW.MagnetShim;
  else
    oldShim = [];
  end
  LoadSystem;
  if ~isempty(oldShim)
    HW.MagnetShim = oldShim;
  end
end
if evalin('base', 'exist(''mySave'', ''var'')')
  mySave = evalin('base', 'mySave');
end
handles.mySave = mySave;
%handles.talker = talker;
handles.HW = HW;
handles.Seq = Seq;
% handles.HW.tRepInit = 0.3;
% handles.HW.ReInit = 1;
handles.Seq.Cal = 0;
handles.Seq.plot = 0;
handles.Seq.plotTR = 0;
handles.Seq.firstTR = 1;
if ~isfield(handles, 'SeqOut'),  handles.SeqOut = [];  end
set(handles.pushbutton_shimX_minus, 'Enable', 'on');
set(handles.pushbutton_shimX_plus, 'Enable', 'on');
set(handles.pushbutton_shimY_minus, 'Enable', 'on');
set(handles.pushbutton_shimY_plus, 'Enable', 'on');
set(handles.pushbutton_shimZ_minus, 'Enable', 'on');
set(handles.pushbutton_shimZ_plus, 'Enable', 'on');
set(handles.pushbutton_save_shim, 'Enable', 'on');

while ((get(handles.checkbox_repetition, 'Value') && get(handles.radiobutton2, 'Value')) || handles.Update)
  handles.Seq.plotAllHandle = [];
  handles.Seq.plotFidHandle = [];
  handles.Seq.plotEchoesHandle = [];
  handles.Seq.plotMaxEchoesHandle = [];
  handles.Update = 0;
  handles.Seq.tOffset = 0;

  handles.Seq.DigitalIOEnable = 0;
  if handles.Seq.DigitalIOEnable
    handles.Seq.DigitalIO.SetTime = [-1e-6;  0;  1e-6 ];  % DigitalIO trigger output at center of excitation pulse
    handles.Seq.DigitalIO.SetValue = [   0;  1;     0 ];  % DigitalIO output 1 raising edge at center of excitation pulse
  end

  if ~isfield(handles, 'SeqOut'),handles.SeqOut.StartSequenceTime=now*24*3600;end
  if (get(handles.checkbox_repetition, 'Value') && ~handles.Seq.firstTR)
    if str2double(get(handles.edit_repetition_time, 'String'))*1e-3 > 2.5 || ...
      (get(handles.checkbox_calibrate, 'Value') && handles.HW.FindFrequencyPause > 2.5)
      % handles.Seq.StartSequenceTime = handles.SeqOut.StartSequenceTime + str2double(get(handles.edit_repetition_time, 'String'))*1e-3;
      handles.Seq.tOffset = handles.HW.tRepInit;
    end
    handles.Seq.StartSequenceTime = [];
    handles.Seq.Reinitialize = 0;
    handles.Seq.TimeFromLastSequence = str2double(get(handles.edit_repetition_time, 'String'))*1e-3 - sum(handles.SeqOut.tRep);
    handles.Seq.TimeToNextSequence = str2double(get(handles.edit_repetition_time, 'String'))*1e-3;
  else
    handles.Seq.StartSequenceTime = [];
    handles.Seq.Reinitialize = 1;
    handles.Seq.TimeFromLastSequence = [];
    handles.Seq.TimeToNextSequence = str2double(get(handles.edit_repetition_time, 'String'))*1e-3;
  end

  handles.Seq.firstTR = 0;
  if ~handles.refreshPlot
    handles.Seq.notCalibrate=double(~get(handles.checkbox_calibrate, 'Value'));
    if isnan(str2double(get(handles.edit_fLarmor, 'String')))
      if ~isa(handles.HW, 'PD.HWClass')
        % update with current settings from base workspace
        handles.HW = evalin('base', 'HW');
      end
      [handles.HW, handles.mySave] = Find_Frequency_Sweep(handles.HW, handles.mySave, 600*handles.Seq.notCalibrate);
      assignin('base', 'mySave', handles.mySave);
      if ~isa(handles.HW, 'PD.HWClass')
        % keep track of potential changes of shim setting during frequency sweep
        HW_old = evalin('base', 'HW');
        handles.HW.MagnetShim = HW_old.MagnetShim;
      end
      if ~isa(HW, 'PD.HWClass'),  assignin('base', 'HW', handles.HW);  end
      figure(handles.figure_SpinEcho);
      set(handles.edit_fLarmor, 'String', ['auto ', num2str(handles.HW.fLarmor/1e6, '%3.6f')]);
    else
      handles.HW.fLarmor = str2double(get(handles.edit_fLarmor, 'String'))*1e6;
      handles.HW.B0 = handles.HW.fLarmor/handles.HW.GammaDef*2*pi;
      if ~handles.Seq.notCalibrate
        if ~isa(handles.HW, 'PD.HWClass')
          % update with current settings from base workspace
          handles.HW = evalin('base', 'HW');
        end
        [handles.HW, handles.mySave] = Find_Frequency_Sweep(handles.HW, handles.mySave, 600*handles.Seq.notCalibrate);
        assignin('base', 'mySave', handles.mySave);
        if ~isa(handles.HW, 'PD.HWClass')
          % keep track of potential changes of shim setting during frequency sweep
          HW_old = evalin('base', 'HW');
          handles.HW.MagnetShim = HW_old.MagnetShim;
        end
        if ~isa(HW, 'PD.HWClass'),  assignin('base', 'HW', handles.HW);  end
        figure(handles.figure_SpinEcho);
        set(handles.edit_fLarmor, 'String', ['auto ', num2str(handles.HW.fLarmor/1e6,'%3.6f')]);
      end
    end
  end

  popup_sel_String = get(handles.pushbutton_update, 'String');
  switch popup_sel_String
    case 'B1+'
      if isnan(str2double(get(handles.edit_p90, 'String')))
        [handles.HW, handles.mySave] = Find_PulseDuration(handles.HW, handles.mySave, 0, 1);
      else
        [handles.HW, handles.mySave] = Find_PulseDuration(handles.HW, handles.mySave, 0, 1, 3, get(handles.edit_p90, 'String')*1e-6);
      end
      if exist('HW', 'var') && isa(HW, 'PD.HWClass')
        ResetStructs;
      else
        LoadSystem;
      end
      handles.HW = HW;
      handles.mySave = mySave;
      assignin('base', 'mySave', handles.mySave);
      if ~isa(HW, 'PD.HWClass'),  assignin('base', 'HW', handles.HW);  end
      set(handles.edit_p90, 'String', ['def ', num2str(handles.HW.tFlip90Def*1e6, '%3.3f')]);
      set(handles.edit_p180, 'String', ['def ', num2str(handles.HW.tFlip180Def*1e6, '%3.3f')]);
      handles.Seq.p90 = handles.HW.tFlip90Def;
      handles.Seq.p180 = handles.HW.tFlip180Def;
      figure(handles.figure_SpinEcho);

    case 'B0'
      t = handles.HW.FindFrequencyPlot;
      handles.HW.FindFrequencyPlot = 1;
      [handles.HW, handles.mySave] = Find_Frequency_Sweep(handles.HW, handles.mySave, 0, [], 1);
      assignin('base', 'mySave', handles.mySave);
      if ~isa(HW, 'PD.HWClass'),  assignin('base', 'HW', handles.HW);  end
      figure(handles.figure_SpinEcho);
      handles.HW.FindFrequencyPlot = t;
      set(handles.edit_fLarmor, 'String', ['auto ', num2str(handles.HW.fLarmor/1e6,'%3.6f')]);

    case 'Shim'
      [handles.HW, handles.mySave] = Find_Shim(handles.HW, handles.mySave, 0, 1);
      assignin('base', 'mySave', handles.mySave);
      if ~isa(HW, 'PD.HWClass'),  assignin('base', 'HW', handles.HW);  end

      set(handles.edit_shimX, 'String', num2str(handles.HW.MagnetShim(1)*1e3, '%.3f'));
      set(handles.edit_shimY, 'String', num2str(handles.HW.MagnetShim(2)*1e3, '%.3f'));
      set(handles.edit_shimZ, 'String', num2str(handles.HW.MagnetShim(3)*1e3, '%.3f'));

    otherwise

      popup_sel_index = get(handles.popupmenu_axes_large, 'Value');
      switch popup_sel_index
        case 1
          handles.Seq.plot = 0;
          handles.Seq.plotAllHandle = handles.axes_large;
          if popup_sel_index == get(handles.popupmenu_axes_small, 'Value')
            set(handles.popupmenu_axes_small, 'Value', 2);
          end

        case 2
          handles.Seq.plotFidHandle = handles.axes_large;
          if popup_sel_index == get(handles.popupmenu_axes_small, 'Value')
            set(handles.popupmenu_axes_small, 'Value', 1);
          end

        case 3
          handles.Seq.plotEchoesHandle = handles.axes_large;
          if popup_sel_index == get(handles.popupmenu_axes_small, 'Value')
            set(handles.popupmenu_axes_small, 'Value', 1);
          end

        case 4
          handles.Seq.plotMaxEchoesHandle = handles.axes_large;
          if popup_sel_index == get(handles.popupmenu_axes_small, 'Value')
            set(handles.popupmenu_axes_small, 'Value', 1);
          end

      end

      popup_sel_index = get(handles.popupmenu_axes_small, 'Value');
      switch popup_sel_index
        case 1
          handles.Seq.plotAllHandle = handles.axes_small;
        case 2
          handles.Seq.plotFidHandle = handles.axes_small;
        case 3
          handles.Seq.plotEchoesHandle = handles.axes_small;
        case 4
          handles.Seq.plotMaxEchoesHandle = handles.axes_small;
      end


      if ~handles.refreshPlot
        if isnan(str2double(get(handles.edit_p90, 'String')))
          set(handles.edit_p90, 'String', ['def ', num2str(handles.HW.tFlip90Def*1e6, '%3.3f')])
          handles.Seq.p90 = handles.HW.tFlip90Def;
        else
          handles.Seq.p90 = str2double(get(handles.edit_p90, 'String'))/1e6;
        end

        if isnan(str2double(get(handles.edit_p180, 'String')))
          set(handles.edit_p180, 'String', ['def ', num2str(handles.HW.tFlip180Def*1e6,'%3.3f')])
          handles.Seq.p180 = handles.HW.tFlip180Def;
        else
          handles.Seq.p180 = str2double(get(handles.edit_p180, 'String'))/1e6;
        end

        if isnan(str2double(get(handles.edit_tEcho, 'String')))
          set(handles.edit_tEcho, 'String', ['def ', num2str(10,'%3.0f')])
          handles.Seq.tEcho = 10e-3;
        else
          handles.Seq.tEcho = str2double(get(handles.edit_tEcho, 'String'))/1e3;
        end

        if isnan(str2double(get(handles.edit_nEchoes, 'String')))
          set(handles.edit_nEchoes, 'String', ['def ', num2str(1,'%3.0f')])
          handles.Seq.nEchos = 1;
        else
          handles.Seq.nEchos = str2double(get(handles.edit_nEchoes, 'String'));
        end

        if get(handles.checkbox_unit_B1, 'Value')
          handles.Seq.B1Amp2Str = 1e6;
          set(handles.text_unit_B1, 'String', [char(181) 'T'])
        else
          handles.Seq.B1Amp2Str = handles.HW.GammaDef/2/pi/1000;
          set(handles.text_unit_B1, 'String', 'kHz')
        end

        if isnan(str2double(get(handles.edit_B1, 'String')))
          set(handles.edit_B1, 'String', ['def ', num2str(handles.HW.TX(handles.iDevice).AmpDef*handles.Seq.B1Amp2Str,'%3.6f')])
          handles.Seq.TXAmp = handles.HW.TX(handles.iDevice).AmpDef;
        else
          handles.Seq.TXAmp = str2double(get(handles.edit_B1, 'String'))/handles.Seq.B1Amp2Str;
        end


        if (handles.Seq.nEchos<=500 && (handles.Seq.nEchos+1)*handles.Seq.tEcho<=200)
          handles.Seq.fast = 1;
        else
          handles.Seq.fast = 0;
          handles.Seq.tOffset = handles.Seq.tOffset + zeros(1, handles.Seq.nEchos);
        end

        handles.Seq.fSample = max((1/(handles.Seq.tEcho*0.1-handles.Seq.p90/2-handles.HW.TX(handles.iDevice).BlankOffset-5e-6)*10), ...
          handles.HW.RX(handles.iDevice).fSample/6250);
        if (handles.Seq.fSample > handles.HW.RX(handles.iDevice).fSample/125*2) && (handles.Seq.nEchos>100)
          handles.Seq.AQFID = 0.5*(handles.HW.RX(handles.iDevice).fSample/125*2)/handles.Seq.fSample;
          handles.Seq.AQEcho = 0.5*(handles.HW.RX(handles.iDevice).fSample/125*2)/handles.Seq.fSample;
          handles.Seq.fSample = max((1/(handles.Seq.tEcho*0.1-handles.Seq.p90/2-handles.HW.TX(handles.iDevice).BlankOffset-5e-6)*10), ...
            handles.HW.RX(handles.iDevice).fSample/6250);
        end

        if handles.Seq.nEchos==0,  handles.Seq.AQFID = 1;  end
        handles.Running = 1;
        if ~isa(handles.HW, 'PD.HWClass')
          % update with current settings before overriding handles
          handles.HW = evalin('base', 'HW');
        end
        guidata(hObject, handles);
        [handles.data, handles.SeqOut] = sequence_EchoStandard(handles.HW, handles.Seq);
        drawnow();
      end

      handles.refreshPlot = 0;

      if ~isfield(handles.SeqOut, 'AQ')
        handles.Seq.plotAllHandle = [];
        handles.Seq.plotFidHandle = [];
        handles.Seq.plotEchoesHandle = [];
        handles.Seq.plotMaxEchoesHandle = [];
      else
        Channel = 1;  % FIXME: Add support for multiple acquisition channels?
        if isfield(handles.SeqOut.AQ(1), 'Device')
          iAQ = find([handles.SeqOut.AQ(:).Channel] == Channel & [handles.SeqOut.AQ(:).Device] == handles.iDevice, 1, 'first');
        else
          iAQ = 1;
        end
      end

      %% plotFidHandle
      if ~isempty(handles.Seq.plotFidHandle)
        FIDdata = handles.data(iAQ).data(1:handles.SeqOut.AQ(iAQ).nSamples(1,1), 1, 1);
        FIDtime = handles.data(iAQ).time_all(1:handles.SeqOut.AQ(iAQ).nSamples(1,1), 1, 1);

        if all(imag(FIDdata)==0)
          plot(handles.Seq.plotFidHandle, FIDtime, real(FIDdata)/HW.RX(handles.iDevice).AmplitudeUnitScale);
        else
          plot(handles.Seq.plotFidHandle, ...
            FIDtime, [abs(FIDdata), real(FIDdata), imag(FIDdata)]/HW.RX(handles.iDevice).AmplitudeUnitScale);
        end
        xlabel(handles.Seq.plotFidHandle,'time / s');
        ylabel(handles.Seq.plotFidHandle, {HW.RX(handles.iDevice).AmplitudeName; HW.RX(handles.iDevice).AmplitudeUnit});

        handles.mySave.plotFidHandleMaxNew = max(abs(FIDdata(~isnan(FIDdata)))) / HW.RX(handles.iDevice).AmplitudeUnitScale;
        handles.mySave.plotFidHandleMinNew = min(min(real(FIDdata(~isnan(FIDdata)))),min(imag(FIDdata(~isnan(FIDdata))))) / HW.RX(handles.iDevice).AmplitudeUnitScale;
        handles.mySave.plotFidHandleMinNew = min(handles.mySave.plotFidHandleMinNew,0);
        [handles.mySave.plotFidHandleMinNew, handles.mySave.plotFidHandleMaxNew] = ...
          get_axisMinMax(handles.mySave.plotFidHandleMinNew, handles.mySave.plotFidHandleMaxNew);
        if~isfield(handles.mySave, 'plotFidHandleMax')
          handles.mySave.plotFidHandleMax = handles.mySave.plotFidHandleMaxNew;
        end
        if~isfield(handles.mySave, 'plotFidHandleMin')
          handles.mySave.plotFidHandleMin = handles.mySave.plotFidHandleMinNew;
        end
        [handles.mySave.plotFidHandleMin, handles.mySave.plotFidHandleMax] = ...
          get_axisMinMax(handles.mySave.plotFidHandleMinNew, handles.mySave.plotFidHandleMaxNew, handles.mySave.plotFidHandleMin, handles.mySave.plotFidHandleMax);
        ylim(handles.Seq.plotFidHandle, [handles.mySave.plotFidHandleMin, handles.mySave.plotFidHandleMax]);
        xlim(handles.Seq.plotFidHandle, [0, +Inf]);
        grid(handles.Seq.plotFidHandle, 'on');
        zoom(handles.Seq.plotFidHandle, 'on');

      end

      %% plotEchoesHandle
      if ~isempty(handles.Seq.plotEchoesHandle)
        ok = 0;
        if size(handles.data(iAQ).data, 2) >= 2
          EchoesData = handles.data(iAQ).data(1:handles.SeqOut.AQ(iAQ).nSamples(2,1),2:end,1);
          EchoesTime = handles.data(iAQ).time_all(1:handles.SeqOut.AQ(iAQ).nSamples(2,1),2:end,1);
          ok = 1;
        elseif size(handles.data(iAQ).data, 3) >= 2
          EchoesData = squeeze(handles.data(iAQ).data(1:handles.SeqOut.AQ(iAQ).nSamples(1,2),1,2:end));
          EchoesTime = squeeze(handles.data(iAQ).time_all(1:handles.SeqOut.AQ(iAQ).nSamples(1,2),1,2:end));
          ok = 1;
        end
        if ok
          if all(imag(EchoesData) == 0)
            plot(handles.Seq.plotEchoesHandle, ...
              EchoesTime, real(EchoesData) / HW.RX(handles.iDevice).AmplitudeUnitScale);
          else
            colorOrder = get(handles.Seq.plotEchoesHandle, 'ColorOrder');
            plot(handles.Seq.plotEchoesHandle, EchoesTime, abs(EchoesData)/HW.RX(handles.iDevice).AmplitudeUnitScale, 'Color', colorOrder(1,:));
            hold(handles.Seq.plotEchoesHandle, 'on');
            plot(handles.Seq.plotEchoesHandle, EchoesTime, real(EchoesData)/HW.RX(handles.iDevice).AmplitudeUnitScale, 'Color', colorOrder(2,:));
            plot(handles.Seq.plotEchoesHandle, EchoesTime, imag(EchoesData)/HW.RX(handles.iDevice).AmplitudeUnitScale, 'Color', colorOrder(3,:));
            hold(handles.Seq.plotEchoesHandle, 'off');
          end
          % set(handles.Seq.plotEchoesHandle, 'ColorOrder', [1 1 1; 0 1 0; 0 0 1]);
          xlabel(handles.Seq.plotEchoesHandle, 'time / s');
          ylabel(handles.Seq.plotEchoesHandle, {HW.RX(handles.iDevice).AmplitudeName; HW.RX(handles.iDevice).AmplitudeUnit});
          % legend(handles.Seq.plotEchoesHandle, 'abs', 'real', 'imag');

          handles.mySave.plotEchoesHandleMaxNew = max(abs(EchoesData(~isnan(EchoesData))))/HW.RX(handles.iDevice).AmplitudeUnitScale;
          handles.mySave.plotEchoesHandleMinNew = min(min(real(EchoesData(~isnan(EchoesData)))), min(imag(EchoesData(~isnan(EchoesData)))))/HW.RX(handles.iDevice).AmplitudeUnitScale;
          handles.mySave.plotEchoesHandleMinNew = min(handles.mySave.plotEchoesHandleMinNew,0);
          [handles.mySave.plotEchoesHandleMinNew, handles.mySave.plotEchoesHandleMaxNew] = ...
            get_axisMinMax(handles.mySave.plotEchoesHandleMinNew, handles.mySave.plotEchoesHandleMaxNew);
          if ~isfield(handles.mySave, 'plotEchoesHandleMax')
            handles.mySave.plotEchoesHandleMax = handles.mySave.plotEchoesHandleMaxNew;
          end
          if ~isfield(handles.mySave, 'plotEchoesHandleMin')
            handles.mySave.plotEchoesHandleMin = handles.mySave.plotEchoesHandleMinNew;
          end
          [handles.mySave.plotEchoesHandleMin, handles.mySave.plotEchoesHandleMax] = ...
            get_axisMinMax(handles.mySave.plotEchoesHandleMinNew, handles.mySave.plotEchoesHandleMaxNew, handles.mySave.plotEchoesHandleMin, handles.mySave.plotEchoesHandleMax);
          ylim(handles.Seq.plotEchoesHandle, [handles.mySave.plotEchoesHandleMin, handles.mySave.plotEchoesHandleMax]);
          xlim(handles.Seq.plotEchoesHandle, [handles.SeqOut.tEcho/2, handles.SeqOut.tEcho*(0.5+handles.SeqOut.nEchos)]);
          grid(handles.Seq.plotEchoesHandle, 'on');
          zoom(handles.Seq.plotEchoesHandle, 'on');
        else
          plot(handles.Seq.plotEchoesHandle, NaN, NaN);
        end

      end


      %% plotAllHandle
      if ~isempty(handles.Seq.plotAllHandle)
        ok = 0;
        if size(handles.data(iAQ).data, 2) >= 2
          EchoesData = handles.data(iAQ).data(1:handles.SeqOut.AQ(iAQ).nSamples(2,1),2:end,1);
          EchoesTime = handles.data(iAQ).time_all(1:handles.SeqOut.AQ(iAQ).nSamples(2,1),2:end,1);
          ok = 1;
        elseif size(handles.data(iAQ).data, 3) >= 2
          EchoesData = squeeze(handles.data(iAQ).data(1:handles.SeqOut.AQ(iAQ).nSamples(1,2),1,2:end));
          EchoesTime = squeeze(handles.data(iAQ).time_all(1:handles.SeqOut.AQ(iAQ).nSamples(1,2),1,2:end));
          ok = 1;
        else
          EchoesData = NaN;
          EchoesTime = NaN;
          ok = 1;
        end

        EchoesData = [EchoesData; NaN(1, size(EchoesData,2))];
        EchoesTime = [EchoesTime; NaN(1, size(EchoesTime,2))];

        FIDdata = handles.data(iAQ).data(1:handles.SeqOut.AQ(iAQ).nSamples(1,1),1,1);
        FIDtime = handles.data(iAQ).time_all(1:handles.SeqOut.AQ(iAQ).nSamples(1,1),1,1);

        AllData = [FIDdata; NaN; EchoesData(:)];
        AllTime = [FIDtime; NaN; EchoesTime(:)];

        if ok
          if all(imag(AllData)==0)
            plot(handles.Seq.plotAllHandle, AllTime, real(AllData)/HW.RX(handles.iDevice).AmplitudeUnitScale);
          else
            plot(handles.Seq.plotAllHandle, ...
              AllTime, abs(AllData)/HW.RX(handles.iDevice).AmplitudeUnitScale, '-', ...
              AllTime, real(AllData)/HW.RX(handles.iDevice).AmplitudeUnitScale, '-', ...
              AllTime, imag(AllData)/HW.RX(handles.iDevice).AmplitudeUnitScale, '-');
          end
          % set(handles.Seq.plotAllHandle,   'ColorOrder', [1 1 1; 0 1 0; 0 0 1]);
          xlabel(handles.Seq.plotAllHandle, 'time / s');
          ylabel(handles.Seq.plotAllHandle, {HW.RX(handles.iDevice).AmplitudeName; HW.RX(handles.iDevice).AmplitudeUnit});
          % legend(handles.Seq.plotAllHandle, {'abs','real','imag'});

          handles.mySave.plotAllHandleMaxNew = max(abs(AllData(~isnan(AllData)))) / HW.RX(handles.iDevice).AmplitudeUnitScale;
          handles.mySave.plotAllHandleMinNew = min(min(real(AllData(~isnan(AllData)))), min(imag(AllData(~isnan(AllData)))))/HW.RX(handles.iDevice).AmplitudeUnitScale;
          handles.mySave.plotAllHandleMinNew = min(handles.mySave.plotAllHandleMinNew, 0);
          [handles.mySave.plotAllHandleMinNew, handles.mySave.plotAllHandleMaxNew] = ...
            get_axisMinMax(handles.mySave.plotAllHandleMinNew, handles.mySave.plotAllHandleMaxNew);
          if ~isfield(handles.mySave, 'plotAllHandleMax')
            handles.mySave.plotAllHandleMax = handles.mySave.plotAllHandleMaxNew;
          end
          if ~isfield(handles.mySave, 'plotAllHandleMin')
            handles.mySave.plotAllHandleMin = handles.mySave.plotAllHandleMinNew; end
          [handles.mySave.plotAllHandleMin, handles.mySave.plotAllHandleMax] = ...
            get_axisMinMax(handles.mySave.plotAllHandleMinNew, handles.mySave.plotAllHandleMaxNew, handles.mySave.plotAllHandleMin, handles.mySave.plotAllHandleMax);
          ylim(handles.Seq.plotAllHandle, [handles.mySave.plotAllHandleMin, handles.mySave.plotAllHandleMax]);
          xlim(handles.Seq.plotAllHandle, [0, +Inf]);
          grid(handles.Seq.plotAllHandle, 'on');

          % zoom(handles.Seq.plotAllHandle,'on');
        end
      end


      %% plotMaxEchoesHandle
      if ~isempty(handles.Seq.plotMaxEchoesHandle)
        ok = 0;
        if size(handles.data(iAQ).data, 2) >= 2
          if round(handles.SeqOut.AQ(iAQ).nSamples(2,1)/2) >= 2
            EchoesData = handles.data(iAQ).data(round(handles.SeqOut.AQ(iAQ).nSamples(2,1)/2),2:end,1);
            EchoesTime = handles.data(iAQ).time_all(round(handles.SeqOut.AQ(iAQ).nSamples(2,1)/2),2:end,1);
            ok = 1;
          end
        elseif size(handles.data(iAQ).data, 3) >= 2
          if round(handles.SeqOut.AQ.nSamples(1,2)/2) >= 2
            EchoesData = squeeze(handles.data(iAQ).data(round(handles.SeqOut.AQ(iAQ).nSamples(1,2)/2),1,2:end));
            EchoesTime = squeeze(handles.data(iAQ).time_all(round(handles.SeqOut.AQ(iAQ).nSamples(1,2)/2),1,2:end));
            ok = 1;
          end
        end
        if ok
          plot(handles.Seq.plotMaxEchoesHandle, ...
            EchoesTime(:), abs(EchoesData(:))/HW.RX(handles.iDevice).AmplitudeUnitScale);
          xlabel(handles.Seq.plotMaxEchoesHandle, 'time / s');
          ylabel(handles.Seq.plotMaxEchoesHandle, {HW.RX(handles.iDevice).AmplitudeName; HW.RX(handles.iDevice).AmplitudeUnit});

          handles.mySave.plotMaxEchoesHandleMaxNew = max(abs(EchoesData(~isnan(EchoesData))))/HW.RX(handles.iDevice).AmplitudeUnitScale;
          handles.mySave.plotMaxEchoesHandleMinNew = 0;
          [handles.mySave.plotMaxEchoesHandleMinNew, handles.mySave.plotMaxEchoesHandleMaxNew] = ...
            get_axisMinMax(handles.mySave.plotMaxEchoesHandleMinNew, handles.mySave.plotMaxEchoesHandleMaxNew);
          if ~isfield(handles.mySave, 'plotMaxEchoesHandleMax')
            handles.mySave.plotMaxEchoesHandleMax = handles.mySave.plotMaxEchoesHandleMaxNew;
          end
          if ~isfield(handles.mySave, 'plotMaxEchoesHandleMin')
            handles.mySave.plotMaxEchoesHandleMin = handles.mySave.plotMaxEchoesHandleMinNew;
          end
          [handles.mySave.plotMaxEchoesHandleMin, handles.mySave.plotMaxEchoesHandleMax] = ...
            get_axisMinMax(handles.mySave.plotMaxEchoesHandleMinNew, handles.mySave.plotMaxEchoesHandleMaxNew, handles.mySave.plotMaxEchoesHandleMin, handles.mySave.plotMaxEchoesHandleMax);
          ylim(handles.Seq.plotMaxEchoesHandle, [handles.mySave.plotMaxEchoesHandleMin, handles.mySave.plotMaxEchoesHandleMax]);
          xlim(handles.Seq.plotMaxEchoesHandle, [0, +Inf]);
          grid(handles.Seq.plotMaxEchoesHandle, 'on');

          zoom(handles.Seq.plotMaxEchoesHandle, 'on');
        else
          plot(handles.Seq.plotMaxEchoesHandle, NaN, NaN);
        end

      end
      %%

      if ~isa(handles.HW, 'PD.HWClass')
        % update with current settings before overriding handles
        handles.HW = evalin('base', 'HW');
      end

      % handles.HW.ReInit = 0;
      % handles.HW.tRepInit = 0.005;
  end

  set(handles.radiobutton2, 'Value', 1);
  handles.Seq.firstTR = 0;
  guidata(hObject, handles);
  uipanel1_SelectionChangeFcn(hObject, eventdata, handles);
end
% handles.HW.tRepInit = 0.3;
% handles.HW.ReInit = 1;
handles.Running = 0;
guidata(hObject, handles);

end


% --- Executes on button press in checkbox_unit_B1.
function checkbox_unit_B1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_unit_B1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_unit_B1
if get(handles.checkbox_unit_B1, 'Value')
  handles.Seq.B1Amp2Str=1e6;
  set(handles.text_unit_B1, 'String', [char(181) 'T'])
else
  handles.Seq.B1Amp2Str=handles.HW.GammaDef/2/pi/1000;
  set(handles.text_unit_B1, 'String', 'kHz')
end

if isfield(handles,'Seq')
  if isfield(handles.Seq,'TXAmp')
    if isfield(handles.Seq,'B1Amp2Str')
      if isnan(str2double(get(handles.edit_B1, 'String')))
        set(handles.edit_B1, 'String', ['def ', num2str(handles.Seq.TXAmp*handles.Seq.B1Amp2Str,'%3.6f')])
      else
        set(handles.edit_B1, 'String', num2str(handles.Seq.TXAmp*handles.Seq.B1Amp2Str,'%3.6f'))
      end
    end
  end
end

end


% --- Executes on button press in checkbox_repetition.
function checkbox_repetition_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_repetition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_repetition
if get(handles.checkbox_repetition, 'Value')
    pushbutton_update_Callback(hObject, eventdata, handles)
end

end


function edit_B1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_B1 as text
%        str2double(get(hObject,'String')) returns contents of edit_B1 as a double

end


% --- Executes during object creation, after setting all properties.
function edit_B1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_B1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

end


% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%  EventName: string 'SelectionChanged' (read only)
%  OldValue: handle of the previously selected object or empty if none was selected
%  NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% handles.Cal.Cal=~get(handles.radiobutton2,'Value');
% handles.Cal.short=get(handles.radiobutton3,'Value');
% handles.Cal.open=get(handles.radiobutton4,'Value');
% handles.Cal.terminated=get(handles.radiobutton5,'Value');

if get(handles.radiobutton2,'Value'),  set(handles.pushbutton_update, 'String', 'Update'); end
if get(handles.radiobutton3,'Value'),  set(handles.pushbutton_update, 'String', 'B1+'); end
if get(handles.radiobutton4,'Value'),  set(handles.pushbutton_update, 'String', 'B0'); end
if get(handles.radiobutton5,'Value'),  set(handles.pushbutton_update, 'String', 'Shim'); end

guidata(hObject, handles);

end


function edit_fLarmor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fLarmor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fLarmor as text
%        str2double(get(hObject,'String')) returns contents of edit_fLarmor as a double

end


% --- Executes during object creation, after setting all properties.
function edit_fLarmor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fLarmor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in checkbox_calibrate.
function checkbox_calibrate_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_calibrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_calibrate

end


function edit_shimX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shimX as text
%        str2double(get(hObject,'String')) returns contents of edit_shimX as a double

shimX = str2double(get(hObject, 'String'));
if isfinite(shimX)
  handles.HW.MagnetShim(1) = shimX/1e3;
else
  set(hObject, 'String', num2str(handles.HW.MagnetShim(1)*1e3, '%.3f'));
end

end


% --- Executes during object creation, after setting all properties.
function edit_shimX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

end


function edit_shimX_inc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shimX_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shimX_inc as text
%        str2double(get(hObject,'String')) returns contents of edit_shimX_inc as a double

shimX = str2double(get(hObject, 'String'));
if ~isfinite(shimX)
  set(hObject, 'String', '0.500');
end

end


% --- Executes during object creation, after setting all properties.
function edit_shimX_inc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shimX_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in pushbutton_shimX_minus.
function pushbutton_shimX_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shimX_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

shimXinc = str2double(get(handles.edit_shimX_inc, 'String'))*1e-3;

if ~isa(handles.HW, 'PD.HWClass'),  handles.HW = evalin('base', 'HW');  end
handles.HW.MagnetShim(1) = handles.HW.MagnetShim(1) - shimXinc;
if ~isa(handles.HW, 'PD.HWClass'),  assignin('base', 'HW', handles.HW);  end

set(handles.edit_shimX, 'String', num2str(handles.HW.MagnetShim(1)*1e3, '%.3f'));

end


% --- Executes on button press in pushbutton_shimX_minus.
function pushbutton_shimX_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shimX_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

shimXinc = str2double(get(handles.edit_shimX_inc, 'String'))*1e-3;

if ~isa(handles.HW, 'PD.HWClass'),  handles.HW = evalin('base', 'HW');  end
handles.HW.MagnetShim(1) = handles.HW.MagnetShim(1) + shimXinc;
if ~isa(handles.HW, 'PD.HWClass'),  assignin('base', 'HW', handles.HW);  end

set(handles.edit_shimX, 'String', num2str(handles.HW.MagnetShim(1)*1e3, '%.3f'));

end


function edit_shimY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shimY as text
%        str2double(get(hObject,'String')) returns contents of edit_shimY as a double

shimY = str2double(get(hObject, 'String'));
if isfinite(shimY)
  handles.HW.MagnetShim(2) = shimY/1e3;
else
  set(hObject, 'String', num2str(handles.HW.MagnetShim(2)*1e3, '%.3f'));
end

end


% --- Executes during object creation, after setting all properties.
function edit_shimY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

end


function edit_shimY_inc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shimY_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shimY_inc as text
%        str2double(get(hObject,'String')) returns contents of edit_shimY_inc as a double

shimY = str2double(get(hObject, 'String'));
if ~isfinite(shimY)
  set(hObject, 'String', '0.500');
end

end


% --- Executes during object creation, after setting all properties.
function edit_shimY_inc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shimY_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in pushbutton_shimX_minus.
function pushbutton_shimY_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shimY_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

shimYinc = str2double(get(handles.edit_shimY_inc, 'String'))*1e-3;

if ~isa(handles.HW, 'PD.HWClass'),  handles.HW = evalin('base', 'HW');  end
handles.HW.MagnetShim(2) = handles.HW.MagnetShim(2) - shimYinc;
if ~isa(handles.HW, 'PD.HWClass'),  assignin('base', 'HW', handles.HW);  end

set(handles.edit_shimY, 'String', num2str(handles.HW.MagnetShim(2)*1e3, '%.3f'));

end


% --- Executes on button press in pushbutton_shimX_minus.
function pushbutton_shimY_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shimY_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

shimYinc = str2double(get(handles.edit_shimY_inc, 'String'))*1e-3;

if ~isa(handles.HW, 'PD.HWClass'),  handles.HW = evalin('base', 'HW');  end
handles.HW.MagnetShim(2) = handles.HW.MagnetShim(2) + shimYinc;
if ~isa(handles.HW, 'PD.HWClass'),  assignin('base', 'HW', handles.HW);  end

set(handles.edit_shimY, 'String', num2str(handles.HW.MagnetShim(2)*1e3, '%.3f'));

end


function edit_shimZ_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shimZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shimY as text
%        str2double(get(hObject,'String')) returns contents of edit_shimY as a double

shimZ = str2double(get(hObject, 'String'));
if isfinite(shimZ)
  handles.HW.MagnetShim(3) = shimZ/1e3;
else
  set(hObject, 'String', num2str(handles.HW.MagnetShim(3)*1e3, '%.3f'));
end

end


% --- Executes during object creation, after setting all properties.
function edit_shimZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shimZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

end


function edit_shimZ_inc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shimZ_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shimY_inc as text
%        str2double(get(hObject,'String')) returns contents of edit_shimY_inc as a double

shimZ = str2double(get(hObject, 'String'));
if ~isfinite(shimZ)
  set(hObject, 'String', '0.500');
end

end


% --- Executes during object creation, after setting all properties.
function edit_shimZ_inc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shimZ_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in pushbutton_shimX_minus.
function pushbutton_shimZ_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shimZ_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

shimZinc = str2double(get(handles.edit_shimZ_inc, 'String'))*1e-3;

if ~isa(handles.HW, 'PD.HWClass'),  handles.HW = evalin('base', 'HW');  end
handles.HW.MagnetShim(3) = handles.HW.MagnetShim(3) - shimZinc;
if ~isa(handles.HW, 'PD.HWClass'),  assignin('base', 'HW', handles.HW);  end

set(handles.edit_shimZ, 'String', num2str(handles.HW.MagnetShim(3)*1e3, '%.3f'));

end


% --- Executes on button press in pushbutton_shimX_minus.
function pushbutton_shimZ_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shimZ_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

shimZinc = str2double(get(handles.edit_shimZ_inc, 'String'))*1e-3;

if ~isa(handles.HW, 'PD.HWClass'),  handles.HW = evalin('base', 'HW');  end
handles.HW.MagnetShim(3) = handles.HW.MagnetShim(3) + shimZinc;
if ~isa(handles.HW, 'PD.HWClass'),  assignin('base', 'HW', handles.HW);  end

set(handles.edit_shimZ, 'String', num2str(handles.HW.MagnetShim(3)*1e3, '%.3f'));

end


% --- Executes on button press in pushbutton_shimX_minus.
function pushbutton_save_shim_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save_shim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isa(handles.HW, 'PD.HWClass'),  handles.HW = evalin('base', 'HW');  end
shimGradIdxStr = sprintf('%d,', find(handles.HW.Grad(handles.iDevice).ShimGradients~=0));
shimElemStr = sprintf('%6.9f, ', handles.HW.Grad(handles.iDevice).AmpOffset(handles.HW.Grad(handles.iDevice).ShimGradients~=0));

ShimSequenceStr = sprintf(' using "%s"', mfilename());

if numel(handles.HW.Grad) > 1
  newCalLine = sprintf(['HW.Grad(%d).AmpOffset([%s]) = [%s]; ', ...
    ' %% %s manually%s, x y z in T/m and B0 in T\n'], ...
    handles.iDevice, shimGradIdxStr(1:end-1), shimElemStr(1:end-2), ...
    datestr(now, 'yyyy-mm-ddTHH:MM:SS'), ShimSequenceStr);
else
  newCalLine = sprintf(['HW.MagnetShim([%s]) = [%s]; ', ...
    ' %% %s manually%s, x y z in T/m and B0 in T\n'], ...
    shimGradIdxStr(1:end-1), shimElemStr(1:end-2), ...
    datestr(now, 'yyyy-mm-ddTHH:MM:SS'), ShimSequenceStr);
end

if ~exist(fileparts(handles.HW.MagnetShimPath), 'dir')
  mkdir(fileparts(handles.HW.MagnetShimPath));
end
fid = fopen(handles.HW.MagnetShimPath, 'a+');
fwrite(fid, newCalLine);
fclose(fid);
disp('A new line was added to the following file: ');
disp(handles.HW.MagnetShimPath);
fprintf('\n');
disp(newCalLine);
fprintf('\n');

end
