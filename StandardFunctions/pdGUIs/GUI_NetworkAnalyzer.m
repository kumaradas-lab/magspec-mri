function varargout = GUI_NetworkAnalyzer(varargin)
%% Program to use Pure Devices driveL as a Network Analyzer
%
%     hf = GUI_NetworkAnalyzer()
%
% Make sure to execute "LoadSystem" before starting the Network Analyzer.
%
% INPUT:
%   none necessary
%
% OUTPUT:
%   hf      handle to the figure of the Network Analyzer.
%
% ------------------------------------------------------------------------
% (C) Copyright 2011-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------

% Last Modified by GUIDE v2.5 05-Mar-2019 11:07:51

%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_NetworkAnalyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_NetworkAnalyzer_OutputFcn, ...
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

% --- Executes just before GUI_NetworkAnalyzer is made visible.
function GUI_NetworkAnalyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_NetworkAnalyzer (see VARARGIN)

% Choose default command line output for GUI_NetworkAnalyzer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using GUI_NetworkAnalyzer.
if strcmp(get(hObject, 'Visible'), 'off')
%     plot(rand(5));
end
if evalin('base', 'exist(''HW'', ''var'')')
    HW = evalin('base', 'HW');
    if ~isempty(HW.Network.fCenter)
        set(handles.edit_fCenter, 'String', num2str(HW.Network.fCenter/1e6,8))
    end
end

handles.popupmenuChoices = ...
  {'Z-Smith chart', 'Y-Smith chart', 'YZ-Smith chart', ...
   'Reflection S11', 'Nyquist plot', 'Reflection S11 raw'};
set([handles.popupmenu_axes1, handles.popupmenu_axes2], ...
    'String', handles.popupmenuChoices);
set(handles.popupmenu_axes1, 'Value', 1);
set(handles.popupmenu_axes2, 'Value', 4);
guidata(hObject, handles);

% UIWAIT makes GUI_NetworkAnalyzer wait for user response (see UIRESUME)
% uiwait(handles.figure_NetworkAnalyzer);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_NetworkAnalyzer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure_NetworkAnalyzer)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure_NetworkAnalyzer,'Name') '?'],...
                     ['Close ' get(handles.figure_NetworkAnalyzer,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure_NetworkAnalyzer)


% --- Executes on selection change in popupmenu_axes1.
function popupmenu_axes1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_axes1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_axes1


% --- Executes during object creation, after setting all properties.
function popupmenu_axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

% set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});



function edit_fCenter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fCenter as text
%        str2double(get(hObject,'String')) returns contents of edit_fCenter as a double


% --- Executes during object creation, after setting all properties.
function edit_fCenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fSpan_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fSpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fSpan as text
%        str2double(get(hObject,'String')) returns contents of edit_fSpan as a double


% --- Executes during object creation, after setting all properties.
function edit_fSpan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fSpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_steps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_steps as text
%        str2double(get(hObject,'String')) returns contents of edit_steps as a double


% --- Executes during object creation, after setting all properties.
function edit_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_power_Callback(hObject, eventdata, handles)
% hObject    handle to edit_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_power as text
%        str2double(get(hObject,'String')) returns contents of edit_power as a double


% --- Executes during object creation, after setting all properties.
function edit_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_axes2.
function popupmenu_axes2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_axes2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_axes2


% --- Executes during object creation, after setting all properties.
function popupmenu_axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_measure.
function pushbutton_measure_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_measure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Update=1;
if ~isfield(handles, 'Running'),  handles.Running = 0; end
if ~isfield(handles, 'Cal'),      handles.Cal.Cal = 0; end
if ~isfield(handles.Cal, 'Cal'),  handles.Cal.Cal = 0; end

mySave = [];
if isfield(handles,'mySave'); mySave=handles.mySave; end
if evalin('base', 'exist(''HW'', ''var'')')
  HW = evalin('base', 'HW');
end
if evalin('base', 'exist(''mySave'', ''var'')')
  mySave = evalin('base', 'mySave');
end
if exist('HW', 'var')
  ResetStructs
else
  mySave_old = mySave;
  LoadSystem
  if ~isempty(mySave_old)
    mySave = mySave_old;
  end
end

handles.mySave = mySave;
% handles.talker = talker;
handles.HW = HW;
handles.Seq = Seq;
handles.HW.tRepInit = 0.3;
handles.Seq.Cal = 0;

% FIXME: Support multiple MMRT devices
if isemptyfield(handles, 'iDevice'),  handles.iDevice = 1;  end

while handles.Update || ...
    (get(handles.checkbox_autorun, 'Value') && get(handles.radiobutton_measure, 'Value'))
  handles.Seq.plotSmith1Handle = [];
  handles.Seq.plotSmith2Handle = [];
  handles.Seq.plotReflectionHandle = [];
  handles.Seq.plotNyquistHandle = [];
  handles.Seq.plotReflectionRawHandle = [];

  doCal = false;

  if ~get(handles.radiobutton_measure, 'Value')
    doCal = true;
    handles.Seq.plotReflectionRawHandle = handles.axes1;
    cla(handles.axes2, 'reset')
    if get(handles.radiobutton_short, 'Value')
      handles.Seq.Cal = 'short';
    end
    if get(handles.radiobutton_open, 'Value')
      handles.Seq.Cal = 'open';
    end
    if get(handles.radiobutton_term, 'Value')
      handles.Seq.Cal = 'terminated';
    end

    handles.Seq.fCenter = handles.HW.MMRT(handles.iDevice).fSystem/2;
    handles.Seq.fSpan   = handles.HW.MMRT(handles.iDevice).fSystem;
    handles.Seq.fSteps  = round(handles.HW.MMRT(handles.iDevice).fSystem/1e6)*10+1;
    handles.Seq.fSample = handles.HW.RX(handles.iDevice).fSample/handles.HW.MMRT(handles.iDevice).fSystem*1e5;
    handles.Seq.tMessung= 100/handles.Seq.fSample;
    % if ~isempty(handles.HW.Network.TXPowerdBm);handles.Seq.TXPowerdBm=0;else handles.Seq.TXPowerdBm=handles.HW.Network.TXPowerdBm;end
  else
    popup_sel_index = get(handles.popupmenu_axes1, 'Value');
    if any(regexpi(handles.popupmenuChoices{popup_sel_index}, 'smith'))
      handles.Seq.plotSmith1Handle = handles.axes1;
      handles.Seq.plotSmith1Type = regexpi(handles.popupmenuChoices{popup_sel_index}, '([yz]*)-', 'tokens');
      handles.Seq.plotSmith1Type = handles.Seq.plotSmith1Type{1}{1};
      if popup_sel_index == get(handles.popupmenu_axes2, 'Value')
        set(handles.popupmenu_axes2, ...
          'Value', find(cellfun(@isempty, regexpi(handles.popupmenuChoices, 'smith')), 1, 'first'));
      end
    elseif any(regexpi(handles.popupmenuChoices{popup_sel_index}, 'S11 raw'))
      handles.Seq.plotReflectionRawHandle = handles.axes1;
      if popup_sel_index==get(handles.popupmenu_axes2, 'Value')
        set(handles.popupmenu_axes2, 'Value', 1);
      end
    elseif any(regexpi(handles.popupmenuChoices{popup_sel_index}, 'S11'))
      handles.Seq.plotReflectionHandle = handles.axes1;
      if popup_sel_index==get(handles.popupmenu_axes2, 'Value')
        set(handles.popupmenu_axes2, 'Value', 1);
      end
    elseif any(regexpi(handles.popupmenuChoices{popup_sel_index}, 'Nyquist'))
      handles.Seq.plotNyquistHandle = handles.axes1;
      if popup_sel_index==get(handles.popupmenu_axes2, 'Value')
        set(handles.popupmenu_axes2, 'Value', 1);
      end
    end

    popup_sel_index = get(handles.popupmenu_axes2, 'Value');
    if any(regexpi(handles.popupmenuChoices{popup_sel_index}, 'smith'))
        handles.Seq.plotSmith2Handle=handles.axes2;
        handles.Seq.plotSmith2Type = regexpi(handles.popupmenuChoices{popup_sel_index}, '([yz]*)-', 'tokens');
        handles.Seq.plotSmith2Type = handles.Seq.plotSmith2Type{1}{1};
    elseif any(regexpi(handles.popupmenuChoices{popup_sel_index}, 'S11 raw'))
        handles.Seq.plotReflectionRawHandle=handles.axes2;
    elseif any(regexpi(handles.popupmenuChoices{popup_sel_index}, 'S11'))
        handles.Seq.plotReflectionHandle=handles.axes2;
    elseif any(regexpi(handles.popupmenuChoices{popup_sel_index}, 'Nyquist'))
        handles.Seq.plotNyquistHandle=handles.axes2;
    end

    handles.Seq.fCenter=str2double(get(handles.edit_fCenter, 'String'))*1e6;
    if isnan(handles.Seq.fCenter)
      if isempty(handles.HW.Network.fCenter)
        t = handles.HW.FindFrequencyPause;
        handles.HW.FindFrequencyPause = 1;

        ChannelDef = handles.HW.TX(handles.iDevice).ChannelDef;
        if (get(handles.checkbox_TRX1, 'Value'))
          handles.HW.TX(handles.iDevice).ChannelDef = 1;
        else
          handles.HW.TX(handles.iDevice).ChannelDef = 2;
        end
        if ~isa(handles.HW, 'PD.HWClass')
          HW = handles.HW;
          LoadCalcHW;
          handles.HW = HW;
        end

        [handles.HW, handles.mySave] = Find_Frequency_Sweep(handles.HW, handles.mySave);

        handles.HW.TX(handles.iDevice).ChannelDef = ChannelDef;
        if ~isa(handles.HW, 'PD.HWClass')
          HW = handles.HW;
          LoadCalcHW;
          handles.HW = HW;
        end

        assignin('base', 'mySave', handles.mySave);
        if ~isa(HW, 'PD.HWClass')
          assignin('base', 'HW', handles.HW);
        end
        figure(handles.figure_NetworkAnalyzer)
        handles.HW.FindFrequencyPause = t;
        handles.Seq.fCenter = handles.HW.fLarmor;
        set(handles.edit_fCenter, 'String', sprintf('fL=%6.6f', handles.Seq.fCenter/1e6));
      else
        handles.Seq.fCenter=handles.HW.Network.fCenter;
        handles.HW.fLarmor=handles.HW.Network.fCenter;
        set(handles.edit_fCenter, 'String', num2str(handles.Seq.fCenter/1e6))
      end
    end
    handles.Seq.fSpan=str2double(get(handles.edit_fSpan, 'String'))*1e6;
    if isnan(handles.Seq.fSpan)
      if isempty(handles.HW.Network.fSpan)
        handles.Seq.fSpan=1e6;
      else
        handles.Seq.fSpan=handles.HW.Network.fSpan;
      end
      set(handles.edit_fSpan, 'String', num2str(handles.Seq.fSpan/1e6))
    end
    handles.Seq.fSteps=str2double(get(handles.edit_steps, 'String'));
    if isnan(handles.Seq.fSteps)
      if isempty(handles.HW.Network.fSteps)
        handles.Seq.fSteps=101;
      else
        handles.Seq.fSteps=handles.HW.Network.fSteps;
      end
      set(handles.edit_steps, 'String', num2str(handles.Seq.fSteps))
    end
    handles.Seq.TXPowerdBm=str2double(get(handles.edit_power, 'String'));
    if isnan(handles.Seq.TXPowerdBm)
      if isempty(handles.HW.Network.TXPowerdBm)
        handles.Seq.TXPowerdBm=0;
      else
        handles.Seq.TXPowerdBm=handles.HW.Network.TXPowerdBm;
      end
      set(handles.edit_power, 'String', num2str(handles.Seq.TXPowerdBm))
    end
    if isnan(str2double(get(handles.edit_cableLength, 'String')))
      set(handles.edit_cableLength, 'String',num2str(handles.HW.AddCableLength))
    end
    handles.HW.AddCableLength=str2double(get(handles.edit_cableLength, 'String'));
    if ~isfield(handles.mySave, 'setBlankCh2')
      handles.mySave.setBlankCh2=1;
      if ~isempty(handles.HW.Network.BlankCh2)
        set(handles.checkbox_blankTx2, 'Value', handles.HW.Network.BlankCh2==1 )
      end
    end
    handles.Seq.BlankCh2 = get(handles.checkbox_blankTx2, 'Value');
    handles.Seq.DampCoil = get(handles.checkbox_DampCoil, 'Value');
  end
  handles.Running=1;
  handles.Update=0;
  [handles.Network, handles.SeqOut, NetworkCal] = sequence_Network(handles.HW, handles.Seq);
  handles.HW.tRepInit=0.05;
  set(handles.radiobutton_measure,'Value',1)
  drawnow
  guidata(hObject, handles);
  uipanel_cal_SelectionChangeFcn(hObject, eventdata, handles)
  if doCal && ~isempty(NetworkCal)
    if evalin('base', 'exist(''HW'', ''var'')') && ~isa(handles.HW, 'PD.HWClass')
      HW = evalin('base', 'HW');
      HW.NetworkCal = NetworkCal;
      assignin('base', 'HW', HW);
    else
      handles.HW.NetworkCal = NetworkCal;
    end
  end
end
handles.HW.tRepInit=0.3;
handles.Running=0;
% handles.talker.Dispose;
guidata(hObject, handles);



% --- Executes on button press in checkbox_blankTx2.
function checkbox_blankTx2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_blankTx2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_blankTx2


% --- Executes on button press in checkbox_autorun.
function checkbox_autorun_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_autorun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_autorun
if get(handles.checkbox_autorun, 'Value')
    pushbutton_measure_Callback(hObject, eventdata, handles)
end



function edit_cableLength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cableLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cableLength as text
%        str2double(get(hObject,'String')) returns contents of edit_cableLength as a double


% --- Executes during object creation, after setting all properties.
function edit_cableLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cableLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_cal.
function uipanel_cal_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_cal
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%  EventName: string 'SelectionChanged' (read only)
%  OldValue: handle of the previously selected object or empty if none was selected
%  NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% handles.Cal.Cal=~get(handles.radiobutton_measure,'Value');
% handles.Cal.short=get(handles.radiobutton_short,'Value');
% handles.Cal.open=get(handles.radiobutton_open,'Value');
% handles.Cal.terminated=get(handles.radiobutton_term,'Value');

if get(handles.radiobutton_measure, 'Value')
    set(handles.pushbutton_measure, 'String', 'Update')
end
if get(handles.radiobutton_short, 'Value')
    set(handles.pushbutton_measure, 'String', 'Short')
end
if get(handles.radiobutton_open, 'Value')
    set(handles.pushbutton_measure, 'String', 'Open')
end
if get(handles.radiobutton_term,'Value')
    set(handles.pushbutton_measure, 'String', 'Term')
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function checkbox_autorun_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_autorun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% if get(handles.checkbox_autorun, 'Value')
%     pushbutton4_Callback(hObject, eventdata, handles)
% end


% --- Executes when figure_NetworkAnalyzer is resized.
function figure_NetworkAnalyzer_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure_NetworkAnalyzer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% user input area
newSize = get(handles.figure_NetworkAnalyzer, 'Position');
inputPosition = get(handles.uipanel_input, 'Position');
inputPosition = [newSize(3)-inputPosition(3)-5, newSize(4)-inputPosition(4)-5, inputPosition(3:4)];
set(handles.uipanel_input, 'Position', inputPosition)
% axes1
panelAx1Position = [0, 0, inputPosition(1)-1, newSize(4)];
set(handles.uipanel_ax1, 'Position', panelAx1Position)
set(handles.axes1, 'OuterPosition', [-60, -30, panelAx1Position(3)+80, panelAx1Position(4)-14])
popup1Position = get(handles.popupmenu_axes1, 'Position');
popup1Position(2) = panelAx1Position(4)-29;
if panelAx1Position(3) < 200+popup1Position(1) +10
  popup1Position(3) = panelAx1Position(3) - popup1Position(1) - 10;
else
  popup1Position(3) = 200;
end
set(handles.popupmenu_axes1, 'Position', [popup1Position(1), panelAx1Position(4)-29, popup1Position(3:4)]);
% axes 2
panelAx2Position = [newSize(3)-inputPosition(3), 0, inputPosition(3), newSize(4)-inputPosition(4)];
set(handles.uipanel_ax2, 'Position', panelAx2Position)
set(handles.axes2, 'OuterPosition', [-75, -30, panelAx2Position(3)+140, panelAx2Position(4)+50])


% --- Executes on button press in checkbox_TRX1.
function checkbox_TRX1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_TRX1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_TRX1
