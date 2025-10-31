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
% (C) Copyright 2012-2021 Pure Devices GmbH, Wuerzburg, Germany
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

% UIWAIT makes GUI_SpinEcho wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_SpinEcho_OutputFcn(hObject, eventdata, handles)
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
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes', 'No', 'Yes');
if ~strcmp(selection, 'Yes')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
if ~handles.Runing
    handles.refreshPlot=1;
    guidata(hObject, handles);
    pushbutton4_Callback(hObject, eventdata, handles)
end


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

% set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
if ~handles.Runing
    handles.refreshPlot=1;
    guidata(hObject, handles);
    pushbutton4_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'Update'),  handles.Update = 1;  end
handles.Update = 1;
if ~isfield(handles, 'Runing'),  handles.Runing = 0;  end  % FIXME: Fix typo "Runing" -> "Running"
if ~isfield(handles, 'Cal'),     handles.Cal.Cal = 0;  end
if ~isfield(handles.Cal, 'Cal'), handles.Cal.Cal = 0;  end
%if isfield(handles, 'talker'), talker = handles.talker;  end
mySave = [];
if isfield(handles, 'mySave'),  mySave = handles.mySave;  end
if ~isfield(handles, 'refreshPlot'),  handles.refreshPlot = 0;  end
handles.Runing = 1;
guidata(hObject, handles);
% FIXME: Support multiple MMRT devices
if isemptyfield(handles, 'iDevice'),  handles.iDevice = 1;  end

if evalin('base', 'exist(''HW'', ''var'')')
  HW = evalin('base', 'HW');
end
if exist('HW', 'var') && isa(HW, 'PD.HW')
  ResetStructs;
else
  LoadSystem;
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
handles.Seq.firtstTR = 1;
if ~isfield(handles, 'SeqOut'),  handles.SeqOut = [];  end
while ((get(handles.checkbox2, 'Value') && get(handles.radiobutton2,'Value')) || handles.Update)
  handles.Seq.plotAllHandle = [];
  handles.Seq.plotFidHandle = [];
  handles.Seq.plotEchoesHandle = [];
  handles.Seq.plotMaxEchoesHandle = [];
  handles.Update = 0;
  handles.Seq.tOffset = 0;

  handles.Seq.DigitalIOEnable = 1;
  if handles.Seq.DigitalIOEnable
    handles.Seq.DigitalIO.SetTime = [-1e-6;  0;  1e-6 ];  % DigitalIO trigger output at center of excitation pulse
    handles.Seq.DigitalIO.SetValue = [   0;  1;     0 ];  % DigitalIO output 1 raising edge at center of excitation pulse
  end

  if ~isfield(handles, 'SeqOut'),handles.SeqOut.StartSequenceTime=now*24*3600;end
  if (get(handles.checkbox2, 'Value') && ~handles.Seq.firtstTR)
    if str2double(get(handles.edit7, 'String'))*1e-3>2.5;
      % handles.Seq.StartSequenceTime = handles.SeqOut.StartSequenceTime + str2double(get(handles.edit7, 'String'))*1e-3;
      handles.Seq.tOffset = handles.HW.tRepInit;
    end
    handles.Seq.StartSequenceTime = [];
    handles.Seq.Reinitialize = 0;
    handles.Seq.TimeFromLastSequence = str2double(get(handles.edit7, 'String'))*1e-3 - sum(handles.SeqOut.tRep);
    handles.Seq.TimeToNextSequence = str2double(get(handles.edit7, 'String'))*1e-3;
  else
    handles.Seq.StartSequenceTime = [];
    handles.Seq.Reinitialize = 1;
    handles.Seq.TimeFromLastSequence = [];
    handles.Seq.TimeToNextSequence = str2double(get(handles.edit7, 'String'))*1e-3;
  end

  handles.Seq.firtstTR = 0;
  if ~handles.refreshPlot
    handles.Seq.notCalibrate=double(~get(handles.checkbox3, 'Value'));
    if isnan(str2double(get(handles.edit6, 'String'))*1e6);
      [handles.HW, handles.mySave] = Find_Frequency_Sweep(handles.HW, handles.mySave, 600*handles.Seq.notCalibrate);
      assignin('base', 'mySave', handles.mySave);
      if ~isa(HW, 'PD.HW'),  assignin('base', 'HW', handles.HW);  end
      figure(handles.figure1);
      set(handles.edit6, 'String', ['auto ', num2str(handles.HW.fLarmor/1e6, '%3.6f')]);
    else
      handles.HW.fLarmor = str2double(get(handles.edit6, 'String'))*1e6;
      handles.HW.B0 = handles.HW.fLarmor/handles.HW.GammaDef*2*pi;
      if ~handles.Seq.notCalibrate
        [handles.HW, handles.mySave] = Find_Frequency_Sweep(handles.HW, handles.mySave, 600*handles.Seq.notCalibrate);
        assignin('base', 'mySave', handles.mySave);
        if ~isa(HW, 'PD.HW'),  assignin('base', 'HW', handles.HW);  end
        figure(handles.figure1);
        set(handles.edit6, 'String', ['auto ', num2str(handles.HW.fLarmor/1e6,'%3.6f')]);
      end
    end
  end

  popup_sel_String = get(handles.pushbutton4, 'String');
  switch popup_sel_String
    case 'B1+'
      if isnan(str2double(get(handles.edit1, 'String'))*1e-6);
        [handles.HW, handles.mySave] = Find_PulseDuration(handles.HW, handles.mySave, 0, 1);
      else
        [handles.HW, handles.mySave] = Find_PulseDuration(handles.HW, handles.mySave, 0, 1, 3, get(handles.edit1, 'String')*1e-6);
      end
      if exist('HW', 'var') && isa(HW, 'PD.HW')
        ResetStructs;
      else
        LoadSystem;
      end
      handles.HW = HW;
      handles.mySave = mySave;
      assignin('base', 'mySave', handles.mySave);
      if ~isa(HW, 'PD.HW'),  assignin('base', 'HW', handles.HW);  end
      set(handles.edit1, 'String', ['def ', num2str(handles.HW.tFlip90Def*1e6, '%3.3f')]);
      set(handles.edit2, 'String', ['def ', num2str(handles.HW.tFlip180Def*1e6, '%3.3f')]);
      handles.Seq.p90 = handles.HW.tFlip90Def;
      handles.Seq.p180 = handles.HW.tFlip180Def;
      figure(handles.figure1);

    case 'B0'
      t = handles.HW.FindFrequencyPlot;
      handles.HW.FindFrequencyPlot = 1;
      [handles.HW, handles.mySave] = Find_Frequency_Sweep(handles.HW, handles.mySave, 0, [], 1);
      assignin('base', 'mySave', handles.mySave);
      if ~isa(HW, 'PD.HW'),  assignin('base', 'HW', handles.HW);  end
      figure(handles.figure1);
      handles.HW.FindFrequencyPlot = t;
      set(handles.edit6, 'String', ['auto ', num2str(handles.HW.fLarmor/1e6,'%3.6f')]);

    case 'Shim'
      [handles.HW, handles.mySave] = Find_Shim(handles.HW, handles.mySave, 0, 1);
      assignin('base', 'mySave', handles.mySave);
      if ~isa(HW, 'PD.HW'),  assignin('base', 'HW', handles.HW);  end

    otherwise

      popup_sel_index = get(handles.popupmenu1, 'Value');
      switch popup_sel_index
        case 1
            handles.Seq.plot = 0;
          handles.Seq.plotAllHandle = handles.axes3;
          if popup_sel_index == get(handles.popupmenu2, 'Value')
            set(handles.popupmenu2, 'Value', 2);
          end

        case 2
          handles.Seq.plotFidHandle = handles.axes3;
          if popup_sel_index == get(handles.popupmenu2, 'Value')
            set(handles.popupmenu2, 'Value', 1);
          end

        case 3
          handles.Seq.plotEchoesHandle = handles.axes3;
          if popup_sel_index == get(handles.popupmenu2, 'Value')
            set(handles.popupmenu2, 'Value', 1);
          end

        case 4
          handles.Seq.plotMaxEchoesHandle = handles.axes3;
          if popup_sel_index == get(handles.popupmenu2, 'Value')
            set(handles.popupmenu2, 'Value', 1);
          end

      end

      popup_sel_index = get(handles.popupmenu2, 'Value');
      switch popup_sel_index
        case 1
          handles.Seq.plotAllHandle = handles.axes2;
        case 2
          handles.Seq.plotFidHandle = handles.axes2;
        case 3
          handles.Seq.plotEchoesHandle = handles.axes2;
        case 4
          handles.Seq.plotMaxEchoesHandle = handles.axes2;
      end


      if ~handles.refreshPlot
        if isnan(str2double(get(handles.edit1, 'String')))
          set(handles.edit1, 'String', ['def ', num2str(handles.HW.tFlip90Def*1e6, '%3.3f')])
          handles.Seq.p90 = handles.HW.tFlip90Def;
        else
          handles.Seq.p90 = str2double(get(handles.edit1, 'String'))/1e6;
        end

        if isnan(str2double(get(handles.edit2, 'String')))
          set(handles.edit2, 'String', ['def ', num2str(handles.HW.tFlip180Def*1e6,'%3.3f')])
          handles.Seq.p180 = handles.HW.tFlip180Def;
        else
          handles.Seq.p180 = str2double(get(handles.edit2, 'String'))/1e6;
        end

        if isnan(str2double(get(handles.edit3, 'String')))
          set(handles.edit3, 'String', ['def ', num2str(10,'%3.0f')])
          handles.Seq.tEcho = 10e-3;
        else
          handles.Seq.tEcho = str2double(get(handles.edit3, 'String'))/1e3;
        end

        if isnan(str2double(get(handles.edit4, 'String')))
          set(handles.edit4, 'String', ['def ', num2str(1,'%3.0f')])
          handles.Seq.nEchos = 1;
        else
          handles.Seq.nEchos = str2double(get(handles.edit4, 'String'));
        end

        if get(handles.checkbox1, 'Value')
          handles.Seq.B1Amp2Str = 1e6;
          set(handles.text10, 'String', [char(181) 'T'])
        else
          handles.Seq.B1Amp2Str = handles.HW.GammaDef/2/pi/1000;
          set(handles.text10, 'String', 'kHz')
        end

        if isnan(str2double(get(handles.edit5, 'String')));
          set(handles.edit5, 'String', ['def ', num2str(handles.HW.TX(handles.iDevice).AmpDef*handles.Seq.B1Amp2Str,'%3.6f')])
          handles.Seq.TXAmp = handles.HW.TX(handles.iDevice).AmpDef;
        else
          handles.Seq.TXAmp = str2double(get(handles.edit5, 'String'))/handles.Seq.B1Amp2Str;
        end


        if (handles.Seq.nEchos<=500 && (handles.Seq.nEchos+1)*handles.Seq.tEcho<=200)
          handles.Seq.fast = 1;
        else
          handles.Seq.fast = 0;
          handles.Seq.tOffset = handles.Seq.tOffset + zeros(1, handles.Seq.nEchos);
        end

        handles.Seq.fSample = max((1/(handles.Seq.tEcho*0.1-handles.Seq.p90/2-handles.HW.TX(handles.iDevice).BlankOffset-5e-6)*10), ...
          handles.HW.RX(handles.iDevice).fSample/6250);
        if (handles.Seq.fSample > handles.HW.RX(handles.iDevice).fSample/125*2) && (handles.Seq.nEchos>100);
          handles.Seq.AQFID = 0.5*(handles.HW.RX(handles.iDevice).fSample/125*2)/handles.Seq.fSample;
          handles.Seq.AQEcho = 0.5*(handles.HW.RX(handles.iDevice).fSample/125*2)/handles.Seq.fSample;
          handles.Seq.fSample = max((1/(handles.Seq.tEcho*0.1-handles.Seq.p90/2-handles.HW.TX(handles.iDevice).BlankOffset-5e-6)*10), ...
            handles.HW.RX(handles.iDevice).fSample/6250);
        end

        if handles.Seq.nEchos==0,  handles.Seq.AQFID = 1;  end
        handles.Runing = 1;
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
        iAQ = find([handles.SeqOut.AQ(:).Channel] == Channel & [handles.SeqOut.AQ(:).Device] == handles.iDevice, 1, 'first');
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
            hold(handles.Seq.plotEchoesHandle, 'all');
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
          end;
          if ~isfield(handles.mySave, 'plotMaxEchoesHandleMin')
            handles.mySave.plotMaxEchoesHandleMin = handles.mySave.plotMaxEchoesHandleMinNew;
          end;
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


      % handles.HW.ReInit = 0;
      % handles.HW.tRepInit = 0.005;
  end

  set(handles.radiobutton2, 'Value', 1);
  handles.Seq.firtstTR = 0;
  guidata(hObject, handles);
  uipanel1_SelectionChangeFcn(hObject, eventdata, handles);
end
% handles.HW.tRepInit = 0.3;
% handles.HW.ReInit = 1;
handles.Runing = 0;
guidata(hObject, handles);



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
        if get(handles.checkbox1, 'Value')
            handles.Seq.B1Amp2Str=1e6;
            set(handles.text10, 'String', [char(181) 'T'])
        else
            handles.Seq.B1Amp2Str=handles.HW.GammaDef/2/pi/1000;
            set(handles.text10, 'String', 'kHz')
        end

        if isfield(handles,'Seq')
            if isfield(handles.Seq,'TXAmp')
                if isfield(handles.Seq,'B1Amp2Str')
                    if isnan(str2double(get(handles.edit5, 'String')));
                        set(handles.edit5, 'String', ['def ', num2str(handles.Seq.TXAmp*handles.Seq.B1Amp2Str,'%3.6f')])
                    else
                        set(handles.edit5, 'String', num2str(handles.Seq.TXAmp*handles.Seq.B1Amp2Str,'%3.6f'))
                    end
                end
            end
        end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
if get(handles.checkbox2, 'Value')
    pushbutton4_Callback(hObject, eventdata, handles)
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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

if get(handles.radiobutton2,'Value'); set(handles.pushbutton4, 'String','Update');end
if get(handles.radiobutton3,'Value'); set(handles.pushbutton4, 'String','B1+');end
if get(handles.radiobutton4,'Value'); set(handles.pushbutton4, 'String','B0');end
if get(handles.radiobutton5,'Value'); set(handles.pushbutton4, 'String','Shim');end

guidata(hObject, handles);



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
