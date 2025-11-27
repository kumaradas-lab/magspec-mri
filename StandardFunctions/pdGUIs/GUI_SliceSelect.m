function varargout = GUI_SliceSelect(varargin)
%% Interface for graphical slice selection
%
%       GUI_SliceSelect()
%
% This GUI can be used to define data used in the AQSlice structure of the most
% common sequences. It can also perform basic (Turbo) Spin Echo and Gradient
% Echo (FLASH) imaging measurements.
%
% It allows the export of the slice configuration data to a file. A possible
% work-flow for that task is:
%   * Open the GUI after configuring the MMRT (e.g. with
%     Demo_Auto_ParameterSearch) by executing the following command:
%           GUI_SliceSelect()
%   * If desired, acquire a 3d localizer:
%       - Select either "Spin Echo" or "Grad Echo" in the "Sequence" popupmenu.
%       - Set the dimension to 3d.
%       - Make sure, the orientation angles are each set to "0" and the slice is
%         centered.
%       - Select a reasonable image size and resolution.
%       - Hit the "Localizer 3d" button and wait until the measurement has
%         concluded.
%   * Select the sequence "Slice Selection" in the "Sequence" popupmenu.
%   * Use the interface controls to set your desired slice configuration.
%   * Choose an identifier for the slice (must consist of characters allowed in
%     Matlab function names).
%   * Hit the "Save AQSlice" button.
%   * A message will be printed to the Command Window stating the name and
%     location of the file with the slice configuration data (e.g.
%     ...\AQSlice_test.m).
%   * You can either use this slice selection data
%       1.  as a function in your sequence definition:
%           e.g.:     Seq.AQSlice(1) = AQSlice_test();
%     or
%       2.  by copying the respective lines and pasting them directly in your
%           sequence definition file. In this case, iSlice must be set in that
%           file (e.g. "iSlice = 1;").
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% Last Modified by GUIDE v2.5 16-Apr-2018 10:00:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_SliceSelect_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_SliceSelect_OutputFcn, ...
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


% --- Executes just before GUI_SliceSelect is made visible.
function GUI_SliceSelect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_SliceSelect (see VARARGIN)

% Choose default command line output for GUI_SliceSelect
if isfield(handles,'mySave')
    HW=handles.HW;
    mySave=handles.mySave;
    Seq=handles.Seq;
    AQ=handles.AQ;
    TX=handles.TX;
    Grad=handles.Grad;
    %talker=handles.talker;
elseif evalin('base', 'exist(''HW'', ''var'')')
    HW = evalin('base', 'HW');
end

if exist('HW', 'var') && isa(HW, 'PD.HWClass')
    if evalin('base', 'exist(''mySave'', ''var'')')
        mySave = evalin('base', 'mySave');
    end
    ResetStructs
else
    LoadSystem
end

handles.HW=HW;
handles.mySave=mySave;
handles.Seq=Seq;
handles.AQ=AQ;
handles.TX=TX;
handles.Grad=Grad;
%handles.talker=talker;

handles.output = hObject;
handles.ButtonDown=0;
handles.MotionFcnBusy=0;
handles.axes1Rot3D=0;
handles.ButtonCenter =0;
handles.ButtonCenterXYZ =0;
handles.ButtonZeroapt=0;
handles.ButtonMeasureSlice=0;
handles.alfa=0;
handles.phi=0;
handles.theta=0;
handles.axes1update=1;
% handles.AqWin1_image=peaks(50);

handles.Seq.AQSlice.Center=[0,0,0];
handles.Seq.AQSlice.normX=[1,0,0];
handles.Seq.AQSlice.normY=[0,1,0];
handles.Seq.AQSlice.normZ=[0,0,1];
handles.Seq.AQSlice.ReceivingPoint2Center=[0,0,0];
handles.Seq.AQSlice.ReceivingPoint=[0,0,0];
handles.Seq.AQSlice.R=0;
handles.nRead=1;
handles.nPhase(1)=1;
handles.nPhase(2)=1;
handles.nPhase(3)=1;
handles.sizeRead=1;
handles.sizePhase(1)=1;
handles.sizePhase(2)=1;
handles.sizePhase(3)=1;
handles.refresh3D=1;
handles.CenterOfImageOffset=[0,0,0];
handles.ButtonMeasure3DLocaliser=0;
handles.CenterOffsetLocaliser=handles.HW.Grad.ImageVolOffset;
handles.LocaliserImageVol=handles.HW.Grad.ImageVol+handles.CenterOffsetLocaliser([1,1,2,2,3,3]);
handles.ButtonReadPhase=1;
handles.ButtonDim=1;
set(handles.edit21,'String',[handles.HW.RootPath, '\User\Save']);

     handles.Vpermut=[ 2 3 1];
     [~,handles.VpermutSort]=sort(handles.Vpermut);

% test=cell2mat(varargin(2));
%             handles.Seq.AQSlice.sizeRead=0.01;
%             handles.Seq.AQSlice.sizePhase=0.01;
%             handles.Seq.AQSlice.R=0;
%             handles.Seq.AQSlice.nRead=50;
%             handles.Seq.AQSlice.nPhase=50;

%             handles.Seq.AQSlice.alfa=get(handles.slider_alfa,'Value')*2*pi;
%             handles.Seq.AQSlice.phi=get(handles.slider_phi,'Value')*2*pi;
%             handles.Seq.AQSlice.theta=get(handles.slider_theta,'Value')*2*pi;


%             if ~isempty(whos('global','image3D'))
%                 global image3D
%                 handles.HW=image3D.HW;
%                 %handles.talker=image3D.talker;
%                 handles.mySave=image3D.mySave;
%                 handles.nx=image3D.AQSlice.nRead;
%                 handles.ny=image3D.AQSlice.nPhase;
%                 handles.nz=image3D.AQSlice.nPhase3D;
%                 [handles.x,handles.y,handles.z] =meshgrid(image3D.xV,image3D.yV,image3D.zV);
%                 handles.v = abs(image3D.dataYXZ);
%             else
%                 handles.LocaliserImageVol=[-0.01,0.01,-0.01,0.01,-0.01,0.01];%[xmin xmax ymin ymax zmin zmax]

            handles.nx=64;
            handles.ny=64;
            handles.nz=64;
            [handles.x,handles.y,handles.z] = meshgrid( linspace(handles.LocaliserImageVol(1),handles.LocaliserImageVol(2),handles.nx),...
                                                        linspace(handles.LocaliserImageVol(3),handles.LocaliserImageVol(4),handles.ny),...
                                                        linspace(handles.LocaliserImageVol(5),handles.LocaliserImageVol(6),handles.nz));
            handles.v=ones(size(handles.x));
            handles.v(((handles.x-handles.CenterOffsetLocaliser(1)).^2+(handles.z-handles.CenterOffsetLocaliser(3)).^2).^0.5>handles.HW.Grad.TubeDiameter/2)=0;
            handles.v(((handles.y-handles.CenterOffsetLocaliser(2))<0)&((handles.x-handles.CenterOffsetLocaliser(1)).^2+(handles.y-handles.CenterOffsetLocaliser(2)).^2+(handles.z-handles.CenterOffsetLocaliser(3)).^2).^0.5>handles.HW.Grad.TubeDiameter/2)=0;
            handles.lim=[min(handles.v(:)) max(handles.v(:))];
%
%             pushbutton1_Callback(hObject, eventdata, handles)

            figure(handles.figureSliceSelect)
            cla(handles.axes1)
%             set(handles.figureSliceSelect,'Renderer','OpenGL')
%             opengl hardware
            view(handles.axes1,3)
            xlabel(handles.axes1,'z')
            ylabel(handles.axes1,'x')
            hzl = zlabel(handles.axes1,'y');
            set(hzl, 'Visible', 'on');
            grid(handles.axes1,'on')
           guidata(hObject, handles);

%       set(handles.axes1,'Projection','perspective')

% %   activelabelhandle(handles.axes1,'xlabel', 'X');
% %   activelabelhandle(handles.axes1,'ylabel', 'Y');
% %   activelabelhandle(handles.axes1,'zlabel', 'Z');

            colormap(handles.axes1,'gray')
%             shading(handles.axes1, 'flat')
            % shading interp
            axis(handles.axes1,'equal')
            axis(handles.axes1,handles.LocaliserImageVol([5,6,1,2,3,4]))
            if handles.axes1Rot3D==0
                rotate3d(handles.axes1,'on');
                handles.axes1Rot3D=1;
            end


% Update handles structure
% handles.axes1update=1;
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using GUI_SliceSelect.
% if strcmp(get(hObject,'Visible'),'off')
%     figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles);
% end
% popupmenu1_Callback(hObject, eventdata, handles)
popupmenu_sequenceSelect_Callback(hObject, eventdata, handles);
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles);

% UIWAIT makes GUI_SliceSelect wait for user response (see UIRESUME)
% uiwait(handles.figureSliceSelect);

end


% --- Outputs from this function are returned to the command line.
function varargout = GUI_SliceSelect_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles;
varargout{1} = handles.output;

end


% --- Executes on button press in pushbutton_reset3d.
function pushbutton_reset3d_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%             if ~isempty(whos('global','image3D'))
%                 global image3D
%                 handles.HW=image3D.HW;
%                 handles.nx=image3D.AQSlice.nRead;
%                 handles.ny=image3D.AQSlice.nPhase;
%                 handles.nz=image3D.AQSlice.nPhase3D;
%                 [handles.x,handles.y,handles.z] =meshgrid(image3D.xV,image3D.yV,image3D.zV);
%                 handles.v = abs(image3D.dataYXZ);
%             else
            handles.CenterOffsetLocaliser=handles.HW.Grad.ImageVolOffset;
            handles.LocaliserImageVol=handles.HW.Grad.ImageVol+handles.CenterOffsetLocaliser([1,1,2,2,3,3]);

            handles.nx=64;
            handles.ny=64;
            handles.nz=64;
            [handles.x,handles.y,handles.z] = meshgrid( linspace(handles.LocaliserImageVol(1),handles.LocaliserImageVol(2),handles.nx)+handles.CenterOffsetLocaliser(1),...
                                                        linspace(handles.LocaliserImageVol(3),handles.LocaliserImageVol(4),handles.ny)+handles.CenterOffsetLocaliser(2),...
                                                        linspace(handles.LocaliserImageVol(5),handles.LocaliserImageVol(6),handles.nz)+handles.CenterOffsetLocaliser(3));
            handles.v=ones(size(handles.x));
            handles.v(((handles.x-handles.CenterOffsetLocaliser(1)).^2+(handles.z-handles.CenterOffsetLocaliser(3)).^2).^0.5>handles.HW.Grad.TubeDiameter/2)=0;
            handles.v(((handles.y-handles.CenterOffsetLocaliser(2))<0)&((handles.x-handles.CenterOffsetLocaliser(1)).^2+(handles.y-handles.CenterOffsetLocaliser(2)).^2+(handles.z-handles.CenterOffsetLocaliser(3)).^2).^0.5>handles.HW.Grad.TubeDiameter/2)=0;
            handles.lim=[min(handles.v(:)) max(handles.v(:))];

%             end
%             handles.axes1update=1;
            handles.refresh3D=1;
            axis(handles.axes1,handles.LocaliserImageVol([5,6,1,2,3,4]))
            guidata(hObject, handles);
%             popupmenu1_Callback(hObject, eventdata, handles)
             figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles);


% if handles.axes1Rot3D==0
%     rotate3d(handles.axes1,'on');
%     handles.axes1Rot3D=1;
% else
%     rotate3d(handles.axes1,'off');
%     handles.axes1Rot3D=0;
% end
% guidata(hObject,handles)

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
printdlg(handles.figureSliceSelect)

end


% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% selection = questdlg(['Close ' get(handles.figureSliceSelect,'Name') '?'],...
%                      ['Close ' get(handles.figureSliceSelect,'Name') '...'],...
%                      'Yes','No','Yes');
% if strcmp(selection,'No')
%     return;
% end
%
% delete(handles.figureSliceSelect)

end


% --- Executes on selection change in popupmenu_sequenceSelect.
function popupmenu_sequenceSelect_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_sequenceSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_sequenceSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_sequenceSelect

popup_sel_index = get(handles.popupmenu_sequenceSelect, 'Value');

% common settings (might be over-written below)
set(handles.pushbutton_Localizer, 'Visible', 'on');
% averages
set(handles.text30, 'Visible', 'on');
set(handles.text31, 'Visible', 'on');
set(handles.text37, 'Visible', 'on');
set(handles.edit13, 'Visible', 'on');
set(handles.text33, 'Visible', 'on');
set(handles.edit14, 'Visible', 'on');
set(handles.text34, 'Visible', 'on');
set(handles.text35, 'Visible', 'on');
% tEcho
set(handles.text17, 'Visible', 'on');
set(handles.edit17, 'Visible', 'on');
set(handles.text39, 'Visible', 'on');
% tRep
set(handles.text18, 'Visible', 'on');
set(handles.edit18, 'Visible', 'on');
set(handles.text40, 'Visible', 'on');
% Spoil
set(handles.text19, 'Visible', 'on');
set(handles.edit19, 'Visible', 'on');
set(handles.text41, 'Visible', 'on');
% TurboFactor
set(handles.text20, 'Visible', 'on');
set(handles.edit20, 'Visible', 'on');
set(handles.text42, 'Visible', 'on');
% Save Path
set(handles.text21, 'Visible', 'on');
set(handles.edit21, 'Visible', 'on');
set(handles.pushbutton_Pathselect, 'Visible', 'on');
% Filename
set(handles.text22, 'String', 'Filename');
set(handles.checkbox_addDate, 'Visible', 'on');
set(handles.checkbox_Save, 'Visible', 'on');
set(handles.pushbutton_startMeasurement, 'String', 'Start Measurement');
switch popup_sel_index
  case PD.SliceSelectType.SpinEcho  % Spin-Echo
    set(handles.text20,'String','Turbo Factor');
    set(handles.text42,'String','x');
    set(handles.edit20,'String','1');
    set(handles.checkbox_Extra1,'Visible','off')
    set(handles.checkbox_Extra2,'Visible','off')
    set(handles.popupmenu_InversionPulse,'Visible','on')
    set(handles.text38,'Visible','on')
  case PD.SliceSelectType.GradEcho % Grad-Echo
    set(handles.text20,'String','T1');
    set(handles.text42,'String','s');
    set(handles.edit20,'String','0.1');
    set(handles.text18,'String','T Rep');
    set(handles.checkbox_Extra1,'Visible','on')
    set(handles.checkbox_Extra2,'Visible','on')
    set(handles.popupmenu_InversionPulse,'Visible','off')
    set(handles.text38,'Visible','off')
  case PD.SliceSelectType.SliceSelect % Slice Selection
    set(handles.pushbutton_Localizer, 'Visible', 'off');
    % averages
    set(handles.text30, 'Visible', 'off');
    set(handles.text31, 'Visible', 'off');
    set(handles.text37, 'Visible', 'off');
    set(handles.edit13, 'Visible', 'off');
    set(handles.text33, 'Visible', 'off');
    set(handles.edit14, 'Visible', 'off');
    set(handles.text34, 'Visible', 'off');
    set(handles.text35, 'Visible', 'off');
    % tEcho
    set(handles.text17, 'Visible', 'off');
    set(handles.edit17, 'Visible', 'off');
    set(handles.text39, 'Visible', 'off');
    % tRep
    set(handles.text18, 'Visible', 'off');
    set(handles.edit18, 'Visible', 'off');
    set(handles.text40, 'Visible', 'off');
    % Spoil
    set(handles.text19, 'Visible', 'off');
    set(handles.edit19, 'Visible', 'off');
    set(handles.text41, 'Visible', 'off');
    % TurboFactor
    set(handles.text20, 'Visible', 'off');
    set(handles.edit20, 'Visible', 'off');
    set(handles.text42, 'Visible', 'off');
    % Save Path
    set(handles.text21, 'Visible', 'off');
    set(handles.edit21, 'Visible', 'off');
    set(handles.pushbutton_Pathselect, 'Visible', 'off');
    % Filename
    set(handles.text22, 'String', 'Slice Identifier');
    set(handles.checkbox_addDate, 'Visible', 'off');
    set(handles.checkbox_Save, 'Visible', 'off');
    set(handles.pushbutton_startMeasurement, 'String', 'Save AQSlice');

end

end


% --- Executes during object creation, after setting all properties.
function popupmenu_sequenceSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_sequenceSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', char(sort(enumeration('PD.SliceSelectType'))));       %%% Sequenzen in SliceSelectType.m eintragen %%%%%%%%%%%%%%%%%%%%%%%%

end


% --- Executes on slider movement.
function slider_alfa_Callback(hObject, eventdata, handles)
% hObject    handle to slider_alfa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function slider_alfa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_alfa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over slider_alfa.
function slider_alfa_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to slider_alfa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes1
%             handles.LocaliserImageVol=[-0.01,0.01,-0.01,0.01,-0.01,0.01];%[xmin xmax ymin ymax zmin zmax]
%
%             handles.nx=64;
%             handles.ny=64;
%             handles.nz=64;
%             [handles.x,handles.y,handles.z] = meshgrid(linspace(handles.LocaliserImageVol(1),handles.LocaliserImageVol(2),handles.nx),linspace(handles.LocaliserImageVol(3),handles.LocaliserImageVol(4),handles.ny),linspace(handles.LocaliserImageVol(5),handles.LocaliserImageVol(6),handles.nz));
%             handles.tubeDiameter=8e-3;
%             handles.v=ones(size(handles.x));
%             handles.v((handles.x.^2+handles.y.^2).^0.5>handles.tubeDiameter/2)=0;
%             handles.v((handles.z<0)&(handles.x.^2+handles.y.^2+handles.z.^2).^0.5>handles.tubeDiameter/2)=0;
%             handles.lim=[min(handles.v(:)) max(handles.v(:))];
%
%             figure(handles.figureSliceSelect)
% %       set(handles.figureSliceSelect,'Renderer','OpenGL')
% %             cla(handles.axes1)
% %             axes(handles.axes1)
%            handles.axes1=axes('units','normal','pos',[.18  .24 .31 .53],'box','on',...
%           'ylim',[handles.y(1) handles.y(end)],...
%           'xlim',[handles.x(1) handles.x(end)],...
%           'zlim',[handles.z(1) handles.z(end)],...
%           'clim',handles.lim);
%             view(handles.axes1,3)
%             view(3)
% %             view(handles.axes1,[0)
% %             colorbar
% %             xlabel(handles.axes1,'z')
% %             ylabel(handles.axes1,'x')
% %             zlabel(handles.axes1,'y')
%             xlabel('z')
%             ylabel('x')
%             zlabel('y')

end


% --- Executes on mouse motion over figure - except title and menu.
function figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figureSliceSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mousexy=get(handles.figureSliceSelect, 'CurrentPoint');
if (mousexy(1)<=101 ) && handles.MotionFcnBusy==0
  handles.MotionFcnBusy = 1;
  guidata(hObject, handles);
%   pause(0.001)
  try
    if or(handles.ButtonMeasureSlice,handles.ButtonMeasure3DLocaliser)
      if isfield(handles,'mySave')
        HW=handles.HW;
        mySave=handles.mySave;
        Seq=handles.Seq;
        AQ=handles.AQ;
        TX=handles.TX;
        Grad=handles.Grad;
        %talker=handles.talker;
      end
      if exist('HW', 'var') && isa(HW, 'PD.HWClass')
        if evalin('base', 'exist(''mySave'', ''var'')')
          mySave = evalin('base', 'mySave');
        end
        ResetStructs
      else
        LoadSystem
      end
      handles.HW=HW;
      handles.mySave=mySave;
      handles.Seq=Seq;
      handles.AQ=AQ;
      handles.TX=TX;
      handles.Grad=Grad;
      %handles.talker=talker;
    end


    if handles.ButtonDim==1
      handles.ButtonDim=0;
%       if get(handles.radiobutton_2D,'Value')
%         set(handles.checkbox_Smooth,'Visible','on')
%       else
%         set(handles.checkbox_Smooth,'Visible','off')
%       end
      if get(handles.radiobutton_0D,'Value')
        set(handles.edit1,'Enable','off');
        set(handles.edit1,'String','1');
        set(handles.edit2,'Enable','off');
        set(handles.edit2,'String','1');
        set(handles.popROS,'Enable','off');
        set(handles.popPOS,'Enable','off');
        set(handles.popupmenu_PhaseOS2,'Enable','off');
        set(handles.popupmenu8,'Enable','off');
        set(handles.edit25,'Enable','off');
        set(handles.edit25,'String','1');
        set(handles.edit26,'Enable','off');
        set(handles.edit26,'String','1');

      elseif get(handles.radiobutton_1D,'Value')
        handles.ButtonReadPhase=1;
        set(handles.edit1,'Enable','off');
        set(handles.edit2,'Enable','off');
        set(handles.edit2,'String','1');
        set(handles.popROS,'Enable','off');
        set(handles.popPOS,'Enable','off');
        set(handles.popupmenu_PhaseOS2,'Enable','off');
        set(handles.popupmenu8,'Enable','off');
        set(handles.edit25,'Enable','off');
        set(handles.edit25,'String','1');
        set(handles.edit26,'Enable','off');
        set(handles.edit26,'String','1');

      elseif get(handles.radiobutton_2D,'Value')
        handles.ButtonReadPhase=1;
        set(handles.edit1,'Enable','off');
        set(handles.edit2,'Enable','off');
        set(handles.edit2,'String','1');
        set(handles.popROS,'Enable','off');
        set(handles.popPOS,'Enable','off');
        set(handles.popupmenu_PhaseOS2,'Enable','on');
        set(handles.popupmenu_PhaseOS2,'Value',2);
        set(handles.popupmenu8,'Enable','off');
%         set(handles.edit6,'String','5');
        set(handles.edit25,'Enable','on');
        set(handles.edit25,'String','32');
        set(handles.edit26,'Enable','off');

      elseif get(handles.radiobutton_3D,'Value')
        handles.ButtonReadPhase=1;
        set(handles.edit1,'Enable','off');
        set(handles.edit2,'Enable','on');
        set(handles.edit2,'String','16');
        set(handles.popROS,'Enable','off');
        set(handles.popPOS,'Enable','on');
        set(handles.popupmenu_PhaseOS2,'Enable','on');
        set(handles.popupmenu_PhaseOS2,'Value',2);
        set(handles.popupmenu8,'Enable','off');
        set(handles.edit25,'Enable','on');
        set(handles.edit25,'String','16');
        set(handles.edit26,'Enable','off');
        set(handles.edit6,'String','Inf');
        set(handles.edit20,'String','8');
        if get(handles.popupmenu_sequenceSelect, 'Value') == PD.SliceSelectType.SpinEcho  % Spin-Echo
          set(handles.text18,'String','Turbo Break')
        end

      elseif get(handles.radiobutton_3DCSI,'Value')
        set(handles.edit1,'Enable','off');
        set(handles.edit1,'String','1');
        set(handles.edit2,'Enable','on');
        set(handles.edit2,'String','1');
        set(handles.popROS,'Enable','off');
        set(handles.popPOS,'Enable','on');
        set(handles.popupmenu_PhaseOS2,'Enable','on');
        set(handles.popupmenu8,'Enable','on');
        set(handles.edit25,'Enable','on');
        set(handles.edit25,'String','1');
        set(handles.edit26,'Enable','on');
        set(handles.edit26,'String','1');

      end

      if str2double(get(handles.edit1,'String'))==1
        set(handles.popROS,'Value',14)
      end
      if str2double(get(handles.edit2,'String'))==1
        set(handles.popPOS,'Value',1)
      end
      if str2double(get(handles.edit25,'String'))==1
        set(handles.popupmenu_PhaseOS2,'Value',1)
      end
      if str2double(get(handles.edit26,'String'))==1
        set(handles.popupmenu8,'Value',1)
      end
    end

    if handles.ButtonReadPhase==1
      handles.ButtonReadPhase=0;
      if get(handles.radiobutton_readout,'Value')
        set(handles.edit1,'Enable','on');
        if get(handles.radiobutton_3D,'Value')
          set(handles.edit1,'String','16');
        else
          set(handles.edit1,'String','32');
        end
        set(handles.popROS,'Enable','on');
        set(handles.popROS,'Value',14)
        set(handles.edit26,'Enable','off');
        set(handles.edit26,'String','1');
        set(handles.popupmenu8,'Enable','off');
        set(handles.popupmenu8,'Value',1)

      elseif get(handles.radiobutton_only_phase,'Value')
        set(handles.edit26,'Enable','on');
        if get(handles.radiobutton_3D,'Value')
          set(handles.edit26,'String','16');
        else
          set(handles.edit26,'String','32');
        end
        set(handles.popupmenu8,'Enable','on');
        set(handles.popupmenu8,'Value',1)
        set(handles.edit1,'Enable','off');
        set(handles.edit1,'String','1');
        set(handles.popROS,'Enable','off');
        set(handles.popROS,'Value',14)
        % set(handles.edit15,'String',1e12)
      end
    end

    mousexy=get(handles.figureSliceSelect, 'CurrentPoint');
% Fixme left or right panel
    if (mousexy(1)<=101 ) && handles.axes1update
      handles.axes1update=0;
      guidata(hObject, handles);
      pause(0.001)
      % axes(handles.axes1);

      handles.Seq.AQSlice.sizeRead=get(handles.slider11,'Value')*(max(abs(handles.LocaliserImageVol(2:2:6)-handles.LocaliserImageVol(1:2:6))));
      handles.Seq.AQSlice.sizePhase(1)=get(handles.slider12,'Value')*(max(abs(handles.LocaliserImageVol(2:2:6)-handles.LocaliserImageVol(1:2:6))));
      handles.Seq.AQSlice.sizePhase(2)=get(handles.slider13,'Value')*(max(abs(handles.LocaliserImageVol(2:2:6)-handles.LocaliserImageVol(1:2:6))));
      handles.Seq.AQSlice.sizePhase(3)=get(handles.slider14,'Value')*(max(abs(handles.LocaliserImageVol(2:2:6)-handles.LocaliserImageVol(1:2:6))));
      handles.Seq.AQSlice.thickness=str2double(get(handles.edit6,'String'))/1000; % mm -> m
      if isnan(handles.Seq.AQSlice.thickness); handles.Seq.AQSlice.thickness=inf;end

      if str2double(get(handles.edit1,'String'))==1
        handles.Seq.AQSlice.sizeRead=1e12;
      end
      if str2double(get(handles.edit2,'String'))==1
        handles.Seq.AQSlice.sizePhase(1)=1e12;
      end
      if str2double(get(handles.edit25,'String'))==1
        handles.Seq.AQSlice.sizePhase(2)=1e12;
      end
      if str2double(get(handles.edit26,'String'))==1
        handles.Seq.AQSlice.sizePhase(3)=1e12;
      end

      set(handles.edit15,'String',num2str(handles.Seq.AQSlice.sizeRead));
      set(handles.edit22,'String',num2str(handles.Seq.AQSlice.sizePhase(1)));
      set(handles.edit23,'String',num2str(handles.Seq.AQSlice.sizePhase(2)));
      set(handles.edit24,'String',num2str(handles.Seq.AQSlice.sizePhase(3)));

      handles.Seq.Loops=str2double(get(handles.edit14,'String'));
      handles.Seq.average=str2double(get(handles.edit13,'String'));
      set(handles.text35,'String',num2str(handles.Seq.average*handles.Seq.Loops));

      % handles.Seq.AQSlice.R=0;
      handles.Seq.AQSlice.nRead=str2double(get(handles.edit1,'String'));
      handles.Seq.AQSlice.nPhase(1)=str2double(get(handles.edit2,'String'));
      handles.Seq.AQSlice.nPhase(2)=str2double(get(handles.edit25,'String'));
      handles.Seq.AQSlice.nPhase(3)=str2double(get(handles.edit26,'String'));
      % handles.LocaliserImageVol=[-0.01,0.01,-0.01,0.01,-0.01,0.01];%[xmin xmax ymin ymax zmin zmax]
      if  handles.nRead~=handles.Seq.AQSlice.nRead || ...
          handles.nPhase(1)~=handles.Seq.AQSlice.nPhase(1) || ...
          handles.nPhase(2)~=handles.Seq.AQSlice.nPhase(2) || ...
          handles.nPhase(3)~=handles.Seq.AQSlice.nPhase(3) || ...
          handles.sizeRead~=handles.Seq.AQSlice.sizeRead || ...
          handles.sizePhase(1)~=handles.Seq.AQSlice.sizePhase(1) || ...
          handles.sizePhase(2)~=handles.Seq.AQSlice.sizePhase(2) || ...
          handles.sizePhase(3)~=handles.Seq.AQSlice.sizePhase(3)
        if handles.Seq.AQSlice.nPhase(1)>1
          [X,Y]=meshgrid((1:handles.Seq.AQSlice.nPhase(2))-handles.Seq.AQSlice.nPhase(2)/2,(1:handles.Seq.AQSlice.nPhase(1))-handles.Seq.AQSlice.nPhase(1)/2);
          handles.AqWin1_image=peaks(X/handles.Seq.AQSlice.nPhase(2)*6*get(handles.slider12,'Value'),Y/handles.Seq.AQSlice.nPhase(1)*6*get(handles.slider11,'Value'));
          % handles.AqWin1_image=zeros(handles.Seq.AQSlice.nPhase(1),handles.Seq.AQSlice.nPhase(2));
        else
          [X,Y]=meshgrid((1:handles.Seq.AQSlice.nPhase(2))-handles.Seq.AQSlice.nPhase(2)/2,(1:handles.Seq.AQSlice.nRead)-handles.Seq.AQSlice.nRead/2);
          handles.AqWin1_image=peaks(X/handles.Seq.AQSlice.nPhase(2)*6*get(handles.slider12,'Value'),Y/handles.Seq.AQSlice.nRead*6*get(handles.slider11,'Value'));
          % handles.AqWin1_image=zeros(handles.Seq.AQSlice.nRead,handles.Seq.AQSlice.nPhase(2));
        end
        % handles.AqWin1_image=zeros(size(X));

        handles.sizeRead=handles.Seq.AQSlice.sizeRead;
        handles.sizePhase(1)=handles.Seq.AQSlice.sizePhase(1);
        handles.sizePhase(2)=handles.Seq.AQSlice.sizePhase(2);
        handles.sizePhase(3)=handles.Seq.AQSlice.sizePhase(3);
        handles.nRead=handles.Seq.AQSlice.nRead;
        handles.nPhase(1)=handles.Seq.AQSlice.nPhase(1);
        handles.nPhase(2)=handles.Seq.AQSlice.nPhase(2);
        handles.nPhase(3)=handles.Seq.AQSlice.nPhase(3);
      end
      % cla(handles.axes1);
      % get(handles.axes1);


      handles.Seq.AQSlice.alfa=handles.alfa;
      handles.Seq.AQSlice.phi=handles.phi;
      handles.Seq.AQSlice.theta=handles.theta;

      if handles.ButtonZeroapt == 1
        handles.ButtonZeroapt = 0;
        set(handles.slider_alfa, 'Value', 0);
        set(handles.slider_phi, 'Value', 0);
        set(handles.slider_theta, 'Value', 0);
      end

      if handles.ButtonCenterXYZ == 1
        handles.ButtonCenterXYZ=0;
        set(handles.slider_y, 'Value', 0.5);
        set(handles.slider_z, 'Value', 0.5);
        set(handles.slider_x, 'Value', 0.5);
        set(handles.slider_phase2, 'Value', 0.5);
        set(handles.slider_phase3_read, 'Value', 0.5);
        set(handles.slider_phase1_slice, 'Value', 0.5);
      end

      handles.alfa = get(handles.slider_alfa, 'Value') * 2*pi;
      handles.phi = get(handles.slider_phi, 'Value') * 2*pi;
      handles.theta = get(handles.slider_theta, 'Value') * 2*pi;


      if handles.alfa~=handles.Seq.AQSlice.alfa, handles.ButtonCenter=1; end
      if handles.phi~=handles.Seq.AQSlice.phi, handles.ButtonCenter=1; end
      if handles.theta~=handles.Seq.AQSlice.theta, handles.ButtonCenter=1; end

      handles.Seq.AQSlice.alfa=handles.alfa;
      handles.Seq.AQSlice.phi=handles.phi;
      handles.Seq.AQSlice.theta=handles.theta;



      handles.CenterOffset =[ get(handles.slider_x,'Value')*(handles.LocaliserImageVol(2)-handles.LocaliserImageVol(1))+handles.LocaliserImageVol(1), ...
                              get(handles.slider_y,'Value')*(handles.LocaliserImageVol(4)-handles.LocaliserImageVol(3))+handles.LocaliserImageVol(3), ...
                              get(handles.slider_z,'Value')*(handles.LocaliserImageVol(6)-handles.LocaliserImageVol(5))+handles.LocaliserImageVol(5)];

      CenterOfImageOffset = [ get(handles.slider_phase1_slice,'Value')*(handles.LocaliserImageVol(2)-handles.LocaliserImageVol(1))-(handles.LocaliserImageVol(2)-handles.LocaliserImageVol(1))/2, ...
                              get(handles.slider_phase2,'Value')*(handles.LocaliserImageVol(4)-handles.LocaliserImageVol(3))-(handles.LocaliserImageVol(4)-handles.LocaliserImageVol(3))/2, ...
                              get(handles.slider_phase3_read,'Value')*(handles.LocaliserImageVol(6)-handles.LocaliserImageVol(5))-(handles.LocaliserImageVol(6)-handles.LocaliserImageVol(5))/2];

      handles.Seq.AQSlice.alfa=handles.alfa;
      handles.Seq.AQSlice.phi=handles.phi;
      handles.Seq.AQSlice.theta=handles.theta;

      % offsetData=nRotate(CenterOfImageOffset(:)-handles.CenterOfImageOffset(:),handles.Seq.AQSlice).';
      offsetData=aptRotate((CenterOfImageOffset(:)-handles.CenterOfImageOffset(:)).',handles.Seq.AQSlice(1).alfa,handles.Seq.AQSlice(1).phi,handles.Seq.AQSlice(1).theta).';
      CenterOffsetAll=handles.CenterOffset+[offsetData(1),offsetData(2),offsetData(3)];
      CenterOffsetSlider = [ (CenterOffsetAll(1)-handles.LocaliserImageVol(1))/(handles.LocaliserImageVol(2)-handles.LocaliserImageVol(1)),...
                             (CenterOffsetAll(2)-handles.LocaliserImageVol(3))/(handles.LocaliserImageVol(4)-handles.LocaliserImageVol(3)),...
                             (CenterOffsetAll(3)-handles.LocaliserImageVol(5))/(handles.LocaliserImageVol(6)-handles.LocaliserImageVol(5))];
      CenterOffsetSlider(CenterOffsetSlider>1)=1;
      CenterOffsetSlider(CenterOffsetSlider<0)=0;
      set(handles.slider_x, 'Value', CenterOffsetSlider(1));
      set(handles.slider_y, 'Value', CenterOffsetSlider(2));
      set(handles.slider_z, 'Value', CenterOffsetSlider(3));
      if handles.ButtonCenter == 1
        set(handles.slider_phase1_slice, 'Value', 0.5);
        set(handles.slider_phase2, 'Value', 0.5);
        set(handles.slider_phase3_read, 'Value', 0.5);
        handles.CenterOfImageOffset = [0,0,0];
        handles.ButtonCenter = 0;
      else
        handles.CenterOfImageOffset = CenterOfImageOffset;
      end
      handles.CenterOffset =[get(handles.slider_x,'Value')*(handles.LocaliserImageVol(2)-handles.LocaliserImageVol(1))+handles.LocaliserImageVol(1), ...
                             get(handles.slider_y,'Value')*(handles.LocaliserImageVol(4)-handles.LocaliserImageVol(3))+handles.LocaliserImageVol(3), ...
                             get(handles.slider_z,'Value')*(handles.LocaliserImageVol(6)-handles.LocaliserImageVol(5))+handles.LocaliserImageVol(5)];



      handles.Seq.AQSlice.Center = handles.CenterOffset;  % attention: changed order
      % [handles.Seq.AQSlice.normV(1),handles.Seq.AQSlice.normV(2),handles.Seq.AQSlice.normV(3)]=sph2cart(handles.Seq.AQSlice.theta,-handles.Seq.AQSlice.phi,1);
      % handles.Seq.AQSlice.R=handles.Seq.AQSlice.Center*handles.Seq.AQSlice.normV.';
      % handles.Seq.AQSlice.Rauf=handles.Seq.AQSlice.R*handles.Seq.AQSlice.normV;
      %
      % handles.Seq.AQSlice.RaufCenter=handles.Seq.AQSlice.Center-handles.Seq.AQSlice.Rauf;
      % handles.Seq.AQSlice.CenterRauf=handles.Seq.AQSlice.Rauf-handles.Seq.AQSlice.Center;
      % handles.Seq.AQSlice.CenterRaufImage=tpaRotate( handles.Seq.AQSlice.CenterRauf',handles.Seq.AQSlice.alfa ,handles.Seq.AQSlice.phi, handles.Seq.AQSlice.theta)';
      % handles.Seq.AQSlice.CenterRaufImage=handles.Seq.AQSlice.CenterRaufImage([3,2,1]).*[1,-1,1];
      % handles.Seq.AQSlice.Center2OriginImage=[handles.Seq.AQSlice.CenterRaufImage([1,2]),handles.Seq.AQSlice.R].*[1,1,-1];

      handles.Seq.AQSlice.normXImage=[1,0,0];
      handles.Seq.AQSlice.normYImage=[0,1,0];
      handles.Seq.AQSlice.normZImage=[0,0,1];
      Angle2Deg=1/(2*pi)*360;
      [Rx, Ry, Rz] = get_aptDegRotationMatrix(handles.Seq.AQSlice.alfa*Angle2Deg, handles.Seq.AQSlice.phi*Angle2Deg, handles.Seq.AQSlice.theta*Angle2Deg);
      handles.Seq.AQSlice.normX = (Rz*(Ry*(Rx*handles.Seq.AQSlice.normXImage.'))).';
      handles.Seq.AQSlice.normY = (Rz*(Ry*(Rx*handles.Seq.AQSlice.normYImage.'))).';
      handles.Seq.AQSlice.normZ = (Rz*(Ry*(Rx*handles.Seq.AQSlice.normZImage.'))).';

      handles.Seq.AQSlice.XOriginImage=(-handles.Seq.AQSlice.Center)*handles.Seq.AQSlice.normX.';
      handles.Seq.AQSlice.YOriginImage=(-handles.Seq.AQSlice.Center)*handles.Seq.AQSlice.normY.';
      handles.Seq.AQSlice.ZOriginImage=(-handles.Seq.AQSlice.Center)*handles.Seq.AQSlice.normZ.';
      handles.Seq.AQSlice.ReceivingPoint=(handles.Seq.AQSlice.Center*handles.Seq.AQSlice.normX.').*handles.Seq.AQSlice.normX;
      handles.Seq.AQSlice.ReceivingPoint2Center=handles.Seq.AQSlice.Center-handles.Seq.AQSlice.ReceivingPoint;

      % handles.Seq.AQSlice.Center2ReceivingPoint=handles.Seq.AQSlice.ReceivingPoint-handles.Seq.AQSlice.Center;
      handles.Seq.AQSlice.Center2OriginImage = [ handles.Seq.AQSlice.XOriginImage,...
                                                 handles.Seq.AQSlice.YOriginImage,...
                                                 handles.Seq.AQSlice.ZOriginImage];
      % disp(['normX              ' num2str(handles.Seq.AQSlice.normX)])
      % disp(['Center2OriginImage ' num2str(handles.Seq.AQSlice.Center2OriginImage)])
      % disp(['Center             ' num2str(handles.Seq.AQSlice.Center)])
      % disp(['normV = ', num2str(handles.Seq.AQSlice.normV)])
      % disp(['R     = ', num2str(handles.Seq.AQSlice.R)])
      % disp(['Rauf  = ', num2str(handles.Seq.AQSlice.Rauf)])
      % disp(['alfa  = ', num2str(handles.Seq.AQSlice.alfa)])
      % disp(['phi   = ', num2str(handles.Seq.AQSlice.phi)])
      % disp(['theta = ', num2str(handles.Seq.AQSlice.theta)])
      %
      % Seq.AmpSlice=[1;0;0];
      % Seq.AmpRead=[0;1;0];
      % Seq.AmpPhase=[0;0;1];
      %
      % [Seq] = GradRotate(Seq, handles.Seq.AQSlice);
      %
      % disp(['AmpSlice = ', num2str(Seq.AmpSlice.','% 3.3f')])
      % disp(['AmpRead  = ', num2str(Seq.AmpRead.','% 3.3f')])
      % disp(['AmpPhase = ', num2str(Seq.AmpPhase.','% 3.3f')])
      guidata(hObject, handles);
      handles = axes1update(handles);


      if handles.ButtonMeasure3DLocaliser==1
        if xor(handles.Seq.AQSlice.nRead>1,handles.Seq.AQSlice.nPhase(3)>1)

        else
          handles.ButtonMeasure3DLocaliser=0;
          warning('Please set nRead>1 xor nPhase(3)>1 to measure a localiser 3D image')
          beep;
        end
        if ~(handles.Seq.AQSlice.nPhase(1)>1)
          handles.ButtonMeasure3DLocaliser=0;
          warning('Please set nPhase(1)>1 to measure a localiser 3D image')
          beep;
        end
        if ~(handles.Seq.AQSlice.nPhase(2)>1)
          handles.ButtonMeasure3DLocaliser=0;
          warning('Please set nPhase(2)>1 to measure a localiser 3D image')
          beep;
        end
        if ~isinf(handles.Seq.AQSlice.thickness)
          handles.ButtonMeasure3DLocaliser=0;
          warning('Please set thickness to Inf to measure a localiser 3D image')
          beep;
        end
        if (handles.Seq.AQSlice.alfa)
          handles.ButtonMeasure3DLocaliser=0;
          warning('Please set alfa=0 to measure a localiser 3D image')
          beep
        end
        if (handles.Seq.AQSlice.phi)
          handles.ButtonMeasure3DLocaliser=0;
          warning('Please set phi=0 to measure a localiser 3D image')
          beep;
        end
        if (handles.Seq.AQSlice.theta)
          handles.ButtonMeasure3DLocaliser=0;
          warning('Please set theta=0 to measure a localiser 3D image')
          beep;
        end
        % if (sum(abs(handles.Seq.AQSlice.Center)))
        %     handles.ButtonMeasure3DLocaliser=0;
        %     warning('Please set Center of XYZ to [0,0,0} to measure a localiser 3D image')
        %     beep;
        % end

      end
      handles.axes1update=1;
      guidata(hObject, handles);
    end

    if handles.ButtonMeasureSlice || handles.ButtonMeasure3DLocaliser

      handles.ButtonMeasureSlice=0;

      % LoadSystem
      %----------- Parameter ------------------------------------------------
      % clear global SliceSelect
      % mySave=[]
      % handles.HW.Grad.TimeDelay=40e-6;
      % handles.HW.Grad.tRamp=handles.HW.Grad.tRamp/3;
      % handles.HW.Grad.tEC=handles.HW.Grad.tEC/3;

      handles.Seq.Loops=str2double(get(handles.edit14,'String'));
      handles.Seq.average=str2double(get(handles.edit13,'String'));
      handles.Seq.averageBreak=str2double(get(handles.edit28,'String'));
      handles.Seq.LoopsBreak=str2double(get(handles.edit28,'String'));
      set(handles.text35,'String',num2str(handles.Seq.average*handles.Seq.Loops));



      handles.Seq.tEcho=str2double(get(handles.edit17,'String'));
      if isnan(str2double(get(handles.edit19,'String')))
        handles.Seq.AQSlice(1).sizePhaseSpoil=[];         %!!!!
      else
        handles.Seq.AQSlice(1).sizePhaseSpoil=str2double(get(handles.edit19,'String'))/1000;   % mm -> m      %!!!!
      end
      handles.Seq.AQSlice.HzPerPixMin=str2double(get(handles.edit29,'String'));
      handles.Seq.AQSlice.thickness=str2double(get(handles.edit6,'String'))/1000; % mm -> m
      if isnan(handles.Seq.AQSlice.thickness); handles.Seq.AQSlice.thickness=inf;end

      if isnan(str2double(get(handles.edit18,'String')))
        handles.Seq.tRep=[];
        handles.Seq.AQSlice(1).TurboBreak=[];
      else
        if strcmp(get(handles.text18,'String'),'T Rep')
          handles.Seq.RepetitionTime=str2double(get(handles.edit18,'String'));
          handles.Seq.AQSlice(1).TurboBreak=[];
        else
          handles.Seq.AQSlice(1).TurboBreak=str2double(get(handles.edit18,'String'));
          handles.Seq.tRep=[];
        end
      end
      t=get(handles.popupmenu_slicePulse, 'String');
      if strcmp(char(t(get(handles.popupmenu_slicePulse, 'Value'))),'auto')
        handles.Seq.AQSlice(1).excitationPulse=[];
      else
        handles.Seq.AQSlice(1).excitationPulse=str2func(char(t(get(handles.popupmenu_slicePulse, 'Value'))));
      end
      t=get(handles.popupmenu_InversionPulse, 'String');
      if strcmp(char(t(get(handles.popupmenu_InversionPulse, 'Value'))),'auto')
        handles.Seq.AQSlice(1).inversionPulse=[];
      else
        handles.Seq.AQSlice(1).inversionPulse=str2func(char(t(get(handles.popupmenu_InversionPulse, 'Value'))));
      end

      t=get(handles.popROS, 'String');
      tv=get(handles.popROS, 'Value');
      if tv==14
        handles.Seq.AQSlice.ReadOS=[];
      else
        handles.Seq.AQSlice.ReadOS=str2double(char(t(tv)));
      end
      t=get(handles.popPOS, 'String');
      handles.Seq.AQSlice.PhaseOS(1)=str2double(char(t(get(handles.popPOS, 'Value'))));
      t=get(handles.popupmenu_PhaseOS2, 'String');
      handles.Seq.AQSlice.PhaseOS(2)=str2double(char(t(get(handles.popupmenu_PhaseOS2, 'Value'))));
      t=get(handles.popupmenu8, 'String');
      handles.Seq.AQSlice.PhaseOS(3)=str2double(char(t(get(handles.popupmenu8, 'Value'))));

      if isnan(str2double(get(handles.edit20,'String')))
        handles.Seq.T1=[];
        handles.Seq.AQSlice(1).TurboBreak=[];
        if strcmp(get(handles.text20,'String'),'T1')
          set(handles.edit20,'String','0.1');
        else
          set(handles.edit20,'String','max');
        end
      end

      if strcmp(get(handles.text20,'String'),'T1')
        handles.Seq.T1=str2double(get(handles.edit20,'String'));         % T1 der Probe Flipwinkel ist acos(exp(-Seq.tRep/Seq.T1))/pi*180
        handles.Seq.AQSlice(1).TurboFactor=1;
      else
        if isnan(str2double(get(handles.edit20,'String')))
          handles.Seq.AQSlice(1).TurboFactor=1e12;
        else
          handles.Seq.AQSlice(1).TurboFactor=str2double(get(handles.edit20,'String'));
        end
        handles.Seq.T1=[];

        factors=factor(prod(handles.Seq.AQSlice.PhaseOS)*prod(handles.Seq.AQSlice.nPhase));
        fsel=zeros(2^(length(factors)),(length(factors)));
        for t=1:length(factors)
          fsel(2^(t-1)+1:2^(t-1):end,t)=1;
        end
        fsel=mod(cumsum(fsel,1),2);
        factors=repmat(factors,size((fsel),1),1);
        factors(fsel==0)=1;
        factors=unique(prod(factors,2));
        handles.Seq.AQSlice(1).TurboFactor=factors(find(factors<=handles.Seq.AQSlice(1).TurboFactor,1,'last'));
        if ~isnan(str2double(get(handles.edit20,'String')))
          set(handles.edit20,'String',num2str(handles.Seq.AQSlice(1).TurboFactor));
        end
        clear factors fsel
      end

      handles.Seq.AQSlice.plot_k_image=1;
      handles.Seq.AQSlice.plotPhase=1;
      handles.Seq.AQSlice.plotB0ppm=0;


      if get(handles.checkbox_plotTR,'Value')
        handles.Seq.plotSeqTR=1:3;
      else
        handles.Seq.plotSeqTR=[];
      end
      if get(handles.checkbox_PlotSequence,'Value')
        handles.Seq.plotSeq=1:3;
      else
        handles.Seq.plotSeq=[];
      end
      if get(handles.checkbox_Plot_kSpace,'Value')
        handles.Seq.AQSlice(1).plotkSpace=1;
      else
        handles.Seq.AQSlice(1).plotkSpace=0;
      end
      if get(handles.checkbox_image,'Value')
        handles.Seq.AQSlice(1).plotImage=1;
      else
        handles.Seq.AQSlice(1).plotImage=0;
      end

      if get(handles.checkbox_Smooth,'Value')
        handles.Seq.AQSlice(1).ZeroFillWindowSize=1.4;                     % Zero Fill Window Size (k-Space)
        handles.Seq.AQSlice(1).ZeroFillFactor=2;                           % Zero fill resolution factor
      else
        handles.Seq.AQSlice(1).ZeroFillWindowSize=[];                     % Zero Fill Window Size (k-Space)
        handles.Seq.AQSlice(1).ZeroFillFactor=[];                         % Zero fill resolution factor
      end

      if get(handles.checkbox_PlotPhase,'Value')
        handles.Seq.AQSlice(1).plotPhase=1;
      else
        handles.Seq.AQSlice(1).plotPhase=0;
      end
      if get(handles.checkbox_Extra1,'Value')
        handles.Seq.AQSlice(1).plotB0ppm=1;
      else
        handles.Seq.AQSlice(1).plotB0ppm=0;
      end
      if get(handles.checkbox_Extra2,'Value')
        handles.Seq.AQSlice(1).plotB0Hz=1;
      else
        handles.Seq.AQSlice(1).plotB0Hz=0;
      end
      % AQSlice.plotB0ppm
      % AQSlice.plotB0Hz
      % AQSlice.plotImage
      % Seq.plotSeq=1:3;
      % handles.HW.Function_While_Measurement=@Function_While_Measurement_test;   % handles.HW=handles.HW.Function_While_Measurement( handles.HW);
      % Seq.showSliceRead=1;
      % handles.Seq.showSlicePhase=(str2double(get(handles.edit27,'String'))==1); %!!!

      popup_sel_index = get(handles.popupmenu_sequenceSelect, 'Value');

      if popup_sel_index ~= PD.SliceSelectType.SliceSelect
        [handles.HW, handles.mySave] = Find_Frequency_Sweep(handles.HW, handles.mySave, 1);  % Find magnet frequency
      end

      handles.Seq.AQSlice(1).PhaseCoordinate=[1 2 3];
      handles.Seq.AQSlice(1).ReadCoordinate=3;
      handles.Seq.AQSlice(1).SliceCoordinate=1;

      try
        switch popup_sel_index
          case PD.SliceSelectType.SpinEcho  % Spin-Echo
            [SeqLoop] = sequence_Spin_Echo(handles.HW, handles.Seq, handles.AQ, handles.TX, handles.Grad,  handles.mySave);
%             handles.HW = SeqLoop.HW;


          case PD.SliceSelectType.GradEcho % Grad-Echo
            [SeqLoop]= sequence_Flash(handles.HW, handles.Seq, handles.AQ, handles.TX, handles.Grad,  handles.mySave);
%             handles.HW = SeqLoop.HW;

          case PD.SliceSelectType.SliceSelect  % Slice Selection
            % write Seq.AQSlice data to file
           err = saveAQSliceToFile(handles, sprintf('AQSlice_%s', get(handles.edit27, 'String')));
           if err < 0
             errordlg(sprintf('Failed to save slice settings to file:\n%s', lastwarn()), 'Error Saving File', 'modal');
           end
           handles.MotionFcnBusy = 0;
           guidata(hObject, handles);
           return;

        end

        if isfield(handles, 'light_handle') && ishghandle(handles.light_handle, 'light'), delete(handles.light_handle); end
        handles.light_handle = camlight();
        set(handles.light_handle, 'Parent', handles.axes1);


        if get(handles.checkbox_Save,'Value')
          if isempty(get(handles.edit27, 'String')) && ~get(handles.checkbox_addDate, 'Value')
            set(handles.checkbox_addDate, 'Value', 1);
          end
          if get(handles.checkbox_addDate,'Value')
            SeqLoop.saveFilename=[get(handles.edit21,'String') '\' get(handles.edit27,'String') ' ' datestr(now, 'dd.mm.yyyy HH-MM-SS.FFF') '.mat'];
          else
            SeqLoop.saveFilename=[get(handles.edit21,'String') '\' get(handles.edit27,'String') '.mat'];
          end
          if ~isfolder(get(handles.edit21,'String')); mkdir(get(handles.edit21,'String'));end
          save(SeqLoop.saveFilename,'SeqLoop')
        end

        if handles.ButtonMeasure3DLocaliser==1
          handles.ButtonMeasure3DLocaliser=0;

          set(handles.slider11,'Value',1)
          set(handles.slider12,'Value',1)
          set(handles.slider13,'Value',1)
          set(handles.slider14,'Value',1)
          % pt=SeqLoop.AQSlice(1).PhaseCoordinate;
          handles.CenterOffsetLocaliser=handles.CenterOffset;
          handles.ButtonCenterXYZ=1;
          % handles.CenterOffset=[0,0,0];
          handles.LocaliserImageVol=[...
                      -SeqLoop.AQSlice(1).sizePhase(1)/2+handles.CenterOffsetLocaliser(1),SeqLoop.AQSlice(1).sizePhase(1)/2+handles.CenterOffsetLocaliser(1),...
                      -SeqLoop.AQSlice(1).sizePhase(2)/2+handles.CenterOffsetLocaliser(2),SeqLoop.AQSlice(1).sizePhase(2)/2+handles.CenterOffsetLocaliser(2),...
                      -SeqLoop.AQSlice(1).sizePhase(3)/2+handles.CenterOffsetLocaliser(3),SeqLoop.AQSlice(1).sizePhase(3)/2+handles.CenterOffsetLocaliser(3)];

          if SeqLoop.AQSlice(1).nRead>1
            % pt(1)=SeqLoop.AQSlice(1).ReadCoordinate;
            handles.LocaliserImageVol(5:6)=[-SeqLoop.AQSlice(1).sizeRead/2+handles.CenterOffsetLocaliser(3),SeqLoop.AQSlice(1).sizeRead/2+handles.CenterOffsetLocaliser(3)];
          end
          % [~,pts]=sort(pt);
          handles.v=abs(permute(squeeze(SeqLoop.data.Image),[3 2 1]));
          handles.nx=size(handles.v,2);
          handles.ny=size(handles.v,1);
          handles.nz=size(handles.v,3);


          [handles.x,handles.y,handles.z] = meshgrid( linspace(handles.LocaliserImageVol(1),handles.LocaliserImageVol(2),handles.nx),...
                                                      linspace(handles.LocaliserImageVol(3),handles.LocaliserImageVol(4),handles.ny),...
                                                      linspace(handles.LocaliserImageVol(5),handles.LocaliserImageVol(6),handles.nz));
          handles.lim=[min(handles.v(:)) max(handles.v(:))];
          axis(handles.axes1,handles.LocaliserImageVol([5,6,1,2,3,4]))
%           handles.axes1update=1;
          handles.refresh3D=1;
          guidata(hObject, handles);
          handles = axes1update(handles);
        end


      catch ME
        warning(ME.message);
        warning(getReport(ME));
        beep();
      end
      if isfield(handles,'mySave')
        HW=handles.HW;
        mySave=handles.mySave;
        Seq=handles.Seq;
        AQ=handles.AQ;
        TX=handles.TX;
        Grad=handles.Grad;
        %talker=handles.talker;
      end
      if exist('HW', 'var') && isa(HW, 'PD.HWClass')
        if evalin('base', 'exist(''mySave'', ''var'')')
          mySave = evalin('base', 'mySave');
        end
        ResetStructs
      else
        LoadSystem
      end
      handles.HW=HW;
      handles.mySave=mySave;
      handles.Seq=Seq;
      handles.AQ=AQ;
      handles.TX=TX;
      handles.Grad=Grad;
      %handles.talker=talker;

    end
  catch ME
    handles.MotionFcnBusy = 0;
    handles.axes1update=1;
    guidata(hObject, handles);
    rethrow(ME);
  end

  handles.MotionFcnBusy = 0;
  guidata(hObject, handles);
end


end


function handles = axes1update(handles)

%  handles.Seq.AQSlice.Center
%  handles.Seq.AQSlice.Center2OriginImage
%              disp('--------------------------------------------------------------------------')
hold(handles.axes1, 'on');
% handles.Seq.AQSlice.CenterRaufImage
if isfield(handles, 'h_Center') && ishandle(handles.h_Center)
  delete(handles.h_Center);
end
handles.h_Center = plot3(handles.axes1,[0;handles.Seq.AQSlice.Center(3)],[0;handles.Seq.AQSlice.Center(1)],[0;handles.Seq.AQSlice.Center(2)],'k','LineWidth',4);
% if isfield(handles, 'h_Center2OriginImageY') && ishandle(handles.h_Center2OriginImageY)
%   delete(handles.h_Center2OriginImageY);
% end
% handles.h_Center2OriginImageY=plot3(handles.axes1,[0;-handles.Seq.AQSlice.Center2OriginImage(3)]-handles.Seq.AQSlice.Center2OriginImage(3),[0;-handles.Seq.AQSlice.Center2OriginImage(1)],[0;-handles.Seq.AQSlice.Center2OriginImage(2)],'--k','LineWidth',4);
if isfield(handles, 'h_ReceivingPoint') && ishandle(handles.h_ReceivingPoint)
  delete(handles.h_ReceivingPoint);
end
handles.h_ReceivingPoint = plot3(handles.axes1,[0;handles.Seq.AQSlice.ReceivingPoint(3)],[0;handles.Seq.AQSlice.ReceivingPoint(1)],[0;handles.Seq.AQSlice.ReceivingPoint(2)],'-r.','LineWidth',1);

if isfield(handles, 'h_ReceivingPoint2Center') && ishandle(handles.h_ReceivingPoint2Center)
  delete(handles.h_ReceivingPoint2Center);
end
handles.h_ReceivingPoint2Center = plot3(handles.axes1,[0;handles.Seq.AQSlice.ReceivingPoint2Center(3)]+handles.Seq.AQSlice.ReceivingPoint(3),[0;handles.Seq.AQSlice.ReceivingPoint2Center(1)]+handles.Seq.AQSlice.ReceivingPoint(1),[0;handles.Seq.AQSlice.ReceivingPoint2Center(2)]+handles.Seq.AQSlice.ReceivingPoint(2),':k','LineWidth',4);

tempLength = sqrt(sum((handles.LocaliserImageVol(2:2:6)-handles.LocaliserImageVol(1:2:5)).^2))/10;
if isfield(handles, 'h_normX') && ishandle(handles.h_normX)
  delete(handles.h_normX);
end
handles.h_normX=plot3(handles.axes1,  [0;handles.Seq.AQSlice.normX(3)*tempLength]+handles.Seq.AQSlice.Center(3),...
                                      [0;handles.Seq.AQSlice.normX(1)*tempLength]+handles.Seq.AQSlice.Center(1),...
                                      [0;handles.Seq.AQSlice.normX(2)*tempLength]+handles.Seq.AQSlice.Center(2),'r','LineWidth',4);
if isfield(handles, 'h_normY') && ishandle(handles.h_normY)
  delete(handles.h_normY);
end
handles.h_normY=plot3(handles.axes1,  [0;handles.Seq.AQSlice.normY(3)*tempLength]+handles.Seq.AQSlice.Center(3),...
                                      [0;handles.Seq.AQSlice.normY(1)*tempLength]+handles.Seq.AQSlice.Center(1),...
                                      [0;handles.Seq.AQSlice.normY(2)*tempLength]+handles.Seq.AQSlice.Center(2),'g','LineWidth',4);
if isfield(handles, 'h_normZ') && ishandle(handles.h_normZ)
  delete(handles.h_normZ);
end
handles.h_normZ=plot3(handles.axes1,  [0;handles.Seq.AQSlice.normZ(3)*tempLength]+handles.Seq.AQSlice.Center(3),...
                                      [0;handles.Seq.AQSlice.normZ(1)*tempLength]+handles.Seq.AQSlice.Center(1),...
                                      [0;handles.Seq.AQSlice.normZ(2)*tempLength]+handles.Seq.AQSlice.Center(2),'b','LineWidth',4);

if isfield(handles, 'h_normXImage') && ishandle(handles.h_normXImage)
  delete(handles.h_normXImage);
end
handles.h_normXImage=plot3(handles.axes1,  [0;handles.Seq.AQSlice.normXImage(3)*tempLength],...
                                      [0;handles.Seq.AQSlice.normXImage(1)*tempLength],...
                                      [0;handles.Seq.AQSlice.normXImage(2)*tempLength],':r','LineWidth',2);
if isfield(handles, 'h_normYImage') && ishandle(handles.h_normYImage)
  delete(handles.h_normYImage);
end
handles.h_normYImage=plot3(handles.axes1,  [0;handles.Seq.AQSlice.normYImage(3)*tempLength],...
                                      [0;handles.Seq.AQSlice.normYImage(1)*tempLength],...
                                      [0;handles.Seq.AQSlice.normYImage(2)*tempLength],':g','LineWidth',2);
if isfield(handles, 'h_normZImage') && ishandle(handles.h_normZImage)
  delete(handles.h_normZImage);
end
handles.h_normZImage=plot3(handles.axes1,  [0;handles.Seq.AQSlice.normZImage(3)*tempLength],...
                                      [0;handles.Seq.AQSlice.normZImage(1)*tempLength],...
                                      [0;handles.Seq.AQSlice.normZImage(2)*tempLength],':b','LineWidth',2);


if isfield(handles, 'h_imageBox') && ishandle(handles.h_imageBox)
  delete(handles.h_imageBox);
end

if handles.Seq.AQSlice.nPhase(1)==1
  box.x=handles.Seq.AQSlice.thickness;
else
  box.x=handles.Seq.AQSlice.sizePhase(1);
end

box.y=handles.Seq.AQSlice.sizePhase(2);
if handles.Seq.AQSlice.nPhase(3)~=1
  box.z=handles.Seq.AQSlice.sizePhase(3);
else
  box.z=handles.Seq.AQSlice.sizeRead;
end
if ~isempty(box.x)
  box.line=[...
     -handles.Seq.AQSlice.normX*box.x/2-handles.Seq.AQSlice.normY*box.y/2-handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    1 bottom
     +handles.Seq.AQSlice.normX*box.x/2-handles.Seq.AQSlice.normY*box.y/2-handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    2
     +handles.Seq.AQSlice.normX*box.x/2+handles.Seq.AQSlice.normY*box.y/2-handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    3
     -handles.Seq.AQSlice.normX*box.x/2+handles.Seq.AQSlice.normY*box.y/2-handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    4
     -handles.Seq.AQSlice.normX*box.x/2-handles.Seq.AQSlice.normY*box.y/2-handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    1
     -handles.Seq.AQSlice.normX*box.x/2-handles.Seq.AQSlice.normY*box.y/2+handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    5 top
     +handles.Seq.AQSlice.normX*box.x/2-handles.Seq.AQSlice.normY*box.y/2+handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    6
     +handles.Seq.AQSlice.normX*box.x/2-handles.Seq.AQSlice.normY*box.y/2-handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    2
     +handles.Seq.AQSlice.normX*box.x/2-handles.Seq.AQSlice.normY*box.y/2+handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    6
     +handles.Seq.AQSlice.normX*box.x/2+handles.Seq.AQSlice.normY*box.y/2+handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    7
     +handles.Seq.AQSlice.normX*box.x/2+handles.Seq.AQSlice.normY*box.y/2-handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    3
     +handles.Seq.AQSlice.normX*box.x/2+handles.Seq.AQSlice.normY*box.y/2+handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    7
     -handles.Seq.AQSlice.normX*box.x/2+handles.Seq.AQSlice.normY*box.y/2+handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    7
     -handles.Seq.AQSlice.normX*box.x/2+handles.Seq.AQSlice.normY*box.y/2-handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    4
     -handles.Seq.AQSlice.normX*box.x/2+handles.Seq.AQSlice.normY*box.y/2+handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center;...    7
     -handles.Seq.AQSlice.normX*box.x/2-handles.Seq.AQSlice.normY*box.y/2+handles.Seq.AQSlice.normZ*box.z/2+handles.Seq.AQSlice.Center       %5 top
     ];

     handles.h_imageBox=plot3(handles.axes1, box.line(:,3),...
                                             box.line(:,1),...
                                             box.line(:,2),'m','LineWidth',1);
end

% plot3([0;handles.Seq.AQSlice.CenterRaufImage(1)],[0;handles.Seq.AQSlice.CenterRaufImage(2)],[0;handles.Seq.AQSlice.CenterRaufImage(3)],'LineWidth',4,'Color','c')
% plot3([0;handles.Seq.AQSlice.RaufCenter(1)]+handles.Seq.AQSlice.Rauf(1),[0;handles.Seq.AQSlice.RaufCenter(2)]+handles.Seq.AQSlice.Rauf(2),[0;handles.Seq.AQSlice.RaufCenter(3)]+handles.Seq.AQSlice.Rauf(3),'LineWidth',4)

% handles.Seq.AQSlice.normV
% handles.Seq.AQSlice.Rauf
% handles.Seq.AQSlice.Center
% handles.Seq.AQSlice.RaufCenter
% handles.Seq.AQSlice.CenterRauf


% handles.Seq.AQSlice.normV*handles.Seq.AQSlice.CenterRauf'

% Y=linspace(-0.5,0.5,handles.Seq.AQSlice.nRead)*handles.Seq.AQSlice.sizeRead+handles.CenterOffset(1);
% X=linspace(-0.5,0.5,handles.Seq.AQSlice.nPhase)*handles.Seq.AQSlice.sizePhase+handles.CenterOffset(2);
if handles.Seq.AQSlice.nPhase(3)>1
  x=handles.CenterOffset(1);
  % y=linspace(-0.5,0.5,handles.Seq.AQSlice.nPhase(1))*handles.Seq.AQSlice.sizePhase(1)+handles.CenterOffset(2);
  % z=linspace(-0.5,0.5,handles.Seq.AQSlice.nPhase(2))*handles.Seq.AQSlice.sizePhase(2)+handles.CenterOffset(3);
  if mod(handles.Seq.AQSlice.nPhase(3),2)
    z=linspace(-0.5+0.5/handles.Seq.AQSlice.nPhase(3),0.5-0.5/handles.Seq.AQSlice.nPhase(3),handles.Seq.AQSlice.nPhase(3))*handles.Seq.AQSlice.sizePhase(3)+handles.CenterOffset(3);
  else
    z=linspace(-0.5,0.5-1/handles.Seq.AQSlice.nPhase(3),handles.Seq.AQSlice.nPhase(3))*handles.Seq.AQSlice.sizePhase(3)+handles.CenterOffset(3);
  end

  if mod(handles.Seq.AQSlice.nPhase(2),2)
    y=linspace(-0.5+0.5/handles.Seq.AQSlice.nPhase(2),0.5-0.5/handles.Seq.AQSlice.nPhase(2),handles.Seq.AQSlice.nPhase(2))*handles.Seq.AQSlice.sizePhase(2)+handles.CenterOffset(2);
  else
    y=linspace(-0.5,0.5-1/handles.Seq.AQSlice.nPhase(2),handles.Seq.AQSlice.nPhase(2))*handles.Seq.AQSlice.sizePhase(2)+handles.CenterOffset(2);
  end

  % z=zeros(handles.Seq.AQSlice.nPhase(1),handles.Seq.AQSlice.nPhase(2))+handles.CenterOffset(3);
else
  x=handles.CenterOffset(1);
  if mod(handles.Seq.AQSlice.nRead,2)
    z=linspace(-0.5+0.5/handles.Seq.AQSlice.nRead,0.5-0.5/handles.Seq.AQSlice.nRead,handles.Seq.AQSlice.nRead)*handles.Seq.AQSlice.sizeRead+handles.CenterOffset(3);
  else
    z=linspace(-0.5,0.5-1/handles.Seq.AQSlice.nRead,handles.Seq.AQSlice.nRead)*handles.Seq.AQSlice.sizeRead+handles.CenterOffset(3);
  end
  if mod(handles.Seq.AQSlice.nPhase(2),2)
    y=linspace(-0.5+0.5/handles.Seq.AQSlice.nPhase(2),0.5-0.5/handles.Seq.AQSlice.nPhase(2),handles.Seq.AQSlice.nPhase(2))*handles.Seq.AQSlice.sizePhase(2)+handles.CenterOffset(2);
  else
    y=linspace(-0.5,0.5-1/handles.Seq.AQSlice.nPhase(2),handles.Seq.AQSlice.nPhase(2))*handles.Seq.AQSlice.sizePhase(2)+handles.CenterOffset(2);
  end

  % z=linspace(-0.5,0.5,handles.Seq.AQSlice.nPhase(2))*handles.Seq.AQSlice.sizePhase(2)+handles.CenterOffset(3);
  % z=zeros(handles.Seq.AQSlice.nRead,handles.Seq.AQSlice.nPhase(2))+handles.CenterOffset(3);
end

if length(y) == 1
  y = [y-(abs(z(end)-z(1)))/10, y+(abs(z(end)-z(1)))/10];
end
if length(z) == 1
  z = [z-(abs(y(end)-y(1)))/10, z+(abs(y(end)-y(1)))/10];
end

[X,Y,Z] = meshgrid(x,y,z);
% h_AqWin1=pcolor(X,Y,abs(AqWin1_image));

% h_AqWin1=surf(handles.axes1,X,Y,Z,abs(handles.AqWin1_image));
h_AqWin1=surf(handles.axes1,squeeze(X),squeeze(Y),squeeze(Z),abs(handles.AqWin1_image),'Clipping','off');
% shading(h_AqWin1 flat
% colorbar
% view(3)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% colormap(gray)
% shading flat
% % shading interp
% axis equal
% axis(handles.LocaliserImageVol)



% set(h_AqWin1,'XData',get(h_AqWin1,'XData')+0.005);

rx=0.0*pi;
ry=0.0*pi;
rz=0.0*pi;


if rx~=0
  rotate(h_AqWin1,[1,0,0],rx*360/2/pi,[0,0,0]);
end
if ry~=0
  rotate(h_AqWin1,[0,1,0],ry*360/2/pi,[0,0,0]);
end
if rz~=0
  rotate(h_AqWin1,[0,0,1],rz*360/2/pi,[0,0,0]);
end
% if 1
%   rotate(h_AqWin1,[1,0,0],0.5*pi*360/2/pi,[0,0,0]);
% end
% handles.Seq.AQSlice.normX=tpaRotate( handles.Seq.AQSlice.normXImage',handles.Seq.AQSlice.alfa ,handles.Seq.AQSlice.phi, handles.Seq.AQSlice.theta)';
% handles.Seq.AQSlice.normY=tpaRotate( handles.Seq.AQSlice.normYImage',handles.Seq.AQSlice.alfa ,handles.Seq.AQSlice.phi, handles.Seq.AQSlice.theta)';
% handles.Seq.AQSlice.normZ=tpaRotate( handles.Seq.AQSlice.normZImage',handles.Seq.AQSlice.alfa ,handles.Seq.AQSlice.phi, handles.Seq.AQSlice.theta)';


if handles.Seq.AQSlice.alfa~=0
  rotate(h_AqWin1,handles.Seq.AQSlice.normXImage,handles.Seq.AQSlice.alfa*360/2/pi,[handles.CenterOffset(1),handles.CenterOffset(2),handles.CenterOffset(3)]);
end
if handles.Seq.AQSlice.phi~=0
  % normYAlfa=aptRotate( handles.Seq.AQSlice.normYImage,handles.Seq.AQSlice.alfa ,0, 0);
  % rotate(h_AqWin1,normYAlfa,handles.Seq.AQSlice.phi*360/2/pi,[handles.CenterOffset(1),handles.CenterOffset(2),handles.CenterOffset(3)]);
  rotate(h_AqWin1,handles.Seq.AQSlice.normYImage,handles.Seq.AQSlice.phi*360/2/pi,[handles.CenterOffset(1),handles.CenterOffset(2),handles.CenterOffset(3)]);
end
if handles.Seq.AQSlice.theta~=0
  % normZAlfaPhi=aptRotate( handles.Seq.AQSlice.normZImage,handles.Seq.AQSlice.alfa ,handles.Seq.AQSlice.phi, 0);
  % rotate(h_AqWin1,normZAlfaPhi,handles.Seq.AQSlice.theta*360/2/pi,[handles.CenterOffset(1),handles.CenterOffset(2),handles.CenterOffset(3)]);
  rotate(h_AqWin1,handles.Seq.AQSlice.normZImage,handles.Seq.AQSlice.theta*360/2/pi,[handles.CenterOffset(1),handles.CenterOffset(2),handles.CenterOffset(3)]);
end
% offsetData=nRotate([-CenterOfImageOffset(2),-CenterOfImageOffset(1),CenterOfImageOffset(3)]',handles.Seq.AQSlice);
%
% set(h_AqWin1,'XData',get(h_AqWin1,'XData')+offsetData(1));
% set(h_AqWin1,'YData',get(h_AqWin1,'YData')+offsetData(2));
% set(h_AqWin1,'ZData',get(h_AqWin1,'ZData')+offsetData(3));


if handles.refresh3D
  if  min(size(handles.v))>=64
    handles.vs=smooth3(abs(handles.v),'gaussian',5);
  elseif min(size(handles.v))>=32
    handles.vs=smooth3(abs(handles.v),'gaussian',3);
  else
    handles.vs=abs(handles.v);
  end
  handles.vsp=permute(handles.vs,handles.Vpermut);
  handles.xyz(1).p=permute(handles.x,handles.Vpermut);
  handles.xyz(2).p=permute(handles.y,handles.Vpermut);
  handles.xyz(3).p=permute(handles.z,handles.Vpermut);

  handles.xp=handles.xyz(handles.VpermutSort(1)).p;
  handles.yp=handles.xyz(handles.VpermutSort(2)).p;
  handles.zp=handles.xyz(handles.VpermutSort(3)).p;

end

handles.xyz(1).s=get(h_AqWin1,'XData');
handles.xyz(2).s=get(h_AqWin1,'YData');
handles.xyz(3).s=get(h_AqWin1,'ZData');

delete(h_AqWin1)

handles.xsp=squeeze( handles.xyz(handles.VpermutSort(1)).s);
handles.ysp=squeeze( handles.xyz(handles.VpermutSort(2)).s);
handles.zsp=squeeze( handles.xyz(handles.VpermutSort(3)).s);

if isfield(handles, 'h_slice') && ishandle(handles.h_slice)
  delete(handles.h_slice);
end
handles.h_slice=slice(handles.axes1,...
    handles.xp,...
    handles.yp,...
    handles.zp,...
    abs(handles.vsp),...
    handles.xsp,...
    handles.ysp,...
    handles.zsp);
set(handles.h_slice,'Clipping','off');
set(handles.h_slice,'FaceColor','flat');
set(handles.h_slice,'EdgeColor',[.5,.5,.5]);
set(handles.h_slice,'EdgeAlpha',0.3);
set(handles.h_slice,'FaceAlpha',1);
set(handles.h_slice,'FaceLighting','none');

% set(handles.figureSliceSelect,'WVisual','RGB 16 bits(05 06 05 00) zdepth 24, Hardware Accelerated, OpenGL, Double Buffered, Window')
% set(handles.figureSliceSelect,'WVisualMode','manual')
if handles.refresh3D
  if isfield(handles,'p') && ishandle(handles.p)
    delete(handles.p);
  end
  if isfield(handles,'pc') && ishandle(handles.pc)
    delete(handles.pc);
  end
  isosurfval = 0.5*max(reshape(convn(abs(handles.vsp),ones(ceil(size(handles.vsp)./4))./numel(ones(ceil(size(handles.vsp)./4))),'same'),numel(handles.vsp),1));
  handles.p = patch(isosurface(...
          handles.xp,...
          handles.yp,...
          handles.zp,...
          abs(handles.vsp),...
          isosurfval),'Parent',handles.axes1);
  handles.pc = patch(isocaps(...
          handles.xp,...
          handles.yp,...
          handles.zp,...
          abs(handles.vsp),...
          isosurfval),'Parent',handles.axes1);
  isonormals(...
          handles.xp,...
          handles.yp,...
          handles.zp,...
          abs(handles.vsp),...
          handles.p);
  set(handles.p,'FaceColor','red','EdgeColor','none','FaceLighting','gouraud');
  set(handles.p,'FaceAlpha',0.3);
  set(handles.pc,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
  set(handles.pc,'FaceAlpha',0.3);
  if isfield(handles, 'light_handle') && ishandle(handles.light_handle)
    delete(handles.light_handle);
  end
  handles.light_handle=camlight;
  set(handles.light_handle,'Parent',handles.axes1)
  handles.refresh3D=0;
end
if isfield(handles, 'ht1') && ishandle(handles.ht1)
  delete(handles.ht1);
end
if isfield(handles, 'ht2') && ishandle(handles.ht2)
  delete(handles.ht2);
end
if isfield(handles, 'ht3') && ishandle(handles.ht3)
  delete(handles.ht3);    
end
if isfield(handles, 'ht4') && ishandle(handles.ht4)
  delete(handles.ht4);
end
handles.ht1=text(handles.xsp(round(end/2),1),handles.ysp(round(end/2),1),handles.zsp(round(end/2),1),'  phase(2) left','Color',[0.5,0.5,0.5]);
handles.ht2=text(handles.xsp(round(end/2),end),handles.ysp(round(end/2),end),handles.zsp(round(end/2),end),'  phase(2) right','Color',[0.5,0.5,0.5]);
if handles.Seq.AQSlice.nPhase(3)>1
  handles.ht3=text(handles.xsp(1,round(end/2)),handles.ysp(1,round(end/2)),handles.zsp(1,round(end/2)),'  phase(3) bottom','Color',[0.5,0.5,0.5]);
  handles.ht4=text(handles.xsp(end,round(end/2)),handles.ysp(end,round(end/2)),handles.zsp(end,round(end/2)),'  phase(3) top','Color',[0.5,0.5,0.5]);
else
  handles.ht3=text(handles.xsp(1,round(end/2)),handles.ysp(1,round(end/2)),handles.zsp(1,round(end/2)),'  read bottom','Color',[0.5,0.5,0.5]);
  handles.ht4=text(handles.xsp(end,round(end/2)),handles.ysp(end,round(end/2)),handles.zsp(end,round(end/2)),'  read top','Color',[0.5,0.5,0.5]);
end
set(handles.ht1,'Parent',handles.axes1)
set(handles.ht2,'Parent',handles.axes1)
set(handles.ht3,'Parent',handles.axes1)
set(handles.ht4,'Parent',handles.axes1)
% handles.Seq.AQSlice.Vphase=([xd(round(end/2),1),yd(round(end/2),1),zd(round(end/2),1)]-[xd(round(end/2),end),yd(round(end/2),end),zd(round(end/2),end)])/sum(([xd(round(end/2),1),yd(round(end/2),1),zd(round(end/2),1)]-[xd(round(end/2),end),yd(round(end/2),end),zd(round(end/2),end)]).^2).^0.5;

% disp(['normV  = ', num2str(handles.Seq.AQSlice.normV,'% 3.3f')])
% disp(['Vread  = ', num2str(handles.Seq.AQSlice.Vread,'% 3.3f')])
% disp(['Vphase = ', num2str(handles.Seq.AQSlice.Vphase,'% 3.3f')])
% 
% disp(['AmpSlice dif = ', num2str(Seq.AmpSlice.'-handles.Seq.AQSlice.normV,'% 3.3f')])
% disp(['AmpRead dif = ', num2str(Seq.AmpRead.'-handles.Seq.AQSlice.Vread,'% 3.3f')])
% disp(['AmpPhase dif = ', num2str(Seq.AmpPhase.'-handles.Seq.AQSlice.Vphase,'% 3.3f')])




% drawnow


hold(handles.axes1, 'off')
% global SliceSelect
% SliceSelect = handles.Seq.AQSlice;
% guidata(hObject,handles)

end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figureSliceSelect_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figureSliceSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ButtonDown=1;
handles.ButtonDown
guidata(hObject,handles)

end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figureSliceSelect_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figureSliceSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ButtonDown=0;
handles.ButtonDown
guidata(hObject,handles)

end

% --- Executes on slider movement.
function slider_phi_Callback(hObject, eventdata, handles)
% hObject    handle to slider_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function slider_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes on slider movement.
function slider_theta_Callback(hObject, eventdata, handles)
% hObject    handle to slider_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function slider_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes on slider movement.
function slider_x_Callback(hObject, eventdata, handles)
% hObject    handle to slider_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function slider_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes on slider movement.
function slider_y_Callback(hObject, eventdata, handles)
% hObject    handle to slider_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function slider_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes on slider movement.
function slider_z_Callback(hObject, eventdata, handles)
% hObject    handle to slider_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function slider_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes on slider movement.
function slider_phase1_slice_Callback(hObject, eventdata, handles)
% hObject    handle to slider_phase1_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.ButtonCenter =1;
guidata(hObject,handles)
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function slider_phase1_slice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_phase1_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes on slider movement.
function slider_phase2_Callback(hObject, eventdata, handles)
% hObject    handle to slider_phase2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.ButtonCenter =1;
guidata(hObject,handles)
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function slider_phase2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_phase2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes on slider movement.
function slider_phase3_read_Callback(hObject, eventdata, handles)
% hObject    handle to slider_phase3_read (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.ButtonCenter =1;
guidata(hObject,handles)
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function slider_phase3_read_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_phase3_read (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes on mouse press over figure background.
function figureSliceSelect_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figureSliceSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end


% --- Executes on button press in pushbutton_center.
function pushbutton_center_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    handles.ButtonCenterXYZ=1;
    guidata(hObject,handles)
    figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

end


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

end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

end


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

end


% --- Executes on slider movement.
function slider11_Callback(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function slider11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes on slider movement.
function slider12_Callback(hObject, eventdata, handles)
% hObject    handle to slider12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function slider12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

end


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

end


function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double

end


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double

end


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double

end


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double
popup_sel_index = get(handles.popupmenu_sequenceSelect, 'Value');
switch popup_sel_index
  case PD.SliceSelectType.SpinEcho  % Spin-Echo
    if str2double(get(handles.edit20,'String'))~=1
      set(handles.text18,'String','Turbo Break');
    else
      set(handles.text18,'String','T Rep');
    end
  case PD.SliceSelectType.GradEcho % Grad-Echo
    % set(handles.text20,'String','T1');
    % set(handles.text42,'String','s');
    % set(handles.edit20,'String','0.1');
end
guidata(hObject, handles);

end


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double

end


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit278 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double

end


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double

end


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on selection change in popROS.
function popROS_Callback(hObject, eventdata, handles)
% hObject    handle to popROS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popROS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popROS

end


% --- Executes during object creation, after setting all properties.
function popROS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popROS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on selection change in popPOS.
function popPOS_Callback(hObject, eventdata, handles)
% hObject    handle to popPOS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popPOS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popPOS

end


% --- Executes during object creation, after setting all properties.
function popPOS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popPOS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on selection change in popupmenu_slicePulse.
function popupmenu_slicePulse_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_slicePulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_slicePulse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_slicePulse

end


% --- Executes during object creation, after setting all properties.
function popupmenu_slicePulse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_slicePulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

pulsesmpd=dir('StandardFunctions/pdFunctions/pulse_*.m');
if isdir('Functions')
    pulsesmf=dir('Functions/pulse_*.m');
    pulsesm={'auto',pulsesmpd.name,pulsesmf.name};
else
    pulsesm={'auto',pulsesmpd.name};
end

pulses=strrep(pulsesm,'.m','');
set(hObject, 'String', pulses);

end


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double

end


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double

end


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on selection change in popupmenu_InversionPulse.
function popupmenu_InversionPulse_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_InversionPulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_InversionPulse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_InversionPulse

end


% --- Executes during object creation, after setting all properties.
function popupmenu_InversionPulse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_InversionPulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
pulsesmpd=dir('StandardFunctions/pdFunctions/pulse_*.m');
if isdir('Functions')
    pulsesmf=dir('Functions/pulse_*.m');
    pulsesm={'auto',pulsesmpd.name,pulsesmf.name};
else
    pulsesm={'auto',pulsesmpd.name};
end

pulses=strrep(pulsesm,'.m','');
set(hObject, 'String', pulses);

end


% --- Executes on button press in pushbutton_startMeasurement.
function pushbutton_startMeasurement_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_startMeasurement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.ButtonMeasureSlice=1;
    guidata(hObject,handles)
    figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes on key press with focus on pushbutton_startMeasurement and none of its controls.
function pushbutton_startMeasurement_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_startMeasurement (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%  Key: name of the key that was pressed, in lower case
%  Character: character interpretation of the key(s) that was pressed
%  Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton_startMeasurement.
function pushbutton_startMeasurement_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_startMeasurement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end


% --- Executes during object creation, after setting all properties.
function slider13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes during object creation, after setting all properties.
function slider14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes during object creation, after setting all properties.
function popupmenu_PhaseOS2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_PhaseOS2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end


function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double
set(handles.slider13,'Value',str2double(get(hObject,'String'))/max(abs(handles.LocaliserImageVol(2:2:6)-handles.LocaliserImageVol(1:2:6))));
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double
set(handles.slider14,'Value',str2double(get(hObject,'String'))/max(abs(handles.LocaliserImageVol(2:2:6)-handles.LocaliserImageVol(1:2:6))));
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double

end


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double

end


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on slider movement.
function slider13_Callback(hObject, eventdata, handles)
% hObject    handle to slider13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes on slider movement.
function slider14_Callback(hObject, eventdata, handles)
% hObject    handle to slider14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double
set(handles.slider12,'Value',str2double(get(hObject,'String'))/max(abs(handles.LocaliserImageVol(2:2:6)-handles.LocaliserImageVol(1:2:6))));
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double
set(handles.slider11,'Value',str2double(get(hObject,'String'))/max(abs(handles.LocaliserImageVol(2:2:6)-handles.LocaliserImageVol(1:2:6))));

end


function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double

end


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on selection change in popupmenu_PhaseOS2.
function popupmenu_PhaseOS2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_PhaseOS2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_PhaseOS2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_PhaseOS2
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8
figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes during object creation, after setting all properties.
function text35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

end


% --- Executes on button press in pushbutton_Localizer.
function pushbutton_Localizer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Localizer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.ButtonMeasure3DLocaliser=1;
    guidata(hObject,handles)
    figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes when uipanel2 is resized.
function uipanel2_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end


% --- Executes when user attempts to close figureSliceSelect.
function figureSliceSelect_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figureSliceSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% %talker=handles.talker;
% finish
delete(hObject);

end


% --- Executes on button press in pushbutton_anglesZero.
function pushbutton_anglesZero_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_anglesZero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.ButtonZeroapt=1;
    guidata(hObject,handles)
    figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --------------------------------------------------------------------
function uipanel_Dim_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_Dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.ButtonDim=1;
    guidata(hObject,handles)
    figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --------------------------------------------------------------------
function uipanel_readphase_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_readphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.ButtonReadPhase=1;
    guidata(hObject,handles)
    figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


function uipanel_readphase_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_readphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.ButtonReadPhase=1;
    guidata(hObject,handles)
    figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


function uipanel_Dim_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_Dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.ButtonDim=1;
    guidata(hObject,handles)
    figureSliceSelect_WindowButtonMotionFcn(hObject, eventdata, handles)

end


% --- Executes on button press in checkbox_plotTR.
function checkbox_plotTR_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plotTR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plotTR

end


% --- Executes on button press in checkbox_PlotSequence.
function checkbox_PlotSequence_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PlotSequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PlotSequence

end


% --- Executes on button press in checkbox_Plot_kSpace.
function checkbox_Plot_kSpace_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Plot_kSpace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Plot_kSpace

end


% --- Executes on button press in checkbox_image.
function checkbox_image_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_image

end


% --- Executes on button press in checkbox_Extra1.
function checkbox_Extra1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Extra1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Extra1

end


% --- Executes on button press in checkbox_Extra2.
function checkbox_Extra2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Extra2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Extra2

end


% --- Executes on button press in checkbox_PlotPhase.
function checkbox_PlotPhase_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PlotPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PlotPhase

end


% --- Executes on button press in pushbutton_Pathselect.
function pushbutton_Pathselect_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Pathselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isdir(get(handles.edit21,'String'))
    handles.save_FolderName = uigetdir(get(handles.edit21,'String'),'Please select save path');
else
    handles.save_FolderName = uigetdir([handles.HW.RootPath, '\User\Save'],'Please select save path');
end
set(handles.edit21,'String',handles.save_FolderName)

end


% --- Executes on button press in checkbox_addDate.
function checkbox_addDate_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_addDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_addDate

end


% --- Executes on button press in checkbox_Save.
function checkbox_Save_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Save

end


% --- Executes on button press in checkbox_Smooth.
function checkbox_Smooth_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Smooth

end


% --- Save AQSlice configuration to .m file
function err = saveAQSliceToFile(handles, filename)

err = 0;
% add file extension if necessary
if numel(filename) > 2 && ~strcmpi(filename(end-1:end), '.m')
  filename(end+1:end+2) = '.m';
end


[~, funcname] = fileparts(filename);

if ~isvarname(funcname)
  warning('PD:AQSlice:FileInvalid', 'FILENAME "%s" is not a valid Matlab function name.', filename);
  err = -1;
  return;
end

fid = fopen(filename, 'w');

if fid < 0
  warning('PD:AQSlice:FileOpen', 'FILENAME "%s" could not be opened.', filename);
  err = -2;
  return;Origin
end

try
  fprintf(fid, '%%%% This file was automatically created by GUI_SliceSelect.\n');
  fprintf(fid, 'function AQSlice = %s()\n\n', funcname);
  fprintf(fid, '%% You can either call this file as a function or copy the data from the next\n');
  fprintf(fid, '%% line until to the line marked with "%%!%%!%%!%%!%%!%%!%%!%%!%%!%%!%%!%%!%%"\n');
  fprintf(fid, '\n\n%% %% AQSlice settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid, 'iSlice = 1;                                                       %% change the value if necessary\n');
  fprintf(fid, '\n\n%% %% Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid, 'Seq.AQSlice(iSlice).nRead = %3d;                                  %% Number of Pixels in read, if nRead>1 nPhase(1)=1\n', handles.Seq.AQSlice(1).nRead);
  fprintf(fid, 'Seq.AQSlice(iSlice).nPhase(1) = %3d;                              %% Number of Pixels in phase(1)\n', handles.Seq.AQSlice(1).nPhase(1));
  fprintf(fid, 'Seq.AQSlice(iSlice).nPhase(2) = %3d;                              %% Number of Pixels in phase(2)\n', handles.Seq.AQSlice(1).nPhase(2));
  fprintf(fid, 'Seq.AQSlice(iSlice).nPhase(3) = %3d;                              %% Number of Pixels in phase(3) for CSI. If nPhase(3)>1, nRead=1 and sizeRead=1e12\n', handles.Seq.AQSlice(1).nPhase(3));
  fprintf(fid, 'Seq.AQSlice(iSlice).HzPerPixMin = %5d;                          %% Bandwidth per pixel in Hz (1/HzPerPixMin= duration of AQ)\n', handles.Seq.AQSlice(1).HzPerPixMin);
  fprintf(fid, 'Seq.AQSlice(iSlice).sizeRead = %.16e;            %% Image size in read in meter (for CSI set to 1e12)\n', handles.Seq.AQSlice(1).sizeRead);
  fprintf(fid, 'Seq.AQSlice(iSlice).sizePhase(1) = %.16e;        %% Image size in phase(1) in meter\n', handles.Seq.AQSlice(1).sizePhase(1));
  fprintf(fid, 'Seq.AQSlice(iSlice).sizePhase(2) = %.16e;        %% Image size in phase(2) in meter\n', handles.Seq.AQSlice(1).sizePhase(2));
  fprintf(fid, 'Seq.AQSlice(iSlice).sizePhase(3) = %.16e;        %% Image size in phase(3) in meter\n', handles.Seq.AQSlice(1).sizePhase(3));
  fprintf(fid, 'Seq.AQSlice(iSlice).thickness = %30s;   %% Slice thickness in meter, used for 2D and 3D! ([] for no slice)\n', matlab.unittest.diagnostics.ConstraintDiagnostic.getDisplayableString(handles.Seq.AQSlice(1).thickness));
  if ~isempty(handles.Seq.AQSlice(1).excitationPulse)
    fprintf(fid, 'Seq.AQSlice(iSlice).excitationPulse = @%-30s;  %% excitation pulse function (type "Pulse_" than press tab for selection of pulses)\n', char(handles.Seq.AQSlice(1).excitationPulse));
  end
  if ~isempty(handles.Seq.AQSlice(1).inversionPulse)
    fprintf(fid, 'Seq.AQSlice(iSlice).inversionPulse = @%-30s;   %% inversion pulse function\n', char(handles.Seq.AQSlice(1).inversionPulse));
  end

  fprintf(fid, '\n\n%% %% Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid, 'Seq.AQSlice(iSlice).ReadOS = %30s;      %% Oversampling read ([] for automatic, recommended >=2) 1...\n', matlab.unittest.diagnostics.ConstraintDiagnostic.getDisplayableString(handles.Seq.AQSlice(1).ReadOS));
  fprintf(fid, 'Seq.AQSlice(iSlice).PhaseOS(1) = %3d;                             %% Oversampling phase(1)  1...\n', handles.Seq.AQSlice(1).PhaseOS(1));
  fprintf(fid, 'Seq.AQSlice(iSlice).PhaseOS(2) = %3d;                             %% Oversampling phase(2)  1...\n', handles.Seq.AQSlice(1).PhaseOS(2));
  fprintf(fid, 'Seq.AQSlice(iSlice).PhaseOS(3) = %3d;                             %% Oversampling phase(3)  1... (set to 1 if nPhase(3)=1;\n', handles.Seq.AQSlice(1).PhaseOS(3));
  fprintf(fid, 'Seq.AQSlice(iSlice).PhaseOS(Seq.AQSlice(iSlice).nPhase==1) = 1;   %% Reset PhaseOS(x) to 1 if nPhase(x) == 1\n');

  fprintf(fid, '\n\n%% %% Orientation in Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  % 18 decimals should be of sufficient precission to represent an angle in double float
  fprintf(fid, 'Seq.AQSlice(iSlice).alfa  = %+.18f;                %% first rotation around x axis in RAD\n', handles.Seq.AQSlice(1).alfa);
  fprintf(fid, 'Seq.AQSlice(iSlice).phi   = %+.18f;                %% second rotation around x axis in RAD\n', handles.Seq.AQSlice(1).phi);
  fprintf(fid, 'Seq.AQSlice(iSlice).theta = %+.18f;                %% third rotation around x axis in RAD\n', handles.Seq.AQSlice(1).theta);
  fprintf(fid, 'Seq.AQSlice(iSlice).Center2OriginImage = [%.16e, %.16e, %.16e]; %% vector from center of the image to the origin in image coordinate system in meter\n', handles.Seq.AQSlice(1).Center2OriginImage);

  fprintf(fid, '\n\n\n%%!%%!%%!%%!%%!%%!%%!%%!%%!%%!%%!%%!%%\n\n\n');
  fprintf(fid, 'AQSlice = Seq.AQSlice(iSlice);\n\nend\n');

  fprintf('AQSlice data saved in file "%s".\n', which(fopen(fid))); % print full path to file in Command Window
catch ME
  warning(getReport(ME));
  err = -3;
end

fclose(fid);

end
