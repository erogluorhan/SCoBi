function varargout = gui_SCoBiMain(varargin)
% GUI_SCOBIMAIN MATLAB code for gui_SCoBiMain.fig
%      GUI_SCOBIMAIN, by itself, creates a new GUI_SCOBIMAIN or raises the existing
%      singleton*.
%
%      H = GUI_SCOBIMAIN returns the handle to a new GUI_SCOBIMAIN or the handle to
%      the existing singleton*.
%
%      GUI_SCOBIMAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SCOBIMAIN.M with the given input arguments.
%
%      GUI_SCOBIMAIN('Property','Value',...) creates a new GUI_SCOBIMAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_SCoBi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_SCoBi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_SCoBiMain

% Last Modified by GUIDE v2.5 26-Oct-2018 12:50:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_SCoBi_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_SCoBi_OutputFcn, ...
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


% --- Executes just before gui_SCoBiMain is made visible.
function gui_SCoBi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_SCoBiMain (see VARARGIN)

% Choose default command line output for gui_SCoBiMain
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

clearvars -global scobiMainGUI
global scobiMainGUI
scobiMainGUI = gui_SCoBiMain_Manager(handles);

% Wait for thee user response
uiwait( hObject );


% %% GET GLOBAL DIRECTORIES
% dir_gui = Directories.getInstance.common_gui;
% 
% 
% set(hObject, 'Units', 'pixels');
% 
% imgfile = [dir_gui filesep 'SCoBI_Illustration.PNG'];
% handles.banner = imread(imgfile); % Read the image file banner.jpg
% info = imfinfo(imgfile); % Determine the size of the image file
% position = get(hObject, 'Position');
% set(hObject, 'Position', [position(1:2) info.Width + 110 info.Height + 105]);
% axes(handles.axes1);
% image(handles.banner)
% set(handles.axes1, ...
%     'Visible', 'off', ...
%     'Units', 'pixels', ...
%     'Position', [50 50 info.Width info.Height]);


% --- Outputs from this function are returned to the command line.
function varargout = gui_SCoBi_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
    varargout{1} = handles.output;
end

global scobiMainGUI
    
% If clicked the exit button
if(~isstruct(handles))
    varargout = cell(30,1);
    return
end
    
varargout = scobiMainGUI.outputFun();

%close main panel
delete(gcf);


% --- Executes on button press in pb_SCoBi_main.
function pb_SCoBi_main_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_SCoBi_Permafrost.
function pb_SCoBi_Permafrost_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi_Permafrost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global scobiMainGUI
scobiMainGUI.syncFromGUI( scobiMainGUI.uiIDs.pb_SCoBi_Permafrost );


% --- Executes on button press in pb_SCoBi_Agriculture.
function pb_SCoBi_Agriculture_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi_Agriculture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global scobiMainGUI
scobiMainGUI.syncFromGUI( scobiMainGUI.uiIDs.pb_SCoBi_Agriculture );


% --- Executes on button press in pb_SCoBi_Forest.
function pb_SCoBi_Forest_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi_Forest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global scobiMainGUI
scobiMainGUI.syncFromGUI( scobiMainGUI.uiIDs.pb_SCoBi_Forest );


% --- Executes on button press in pb_SCoBi_Topography.
function pb_SCoBi_Topography_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi_Topography (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global scobiMainGUI
scobiMainGUI.syncFromGUI( scobiMainGUI.uiIDs.pb_SCoBi_Topography );


% --- Executes on button press in pb_SCoBi_Soil.
function pb_SCoBi_Root_Zone_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi_Soil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global scobiMainGUI
scobiMainGUI.syncFromGUI( scobiMainGUI.uiIDs.pb_SCoBi_Root_Zone );


% --- Executes on button press in pb_SCoBi_Snow.
function pb_SCoBi_Snow_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi_Snow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global scobiMainGUI
scobiMainGUI.syncFromGUI( scobiMainGUI.uiIDs.pb_SCoBi_Snow );


% --- Executes on button press in pb_SCoBi_Soil.
function pb_SCoBi_Soil_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi_Soil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global scobiMainGUI
scobiMainGUI.syncFromGUI( scobiMainGUI.uiIDs.pb_SCoBi_Soil );


% --- Executes on button press in pb_SCoBi_Wetland.
function pb_SCoBi_Wetland_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi_Wetland (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global scobiMainGUI
scobiMainGUI.syncFromGUI( scobiMainGUI.uiIDs.pb_SCoBi_Wetland );


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in pb_about.
function pb_about_Callback(hObject, eventdata, handles)
% hObject    handle to pb_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global scobiMainGUI
scobiMainGUI.syncFromGUI( scobiMainGUI.uiIDs.pb_about );


% --- Executes on button press in pb_documents.
function pb_documents_Callback(hObject, eventdata, handles)
% hObject    handle to pb_documents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global scobiMainGUI
scobiMainGUI.syncFromGUI( scobiMainGUI.uiIDs.pb_documents );
