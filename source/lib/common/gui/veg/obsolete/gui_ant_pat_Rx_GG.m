function varargout = gui_ant_pat_Rx_GG(varargin)
% GUI_ANT_PAT_RX_GG MATLAB code for gui_ant_pat_Rx_GG.fig
%      GUI_ANT_PAT_RX_GG, by itself, creates a new GUI_ANT_PAT_RX_GG or raises the existing
%      singleton*.
%
%      H = GUI_ANT_PAT_RX_GG returns the handle to a new GUI_ANT_PAT_RX_GG or the handle to
%      the existing singleton*.
%
%      GUI_ANT_PAT_RX_GG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ANT_PAT_RX_GG.M with the given input arguments.
%
%      GUI_ANT_PAT_RX_GG('Property','Value',...) creates a new GUI_ANT_PAT_RX_GG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_ant_pat_Rx_GG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_ant_pat_Rx_GG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_ant_pat_Rx_GG

% Last Modified by GUIDE v2.5 24-Jul-2018 20:25:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_ant_pat_Rx_GG_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_ant_pat_Rx_GG_OutputFcn, ...
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


% --- Executes just before gui_ant_pat_Rx_GG is made visible.
function gui_ant_pat_Rx_GG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_ant_pat_Rx_GG (see VARARGIN)

% % Choose default command line output for gui_ant_pat_Rx_GG
% handles.output = hObject;
% 
% % Update handles structure
% guidata(hObject, handles);

clearvars -global antPatRxGgGUI
global antPatRxGgGUI
antPatRxGgGUI = AntPatRxGgGUIManager(handles);

% Wait for thee user response
uiwait( hObject );


% --- Outputs from this function are returned to the command line.
function varargout = gui_ant_pat_Rx_GG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;

global antPatRxGgGUI
    
% If clicked the exit button
if(~isstruct(handles))
    varargout = cell(3,1);
    return
end
    
varargout = antPatRxGgGUI.outputFun();

%close main panel
delete(gcf);



function edit_hpbw_deg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hpbw_deg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hpbw_deg as text
%        str2double(get(hObject,'String')) returns contents of edit_hpbw_deg as a double

global antPatRxGgGUI
antPatRxGgGUI.syncFromGUI( antPatRxGgGUI.uiIDs.edit_hpbw_deg );


% --- Executes during object creation, after setting all properties.
function edit_hpbw_deg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hpbw_deg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SLL_dB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SLL_dB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SLL_dB as text
%        str2double(get(hObject,'String')) returns contents of edit_SLL_dB as a double

global antPatRxGgGUI
antPatRxGgGUI.syncFromGUI( antPatRxGgGUI.uiIDs.edit_SLL_dB );


% --- Executes during object creation, after setting all properties.
function edit_SLL_dB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SLL_dB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_ok.
function pb_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global antPatRxGgGUI
antPatRxGgGUI.syncFromGUI( antPatRxGgGUI.uiIDs.pb_ok );


% --- Executes on button press in pb_cancel.
function pb_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global antPatRxGgGUI
antPatRxGgGUI.syncFromGUI( antPatRxGgGUI.uiIDs.pb_cancel );



function edit_XPL_dB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_XPL_dB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_XPL_dB as text
%        str2double(get(hObject,'String')) returns contents of edit_XPL_dB as a double

global antPatRxGgGUI
antPatRxGgGUI.syncFromGUI( antPatRxGgGUI.uiIDs.edit_XPL_dB );


% --- Executes during object creation, after setting all properties.
function edit_XPL_dB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_XPL_dB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
