function varargout = gui_SCoBi(varargin)
% GUI_SCOBI MATLAB code for gui_SCoBi.fig
%      GUI_SCOBI, by itself, creates a new GUI_SCOBI or raises the existing
%      singleton*.
%
%      H = GUI_SCOBI returns the handle to a new GUI_SCOBI or the handle to
%      the existing singleton*.
%
%      GUI_SCOBI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SCOBI.M with the given input arguments.
%
%      GUI_SCOBI('Property','Value',...) creates a new GUI_SCOBI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_SCoBi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_SCoBi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_SCoBi

% Last Modified by GUIDE v2.5 02-Oct-2018 13:14:40

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


% --- Executes just before gui_SCoBi is made visible.
function gui_SCoBi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_SCoBi (see VARARGIN)

% % Choose default command line output for gui_SCoBi
% handles.output = hObject;
% 
% % Update handles structure
% guidata(hObject, handles);

clearvars -global guiSCoBiManager
global guiSCoBiManager
simulator_id = varargin{1};
guiSCoBiManager = gui_SCoBi_Manager( handles, simulator_id );

% Wait for thee user response
uiwait( hObject );



% --- Outputs from this function are returned to the command line.
function varargout = gui_SCoBi_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % Get default command line output from handles structure
% varargout{1} = handles.output;
global guiSCoBiManager
    
% If clicked the exit button
if(~isstruct(handles))
    varargout = cell(1,2);
    return
end
    
varargout = guiSCoBiManager.outputFun();

%close main panel
delete(gcf);




% --- Executes on selection change in popup_sim_mode.
function popup_sim_mode_Callback(hObject, eventdata, handles)
% hObject    handle to popup_sim_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_sim_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_sim_mode

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.popup_sim_mode );


% --- Executes during object creation, after setting all properties.
function popup_sim_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_sim_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_campaign_Callback(hObject, eventdata, handles)
% hObject    handle to edit_campaign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_campaign as text
%        str2double(get(hObject,'String')) returns contents of edit_campaign as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_campaign );


% --- Executes during object creation, after setting all properties.
function edit_campaign_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_campaign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_gnd_cover.
function popup_gnd_cover_Callback(hObject, eventdata, handles)
% hObject    handle to popup_gnd_cover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_gnd_cover contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_gnd_cover

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.popup_gnd_cover );


% --- Executes during object creation, after setting all properties.
function popup_gnd_cover_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_gnd_cover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_write_attenuation.
function cb_write_attenuation_Callback(hObject, eventdata, handles)
% hObject    handle to cb_write_attenuation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_write_attenuation

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.cb_write_attenuation );


% --- Executes on button press in cb_include_in_master_sim_file.
function cb_include_in_master_sim_file_Callback(hObject, eventdata, handles)
% hObject    handle to cb_include_in_master_sim_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_include_in_master_sim_file

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.cb_include_in_master_sim_file );



function edit_f_MHz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_f_MHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_f_MHz as text
%        str2double(get(hObject,'String')) returns contents of edit_f_MHz as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_f_MHz );


% --- Executes during object creation, after setting all properties.
function edit_f_MHz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_f_MHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_r_Tx_km_Callback(hObject, eventdata, handles)
% hObject    handle to edit_r_Tx_km (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_r_Tx_km as text
%        str2double(get(hObject,'String')) returns contents of edit_r_Tx_km as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_r_Tx_km );


% --- Executes during object creation, after setting all properties.
function edit_r_Tx_km_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_r_Tx_km (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_EIRP_dB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_EIRP_dB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_EIRP_dB as text
%        str2double(get(hObject,'String')) returns contents of edit_EIRP_dB as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_EIRP_dB );


% --- Executes during object creation, after setting all properties.
function edit_EIRP_dB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_EIRP_dB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_pol_Tx.
function popup_pol_Tx_Callback(hObject, eventdata, handles)
% hObject    handle to popup_pol_Tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_pol_Tx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_pol_Tx

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.popup_pol_Tx );


% --- Executes during object creation, after setting all properties.
function popup_pol_Tx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_pol_Tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_hr_m_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hr_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hr_m as text
%        str2double(get(hObject,'String')) returns contents of edit_hr_m as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_hr_m );


% --- Executes during object creation, after setting all properties.
function edit_hr_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hr_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_G0r_dB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_G0r_dB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_G0r_dB as text
%        str2double(get(hObject,'String')) returns contents of edit_G0r_dB as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_G0r_dB );


% --- Executes during object creation, after setting all properties.
function edit_G0r_dB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_G0r_dB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_hpbw_deg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hpbw_deg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hpbw_deg as text
%        str2double(get(hObject,'String')) returns contents of edit_hpbw_deg as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_hpbw_deg );


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

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_SLL_dB );


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



function edit_XPL_dB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_XPL_dB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_XPL_dB as text
%        str2double(get(hObject,'String')) returns contents of edit_XPL_dB as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_XPL_dB );


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


% --- Executes on selection change in popup_pol_Rx.
function popup_pol_Rx_Callback(hObject, eventdata, handles)
% hObject    handle to popup_pol_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_pol_Rx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_pol_Rx

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.popup_pol_Rx );


% --- Executes during object creation, after setting all properties.
function popup_pol_Rx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_pol_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_load_inputs.
function pb_load_inputs_Callback(hObject, eventdata, handles)
% hObject    handle to pb_load_inputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_load_inputs );


% --- Executes on button press in pb_save_inputs.
function pb_save_inputs_Callback(hObject, eventdata, handles)
% hObject    handle to pb_save_inputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_save_inputs );


% --- Executes on button press in pb_SCoBi.
function pb_SCoBi_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_SCoBi );


% --- Executes on button press in pb_exit.
function pb_exit_Callback(hObject, eventdata, handles)
% hObject    handle to pb_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_exit );


        


% --- Executes on selection change in popup_ant_pat_Rx.
function popup_ant_pat_Rx_Callback(hObject, eventdata, handles)
% hObject    handle to popup_ant_pat_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_ant_pat_Rx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_ant_pat_Rx

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.popup_ant_pat_Rx );


% --- Executes during object creation, after setting all properties.
function popup_ant_pat_Rx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_ant_pat_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_diel_model.
function popup_diel_model_Callback(hObject, eventdata, handles)
% hObject    handle to popup_diel_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_diel_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_diel_model

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.popup_diel_model );


% --- Executes during object creation, after setting all properties.
function popup_diel_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_diel_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_ant_pat_Rx_file.
function pb_ant_pat_Rx_file_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ant_pat_Rx_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_ant_pat_Rx_file );


% --- Executes on button press in pb_config_inputs_file.
function pb_config_inputs_file_Callback(hObject, eventdata, handles)
% hObject    handle to pb_config_inputs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_config_inputs_file );


% --- Executes on button press in pb_veg_inputs_file.
function pb_veg_inputs_file_Callback(hObject, eventdata, handles)
% hObject    handle to pb_veg_inputs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_veg_inputs_file );



function edit_config_inputs_file_Callback(hObject, eventdata, handles)
% hObject    handle to edit_config_inputs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_config_inputs_file as text
%        str2double(get(hObject,'String')) returns contents of edit_config_inputs_file as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_config_inputs_file );


% --- Executes during object creation, after setting all properties.
function edit_config_inputs_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_config_inputs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ant_pat_Rx_file_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ant_pat_Rx_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ant_pat_Rx_file as text
%        str2double(get(hObject,'String')) returns contents of edit_ant_pat_Rx_file as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_ant_pat_Rx_file );


% --- Executes during object creation, after setting all properties.
function edit_ant_pat_Rx_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ant_pat_Rx_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_veg_inputs_file_Callback(hObject, eventdata, handles)
% hObject    handle to edit_veg_inputs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_veg_inputs_file as text
%        str2double(get(hObject,'String')) returns contents of edit_veg_inputs_file as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_veg_inputs_file );


% --- Executes during object creation, after setting all properties.
function edit_veg_inputs_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_veg_inputs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on edit_config_inputs_file and none of its controls.
function edit_config_inputs_file_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_config_inputs_file (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popup_orientation_Rx.
function popup_orientation_Rx_Callback(hObject, eventdata, handles)
% hObject    handle to popup_orientation_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_orientation_Rx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_orientation_Rx

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.popup_orientation_Rx );


% --- Executes during object creation, after setting all properties.
function popup_orientation_Rx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_orientation_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_th0_Rx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_th0_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_th0_Rx as text
%        str2double(get(hObject,'String')) returns contents of edit_th0_Rx as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_th0_Rx );


% --- Executes during object creation, after setting all properties.
function edit_th0_Rx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_th0_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ph0_Rx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ph0_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ph0_Rx as text
%        str2double(get(hObject,'String')) returns contents of edit_ph0_Rx as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_ph0_Rx );


% --- Executes during object creation, after setting all properties.
function edit_ph0_Rx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ph0_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ant_pat_res_Rx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ant_pat_res_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ant_pat_res_Rx as text
%        str2double(get(hObject,'String')) returns contents of edit_ant_pat_res_Rx as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_ant_pat_res_Rx );


% --- Executes during object creation, after setting all properties.
function edit_ant_pat_res_Rx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ant_pat_res_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_ant_pat_Rx_GG_input.
function pb_ant_pat_Rx_GG_input_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ant_pat_Rx_GG_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_ant_pat_Rx_GG_input );


% --- Executes during object creation, after setting all properties.
function panel_main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to panel_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popup_orientation_Tx.
function popup_orientation_Tx_Callback(hObject, eventdata, handles)
% hObject    handle to popup_orientation_Tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_orientation_Tx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_orientation_Tx

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.popup_orientation_Tx );


% --- Executes during object creation, after setting all properties.
function popup_orientation_Tx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_orientation_Tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_th0_Tx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_th0_Tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_th0_Tx as text
%        str2double(get(hObject,'String')) returns contents of edit_th0_Tx as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_th0_Tx );


% --- Executes during object creation, after setting all properties.
function edit_th0_Tx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_th0_Tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ph0_Tx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ph0_Tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ph0_Tx as text
%        str2double(get(hObject,'String')) returns contents of edit_ph0_Tx as a double

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.edit_ph0_Tx );


% --- Executes during object creation, after setting all properties.
function edit_ph0_Tx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ph0_Tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_SCoBi_Illustration.
function pb_SCoBi_Illustration_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi_Illustration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_SCoBi_Illustration );


% --- Executes on button press in pb_Forest.
function pb_Forest_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Forest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_Forest );


% --- Executes on button press in pb_Snow.
function pb_Snow_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Snow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_Snow );


% --- Executes on button press in pb_Soil.
function pb_Soil_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Soil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_Soil );


% --- Executes on button press in pb_Topography.
function pb_Topography_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Topography (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_Topography );


% --- Executes on button press in pb_Root_zone.
function pb_Root_zone_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Root_zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_Root_zone );


% --- Executes on button press in pb_Permafrost.
function pb_Permafrost_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Permafrost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_Permafrost );


% --- Executes on button press in pb_Agriculture.
function pb_Agriculture_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Agriculture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.pb_Agriculture );


% --- Executes on selection change in popup_gnd_structure.
function popup_gnd_structure_Callback(hObject, eventdata, handles)
% hObject    handle to popup_gnd_structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_gnd_structure contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_gnd_structure

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.popup_gnd_structure );


% --- Executes during object creation, after setting all properties.
function popup_gnd_structure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_gnd_structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_3rd_order.
function cb_3rd_order_Callback(hObject, eventdata, handles)
% hObject    handle to cb_3rd_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_3rd_order

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.cb_3rd_order );


% --- Executes on button press in cb_2nd_order.
function cb_2nd_order_Callback(hObject, eventdata, handles)
% hObject    handle to cb_2nd_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_2nd_order

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.cb_2nd_order );


% --- Executes on button press in cb_logistic_regression.
function cb_logistic_regression_Callback(hObject, eventdata, handles)
% hObject    handle to cb_logistic_regression (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_logistic_regression

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.cb_logistic_regression );


% --- Executes on button press in cb_discrete_slab.
function cb_discrete_slab_Callback(hObject, eventdata, handles)
% hObject    handle to cb_discrete_slab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_discrete_slab

global guiSCoBiManager
guiSCoBiManager.syncFromGUI( guiSCoBiManager.uiIDs.cb_discrete_slab );


% --- Executes on button press in pb_Wetland.
function pb_Wetland_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Wetland (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
