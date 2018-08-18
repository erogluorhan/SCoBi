function varargout = gui_SCoBi_Veg(varargin)
% GUI_SCOBI_VEG MATLAB code for gui_SCoBi_Veg.fig
%      GUI_SCOBI_VEG, by itself, creates a new GUI_SCOBI_VEG or raises the existing
%      singleton*.
%
%      H = GUI_SCOBI_VEG returns the handle to a new GUI_SCOBI_VEG or the handle to
%      the existing singleton*.
%
%      GUI_SCOBI_VEG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SCOBI_VEG.M with the given input arguments.
%
%      GUI_SCOBI_VEG('Property','Value',...) creates a new GUI_SCOBI_VEG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_SCoBi_Veg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_SCoBi_Veg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_SCoBi_Veg

% Last Modified by GUIDE v2.5 16-Aug-2018 13:46:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_SCoBi_Veg_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_SCoBi_Veg_OutputFcn, ...
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


% --- Executes just before gui_SCoBi_Veg is made visible.
function gui_SCoBi_Veg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_SCoBi_Veg (see VARARGIN)

% % Choose default command line output for gui_SCoBi_Veg
% handles.output = hObject;
% 
% % Update handles structure
% guidata(hObject, handles);

clearvars -global guiSCoBiVegManager
global guiSCoBiVegManager
simulator_id = varargin{1};
guiSCoBiVegManager = gui_SCoBi_Veg_Manager( handles, simulator_id );

% Wait for thee user response
uiwait( hObject );



% --- Outputs from this function are returned to the command line.
function varargout = gui_SCoBi_Veg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % Get default command line output from handles structure
% varargout{1} = handles.output;
global guiSCoBiVegManager
    
% If clicked the exit button
if(~isstruct(handles))
    varargout = cell(1,1);
    return
end
    
varargout = guiSCoBiVegManager.outputFun();

%close main panel
delete(gcf);




% --- Executes on selection change in popup_sim_mode.
function popup_sim_mode_Callback(hObject, eventdata, handles)
% hObject    handle to popup_sim_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_sim_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_sim_mode

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.popup_sim_mode );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_campaign );


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



function edit_campaign_date_Callback(hObject, eventdata, handles)
% hObject    handle to edit_campaign_date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_campaign_date as text
%        str2double(get(hObject,'String')) returns contents of edit_campaign_date as a double

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_campaign_date );


% --- Executes during object creation, after setting all properties.
function edit_campaign_date_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_campaign_date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_plot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_plot as text
%        str2double(get(hObject,'String')) returns contents of edit_plot as a double

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_plot );


% --- Executes during object creation, after setting all properties.
function edit_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_plot (see GCBO)
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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.popup_gnd_cover );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.cb_write_attenuation );


% --- Executes on button press in cb_include_in_master_sim_file.
function cb_include_in_master_sim_file_Callback(hObject, eventdata, handles)
% hObject    handle to cb_include_in_master_sim_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_include_in_master_sim_file

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.cb_include_in_master_sim_file );



function edit_veg_plant_Callback(hObject, eventdata, handles)
% hObject    handle to edit_veg_plant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_veg_plant as text
%        str2double(get(hObject,'String')) returns contents of edit_veg_plant as a double

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_veg_plant );


% --- Executes during object creation, after setting all properties.
function edit_veg_plant_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_veg_plant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Nr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Nr as text
%        str2double(get(hObject,'String')) returns contents of edit_Nr as a double

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_Nr );


% --- Executes during object creation, after setting all properties.
function edit_Nr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Nr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Nfz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Nfz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Nfz as text
%        str2double(get(hObject,'String')) returns contents of edit_Nfz as a double

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_Nfz );


% --- Executes during object creation, after setting all properties.
function edit_Nfz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Nfz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_f_MHz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_f_MHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_f_MHz as text
%        str2double(get(hObject,'String')) returns contents of edit_f_MHz as a double

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_f_MHz );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_r_Tx_km );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_EIRP_dB );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.popup_pol_Tx );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_hr_m );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_G0r_dB );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_hpbw_deg );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_SLL_dB );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_XPL_dB );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.popup_pol_Rx );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.pb_load_inputs );


% --- Executes on button press in pb_save_inputs.
function pb_save_inputs_Callback(hObject, eventdata, handles)
% hObject    handle to pb_save_inputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.pb_save_inputs );


% --- Executes on button press in pb_SCoBi.
function pb_SCoBi_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SCoBi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.pb_SCoBi );


% --- Executes on button press in pb_exit.
function pb_exit_Callback(hObject, eventdata, handles)
% hObject    handle to pb_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.pb_exit );


        


% --- Executes on selection change in popup_ant_pat_Rx.
function popup_ant_pat_Rx_Callback(hObject, eventdata, handles)
% hObject    handle to popup_ant_pat_Rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_ant_pat_Rx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_ant_pat_Rx

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.popup_ant_pat_Rx );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.popup_diel_model );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.pb_ant_pat_Rx_file );


% --- Executes on button press in pb_config_inputs_file.
function pb_config_inputs_file_Callback(hObject, eventdata, handles)
% hObject    handle to pb_config_inputs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.pb_config_inputs_file );


% --- Executes on button press in pb_veg_inputs_file.
function pb_veg_inputs_file_Callback(hObject, eventdata, handles)
% hObject    handle to pb_veg_inputs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.pb_veg_inputs_file );



function edit_config_inputs_file_Callback(hObject, eventdata, handles)
% hObject    handle to edit_config_inputs_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_config_inputs_file as text
%        str2double(get(hObject,'String')) returns contents of edit_config_inputs_file as a double

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_config_inputs_file );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_ant_pat_Rx_file );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_veg_inputs_file );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.popup_orientation_Rx );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_th0_Rx );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_ph0_Rx );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.edit_ant_pat_res_Rx );


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

global guiSCoBiVegManager
guiSCoBiVegManager.syncFromGUI( guiSCoBiVegManager.uiIDs.pb_ant_pat_Rx_GG_input );
