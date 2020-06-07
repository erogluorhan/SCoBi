
classdef gui_SCoBi_Manager < SCoBiGUIManagers
% class gui_SCoBi_Manager 
%
%   Implements handles for gui_SCoBi.fig and gui_SCoBi.m.
%   
%   [simulator_id, inputStruct ] = gui_SCoBi(simulator_id);
%
%   INPUT:
%   simulator_id:   An integer that should be given considering 
%   the Constants.SIMULATORS cell array.
%
%   See also SCoBiGUIManagers, gui_SCoBiMain_Manager.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd
%   Adapted from gui_goGPS.m class of goGPS v0.4.3 software

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0

    %GUI_SCOBI_MANAGER This class implements handles for GUI elements, and 
    % performs the GUI-related functionalities
    
    
    properties(GetAccess = 'private', SetAccess = 'private')
        
        % Struct data to provide input parameters as output of GUI
        inputStruct
           
    end
    
    methods
        
        % Constructor
        function obj = gui_SCoBi_Manager( handles, simulator_id )           
            
            % Call superclass constructor
            obj = obj@SCoBiGUIManagers( handles, simulator_id );
            
        end


        % Check if ant_pat_Rx value is Cosine to the power n
        function result = is_popup_ant_pat_Rx_cos_pow_n(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_ant_pat_Rx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_ant_pat_Rx) == Constants.ID_RX_COS_POW_N );
            
        end


        % Check if ant_pat_Rx value is Generalized-Gaussian
        function result = is_popup_ant_pat_Rx_GG(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_ant_pat_Rx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_ant_pat_Rx) == Constants.ID_RX_GG );
            
        end


        % Check if ant_pat_Rx value is User-define
        function result = is_popup_ant_pat_Rx_user_defined(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_ant_pat_Rx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_ant_pat_Rx) == Constants.ID_RX_USER_DEFINED );
            
        end


        % Check if diel_model value is Dobson
        function result = is_popup_diel_model_dobson(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_diel_model );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_diel_model) == Constants.ID_DIEL_DOBSON );
            
        end


        % Check if diel_model value is Mironov
        function result = is_popup_diel_model_mironov(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_diel_model );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_diel_model) == Constants.ID_DIEL_MIRONOV );
            
        end


        % Check if diel_model value is Wang
        function result = is_popup_diel_model_wang(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_diel_model );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_diel_model) == Constants.ID_DIEL_WANG );
            
        end


        % Check if gnd_cover value is Bare-soil
        function result = is_popup_gnd_cover_bare_soil(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_gnd_cover );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_gnd_cover) == Constants.ID_BARE_SOIL );
            
        end


        % Check if gnd_cover value is Vegetation
        function result = is_popup_gnd_cover_vegetation(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_gnd_cover );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_gnd_cover) == Constants.ID_VEG_COVER );
            
        end


        % Check if gnd_structure value is Single-layered
        function result = is_popup_gnd_structure_single_layered(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_gnd_structure );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_gnd_structure) == Constants.ID_GND_SINGLE_LAYERED );
            
        end


        % Check if gnd_structure value is Multi-layered
        function result = is_popup_gnd_structure_multi_layered(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_gnd_structure );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_gnd_structure) == Constants.ID_GND_MULTI_LAYERED );
            
        end


        % Check if orientation_Rx value is Fixed
        function result = is_popup_orientation_Rx_fixed(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_orientation_Rx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_orientation_Rx) == Constants.ID_RX_FIXED );
            
        end


        % Check if orientation_Rx value is Specular-facing
        function result = is_popup_orientation_Rx_specular_facing(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_orientation_Rx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_orientation_Rx) == Constants.ID_RX_SPECULAR_FACING );
            
        end


        % Check if orientation_Tx value is Fixed
        function result = is_popup_orientation_Tx_geostationary(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_orientation_Tx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_orientation_Tx) == Constants.ID_TX_GEOSTATIONARY );
            
        end


        % Check if orientation_Tx value is Specular-facing
        function result = is_popup_orientation_Tx_variable(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_orientation_Tx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_orientation_Tx) == Constants.ID_TX_VARIABLE );
            
        end


        % Check if sim_mode value is Snapshot
        function result = is_popup_sim_mode_snapshot(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_sim_mode );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_sim_mode) == Constants.ID_SNAPSHOT );
            
        end


        % Check if sim_mode value is Time-series
        function result = is_popup_sim_mode_time_series(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_sim_mode );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_sim_mode) == Constants.ID_TIME_SERIES );
            
        end
        
        
        function funout = outputFun(obj)
            
            funout{1} = obj.simulator_id;
            
            funout{2} = obj.inputStruct;
            
        end
        
        
        % EVENT MANAGER
        % When an element is modified in the GUI, this must be called!
        % Sync the internal copy of the interface with the GUI
        function syncFromGUI(obj, idEl)   
            
            
            %% GET GLOBAL DIRECTORIES
            dir_input_sys = Directories.getInstance.input_sys;
            dir_input_config = Directories.getInstance.input_config;
            dir_input_ant_pat_Rx = Directories.getInstance.input_ant_pat_Rx;
            dir_input_veg = Directories.getInstance.input_veg;
            
            
            if (nargin == 1)
                idEl = obj.uiIDs.popup_sim_mode;
            end

            % First check if the value of every GUI element is valid
            isGuiValuesValid = obj.checkUIvaluesValidity(idEl);
            
            if isGuiValuesValid
                
                obj.getFlag(idEl) = true;

                % Read all the values of the elements
                obj.getAllElContent();
                obj.setAllElContent();
            
          
                %% SIMULATION SETTINGS
                % sim_mode
                % If sim_mode popup value is changed            
                if sum(intersect(idEl, obj.uiIDs.popup_sim_mode)) > 0  

                  obj.updateGUI();

                end

                % If gnd_cover popup value is changed            
                if sum(intersect(idEl, obj.uiIDs.popup_gnd_cover)) > 0                         

                  % If selected gnd_cover is Bare-soil
                  if obj.is_popup_gnd_cover_bare_soil()

                      obj.setElStatus(obj.uiGroups.on_popup_gnd_cover, 0, 0);                  

                  % Else if selected gnd_cover is Vegetation
                  elseif obj.is_popup_gnd_cover_vegetation()

                      obj.setElStatus([obj.uiGroups.on_popup_gnd_cover], 1, 0);

                  end

                  obj.updateGUI();

                end


                %% TRANSMITTER INPUTS
                % If orientation_Tx popup value is changed            
                if sum(intersect(idEl, obj.uiIDs.popup_orientation_Tx)) > 0                

                  % If selected orientation_Tx is Geo-stationary
                  if obj.is_popup_orientation_Tx_geostationary()

                      obj.setElStatus(obj.uiGroups.on_popup_orientation_Tx_geostationary, 1, 0);  
                      obj.setElStatus(obj.uiGroups.on_popup_orientation_Tx_variable, 0, 0);

                      obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Tx_geostationary, 'on' );
                      obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Tx_variable, 'off' );

                  % Else if selected orientation_Tx is Variable
                  elseif obj.is_popup_orientation_Tx_variable()

                      obj.setElStatus([obj.uiGroups.on_popup_orientation_Tx_geostationary], 0, 0);
                      obj.setElStatus(obj.uiGroups.on_popup_orientation_Tx_variable, 1, 0);

                      obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Tx_geostationary, 'off' );
                      obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Tx_variable, 'on' );

                  end

                  obj.updateGUI();

                end


                %% RECEIVER INPUTS
                % If orientation_Rx popup value is changed            
                if sum(intersect(idEl, obj.uiIDs.popup_orientation_Rx)) > 0                

                  % If selected orientation_Rx is Fixed
                  if obj.is_popup_orientation_Rx_fixed()

                      obj.setElStatus(obj.uiGroups.on_popup_orientation_Rx_fixed, 1, 0);  
                      obj.setElStatus(obj.uiGroups.on_popup_orientation_Rx_specular_facing, 0, 0); 

                      obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Rx_fixed, 'on' );
                      obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Rx_specular_facing, 'off' );                 

                  % Else if selected orientation_Rx is Specular-facing
                  elseif obj.is_popup_orientation_Rx_specular_facing()

                      obj.setElStatus(obj.uiGroups.on_popup_orientation_Rx_fixed, 0, 0);  
                      obj.setElStatus(obj.uiGroups.on_popup_orientation_Rx_specular_facing, 1, 0); 

                      obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Rx_fixed, 'off' );
                      obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Rx_specular_facing, 'on' );

                  end

                  obj.updateGUI();

                end

                % If ant_pat_Rx popup value is changed            
                if sum(intersect(idEl, obj.uiIDs.popup_ant_pat_Rx)) > 0  

                  persistent last_ant_pat_res_Rx_val               

                  % If selected ant_pat_Rx is Generalized-Gaussian
                  if obj.is_popup_ant_pat_Rx_GG()                  

                      if isnumeric( last_ant_pat_res_Rx_val ) && ~isempty(last_ant_pat_res_Rx_val)  

                          obj.setElVal( obj.uiIDs.edit_ant_pat_res_Rx, num2str(last_ant_pat_res_Rx_val) );

                      end

                      obj.setElStatus([obj.uiGroups.on_popup_ant_pat_Rx_GG], 1, 0); 
                      obj.setElStatus([obj.uiGroups.on_popup_ant_pat_Rx_user_defined], 0, 0); 

                      obj.setGuiElVisibility( obj.uiGroups.on_popup_ant_pat_Rx_GG, 'on' );
                      obj.setGuiElVisibility( obj.uiGroups.visibility_on_popup_ant_pat_Rx_user_defined, 'off' );                  


                  % Else if selected ant_pat_Rx is User-defined
                  elseif obj.is_popup_ant_pat_Rx_user_defined()

                      last_ant_pat_res_Rx_val = str2double( obj.getElVal( obj.uiIDs.edit_ant_pat_res_Rx ) );

                      obj.setElStatus([obj.uiGroups.on_popup_ant_pat_Rx_GG], 0, 0); 
                      obj.setElStatus([obj.uiGroups.on_popup_ant_pat_Rx_user_defined], 1, 0); 

                      obj.setGuiElVisibility( obj.uiGroups.on_popup_ant_pat_Rx_GG, 'off' );
                      obj.setGuiElVisibility( obj.uiGroups.visibility_on_popup_ant_pat_Rx_user_defined, 'on' );                

                  % Else if selected ant_pat_Rx is Cosine-to-the-power-n
                  elseif obj.is_popup_ant_pat_Rx_cos_pow_n()

                      % Display a warning that the method is not
                      % implemented yet
                      waitfor(msgbox('WARNING: This method is not yet implemented!'));

                      % Change the popup value to default GG
                      obj.setElVal(obj.uiIDs.popup_ant_pat_Rx, Constants.ID_RX_GG); 

                      obj.syncFromGUI(obj.uiIDs.popup_ant_pat_Rx);

                  end

                  obj.updateGUI();

                end


                %% GROUND INPUTS
                % If diel_model popup value is changed            
                if sum(intersect(idEl, obj.uiIDs.popup_diel_model)) > 0 

                  obj.updateGUI();

                end

                % If gnd_structure popup value is changed            
                if sum(intersect(idEl, obj.uiIDs.popup_gnd_structure)) > 0                

                  % If selected gnd_structure is Single-layered
                  if obj.is_popup_gnd_structure_single_layered()

                      obj.setElStatus(obj.uiGroups.on_popup_gnd_structure_single_layered, 1, 0);  
                      obj.setElStatus(obj.uiGroups.on_popup_gnd_structure_multi_layered, 0, 0); 

                      obj.setGuiElVisibility( obj.uiGroups.on_popup_gnd_structure_single_layered, 'on' );
                      obj.setGuiElVisibility( obj.uiGroups.on_popup_gnd_structure_multi_layered, 'off' );                 

                  % Else if selected gnd_structure is Multi-layered
                  elseif obj.is_popup_gnd_structure_multi_layered()

                      obj.setElStatus(obj.uiGroups.on_popup_gnd_structure_single_layered, 0, 0);  
                      obj.setElStatus(obj.uiGroups.on_popup_gnd_structure_multi_layered, 1, 0); 

                      obj.setGuiElVisibility( obj.uiGroups.on_popup_gnd_structure_single_layered, 'off' );
                      obj.setGuiElVisibility( obj.uiGroups.on_popup_gnd_structure_multi_layered, 'on' );

                  end

                  obj.updateGUI();

                end


                %% PUSH BUTTONS

                % on the sim: About
                if sum(intersect(idEl, obj.uiIDs.pb_about)) > 0

                    obj.openAboutDialog();

                end

                % on the sim: Documents
                if sum(intersect(idEl, obj.uiIDs.pb_documents)) > 0                

                    obj.openDocumentsDialog();

                end

                % on Forest button
                if sum(intersect(idEl, obj.uiIDs.pb_Forest)) > 0

                    obj.init( obj.handles, Constants.ID_SIM_FOREST );

                end

                % on Agriculture button
                if sum(intersect(idEl, obj.uiIDs.pb_Agriculture)) > 0

                    obj.init( obj.handles, Constants.ID_SIM_AGRICULTURE );

                end

                % on Root-zone button
                if sum(intersect(idEl, obj.uiIDs.pb_Root_zone)) > 0

                    obj.init( obj.handles, Constants.ID_SIM_ROOT_ZONE );

                end

                % on Soil button
                if sum(intersect(idEl, obj.uiIDs.pb_Soil)) > 0

                    obj.init( obj.handles, Constants.ID_SIM_SOIL );

                end

                % on the other simulation buttons
                if sum(intersect(idEl, obj.uiIDs.pb_Snow)) > 0 ...
                       || sum(intersect(idEl, obj.uiIDs.pb_Topography)) > 0 ...
                       || sum(intersect(idEl, obj.uiIDs.pb_Permafrost)) > 0

                    waitfor(msgbox('WARNING: This simulation mode has not yet been implemented!'));

                end

                % on Load Settings
                if sum(intersect(idEl, obj.uiIDs.pb_load_inputs)) > 0

                    filter = {strcat(dir_input_sys, '\*.mat')};
                    [file, path] = uigetfile(filter, 'Load GUI From Input File');


                    % If file and path are not empty
                    if length(path)>1 && length(file)>1

                        obj.inputStruct = obj.loadGUIFromInputFile( file, path );

                        obj.saveLastInputFile(path, file);

                    end

                end

                % on Save Settings
                if sum(intersect(idEl, obj.uiIDs.pb_save_inputs)) > 0

                    obj.saveGUIToInputFile( Constants.ID_GUI_SAVE_AS );

                end

                % on SCoBi
                if sum(intersect(idEl, obj.uiIDs.pb_SCoBi)) > 0

                    % Try to save GUI values to an input file first
                    % It will perform a check for any changes to the loaded or recently 
                    % saved inputs   
                    savingResult = obj.saveGUIToInputFile( Constants.ID_GUI_SAVE );

                    if savingResult == 1 || savingResult == 2

                        uiresume(obj.handles.panel_main);

                    end

                end   

                % on Browse Configuration Inputs File
                if sum(intersect(idEl, obj.uiIDs.pb_config_inputs_file)) > 0

                    filter = {strcat(dir_input_config, '\*.xlsx')};
                    [file, path] = uigetfile(filter, 'Load Configuration Inputs File');


                    % If file and path are not empty
                    if length(path)>1 && length(file)>1

                       filename = strcat( path, '\',file );

                       obj.setElVal(obj.uiIDs.edit_config_inputs_file, filename, 0);  

                    end

                end  

                % on Browse Antenna Pattern Inputs File
                if sum(intersect(idEl, obj.uiIDs.pb_ant_pat_Rx_file)) > 0

                    filter = {strcat(dir_input_ant_pat_Rx, '\*.xlsx')};
                    [file, path] = uigetfile(filter, 'Load Antenna Pattern Inputs File');


                    % If file and path are not empty
                    if length(path)>1 && length(file)>1

                       filename = strcat( path, '\',file );

                       obj.setElVal(obj.uiIDs.edit_ant_pat_Rx_file, filename, 0);  

                    end

                end

                % on Browse Antenna Pattern Inputs File
                if sum(intersect(idEl, obj.uiIDs.pb_veg_inputs_file)) > 0

                    filter = {strcat(dir_input_veg, '\*.xlsx')};
                    [file, path] = uigetfile(filter, 'Load Vegetation Inputs File');


                    % If file and path are not empty
                    if length(path)>1 && length(file)>1

                       filename = strcat( path, '\',file );

                       obj.setElVal(obj.uiIDs.edit_veg_inputs_file, filename, 0);  

                    end

                end

                obj.onoffUIEl();
                obj.checkUIdependencies();

                % on Exit
                if sum(intersect(idEl, obj.uiIDs.pb_exit)) > 0
                    obj.closeGUI(obj);
                end
                
            end   % Validity check for GUI elements' values
            
        end     % function syncFromGUI
        
        
        
    end
    
    
    
    % Internal initialization functions
    methods(Access = 'protected')
        
        % Initialize the instance
        function init( obj, handles, simulator_id )
            
            tic;
            
            % Call superclass's init() function
            init@SCoBiGUIManagers( obj, handles, simulator_id );
            
            t0 = toc;
            
            fprintf('SCoBi GUI initialization completed in %.2f seconds\n', t0);
        
        end
        
        function initPopupMenus(obj)
            obj.init_popup_sim_mode();
            obj.init_popup_gnd_cover();
            obj.init_popup_pol_Tx();
            obj.init_popup_pol_Rx();
            obj.init_popup_ant_pat_Rx();
            obj.init_popup_diel_model();
            obj.init_popup_gnd_structure();
        end

        
        % Fill the popup_ant_pat_Rx
        function init_popup_ant_pat_Rx(obj, str)
            
            if nargin < 2
                str = Constants.RX_ANT_PATTERNS;
            end
            
            value = get(obj.handles.popup_ant_pat_Rx,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_ant_pat_Rx,'Value', value);
            set(obj.handles.popup_ant_pat_Rx,'String', str);
        end

        
        % Fill the popup_diel_model
        function init_popup_diel_model(obj, str)
            
            if nargin < 2
                str = Constants.DIEL_MODELS;
            end
            
            value = get(obj.handles.popup_diel_model,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_diel_model,'Value', value);
            set(obj.handles.popup_diel_model,'String', str);
        end

        
        % Fill the popup_gnd_cover
        function init_popup_gnd_cover(obj, str)
            
            if nargin < 2
                str = Constants.GND_COVERS;
            end
            
            value = get(obj.handles.popup_gnd_cover,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_gnd_cover,'Value', value);
            set(obj.handles.popup_gnd_cover,'String', str);
        end

        
        % Fill the popup_gnd_structure
        function init_popup_gnd_structure(obj, str)
            
            if nargin < 2
                str = Constants.GND_STRUCTURES;
            end
            
            value = get(obj.handles.popup_gnd_structure,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_gnd_structure,'Value', value);
            set(obj.handles.popup_gnd_structure,'String', str);
        end

        
        % Fill the popup_orientation_Tx
        function init_popup_orientation_Tx(obj, str)
            
            if nargin < 2
                str = Constants.TX_ORIENTATIONS;
            end
            
            value = get(obj.handles.popup_orientation_Tx,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_orientation_Tx,'Value', value);
            set(obj.handles.popup_orientation_Tx,'String', str);
        end

        
        % Fill the popup_orientation_Rx
        function init_popup_orientation_Rx(obj, str)
            
            if nargin < 2
                str = Constants.RX_ORIENTATIONS;
            end
            
            value = get(obj.handles.popup_orientation_Rx,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_orientation_Rx,'Value', value);
            set(obj.handles.popup_orientation_Rx,'String', str);
        end

        
        % Fill the popup_pol_Rx
        function init_popup_pol_Rx(obj, str)
            
            if nargin < 2
                % Do not show polarization H and V since X and Y stands for
                % them
                str = Constants.POLARIZATIONS(1, Constants.ID_POL_R : Constants.ID_POL_Y);
            end
            
            value = get(obj.handles.popup_pol_Rx,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_pol_Rx,'Value', value);
            set(obj.handles.popup_pol_Rx,'String', str);
        end

        
        % Fill the popup_pol_Tx
        function init_popup_pol_Tx(obj, str)
            
            if nargin < 2
                % Do not show polarization H and V since X and Y stands for
                % them
                str = Constants.POLARIZATIONS(1, Constants.ID_POL_R : Constants.ID_POL_Y);
            end
            
            value = get(obj.handles.popup_pol_Tx,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_pol_Tx,'Value', value);
            set(obj.handles.popup_pol_Tx,'String', str);
        end

        
        % Fill the popup_sim_mode
        function init_popup_sim_mode(obj, str)
            
            if nargin < 2
                str = Constants.SIM_MODES;
            end
            
            value = get(obj.handles.popup_sim_mode,'Value');
            value = min( 1, max(length(str), value) );
            set(obj.handles.popup_sim_mode,'Value', value);
            set(obj.handles.popup_sim_mode,'String', str);
        end


        %% INIT UI IDs
        % Assign an id (integer) to each UI element
        % Steps to add a new UI element and manage it from this class:
        %  1. Add an id (with incremental i value)
        %      i=i+1;  id.<name> = i;
        %  2. Associate to this id the correct handle
        %      pointers(i) = obj.handles.<handle>;
        %  3. Add this element's id to the corresponding group, if any
        %      groupIDs.<group> = [id.<name 1> ... id.<name n>];         
        function initUIAccess(obj)
            
                        
          %% MAIN PANEL  
          i = 1;          id.panel_main = i;                pointers(i) = obj.handles.panel_main;
            
          i = i+1;        id.text_description = i;          pointers(i) = obj.handles.text_description;
          
          
          i = i+1;        id.pb_SCoBi_Illustration = i;     pointers(i) = obj.handles.pb_SCoBi_Illustration;  
          i = i+1;        id.pb_about = i;                  pointers(i) = obj.handles.pb_about; 
          i = i+1;        id.pb_documents = i;              pointers(i) = obj.handles.pb_documents;   
          i = i+1;        id.pb_Forest = i;                 pointers(i) = obj.handles.pb_Forest; 
          i = i+1;        id.pb_Snow = i;                   pointers(i) = obj.handles.pb_Snow;   
          i = i+1;        id.pb_Soil = i;                   pointers(i) = obj.handles.pb_Soil;   
          i = i+1;        id.pb_Topography = i;             pointers(i) = obj.handles.pb_Topography;   
          i = i+1;        id.pb_Root_zone = i;              pointers(i) = obj.handles.pb_Root_zone;   
          i = i+1;        id.pb_Permafrost = i;             pointers(i) = obj.handles.pb_Permafrost;   
          i = i+1;        id.pb_Agriculture = i;            pointers(i) = obj.handles.pb_Agriculture;     
          i = i+1;        id.pb_load_inputs = i;            pointers(i) = obj.handles.pb_load_inputs;          
          i = i+1;        id.pb_save_inputs = i;            pointers(i) = obj.handles.pb_save_inputs;        
          i = i+1;        id.pb_SCoBi = i;                  pointers(i) = obj.handles.pb_SCoBi;       
          i = i+1;        id.pb_exit = i;                   pointers(i) = obj.handles.pb_exit;
            
            
          %% PANELS            
          i = i+1;        id.panel_sim_settings = i;        pointers(i) = obj.handles.panel_sim_settings;
          i = i+1;        id.panel_Tx_inputs = i;           pointers(i) = obj.handles.panel_Tx_inputs;
          i = i+1;        id.panel_Rx_inputs = i;           pointers(i) = obj.handles.panel_Rx_inputs;
          i = i+1;        id.panel_gnd_inputs = i;          pointers(i) = obj.handles.panel_gnd_inputs;
          i = i+1;        id.panel_input_files = i;         pointers(i) = obj.handles.panel_input_files;
                        
          
          %% SIMULATION SETTINGS PANEL ELEMENTS   
          i = i+1;        id.text_campaign = i;                 pointers(i) = obj.handles.text_campaign;
          i = i+1;        id.edit_campaign = i;                 pointers(i) = obj.handles.edit_campaign;         
          i = i+1;        id.text_sim_mode = i;             pointers(i) = obj.handles.text_sim_mode;  
          i = i+1;        id.popup_sim_mode = i;            pointers(i) = obj.handles.popup_sim_mode;
          i = i+1;        id.text_gnd_cover = i;            pointers(i) = obj.handles.text_gnd_cover;
          i = i+1;        id.popup_gnd_cover = i;           pointers(i) = obj.handles.popup_gnd_cover;
          
          i = i+1;        id.panel_preferences = i;                 pointers(i) = obj.handles.panel_preferences;   % Another panel inside a panel
          i = i+1;        id.cb_write_attenuation = i;              pointers(i) = obj.handles.cb_write_attenuation;
          i = i+1;        id.cb_include_in_master_sim_file = i;     pointers(i) = obj.handles.cb_include_in_master_sim_file;
          i = i+1;        id.cb_calculate_penetration_depth = i;     pointers(i) = obj.handles.cb_calculate_penetration_depth;
          
            
          groupIDs.sim_settings = [id.panel_sim_settings id.text_campaign : id.cb_include_in_master_sim_file];
          
          
          %% TRANSMITTER (Tx) INPUTS PANEL ELEMENTS            
          i = i+1;        id.text_f_MHz = i;                pointers(i) = obj.handles.text_f_MHz;
          i = i+1;        id.edit_f_MHz = i;                pointers(i) = obj.handles.edit_f_MHz;
          i = i+1;        id.text_MHz = i;                  pointers(i) = obj.handles.text_MHz;
          i = i+1;        id.text_r_Tx_km = i;              pointers(i) = obj.handles.text_r_Tx_km;
          i = i+1;        id.edit_r_Tx_km = i;              pointers(i) = obj.handles.edit_r_Tx_km;
          i = i+1;        id.text_km_r_Tx = i;              pointers(i) = obj.handles.text_km_r_Tx;
          i = i+1;        id.text_EIRP_dB = i;              pointers(i) = obj.handles.text_EIRP_dB;
          i = i+1;        id.edit_EIRP_dB = i;              pointers(i) = obj.handles.edit_EIRP_dB;
          i = i+1;        id.text_dB_EIRP = i;              pointers(i) = obj.handles.text_dB_EIRP;
          i = i+1;        id.text_pol_Tx = i;               pointers(i) = obj.handles.text_pol_Tx;
          i = i+1;        id.popup_pol_Tx = i;              pointers(i) = obj.handles.popup_pol_Tx;
          i = i+1;        id.text_orientation_Tx = i;       pointers(i) = obj.handles.text_orientation_Tx;
          i = i+1;        id.popup_orientation_Tx = i;      pointers(i) = obj.handles.popup_orientation_Tx;
          % Geo-stationary Transmitter Input
          i = i+1;        id.text_el0_Tx = i;               pointers(i) = obj.handles.text_el0_Tx;
          i = i+1;        id.edit_el0_Tx = i;               pointers(i) = obj.handles.edit_el0_Tx;
          i = i+1;        id.text_deg_el0_Tx = i;           pointers(i) = obj.handles.text_deg_el0_Tx;
          i = i+1;        id.text_ph0_Tx = i;               pointers(i) = obj.handles.text_ph0_Tx;
          i = i+1;        id.edit_ph0_Tx = i;               pointers(i) = obj.handles.edit_ph0_Tx;
          i = i+1;        id.text_deg_ph0_Tx = i;           pointers(i) = obj.handles.text_deg_ph0_Tx;
          % Variable-orientation Transmitter UI element
          i = i+1;        id.text_orientation_Tx_variable = i;           pointers(i) = obj.handles.text_orientation_Tx_variable;
            
          groupIDs.Tx_inputs                            = [id.panel_Tx_inputs, id.text_f_MHz : id.text_orientation_Tx_variable];
          
          
          %% RECEIVER (Rx) INPUTS PANEL ELEMENTS           
          i = i+1;        id.text_hr_m = i;                 pointers(i) = obj.handles.text_hr_m;
          i = i+1;        id.edit_hr_m = i;                 pointers(i) = obj.handles.edit_hr_m;
          i = i+1;        id.text_m_hr = i;                 pointers(i) = obj.handles.text_m_hr;
          i = i+1;        id.text_G0r_dB = i;               pointers(i) = obj.handles.text_G0r_dB;
          i = i+1;        id.edit_G0r_dB = i;               pointers(i) = obj.handles.edit_G0r_dB;
          i = i+1;        id.text_dB_G0r = i;               pointers(i) = obj.handles.text_dB_G0r;
          i = i+1;        id.text_pol_Rx = i;               pointers(i) = obj.handles.text_pol_Rx;
          i = i+1;        id.popup_pol_Rx = i;              pointers(i) = obj.handles.popup_pol_Rx;
          i = i+1;        id.text_orientation_Rx = i;       pointers(i) = obj.handles.text_orientation_Rx;
          i = i+1;        id.popup_orientation_Rx = i;      pointers(i) = obj.handles.popup_orientation_Rx;
          % Fixed-orientation Receiver Input
          i = i+1;        id.text_th0_Rx = i;               pointers(i) = obj.handles.text_th0_Rx;
          i = i+1;        id.edit_th0_Rx = i;               pointers(i) = obj.handles.edit_th0_Rx;
          i = i+1;        id.text_deg_th0_Rx = i;           pointers(i) = obj.handles.text_deg_th0_Rx;
          i = i+1;        id.text_ph0_Rx = i;               pointers(i) = obj.handles.text_ph0_Rx;
          i = i+1;        id.edit_ph0_Rx = i;               pointers(i) = obj.handles.edit_ph0_Rx;
          i = i+1;        id.text_deg_ph0_Rx = i;           pointers(i) = obj.handles.text_deg_ph0_Rx;
          % Specular-facing orientation
          i = i+1;        id.text_orientation_Rx_specular_facing = i;    pointers(i) = obj.handles.text_orientation_Rx_specular_facing;
          % Antenna Pattern -related Input
          i = i+1;        id.panel_ant_pat_Rx = i;          pointers(i) = obj.handles.panel_ant_pat_Rx;
          i = i+1;        id.text_ant_pat_Rx = i;           pointers(i) = obj.handles.text_ant_pat_Rx;
          i = i+1;        id.popup_ant_pat_Rx = i;          pointers(i) = obj.handles.popup_ant_pat_Rx;
          % Antenna Pattern GG Input
          i = i+1;        id.text_hpbw_deg = i;             pointers(i) = obj.handles.text_hpbw_deg;
          i = i+1;        id.edit_hpbw_deg = i;             pointers(i) = obj.handles.edit_hpbw_deg;
          i = i+1;        id.text_deg_hpbw = i;             pointers(i) = obj.handles.text_deg_hpbw;
          i = i+1;        id.text_SLL_dB = i;               pointers(i) = obj.handles.text_SLL_dB;
          i = i+1;        id.edit_SLL_dB = i;               pointers(i) = obj.handles.edit_SLL_dB;
          i = i+1;        id.text_dB_SLL = i;               pointers(i) = obj.handles.text_dB_SLL;
          i = i+1;        id.text_XPL_dB = i;               pointers(i) = obj.handles.text_XPL_dB;
          i = i+1;        id.edit_XPL_dB = i;               pointers(i) = obj.handles.edit_XPL_dB;
          i = i+1;        id.text_dB_XPL = i;               pointers(i) = obj.handles.text_dB_XPL;
          i = i+1;        id.text_ant_pat_res_Rx = i;       pointers(i) = obj.handles.text_ant_pat_res_Rx;
          i = i+1;        id.edit_ant_pat_res_Rx = i;       pointers(i) = obj.handles.edit_ant_pat_res_Rx;
          i = i+1;        id.text_deg_ant_pat_res_Rx = i;   pointers(i) = obj.handles.text_deg_ant_pat_res_Rx;  
          % Antenna Pattern  User-defined
          i = i+1;        id.text_ant_pat_Rx_user_defined = i;	pointers(i) = obj.handles.text_ant_pat_Rx_user_defined;
            
          groupIDs.Rx_inputs                               = [id.panel_Rx_inputs, id.text_hr_m : id.popup_ant_pat_Rx];
          groupIDs.ant_pat_Rx_GG_inputs                    = id.text_hpbw_deg : id.text_deg_ant_pat_res_Rx;
          
          
          %% GROUND INPUTS PANEL ELEMENTS               
          i = i+1;        id.text_diel_model = i;           pointers(i) = obj.handles.text_diel_model;  
          i = i+1;        id.popup_diel_model = i;          pointers(i) = obj.handles.popup_diel_model;
          i = i+1;        id.text_gnd_structure = i;        pointers(i) = obj.handles.text_gnd_structure;
          i = i+1;        id.popup_gnd_structure = i;       pointers(i) = obj.handles.popup_gnd_structure;          
          i = i+1;        id.panel_diel_profiles = i;       pointers(i) = obj.handles.panel_diel_profiles;
          i = i+1;        id.cb_discrete_slab = i;          pointers(i) = obj.handles.cb_discrete_slab;
          i = i+1;        id.cb_logistic_regression = i;    pointers(i) = obj.handles.cb_logistic_regression;
          i = i+1;        id.cb_2nd_order = i;              pointers(i) = obj.handles.cb_2nd_order;
          i = i+1;        id.cb_3rd_order = i;              pointers(i) = obj.handles.cb_3rd_order;
          i = i+1;        id.text_gnd_structure_single_layered = i;      pointers(i) = obj.handles.text_gnd_structure_single_layered;
          
            
          groupIDs.gnd_inputs = [id.panel_gnd_inputs id.text_diel_model : id.popup_gnd_structure];
          groupIDs.gnd_multi_layered_inputs = id.panel_diel_profiles : id.cb_3rd_order;
          
          
          %% INPUT FILES PANEL ELEMENTS
          i = i+1;        id.text_config_inputs_file = i;      pointers(i) = obj.handles.text_config_inputs_file;
          i = i+1;        id.edit_config_inputs_file = i;      pointers(i) = obj.handles.edit_config_inputs_file;
          i = i+1;        id.pb_config_inputs_file = i;        pointers(i) = obj.handles.pb_config_inputs_file;
          i = i+1;        id.text_ant_pat_Rx_file = i;      pointers(i) = obj.handles.text_ant_pat_Rx_file;
          i = i+1;        id.edit_ant_pat_Rx_file = i;      pointers(i) = obj.handles.edit_ant_pat_Rx_file;
          i = i+1;        id.pb_ant_pat_Rx_file = i;        pointers(i) = obj.handles.pb_ant_pat_Rx_file;
          i = i+1;        id.text_veg_inputs_file = i;      pointers(i) = obj.handles.text_veg_inputs_file;
          i = i+1;        id.edit_veg_inputs_file = i;      pointers(i) = obj.handles.edit_veg_inputs_file;
          i = i+1;        id.pb_veg_inputs_file = i;        pointers(i) = obj.handles.pb_veg_inputs_file;
            

          %% GUI ENABLE/DISABLE GROUPS                              
          % On gnd_cover change
          % Because there are currently only two ground covers, it can be 
          % adjusted by only on_popup_gnd_cover. If three or more items exist in 
          % the future, each should have its own enable/disable group) 
          groupIDs.on_popup_gnd_cover = [id.text_veg_inputs_file ...
                                   id.edit_veg_inputs_file ...
                                   id.pb_veg_inputs_file ];
                                
          % On Tx_orientation change
          % Because there are currently only two Tx orientations, it can be
          % adjusted by only on_Tx_orientation. If more items exist in the
          % future, each should have its own enable/disable group)             
          groupIDs.on_popup_orientation_Tx_geostationary = [id.text_el0_Tx ...
                                    id.edit_el0_Tx ...    
                                    id.text_deg_el0_Tx ...
                                    id.text_ph0_Tx ...
                                    id.edit_ph0_Tx ...
                                    id.text_deg_ph0_Tx];
          groupIDs.on_popup_orientation_Tx_variable = id.text_orientation_Tx_variable;
                                
          % On Rx_orientation change
          % Because there are currently only two Rx orientations, it can be
          % adjusted by only on_Rx_orientation. If more items exist in the
          % future, each should have its own enable/disable group)             
          groupIDs.on_popup_orientation_Rx_fixed = [id.text_th0_Rx ...
                                                id.edit_th0_Rx ...
                                                id.text_deg_th0_Rx ...
                                                id.text_ph0_Rx ...
                                                id.edit_ph0_Rx ...
                                                id.text_deg_ph0_Rx];
          groupIDs.on_popup_orientation_Rx_specular_facing = id.text_orientation_Rx_specular_facing;
          
          % On ant_pat_Rx change
          groupIDs.on_popup_ant_pat_Rx_GG = [id.text_ant_pat_res_Rx ...
                                             id.edit_ant_pat_res_Rx ...
                                             id.text_deg_ant_pat_res_Rx ...
                                             groupIDs.ant_pat_Rx_GG_inputs ];
                                      
          groupIDs.on_popup_ant_pat_Rx_user_defined = [id.text_ant_pat_Rx_file ...
                                                       id.edit_ant_pat_Rx_file ...
                                                       id.pb_ant_pat_Rx_file ];
          groupIDs.visibility_on_popup_ant_pat_Rx_user_defined = id.text_ant_pat_Rx_user_defined;
                            
          
          % On gnd_structure change        
          groupIDs.on_popup_gnd_structure_multi_layered = groupIDs.gnd_multi_layered_inputs;
          groupIDs.on_popup_gnd_structure_single_layered = id.text_gnd_structure_single_layered;
                                      
          
          groupIDs.input_files = [id.panel_input_files, ...
                               id.text_config_inputs_file : id.text_veg_inputs_file];         

            
          %% Initialize object properties
          [groupIDs.gPanels, groupIDs.strEl, groupIDs.valEl] = obj.autoElClassification(pointers);
            
          % Save in object
          obj.uiIDs = id;
          obj.uiGroups = groupIDs;
          obj.uiPointers = pointers(1:i);
            
          obj.curState = false(i,1);
          obj.newState = false(i,1);
          obj.setFlag = false(i,1);
          obj.getFlag = false(i,1);
            
        end
        
        
        % Function that runs inside onoffUIEl
        % Test every logical dependence in the GUI
        % E.g. a flag that activate other fields
        function checkUIdependencies(obj)              
            
            % gnd_cover
            if obj.is_popup_gnd_cover_bare_soil()

                  obj.setElStatus(obj.uiGroups.on_popup_gnd_cover, 0, 0);
                  
            end
              
            
            % popup_orientation_Tx
            if obj.is_popup_orientation_Tx_geostationary()

                  obj.setElStatus(obj.uiGroups.on_popup_orientation_Tx_geostationary, 1, 0);
                  obj.setElStatus(obj.uiGroups.on_popup_orientation_Tx_variable, 0, 0);
                  
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Tx_geostationary, 'on' );
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Tx_variable, 'off' );
                  
            elseif obj.is_popup_orientation_Tx_variable()

                  obj.setElStatus(obj.uiGroups.on_popup_orientation_Tx_geostationary, 0, 0);
                  obj.setElStatus(obj.uiGroups.on_popup_orientation_Tx_variable, 1, 0);
                  
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Tx_geostationary, 'off' );
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Tx_variable, 'on' );
                  
            end
              
            
            % popup_orientation_Rx
            if obj.is_popup_orientation_Rx_fixed()

                  obj.setElStatus(obj.uiGroups.on_popup_orientation_Rx_fixed, 1, 0);
                  obj.setElStatus(obj.uiGroups.on_popup_orientation_Rx_specular_facing, 0, 0);
                  
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Rx_fixed, 'on' );
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Rx_specular_facing, 'off' );
                  
            elseif obj.is_popup_orientation_Rx_specular_facing()

                  obj.setElStatus(obj.uiGroups.on_popup_orientation_Rx_fixed, 0, 0);
                  obj.setElStatus(obj.uiGroups.on_popup_orientation_Rx_specular_facing, 1, 0);
                  
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Rx_fixed, 'off' );
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_orientation_Rx_specular_facing, 'on' );
                  
            end
              
            
            % popup_ant_pat_Rx
            if obj.is_popup_ant_pat_Rx_user_defined()

                  obj.setElStatus([obj.uiGroups.on_popup_ant_pat_Rx_GG], 0, 0); 
                  obj.setElStatus([obj.uiGroups.on_popup_ant_pat_Rx_user_defined], 1, 0); 
                  
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_ant_pat_Rx_GG, 'off' );
                  obj.setGuiElVisibility( obj.uiGroups.visibility_on_popup_ant_pat_Rx_user_defined, 'on' ); 
                      
                  obj.setElVal( obj.uiIDs.edit_ant_pat_res_Rx, num2str(0) );
                  
            elseif obj.is_popup_ant_pat_Rx_GG()

                  obj.setElStatus([obj.uiGroups.on_popup_ant_pat_Rx_GG], 1, 0); 
                  obj.setElStatus([obj.uiGroups.on_popup_ant_pat_Rx_user_defined], 0, 0); 
                  
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_ant_pat_Rx_GG, 'on' );
                  obj.setGuiElVisibility( obj.uiGroups.visibility_on_popup_ant_pat_Rx_user_defined, 'off' ); 
            end
              
            
            % popup_gnd_structure
            if obj.is_popup_gnd_structure_single_layered()

                  obj.setElStatus(obj.uiGroups.on_popup_gnd_structure_single_layered, 1, 0);
                  obj.setElStatus(obj.uiGroups.on_popup_gnd_structure_multi_layered, 0, 0);
                  
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_gnd_structure_single_layered, 'on' );
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_gnd_structure_multi_layered, 'off' );
                  
            elseif obj.is_popup_gnd_structure_multi_layered()

                  obj.setElStatus(obj.uiGroups.on_popup_gnd_structure_single_layered, 0, 0);
                  obj.setElStatus(obj.uiGroups.on_popup_gnd_structure_multi_layered, 1, 0);
                  
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_gnd_structure_single_layered, 'off' );
                  obj.setGuiElVisibility( obj.uiGroups.on_popup_gnd_structure_multi_layered, 'on' );
                  
            end
            
        end
        
        
        
        % Function that runs inside syncFromGUI
        % Tests the validity of the value of every GUI element
        function isValid = checkUIvaluesValidity(obj, idEl)
           
            
            isValid = true;
            
            if ( sum(intersect(idEl, obj.uiIDs.pb_SCoBi)) > 0 ) || ...
               ( sum(intersect(idEl, obj.uiIDs.pb_save_inputs)) > 0 )
            
                % CAMPAIGN
                hObject =  obj.uiPointers( obj.uiIDs.edit_campaign );
                str = get( hObject, 'String' );

                % Check if it is a numerical value
                if isempty( str )
                    errordlg('Campaign cannot be empty!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return
                end


                % F_MHZ
                hObject = obj.uiPointers( obj.uiIDs.edit_f_MHz );
                str = get( hObject, 'String' );

                % Check if it is a numerical value
                if isempty( str )
                    errordlg('Frequency cannot be empty!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return

                elseif isempty( str2num(str) )
                    errordlg('Frequency must be numerical!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return

                % Check if it is a positive number
                elseif str2num(str) <= 0
                    errordlg('Frequency must be a positive number!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false; 
                    return
                end


                % R_TX_KM
                hObject = obj.uiPointers( obj.uiIDs.edit_r_Tx_km );
                str = get( hObject, 'String' );

                % Check if it is a numerical value
                if isempty( str )
                    errordlg('Range cannot be empty!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return

                elseif isempty( str2num(str) )
                    errordlg('Range must be numerical!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return

                % Check if it is a positive number
                elseif str2num(str) <= ( Constants.R_EARTH / Constants.KM_TO_M )
                    errordlg('Range must be higher than Earth''s radius (', num2str(Constants.R_EARTH / Constants.KM_TO_M), ' km)!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false; 
                    return
                end


                % EIRP_DB
                hObject = obj.uiPointers( obj.uiIDs.edit_EIRP_dB );
                str = get( hObject, 'String' );

                % Check if it is a numerical value
                if isempty( str )
                    errordlg('EIRP cannot be empty!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return

                elseif isempty( str2num(str) )
                    errordlg('EIRP must be numerical!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return

                % Check if it is a positive number
                elseif str2num(str) < 0
                    errordlg('EIRP cannot be negative!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false; 
                    return
                end


                % EL0_TX and PH0_TX (TRANSMITTER ELEVATION AND AZIMUTH ANGLES)          
                id_orientation_Tx = obj.getElVal(obj.uiIDs.popup_orientation_Tx);

                % It is required when the transmitter is a Geo-stationary
                % satellite
                if id_orientation_Tx == Constants.ID_TX_GEOSTATIONARY

                    % el0_Tx
                    hObject = obj.uiPointers( obj.uiIDs.edit_el0_Tx );
                    str = get( hObject, 'String' ); 

                    % Check if it is a numerical value
                    if isempty( str )
                        errordlg('Transmitter Elevation angle cannot be empty!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif isempty( str2num(str) )
                        errordlg('Transmitter Elevation angle must be numerical!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return
                    end

                    % ph0_Tx
                    hObject = obj.uiPointers( obj.uiIDs.edit_ph0_Tx );
                    str = get( hObject, 'String' ); 

                    % Check if it is a numerical value
                    if isempty( str )
                        errordlg('Transmitter Azimuth angle cannot be empty!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif isempty( str2num(str) )
                        errordlg('Transmitter Azimuth angle must be numerical!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return
                    end

                end


                % DIELECTRIC PROFILE FITTING FUNCTIONS           
                id_gnd_structure = obj.getElVal(obj.uiIDs.popup_gnd_structure);

                % It is required when the ground is Multi-layered
                if id_gnd_structure == Constants.ID_GND_MULTI_LAYERED

                    val2nd = obj.getElVal(obj.uiIDs.cb_2nd_order);
                    val3rd = obj.getElVal(obj.uiIDs.cb_3rd_order);
                    valLogistic = obj.getElVal(obj.uiIDs.cb_logistic_regression);
                    valDiscrete = obj.getElVal(obj.uiIDs.cb_discrete_slab);

                    % Check if it is a numerical value
                    if val2nd == 0 && val3rd == 0  && valLogistic == 0  && valDiscrete == 0 
                        errordlg('None of the Dielectric Profile fitting functions is selected. At least check one of them!', 'Invalid Input', 'modal');
                        hObject = obj.uiPointers( obj.uiIDs.cb_2nd_order );
                        uicontrol(hObject)
                        isValid = false;
                        return
                    end

                end


                % HR_M (RECEIVER ALTITUDE)
                hObject = obj.uiPointers( obj.uiIDs.edit_hr_m );
                str = get( hObject, 'String' );

                % Check if it is a numerical value
                if isempty( str )
                    errordlg('Receiver altitude cannot be empty!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return

                elseif isempty( str2num(str) )
                    errordlg('Receiver altitude must be numerical!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return

                % Check if it is a positive number
                elseif str2num(str) <= 0
                    errordlg('Receiver altitude must be positive!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false; 
                    return
                end


                % G0R_DB (RECEIVER GAIN)
                hObject = obj.uiPointers( obj.uiIDs.edit_G0r_dB );
                str = get( hObject, 'String' );

                % Check if it is a numerical value
                if isempty( str )
                    errordlg('Receiver gain cannot be empty!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return

                elseif isempty( str2num(str) )
                    errordlg('Receiver gain must be numerical!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return

                % Check if it is a positive number
                elseif str2num(str) < 0
                    errordlg('Receiver gain cannot be negative!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false; 
                    return
                end


                % TH0_RX and PH0_RX (RECEIVER ZENITH AND AZIMUTH OBSERVATION ANGLES)            
                id_orientation_Rx = obj.getElVal(obj.uiIDs.popup_orientation_Rx);

                % It is required when the receiver has a Fixed orientation
                if id_orientation_Rx == Constants.ID_RX_FIXED

                    % th0_Rx
                    hObject = obj.uiPointers( obj.uiIDs.edit_th0_Rx );
                    str = get( hObject, 'String' ); 

                    % Check if it is a numerical value
                    if isempty( str )
                        errordlg('Receiver''s zenith observation angle cannot be empty!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif isempty( str2num(str) )
                        errordlg('Receiver''s zenith observation angle must be numerical!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return
                    end

                    % ph0_Rx
                    hObject = obj.uiPointers( obj.uiIDs.edit_ph0_Rx );
                    str = get( hObject, 'String' ); 

                    % Check if it is a numerical value
                    if isempty( str )
                        errordlg('Receiver''s azimuth observation angle cannot be empty!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif isempty( str2num(str) )
                        errordlg('Receiver''s azimuth observation angle must be numerical!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return
                    end

                end


                % ANTENNA PATTERN GENERALIZED-GAUSSIAN PARAMETERS            
                id_ant_pat_Rx = obj.getElVal(obj.uiIDs.popup_ant_pat_Rx);

                % It is required when the receiver is chosen to have a 
                % Generalized-Gaussian pattern
                if id_ant_pat_Rx == Constants.ID_RX_GG

                    % Hpbw
                    hObject = obj.uiPointers( obj.uiIDs.edit_hpbw_deg );
                    str = get( hObject, 'String' ); 

                    % Check if it is a numerical value
                    if isempty( str )
                        errordlg('Receiver''s HPBW cannot be empty!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif isempty( str2num(str) )
                        errordlg('Receiver''s HPBW must be numerical!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif ( str2num(str) <= 0 ) || ( str2num(str) > Constants.RX_ANT_PAT_GG_HPBW_MAX )
                        errordlg( strcat( 'Receiver''s HPBW should be in the interval: (0, ', num2str(Constants.RX_ANT_PAT_GG_HPBW_MAX), ']'), 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return
                    end

                    % SLL_dB
                    hObject = obj.uiPointers( obj.uiIDs.edit_SLL_dB );
                    str = get( hObject, 'String' ); 

                    % Check if it is a numerical value
                    if isempty( str )
                        errordlg('Receiver''s side-lobe level cannot be empty!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif isempty( str2num(str) )
                        errordlg('Receiver''s side-lobe level must be numerical!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif ( str2num(str) < Constants.RX_ANT_PAT_GG_SLL_MIN ) || ( str2num(str) > Constants.RX_ANT_PAT_GG_SLL_MAX )
                        errordlg( strcat( 'Receiver''s side-lobe level should be in the interval: [', num2str(Constants.RX_ANT_PAT_GG_SLL_MIN), ',', num2str(Constants.RX_ANT_PAT_GG_SLL_MAX), ']'), 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return
                    end

                    % XPL_dB
                    hObject = obj.uiPointers( obj.uiIDs.edit_XPL_dB );
                    str = get( hObject, 'String' ); 

                    % Check if it is a numerical value
                    if isempty( str )
                        errordlg('Receiver''s cross-polarization level cannot be empty!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif isempty( str2num(str) )
                        errordlg('Receiver''s cross-polarization level must be numerical!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif ( str2num(str) < Constants.RX_ANT_PAT_GG_XPL_MIN ) || ( str2num(str) > Constants.RX_ANT_PAT_GG_XPL_MAX )
                        errordlg( strcat( 'Receiver''s cross-polarization level should be in the interval: [', num2str(Constants.RX_ANT_PAT_GG_XPL_MIN), ',', num2str(Constants.RX_ANT_PAT_GG_XPL_MAX), ']'), 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return
                    end

                    % ant_pat_res_Rx
                    hObject = obj.uiPointers( obj.uiIDs.edit_ant_pat_res_Rx );
                    str = get( hObject, 'String' ); 

                    % Check if it is a numerical value
                    if isempty( str )
                        errordlg('Receiver''s antenna pattern resolution cannot be empty!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif isempty( str2num(str) )
                        errordlg('Receiver''s antenna pattern resolution must be numerical!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return

                    elseif str2num(str) <= 0
                        errordlg('Receiver''s antenna pattern resolution must be positive!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return
                    end

                end


                % CONFIGURATION INPUTS FILE 
                hObject = obj.uiPointers( obj.uiIDs.edit_config_inputs_file );
                str = get( hObject, 'String' ); 
                
                [filepath, name, ext] = fileparts(str);

                % Check if it is empty
                if isempty( str )
                    errordlg('Configuration Inputs File cannot be empty!', 'Invalid Input', 'modal');
                    uicontrol(hObject)
                    isValid = false;
                    return
                    
                elseif ~isempty(filepath) && ~isempty(name) && ~isempty(ext)
                    if ~strcmp( ext, '.xls') && ~strcmp( ext, '.xlsx')
                        errordlg('There is a problem with Configuration Inputs File! It might be empty or might not be an Excel file.', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return
                    end
                end


                % ANTENNA PATTERN FILE   
                id_ant_pat_Rx = obj.getElVal(obj.uiIDs.popup_ant_pat_Rx);

                % It is required when the receiver is chosen to have a 
                % User-defined pattern
                if id_ant_pat_Rx == Constants.ID_RX_USER_DEFINED

                    hObject = obj.uiPointers( obj.uiIDs.edit_ant_pat_Rx_file );
                    str = get( hObject, 'String' ); 
                
                    [filepath, name, ext] = fileparts(str);

                    % Check if it is empty
                    if isempty( str )
                        errordlg('Antenna Pattern File cannot be empty!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return
                    
                    elseif ~isempty(filepath) && ~isempty(name) && ~isempty(ext)
                        if ~strcmp( ext, '.xls') && ~strcmp( ext, '.xlsx')
                            errordlg('There is a problem with Antenna Pattern File! It might be empty or might not be an Excel file.', 'Invalid Input', 'modal');
                            uicontrol(hObject)
                            isValid = false;
                            return
                        end
                    end

                end


                % VEGETATION INPUTS FILE  
                id_gnd_cover = obj.getElVal(obj.uiIDs.popup_gnd_cover);

                % It is required when the ground cover is chosen to be 
                % vegetation 
                if id_gnd_cover == Constants.ID_VEG_COVER

                    hObject = obj.uiPointers( obj.uiIDs.edit_veg_inputs_file );
                    str = get( hObject, 'String' ); 
                
                    [filepath, name, ext] = fileparts(str);

                    % Check if it is empty
                    if isempty( str )
                        errordlg('Vegetation Inputs File cannot be empty!', 'Invalid Input', 'modal');
                        uicontrol(hObject)
                        isValid = false;
                        return
                    
                    elseif ~isempty(filepath) && ~isempty(name) && ~isempty(ext)
                        if ~strcmp( ext, '.xls') && ~strcmp( ext, '.xlsx')
                            errordlg('There is a problem with Vegetation Inputs File! It might be empty or might not be an Excel file.', 'Invalid Input', 'modal');
                            uicontrol(hObject)
                            isValid = false;
                            return
                        end
                    end

                end
                
            end  % End of GUI element check
            
                        
        end     % End of function
        
        
        
        % Get enable / disable status of the element of the interface
        % This function should be called only once, later in the code
        % the status is kept updated
        function getAllElStatus(obj)            
            panels = false( length(obj.uiPointers), 1 );
            panels( obj.uiGroups.gPanels) = true;          % logical indexes of the panels
            
            idEl = 1 : length( obj.uiPointers );
            idEl = idEl(~panels);						% id of elements that are not panels
            idEl = idEl(obj.uiIDs.panel_sim_settings : end);	% The only elements to be considered starts 
            											% from the first panel sim_settings
            
            % For each panel
            for i = 1 : length(obj.uiGroups.gPanels)
            	obj.curState(obj.uiGroups.gPanels(i)) = obj.isGuiPanelOn( obj.uiPointers(obj.uiGroups.gPanels(i)) );
           	end
           	
           	% For all the other elements
           	for i=1:length(idEl)
            	obj.curState(idEl(i)) = obj.isGuiElOn(obj.uiPointers(idEl(i)));
            end
            
            obj.newState = obj.curState;
        end      
        
        
        function getSpecificElContent(obj) 
            
            idEl = 1 : length(obj.uiPointers);              % All the elements

            % Sets of panels
            panels = false(length(obj.uiPointers),1);     % init logical panels group
            panels(obj.uiGroups.gPanels) = true;           % set logical indexes of the panels
            % modified panels
            mPanel = idEl(panels & obj.getFlag);
                
            % Sets of text elements                
            textEl = false(length(obj.uiPointers),1);     % init logical text elements group
            textEl(obj.uiGroups.strEl) = true;              % set logical indexes of the text elements
                
            % Modified text elements
            mTextEl = idEl(textEl & obj.getFlag);

            mIdEl = idEl(~panels & ~textEl & obj.getFlag); % id of elements that are not panels nor text elements that have been modified
                
            % For each modified panel
            for i=1:length(mPanel)
                obj.curVal{mPanel(i)} = obj.getGuiElTitle(obj.uiPointers(mPanel(i)));
            end
                
            % For each modified txt Element
            for i=1:length(mTextEl)
                obj.curVal{mTextEl(i)} = obj.getGuiElStr(obj.uiPointers(mTextEl(i)));
            end

            % For all the other modified elements
            for i=1:length(mIdEl)
                obj.curVal{mIdEl(i)} = obj.getGuiElVal(obj.uiPointers(mIdEl(i)));
            end
            
        end
        
        
        % Initialize specific GUI elements to this class
        function initSpecificGUI(obj)
            
            
            %% GET GLOBAL DIRECTORIES
            dir_gui_last_input = Directories.getInstance.scobi_gui_last_input;
            dir_input_sys = Directories.getInstance.input_sys;
            
            
            % Read the last input file, if any
            lastInputFile = [];
            
            lastInputFileName = Constants.LAST_INPUT_FILENAMES{ 1, obj.simulator_id };
            defaultInputFileName = Constants.DEFAULT_INPUT_FILENAMES{ 1, obj.simulator_id };
            
            if exist( [strcat(dir_gui_last_input, '\') lastInputFileName], 'file' )
                
                filename = strcat( dir_gui_last_input, '\', lastInputFileName );
                
                load(filename);
                
                lastInputFile = lastInput.lastInputFileName;
                
            end
            
            % If there is an input file from the last execution of SCoBi
            if ~isempty( lastInputFile ) && exist( lastInputFile,'file' )
                
                obj.inputStruct = obj.loadGUIFromInputFile( [], lastInputFile );
                
            % Else if the default input exists, load it
            elseif exist([strcat(dir_input_sys, '/') defaultInputFileName],'file')
                
                obj.inputStruct = obj.loadGUIFromInputFile( defaultInputFileName );
                
            % Else, make a warning about it and open an empty GUI
            else
                
                waitfor(msgbox('WARNING: There is no default input, you should add input files or give inputs manually!'));
                
                obj.loadEmptyGUI();
                
            end
            
        end
        
        
        % Load the state of the gui from a matlab file.
        function loadEmptyGUI( obj )

            
            %%   SIMULATION SETTINGS
            obj.setElVal(obj.uiIDs.edit_campaign, '', 0);
            
            % Simulation mode: Snapshot OR Time-series
            obj.init_popup_sim_mode();
            obj.setElVal( obj.uiIDs.popup_sim_mode, Constants.ID_SNAPSHOT, 0 );
            
            % Ground cover: Bare-soil OR Vegetation
            obj.init_popup_gnd_cover();
            obj.setElVal(obj.uiIDs.popup_gnd_cover, Constants.ID_VEG_COVER, 0);
            
            % Flag to write Attenuation to Excel file
            obj.setElVal(obj.uiIDs.cb_write_attenuation, 0, 0);
            
            % Flag to include the simulation in the Master Simulation file
            obj.setElVal(obj.uiIDs.cb_include_in_master_sim_file, 0, 0);
            
            % Flag to calculate and save penetration depth (multilayer)
            obj.setElVal(obj.uiIDs.cb_calculate_penetration_depth, 0, 0);
            
            
            %% TRANSMITTER (Tx) INPUTS
            obj.setElVal(obj.uiIDs.edit_f_MHz, [], 0);
            
            obj.setElVal(obj.uiIDs.edit_r_Tx_km, [], 0);
            
            obj.setElVal(obj.uiIDs.edit_EIRP_dB, [], 0);
            
            obj.init_popup_pol_Tx();
            obj.setElVal(obj.uiIDs.popup_pol_Tx, Constants.ID_POL_R, 0);
            
            obj.init_popup_orientation_Tx();
            obj.setElVal(obj.uiIDs.popup_orientation_Tx, Constants.ID_TX_GEOSTATIONARY, 0);
                
            obj.setElVal(obj.uiIDs.edit_el0_Tx, [], 0);

            obj.setElVal(obj.uiIDs.edit_ph0_Tx, [], 0);
            
            
            %% RECEIVER (Rx) INPUTS
            obj.setElVal(obj.uiIDs.edit_hr_m, [], 0);
            
            obj.setElVal(obj.uiIDs.edit_G0r_dB, [], 0);
            
            obj.init_popup_pol_Rx();
            obj.setElVal(obj.uiIDs.popup_pol_Rx, Constants.ID_POL_R, 0);
            
            obj.init_popup_orientation_Rx();
            obj.setElVal(obj.uiIDs.popup_orientation_Rx, Constants.ID_RX_FIXED, 0);
                
            obj.setElVal(obj.uiIDs.edit_th0_Rx, [], 0);

            obj.setElVal(obj.uiIDs.edit_ph0_Rx, [], 0);
            
            obj.init_popup_ant_pat_Rx();
            obj.setElVal(obj.uiIDs.popup_ant_pat_Rx, Constants.ID_RX_GG, 0);
            
            obj.setElVal(obj.uiIDs.edit_ant_pat_res_Rx, [], 0);

            obj.setElVal(obj.uiIDs.edit_hpbw_deg, [], 0);

            obj.setElVal(obj.uiIDs.edit_SLL_dB, [], 0);

            obj.setElVal(obj.uiIDs.edit_XPL_dB, [], 0);
            
            
            %% GROUND INPUTS
            obj.init_popup_diel_model();
            obj.setElVal(obj.uiIDs.popup_diel_model, Constants.ID_DIEL_DOBSON, 0);
            
            obj.init_popup_gnd_structure();
            obj.setElVal(obj.uiIDs.popup_gnd_structure, Constants.ID_GND_SINGLE_LAYERED, 0);            
            
            % Flag to calculate Discrete-slab for dielectric profile 
            obj.setElVal(obj.uiIDs.cb_discrete_slab, 0, 0);            
            
            % Flag to calculate Logistic regression for dielectric profile 
            obj.setElVal(obj.uiIDs.cb_logistic_regression, 0, 0);            
            
            % Flag to calculate 2nd-order poly-fit for dielectric profile 
            obj.setElVal(obj.uiIDs.cb_2nd_order, 0, 0);            
            
            % Flag to calculate 3rd order poly-fit for dielectric profile 
            obj.setElVal(obj.uiIDs.cb_3rd_order, 0, 0);
            
            %% INPUT FILES  
            obj.setElVal(obj.uiIDs.edit_config_inputs_file, '', 0);     
            
            obj.setElVal(obj.uiIDs.edit_ant_pat_Rx_file, '', 0);
            
            obj.setElVal(obj.uiIDs.edit_veg_inputs_file, '', 1);

            obj.updateGUI();
            
            %% LOAD GUI
            % Check all the dependencies
            obj.syncFromGUI(obj.uiIDs.popup_sim_mode);
            
        end
        
        
        % Load the state of the gui from a matlab file.
        function [inputStruct] = loadGUIFromInputFile( obj, filename, pathname )
            
            
            %% GET GLOBAL DIRECTORIES
            dir_input_sys = Directories.getInstance.input_sys;  
            dir_input_ant_pat_Rx = Directories.getInstance.input_ant_pat_Rx;
            dir_input_config = Directories.getInstance.input_config;   
            dir_input_veg = Directories.getInstance.input_veg;               
            
            
            %% LOAD THE VARIABLE "inputStruct" FROM INPUT FILE
            isDefaultInput = 0;
            
            % If the default input file is used
            if nargin == 2
                
                pathname = dir_input_sys;
                
                isDefaultInput = 1;
                
            end
            
            % TO-DO: Excel Input Files filenames are given with relative 
            % path in default inputs since every user's source code can
            % have a different absolute path. However this is a problem if
            % the user plays with the default inputs and save the changes
            % into one of the default input files rather than creating a
            % new input file.
            
            % If filename is empty, then it means this function is called
            % by ...
            if ~isempty( filename )
                
                % Load the input file
                filename = strcat( pathname, '\',filename );
                load( filename ); % the file contains the variable inputStruct
            
            % If filename is empty, then it means this function is called 
            % by ... and pathname contains all    
            else
                
                load( pathname );
                
            end

            % TO-DO: Excel Input Files filenames are given with relative 
            % path in default inputs since every user's source code can
            % have a different absolute path. However this is a problem if
            % the user plays with the default inputs and save the changes
            % into one of the default input files rather than creating a
            % new input file.
            
            % If the default input is used, make the full input path 
            if isDefaultInput
                
                inputStruct.config_inputs_file = strcat(dir_input_config, '\', inputStruct.config_inputs_file);
                
                if strcmp( inputStruct.ant_pat_Rx, Constants.RX_ANT_PATTERNS{1, Constants.ID_RX_USER_DEFINED } ) 
                
                    inputStruct.ant_pat_Rx_file = strcat(dir_input_ant_pat_Rx, '\', inputStruct.ant_pat_Rx_file);
                    
                end
                
                if strcmp( inputStruct.gnd_cover, Constants.GND_COVERS{1, Constants.ID_VEG_COVER } ) 
                
                    inputStruct.veg_inputs_file = strcat(dir_input_veg, '\', inputStruct.veg_inputs_file);
                    
                end
                
            end
            
            
            %% SET GUI ELEMENTS' VALUES FROM inputStruct
            %%   SIMULATION SETTINGS        
            obj.setElVal(obj.uiIDs.edit_campaign, inputStruct.campaign, 0);
            
            % Simulation Mode: Snapshot OR Time-series
            obj.init_popup_sim_mode();
            sim_mode_id = findElementIdInCell( Constants.SIM_MODES, inputStruct.sim_mode );
            obj.setElVal( obj.uiIDs.popup_sim_mode, sim_mode_id, 0 );
            
            % Ground cover: Bare-soil OR Vegetation
            obj.init_popup_gnd_cover();
            gnd_cover_id = findElementIdInCell( Constants.GND_COVERS, inputStruct.gnd_cover );
            obj.setElVal(obj.uiIDs.popup_gnd_cover, gnd_cover_id, 0);
            
            % Flag to write Attenuation to Excel file
            obj.setElVal(obj.uiIDs.cb_write_attenuation, inputStruct.write_attenuation, 0);
            
            % Flag to include the simulation in the Master Simulation file
            obj.setElVal(obj.uiIDs.cb_include_in_master_sim_file, inputStruct.include_in_master_sim_file, 0);
            
            % Flag to calculate penetration depth (multilayer)
            obj.setElVal(obj.uiIDs.cb_calculate_penetration_depth, inputStruct.calculate_penetration_depth, 0);
            
            
            %% TRANSMITTER (Tx) INPUTS
            obj.setElVal(obj.uiIDs.edit_f_MHz, num2str(inputStruct.f_MHz), 0);
            
            obj.setElVal(obj.uiIDs.edit_r_Tx_km, num2str(inputStruct.r_Tx_km), 0);
            
            obj.setElVal(obj.uiIDs.edit_EIRP_dB, num2str(inputStruct.EIRP_dB), 0);
            
            obj.init_popup_pol_Tx();
            pol_Tx_id = findElementIdInCell( Constants.POLARIZATIONS, inputStruct.pol_Tx );
            obj.setElVal(obj.uiIDs.popup_pol_Tx, pol_Tx_id, 0);
            
            obj.init_popup_orientation_Tx();
            orientation_Tx_id = findElementIdInCell( Constants.TX_ORIENTATIONS, inputStruct.orientation_Tx );
            obj.setElVal(obj.uiIDs.popup_orientation_Tx, orientation_Tx_id, 0);
            
            % If transmitter orientation is Geo-stationary, then load 
            % incidence and azimuth agles
            if orientation_Tx_id == Constants.ID_TX_GEOSTATIONARY
                
                obj.setElVal(obj.uiIDs.edit_el0_Tx, num2str(inputStruct.el0_Tx_deg), 0);
            
                obj.setElVal(obj.uiIDs.edit_ph0_Tx, num2str(inputStruct.ph0_Tx_deg), 0);
                
            end
            
            
            %% RECEIVER (Rx) INPUTS
            obj.setElVal(obj.uiIDs.edit_hr_m, num2str(inputStruct.hr_m), 0);
            
            obj.setElVal(obj.uiIDs.edit_G0r_dB, num2str(inputStruct.G0r_dB), 0);
            
            obj.init_popup_pol_Rx();
            pol_Rx_id = findElementIdInCell( Constants.POLARIZATIONS, inputStruct.pol_Rx );
            obj.setElVal(obj.uiIDs.popup_pol_Rx, pol_Rx_id, 0);
            
            obj.init_popup_orientation_Rx();
            orientation_Rx_id = findElementIdInCell( Constants.RX_ORIENTATIONS, inputStruct.orientation_Rx );
            obj.setElVal(obj.uiIDs.popup_orientation_Rx, orientation_Rx_id, 0);
            
            % If receiver orientation is fixed, then load incidence and
            % azimuth agles
            if orientation_Rx_id == Constants.ID_RX_FIXED
                
                obj.setElVal(obj.uiIDs.edit_th0_Rx, num2str(inputStruct.th0_Rx_deg), 0);
            
                obj.setElVal(obj.uiIDs.edit_ph0_Rx, num2str(inputStruct.ph0_Rx_deg), 0);
                
            end
            
            obj.init_popup_ant_pat_Rx();
            ant_pat_Rx_id = findElementIdInCell( Constants.RX_ANT_PATTERNS, inputStruct.ant_pat_Rx );
            obj.setElVal(obj.uiIDs.popup_ant_pat_Rx, ant_pat_Rx_id, 0);
            
            % If receiver antenna pattern is Generalized-Gaussian, then 
            % load the corresponding inputs
            if ant_pat_Rx_id == Constants.ID_RX_GG
            
                obj.setElVal(obj.uiIDs.edit_ant_pat_res_Rx, num2str(inputStruct.ant_pat_res_deg_Rx), 0);
            
                obj.setElVal(obj.uiIDs.edit_hpbw_deg, num2str(inputStruct.hpbw_deg), 0);
                
                obj.setElVal(obj.uiIDs.edit_SLL_dB, num2str(inputStruct.SLL_dB), 0);
                
                obj.setElVal(obj.uiIDs.edit_XPL_dB, num2str(inputStruct.XPL_dB), 0);
                
            end
            
            
            %% GROUND INPUTS
            obj.init_popup_diel_model();
            diel_model_id = findElementIdInCell( Constants.DIEL_MODELS, inputStruct.diel_model );
            obj.setElVal(obj.uiIDs.popup_diel_model, diel_model_id, 0);
            
            obj.init_popup_gnd_structure();
            gnd_structure_id = findElementIdInCell( Constants.GND_STRUCTURES, inputStruct.gnd_structure );
            obj.setElVal(obj.uiIDs.popup_gnd_structure, gnd_structure_id, 0);            
            
            if gnd_structure_id == Constants.ID_GND_MULTI_LAYERED
                
                % Flag to calculate Discrete-slab for dielectric profile 
                obj.setElVal(obj.uiIDs.cb_discrete_slab, inputStruct.calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_SLAB,1), 0);            

                % Flag to calculate Logistic regression for dielectric profile 
                obj.setElVal(obj.uiIDs.cb_logistic_regression, inputStruct.calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_LOGISTIC,1), 0);            

                % Flag to calculate 2nd-order poly-fit for dielectric profile 
                obj.setElVal(obj.uiIDs.cb_2nd_order, inputStruct.calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_2ND_ORDER,1), 0);            

                % Flag to calculate 3rd order poly-fit for dielectric profile 
                obj.setElVal(obj.uiIDs.cb_3rd_order, inputStruct.calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_3RD_ORDER,1), 0);
                
            end
            
            
            %% INPUT FILES  
            obj.setElVal(obj.uiIDs.edit_config_inputs_file, inputStruct.config_inputs_file, 0);     
            
            if ant_pat_Rx_id == Constants.ID_RX_USER_DEFINED
            
                obj.setElVal(obj.uiIDs.edit_ant_pat_Rx_file, inputStruct.ant_pat_Rx_file, 0);
                
            end
            
            if gnd_cover_id == Constants.ID_VEG_COVER
            
                obj.setElVal(obj.uiIDs.edit_veg_inputs_file, inputStruct.veg_inputs_file, 1);
                
            end

            obj.updateGUI();
            
            %% LOAD GUI
            % Check all the dependencies
            obj.syncFromGUI(obj.uiIDs.popup_sim_mode);
            
        end
        
        % Save the state of the GUI to a matlab file
        function savingResult = saveGUIToInputFile( obj, saveOption )
        % savingResult: 0 - Didn't save. Not good to go
        % savingResult: 1 - Saved. Good to go
        % savingResult: 2 - Didn't save, but good to go
        
        
        %% GET GLOBAL DIRECTORIES
        dir_input_sys = Directories.getInstance.input_sys;
     
            
        %% SIMULATION SETTINGS PANEL ELEMENTS
        inputStruct.campaign        = obj.getElVal(obj.uiIDs.edit_campaign);
        
        sim_mode_id                     = obj.getElVal(obj.uiIDs.popup_sim_mode);            
        inputStruct.sim_mode            = Constants.SIM_MODES{ 1, sim_mode_id };

        gnd_cover_id                    = obj.getElVal(obj.uiIDs.popup_gnd_cover);            
        inputStruct.gnd_cover           = Constants.GND_COVERS{ 1, gnd_cover_id };

        inputStruct.write_attenuation	= obj.getElVal(obj.uiIDs.cb_write_attenuation);

        inputStruct.include_in_master_sim_file	= obj.getElVal(obj.uiIDs.cb_include_in_master_sim_file);
        
        inputStruct.calculate_penetration_depth	= obj.getElVal(obj.uiIDs.cb_calculate_penetration_depth);


        %% TRANSMITTER (Tx) INPUTS PANEL ELEMENTS  
        inputStruct.f_MHz     = str2double(obj.getElVal(obj.uiIDs.edit_f_MHz));

        inputStruct.r_Tx_km	= str2double(obj.getElVal(obj.uiIDs.edit_r_Tx_km));

        inputStruct.EIRP_dB	= str2double(obj.getElVal(obj.uiIDs.edit_EIRP_dB));

        pol_Tx_id       = obj.getElVal(obj.uiIDs.popup_pol_Tx);
        inputStruct.pol_Tx	= Constants.POLARIZATIONS{ 1, pol_Tx_id };

        orientation_Tx_id       = obj.getElVal(obj.uiIDs.popup_orientation_Tx);
        inputStruct.orientation_Tx	= Constants.TX_ORIENTATIONS{ 1, orientation_Tx_id };

        % If transmitter orientation is Geo-stationary, then assign 
        % incidence and azimuth angles
        if orientation_Tx_id == Constants.ID_TX_GEOSTATIONARY
            
            inputStruct.el0_Tx_deg	= str2double(obj.getElVal(obj.uiIDs.edit_el0_Tx));
            inputStruct.ph0_Tx_deg	= str2double(obj.getElVal(obj.uiIDs.edit_ph0_Tx));
            
        end


        %% RECEIVER (Rx) INPUTS PANEL ELEMENTS 
        inputStruct.hr_m      = str2double(obj.getElVal(obj.uiIDs.edit_hr_m));

        inputStruct.G0r_dB	= str2double(obj.getElVal(obj.uiIDs.edit_G0r_dB));

        pol_Rx_id       = obj.getElVal(obj.uiIDs.popup_pol_Rx);
        inputStruct.pol_Rx	= Constants.POLARIZATIONS{ 1, pol_Rx_id };

        orientation_Rx_id       = obj.getElVal(obj.uiIDs.popup_orientation_Rx);
        inputStruct.orientation_Rx	= Constants.RX_ORIENTATIONS{ 1, orientation_Rx_id };

        % If receiver orientation is fixed, then assign incidence and
        % azimuth angles
        if orientation_Rx_id == Constants.ID_RX_FIXED
            
            inputStruct.th0_Rx_deg	= str2double(obj.getElVal(obj.uiIDs.edit_th0_Rx));
            inputStruct.ph0_Rx_deg	= str2double(obj.getElVal(obj.uiIDs.edit_ph0_Rx));
            
        end

        ant_pat_Rx_id       = obj.getElVal(obj.uiIDs.popup_ant_pat_Rx);
        inputStruct.ant_pat_Rx	= Constants.RX_ANT_PATTERNS{ 1, ant_pat_Rx_id };

        % If receiver antenna pattern is Generalized-Gaussian, then 
        % assign GG params
        if ant_pat_Rx_id == Constants.ID_RX_GG

            inputStruct.hpbw_deg = str2double( obj.getElVal(obj.uiIDs.edit_hpbw_deg) );
            inputStruct.SLL_dB	= str2double( obj.getElVal(obj.uiIDs.edit_SLL_dB) );
            inputStruct.XPL_dB	= str2double( obj.getElVal(obj.uiIDs.edit_XPL_dB) );
            inputStruct.ant_pat_res_deg_Rx	= str2double(obj.getElVal(obj.uiIDs.edit_ant_pat_res_Rx));

        end


        %% GROUND INPUTS PANEL ELEMENTS
        diel_model_id       = obj.getElVal(obj.uiIDs.popup_diel_model);
        inputStruct.diel_model	= Constants.DIEL_MODELS{ 1, diel_model_id };

        gnd_structure_id            = obj.getElVal(obj.uiIDs.popup_gnd_structure);
        inputStruct.gnd_structure	= Constants.GND_STRUCTURES{ 1, gnd_structure_id };

        % If receiver antenna pattern is Generalized-Gaussian, then 
        % assign GG params
        if gnd_structure_id == Constants.ID_GND_MULTI_LAYERED

            calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_SLAB,1) = obj.getElVal(obj.uiIDs.cb_discrete_slab);
            calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_LOGISTIC,1) = obj.getElVal(obj.uiIDs.cb_logistic_regression);
            calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_2ND_ORDER,1) = obj.getElVal(obj.uiIDs.cb_2nd_order);
            calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_3RD_ORDER,1) = obj.getElVal(obj.uiIDs.cb_3rd_order);
            
            inputStruct.calc_diel_profile_fit_functions = calc_diel_profile_fit_functions;

        end


        %% INPUT FILES PANEL ELEMENTS
        inputStruct.config_inputs_file	= obj.getElVal(obj.uiIDs.edit_config_inputs_file);

        if ant_pat_Rx_id == Constants.ID_RX_USER_DEFINED

            inputStruct.ant_pat_Rx_file	= obj.getElVal(obj.uiIDs.edit_ant_pat_Rx_file);
            
        end

        if gnd_cover_id == Constants.ID_VEG_COVER
        
            inputStruct.veg_inputs_file = obj.getElVal(obj.uiIDs.edit_veg_inputs_file); 
            
        end


        %% SAVE ALL TO FILE
        % If save option is "save" AND inputStructs should be different, OR
        % save option is "save as", then open a saving dialog
        if ( saveOption == Constants.ID_GUI_SAVE && ~isequaln( obj.inputStruct, inputStruct ) ) ...
            || saveOption == Constants.ID_GUI_SAVE_AS

            filter = {strcat(dir_input_sys, '\*.mat')};
            [file, path] = uiputfile(filter, 'Save Input File');

            % If file and path are not empty
            if length(path)>1 && length(file)>1

                savingResult = 1;

                filename = strcat( path, file );

                save(filename, 'inputStruct');

                obj.inputStruct = inputStruct;

                obj.saveLastInputFile(  path, file );

            else

                % No valid save file info. Do not save and not good to go
                savingResult = 0;

            end

        else
            
            % No change on the input file; so, do not save but good to go.
            savingResult = 2;

        end
            
        end
        
        
        function saveLastInputFile(obj, path, file)
            
            
            %% GET GLOBAL DIRECTORIES
            dir_gui_last_input = Directories.getInstance.scobi_gui_last_input;
            

            lastInputFileName = Constants.LAST_INPUT_FILENAMES{ 1, obj.simulator_id };
            
            filename = strcat( dir_gui_last_input, '/', lastInputFileName );
            
            if ~exist(dir_gui_last_input, 'dir')
                mkdir(dir_gui_last_input)
            end
            
            lastInput.lastInputFileName = strcat( path, file );
            
            save(filename, 'lastInput');

        end
        
        
        function setSpecificElContent(obj)
                
            idEl = 1:length(obj.uiPointers);              % All the elements

            % Sets of panels
            panels = false(length(obj.uiPointers),1);     % init logical panels group
            panels(obj.uiGroups.gPanels) = true;           % set logical indexes of the panels
            % modified panels
            mPanel = idEl(panels & obj.setFlag);         
                
            % Sets of text elements
            textEl = false(length(obj.uiPointers),1);     % init logical text elements group
            textEl(obj.uiGroups.strEl) = true;              % set logical indexes of the text elements
            % Modified text elements
            mTextEl = idEl(textEl & obj.setFlag);

            mIdEl = idEl(~panels & ~textEl & obj.setFlag); % id of elements that are not panels nor text elements that have been modified
                
            % For each modified panel
            for i=1:length(mPanel)
                obj.setGuiElTitle(obj.uiPointers(mPanel(i)), obj.newVal{mPanel(i)});
            end
                
            % For each modified txt element
            for i=1:length(mTextEl)
                obj.setGuiElStr(obj.uiPointers(mTextEl(i)), obj.newVal{mTextEl(i)});
            end

            % For all the other modified elements
            for i=1:length(mIdEl)
                obj.setGuiElVal(obj.uiPointers(mIdEl(i)), obj.newVal{mIdEl(i)});
            end
            
        end 
        
    end
    
end

