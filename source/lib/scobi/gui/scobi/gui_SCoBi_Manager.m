% TO-DO: Check comments, copyrights, etc.
classdef gui_SCoBi_Manager < SCoBiGUIManagers
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
            result = isOn && (obj.getElVal(obj.uiIDs.popup_ant_pat_Rx) == Constants.id_Rx_cos_pow_n );
            
        end


        % Check if ant_pat_Rx value is Generalized-Gaussian
        function result = is_popup_ant_pat_Rx_GG(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_ant_pat_Rx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_ant_pat_Rx) == Constants.id_Rx_GG );
            
        end


        % Check if ant_pat_Rx value is User-define
        function result = is_popup_ant_pat_Rx_user_defined(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_ant_pat_Rx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_ant_pat_Rx) == Constants.id_Rx_user_defined );
            
        end


        % Check if diel_model value is Dobson
        function result = is_popup_diel_model_dobson(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_diel_model );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_diel_model) == Constants.id_diel_dobson );
            
        end


        % Check if diel_model value is Mironov
        function result = is_popup_diel_model_mironov(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_diel_model );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_diel_model) == Constants.id_diel_mironov );
            
        end


        % Check if diel_model value is Wang
        function result = is_popup_diel_model_wang(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_diel_model );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_diel_model) == Constants.id_diel_wang );
            
        end


        % Check if gnd_cover value is Bare-soil
        function result = is_popup_gnd_cover_bare_soil(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_gnd_cover );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_gnd_cover) == Constants.id_bare_soil );
            
        end


        % Check if gnd_cover value is Vegetation
        function result = is_popup_gnd_cover_vegetation(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_gnd_cover );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_gnd_cover) == Constants.id_veg_cover );
            
        end


        % Check if gnd_structure value is Single-layered
        function result = is_popup_gnd_structure_single_layered(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_gnd_structure );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_gnd_structure) == Constants.id_gnd_single_layered );
            
        end


        % Check if gnd_structure value is Multi-layered
        function result = is_popup_gnd_structure_multi_layered(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_gnd_structure );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_gnd_structure) == Constants.id_gnd_multi_layered );
            
        end


        % Check if orientation_Rx value is Fixed
        function result = is_popup_orientation_Rx_fixed(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_orientation_Rx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_orientation_Rx) == Constants.id_Rx_fixed );
            
        end


        % Check if orientation_Rx value is Specular-facing
        function result = is_popup_orientation_Rx_specular_facing(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_orientation_Rx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_orientation_Rx) == Constants.id_Rx_specular_facing );
            
        end


        % Check if orientation_Tx value is Fixed
        function result = is_popup_orientation_Tx_geostationary(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_orientation_Tx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_orientation_Tx) == Constants.id_Tx_geostationary );
            
        end


        % Check if orientation_Tx value is Specular-facing
        function result = is_popup_orientation_Tx_variable(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_orientation_Tx );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_orientation_Tx) == Constants.id_Tx_variable );
            
        end


        % Check if sim_mode value is Snapshot
        function result = is_popup_sim_mode_snapshot(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_sim_mode );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_sim_mode) == Constants.id_snapshot );
            
        end


        % Check if sim_mode value is Time-series
        function result = is_popup_sim_mode_time_series(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_sim_mode );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_sim_mode) == Constants.id_time_series );
            
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

              else

                  % TO-DO: Handle Exception?

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

              else

                  % TO-DO: Handle Exception?

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

              else

                  % TO-DO: Handle Exception?

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
                  obj.setElVal(obj.uiIDs.popup_ant_pat_Rx, Constants.id_Rx_GG); 
        
                  obj.syncFromGUI(obj.uiIDs.popup_ant_pat_Rx);

              else

                  % TO-DO: Handle Exception?

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

              else

                  % TO-DO: Handle Exception?

              end
              
              obj.updateGUI();
              
          end
          
            
          %% PUSH BUTTONS               
            % on Forest button
            if sum(intersect(idEl, obj.uiIDs.pb_Forest)) > 0
                    
                obj.init( obj.handles, Constants.id_sim_forest );
                
            end
            
            % on Agriculture button
            if sum(intersect(idEl, obj.uiIDs.pb_Agriculture)) > 0
                    
                obj.init( obj.handles, Constants.id_sim_agriculture );
                
            end
            
            % on Root-zone button
            if sum(intersect(idEl, obj.uiIDs.pb_Root_zone)) > 0
                    
                obj.init( obj.handles, Constants.id_sim_root_zone );
                
            end
            
            % on Soil button
            if sum(intersect(idEl, obj.uiIDs.pb_Soil)) > 0
                    
                obj.init( obj.handles, Constants.id_sim_soil );
                
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
                    
                obj.saveGUIToInputFile( Constants.id_GUI_save_as );
                
            end
            
            % on SCoBi
            if sum(intersect(idEl, obj.uiIDs.pb_SCoBi)) > 0

                % Try to save GUI values to an input file first
                % It will perform a check for any changes to the loaded or recently 
                % saved inputs   
                savingResult = obj.saveGUIToInputFile( Constants.id_GUI_save );

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
            
        end
        
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
                str = Constants.Rx_ant_pats;
            end
            
            value = get(obj.handles.popup_ant_pat_Rx,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_ant_pat_Rx,'Value', value);
            set(obj.handles.popup_ant_pat_Rx,'String', str);
        end

        
        % Fill the popup_diel_model
        function init_popup_diel_model(obj, str)
            
            if nargin < 2
                str = Constants.diel_models;
            end
            
            value = get(obj.handles.popup_diel_model,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_diel_model,'Value', value);
            set(obj.handles.popup_diel_model,'String', str);
        end

        
        % Fill the popup_gnd_cover
        function init_popup_gnd_cover(obj, str)
            
            if nargin < 2
                str = Constants.gnd_covers;
            end
            
            value = get(obj.handles.popup_gnd_cover,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_gnd_cover,'Value', value);
            set(obj.handles.popup_gnd_cover,'String', str);
        end

        
        % Fill the popup_gnd_structure
        function init_popup_gnd_structure(obj, str)
            
            if nargin < 2
                str = Constants.gnd_structures;
            end
            
            value = get(obj.handles.popup_gnd_structure,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_gnd_structure,'Value', value);
            set(obj.handles.popup_gnd_structure,'String', str);
        end

        
        % Fill the popup_orientation_Tx
        function init_popup_orientation_Tx(obj, str)
            
            if nargin < 2
                str = Constants.Tx_orientations;
            end
            
            value = get(obj.handles.popup_orientation_Tx,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_orientation_Tx,'Value', value);
            set(obj.handles.popup_orientation_Tx,'String', str);
        end

        
        % Fill the popup_orientation_Rx
        function init_popup_orientation_Rx(obj, str)
            
            if nargin < 2
                str = Constants.Rx_orientations;
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
                str = Constants.polarizations(1, Constants.id_pol_R : Constants.id_pol_Y);
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
                str = Constants.polarizations(1, Constants.id_pol_R : Constants.id_pol_Y);
            end
            
            value = get(obj.handles.popup_pol_Tx,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_pol_Tx,'Value', value);
            set(obj.handles.popup_pol_Tx,'String', str);
        end

        
        % Fill the popup_sim_mode
        function init_popup_sim_mode(obj, str)
            
            if nargin < 2
                str = Constants.sim_modes;
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
          i = i+1;        id.text_th0_Tx = i;               pointers(i) = obj.handles.text_th0_Tx;
          i = i+1;        id.edit_th0_Tx = i;               pointers(i) = obj.handles.edit_th0_Tx;
          i = i+1;        id.text_deg_th0_Tx = i;           pointers(i) = obj.handles.text_deg_th0_Tx;
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
          groupIDs.on_popup_orientation_Tx_geostationary = [id.text_th0_Tx ...
                                    id.edit_th0_Tx ...    
                                    id.text_deg_th0_Tx ...
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
            
            lastInputFileName = Constants.lastInputFileNames{ 1, obj.simulator_id };
            defaultInputFileName = Constants.defaultInputFileNames{ 1, obj.simulator_id };
            
            if exist( [strcat(dir_gui_last_input, '\') lastInputFileName], 'file' )
                
                filename = strcat( dir_gui_last_input, '\', lastInputFileName );
                
                load(filename);
                
                lastInputFile = lastInput.lastInputFileName;
                % ConstantNames.lastInputFile contains the lastInputFileName 
                % string that might be saved during the last SCoBi run
                
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
            obj.setElVal( obj.uiIDs.popup_sim_mode, Constants.id_snapshot, 0 );
            
            % Ground cover: Bare-soil OR Vegetation
            obj.init_popup_gnd_cover();
            obj.setElVal(obj.uiIDs.popup_gnd_cover, Constants.id_veg_cover, 0);
            
            % Flag to write Attenuation to Excel file
            obj.setElVal(obj.uiIDs.cb_write_attenuation, 0, 0);
            
            % Flag to include the simulation in the Master Simulation file
            obj.setElVal(obj.uiIDs.cb_include_in_master_sim_file, 0, 0);
            
            
            %% TRANSMITTER (Tx) INPUTS
            obj.setElVal(obj.uiIDs.edit_f_MHz, [], 0);
            
            obj.setElVal(obj.uiIDs.edit_r_Tx_km, [], 0);
            
            obj.setElVal(obj.uiIDs.edit_EIRP_dB, [], 0);
            
            obj.init_popup_pol_Tx();
            obj.setElVal(obj.uiIDs.popup_pol_Tx, Constants.id_pol_R, 0);
            
            obj.init_popup_orientation_Tx();
            obj.setElVal(obj.uiIDs.popup_orientation_Tx, Constants.id_Tx_geostationary, 0);
                
            obj.setElVal(obj.uiIDs.edit_th0_Tx, [], 0);

            obj.setElVal(obj.uiIDs.edit_ph0_Tx, [], 0);
            
            
            %% RECEIVER (Rx) INPUTS
            obj.setElVal(obj.uiIDs.edit_hr_m, [], 0);
            
            obj.setElVal(obj.uiIDs.edit_G0r_dB, [], 0);
            
            obj.init_popup_pol_Rx();
            obj.setElVal(obj.uiIDs.popup_pol_Rx, Constants.id_pol_R, 0);
            
            obj.init_popup_orientation_Rx();
            obj.setElVal(obj.uiIDs.popup_orientation_Rx, Constants.id_Rx_fixed, 0);
                
            obj.setElVal(obj.uiIDs.edit_th0_Rx, [], 0);

            obj.setElVal(obj.uiIDs.edit_ph0_Rx, [], 0);
            
            obj.init_popup_ant_pat_Rx();
            obj.setElVal(obj.uiIDs.popup_ant_pat_Rx, Constants.id_Rx_GG, 0);
            
            obj.setElVal(obj.uiIDs.edit_ant_pat_res_Rx, [], 0);

            obj.setElVal(obj.uiIDs.edit_hpbw_deg, [], 0);

            obj.setElVal(obj.uiIDs.edit_SLL_dB, [], 0);

            obj.setElVal(obj.uiIDs.edit_XPL_dB, [], 0);
            
            
            %% GROUND INPUTS
            obj.init_popup_diel_model();
            obj.setElVal(obj.uiIDs.popup_diel_model, Constants.id_diel_dobson, 0);
            
            obj.init_popup_gnd_structure();
            obj.setElVal(obj.uiIDs.popup_gnd_structure, Constants.id_gnd_single_layered, 0);            
            
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
            
            % TO-DO:
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

            % TO-DO:
            % If the default input is used, make the full input path 
            if isDefaultInput
                
                inputStruct.config_inputs_file = strcat(dir_input_config, '\', inputStruct.config_inputs_file);
                
                if strcmp( inputStruct.ant_pat_Rx, Constants.Rx_ant_pats{1, Constants.id_Rx_user_defined } ) 
                
                    inputStruct.ant_pat_Rx_file = strcat(dir_input_ant_pat_Rx, '\', inputStruct.ant_pat_Rx_file);
                    
                end
                
                if strcmp( inputStruct.gnd_cover, Constants.gnd_covers{1, Constants.id_veg_cover } ) 
                
                    inputStruct.veg_inputs_file = strcat(dir_input_veg, '\', inputStruct.veg_inputs_file);
                    
                end
                
            end
            
            
            %% SET GUI ELEMENTS' VALUES FROM inputStruct
            %%   SIMULATION SETTINGS        
            obj.setElVal(obj.uiIDs.edit_campaign, inputStruct.campaign, 0);
            
            % Simulation Mode: Snapshot OR Time-series
            obj.init_popup_sim_mode();
            sim_mode_id = findElementIdInCell( Constants.sim_modes, inputStruct.sim_mode );
            obj.setElVal( obj.uiIDs.popup_sim_mode, sim_mode_id, 0 );
            
            % Ground cover: Bare-soil OR Vegetation
            obj.init_popup_gnd_cover();
            gnd_cover_id = findElementIdInCell( Constants.gnd_covers, inputStruct.gnd_cover );
            obj.setElVal(obj.uiIDs.popup_gnd_cover, gnd_cover_id, 0);
            
            % Flag to write Attenuation to Excel file
            obj.setElVal(obj.uiIDs.cb_write_attenuation, inputStruct.write_attenuation, 0);
            
            % Flag to include the simulation in the Master Simulation file
            obj.setElVal(obj.uiIDs.cb_include_in_master_sim_file, inputStruct.include_in_master_sim_file, 0);
            
            
            %% TRANSMITTER (Tx) INPUTS
            obj.setElVal(obj.uiIDs.edit_f_MHz, num2str(inputStruct.f_MHz), 0);
            
            obj.setElVal(obj.uiIDs.edit_r_Tx_km, num2str(inputStruct.r_Tx_km), 0);
            
            obj.setElVal(obj.uiIDs.edit_EIRP_dB, num2str(inputStruct.EIRP_dB), 0);
            
            obj.init_popup_pol_Tx();
            pol_Tx_id = findElementIdInCell( Constants.polarizations, inputStruct.pol_Tx );
            obj.setElVal(obj.uiIDs.popup_pol_Tx, pol_Tx_id, 0);
            
            obj.init_popup_orientation_Tx();
            orientation_Tx_id = findElementIdInCell( Constants.Tx_orientations, inputStruct.orientation_Tx );
            obj.setElVal(obj.uiIDs.popup_orientation_Tx, orientation_Tx_id, 0);
            
            % If transmitter orientation is Geo-stationary, then load 
            % incidence and azimuth agles
            if orientation_Tx_id == Constants.id_Tx_geostationary
                
                obj.setElVal(obj.uiIDs.edit_th0_Tx, num2str(inputStruct.th0_Tx_deg), 0);
            
                obj.setElVal(obj.uiIDs.edit_ph0_Tx, num2str(inputStruct.ph0_Tx_deg), 0);
                
            end
            
            
            %% RECEIVER (Rx) INPUTS
            obj.setElVal(obj.uiIDs.edit_hr_m, num2str(inputStruct.hr_m), 0);
            
            obj.setElVal(obj.uiIDs.edit_G0r_dB, num2str(inputStruct.G0r_dB), 0);
            
            obj.init_popup_pol_Rx();
            pol_Rx_id = findElementIdInCell( Constants.polarizations, inputStruct.pol_Rx );
            obj.setElVal(obj.uiIDs.popup_pol_Rx, pol_Rx_id, 0);
            
            obj.init_popup_orientation_Rx();
            orientation_Rx_id = findElementIdInCell( Constants.Rx_orientations, inputStruct.orientation_Rx );
            obj.setElVal(obj.uiIDs.popup_orientation_Rx, orientation_Rx_id, 0);
            
            % If receiver orientation is fixed, then load incidence and
            % azimuth agles
            if orientation_Rx_id == Constants.id_Rx_fixed
                
                obj.setElVal(obj.uiIDs.edit_th0_Rx, num2str(inputStruct.th0_Rx_deg), 0);
            
                obj.setElVal(obj.uiIDs.edit_ph0_Rx, num2str(inputStruct.ph0_Rx_deg), 0);
                
            end
            
            obj.init_popup_ant_pat_Rx();
            ant_pat_Rx_id = findElementIdInCell( Constants.Rx_ant_pats, inputStruct.ant_pat_Rx );
            obj.setElVal(obj.uiIDs.popup_ant_pat_Rx, ant_pat_Rx_id, 0);
            
            % If receiver antenna pattern is Generalized-Gaussian, then 
            % load the corresponding inputs
            if ant_pat_Rx_id == Constants.id_Rx_GG
            
                obj.setElVal(obj.uiIDs.edit_ant_pat_res_Rx, num2str(inputStruct.ant_pat_res_deg_Rx), 0);
            
                obj.setElVal(obj.uiIDs.edit_hpbw_deg, num2str(inputStruct.hpbw_deg), 0);
                
                obj.setElVal(obj.uiIDs.edit_SLL_dB, num2str(inputStruct.SLL_dB), 0);
                
                obj.setElVal(obj.uiIDs.edit_XPL_dB, num2str(inputStruct.XPL_dB), 0);
                
            end
            
            
            %% GROUND INPUTS
            obj.init_popup_diel_model();
            diel_model_id = findElementIdInCell( Constants.diel_models, inputStruct.diel_model );
            obj.setElVal(obj.uiIDs.popup_diel_model, diel_model_id, 0);
            
            obj.init_popup_gnd_structure();
            gnd_structure_id = findElementIdInCell( Constants.gnd_structures, inputStruct.gnd_structure );
            obj.setElVal(obj.uiIDs.popup_gnd_structure, gnd_structure_id, 0);            
            
            if gnd_structure_id == Constants.id_gnd_multi_layered
                
                % Flag to calculate Discrete-slab for dielectric profile 
                obj.setElVal(obj.uiIDs.cb_discrete_slab, inputStruct.calc_diel_profile_fit_functions(Constants.id_diel_slab,1), 0);            

                % Flag to calculate Logistic regression for dielectric profile 
                obj.setElVal(obj.uiIDs.cb_logistic_regression, inputStruct.calc_diel_profile_fit_functions(Constants.id_diel_logistic,1), 0);            

                % Flag to calculate 2nd-order poly-fit for dielectric profile 
                obj.setElVal(obj.uiIDs.cb_2nd_order, inputStruct.calc_diel_profile_fit_functions(Constants.id_diel_2nd_order,1), 0);            

                % Flag to calculate 3rd order poly-fit for dielectric profile 
                obj.setElVal(obj.uiIDs.cb_3rd_order, inputStruct.calc_diel_profile_fit_functions(Constants.id_diel_3rd_order,1), 0);
                
            end
            
            
            %% INPUT FILES  
            obj.setElVal(obj.uiIDs.edit_config_inputs_file, inputStruct.config_inputs_file, 0);     
            
            if ant_pat_Rx_id == Constants.id_Rx_user_defined
            
                obj.setElVal(obj.uiIDs.edit_ant_pat_Rx_file, inputStruct.ant_pat_Rx_file, 0);
                
            end
            
            if gnd_cover_id == Constants.id_veg_cover
            
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
        inputStruct.sim_mode            = Constants.sim_modes{ 1, sim_mode_id };

        gnd_cover_id                    = obj.getElVal(obj.uiIDs.popup_gnd_cover);            
        inputStruct.gnd_cover           = Constants.gnd_covers{ 1, gnd_cover_id };

        inputStruct.write_attenuation	= obj.getElVal(obj.uiIDs.cb_write_attenuation);

        inputStruct.include_in_master_sim_file	= obj.getElVal(obj.uiIDs.cb_include_in_master_sim_file);


        %% TRANSMITTER (Tx) INPUTS PANEL ELEMENTS  
        inputStruct.f_MHz     = str2double(obj.getElVal(obj.uiIDs.edit_f_MHz));

        inputStruct.r_Tx_km	= str2double(obj.getElVal(obj.uiIDs.edit_r_Tx_km));

        inputStruct.EIRP_dB	= str2double(obj.getElVal(obj.uiIDs.edit_EIRP_dB));

        pol_Tx_id       = obj.getElVal(obj.uiIDs.popup_pol_Tx);
        inputStruct.pol_Tx	= Constants.polarizations{ 1, pol_Tx_id };

        orientation_Tx_id       = obj.getElVal(obj.uiIDs.popup_orientation_Tx);
        inputStruct.orientation_Tx	= Constants.Tx_orientations{ 1, orientation_Tx_id };

        % If transmitter orientation is Geo-stationary, then assign 
        % incidence and azimuth angles
        if orientation_Tx_id == Constants.id_Tx_geostationary
            
            inputStruct.th0_Tx_deg	= str2double(obj.getElVal(obj.uiIDs.edit_th0_Tx));
            inputStruct.ph0_Tx_deg	= str2double(obj.getElVal(obj.uiIDs.edit_ph0_Tx));
            
        end


        %% RECEIVER (Rx) INPUTS PANEL ELEMENTS 
        inputStruct.hr_m      = str2double(obj.getElVal(obj.uiIDs.edit_hr_m));

        inputStruct.G0r_dB	= str2double(obj.getElVal(obj.uiIDs.edit_G0r_dB));

        pol_Rx_id       = obj.getElVal(obj.uiIDs.popup_pol_Rx);
        inputStruct.pol_Rx	= Constants.polarizations{ 1, pol_Rx_id };

        orientation_Rx_id       = obj.getElVal(obj.uiIDs.popup_orientation_Rx);
        inputStruct.orientation_Rx	= Constants.Rx_orientations{ 1, orientation_Rx_id };

        % If receiver orientation is fixed, then assign incidence and
        % azimuth angles
        if orientation_Rx_id == Constants.id_Rx_fixed
            
            inputStruct.th0_Rx_deg	= str2double(obj.getElVal(obj.uiIDs.edit_th0_Rx));
            inputStruct.ph0_Rx_deg	= str2double(obj.getElVal(obj.uiIDs.edit_ph0_Rx));
            
        end

        ant_pat_Rx_id       = obj.getElVal(obj.uiIDs.popup_ant_pat_Rx);
        inputStruct.ant_pat_Rx	= Constants.Rx_ant_pats{ 1, ant_pat_Rx_id };

        % If receiver antenna pattern is Generalized-Gaussian, then 
        % assign GG params
        if ant_pat_Rx_id == Constants.id_Rx_GG

            inputStruct.hpbw_deg = str2double( obj.getElVal(obj.uiIDs.edit_hpbw_deg) );
            inputStruct.SLL_dB	= str2double( obj.getElVal(obj.uiIDs.edit_SLL_dB) );
            inputStruct.XPL_dB	= str2double( obj.getElVal(obj.uiIDs.edit_XPL_dB) );
            inputStruct.ant_pat_res_deg_Rx	= str2double(obj.getElVal(obj.uiIDs.edit_ant_pat_res_Rx));

        end


        %% GROUND INPUTS PANEL ELEMENTS
        diel_model_id       = obj.getElVal(obj.uiIDs.popup_diel_model);
        inputStruct.diel_model	= Constants.diel_models{ 1, diel_model_id };

        gnd_structure_id            = obj.getElVal(obj.uiIDs.popup_gnd_structure);
        inputStruct.gnd_structure	= Constants.gnd_structures{ 1, gnd_structure_id };

        % If receiver antenna pattern is Generalized-Gaussian, then 
        % assign GG params
        if gnd_structure_id == Constants.id_gnd_multi_layered

            calc_diel_profile_fit_functions(Constants.id_diel_slab,1) = obj.getElVal(obj.uiIDs.cb_discrete_slab);
            calc_diel_profile_fit_functions(Constants.id_diel_logistic,1) = obj.getElVal(obj.uiIDs.cb_logistic_regression);
            calc_diel_profile_fit_functions(Constants.id_diel_2nd_order,1) = obj.getElVal(obj.uiIDs.cb_2nd_order);
            calc_diel_profile_fit_functions(Constants.id_diel_3rd_order,1) = obj.getElVal(obj.uiIDs.cb_3rd_order);
            
            inputStruct.calc_diel_profile_fit_functions = calc_diel_profile_fit_functions;

        end


        %% INPUT FILES PANEL ELEMENTS
        inputStruct.config_inputs_file	= obj.getElVal(obj.uiIDs.edit_config_inputs_file);

        if ant_pat_Rx_id == Constants.id_Rx_user_defined

            inputStruct.ant_pat_Rx_file	= obj.getElVal(obj.uiIDs.edit_ant_pat_Rx_file);
            
        end

        if gnd_cover_id == Constants.id_veg_cover
        
            inputStruct.veg_inputs_file = obj.getElVal(obj.uiIDs.edit_veg_inputs_file); 
            
        end


        %% SAVE ALL TO FILE
        % If save option is "save" AND inputStructs should be different, OR
        % save option is "save as", then open a saving dialog
        if ( saveOption == Constants.id_GUI_save && ~isequaln( obj.inputStruct, inputStruct ) ) ...
            || saveOption == Constants.id_GUI_save_as

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
            

            lastInputFileName = Constants.lastInputFileNames{ 1, obj.simulator_id };
            
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

