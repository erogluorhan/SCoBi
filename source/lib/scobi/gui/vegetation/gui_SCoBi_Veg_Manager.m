% TO-DO: Check comments, copyrights, etc.
classdef gui_SCoBi_Veg_Manager < SCoBiGUIManagers
    %GUI_SCO_VEG_MANAGER This class implements handles for GUI elements, and 
    % performs the GUI-related functionalities
    
    
    properties(GetAccess = 'private', SetAccess = 'private')
        
        
        %% INPUT FILES  
        % Default input file that should be under \SCoBi\source\input\sys
        defaultInputFileName = 'default_input-scobi_veg.mat';   
        
        % Full path and name of the input file that is studied in the last
        % execution of SCoBi is kept in the following file, if any 
        lastInputFile = 'last_input.mat';
        
        % Struct data to provide input parameters as output of GUI
        inputStruct
           
    end
    
    methods
        
        % Constructor
        function obj = gui_SCoBi_Veg_Manager( handles )            
            
            % Call superclass constructor
            obj = obj@SCoBiGUIManagers( handles );
            
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


        % Check if veg_method value is Homogenous
        function result = is_popup_veg_method_homogenous(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_veg_method );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_veg_method) == Constants.id_veg_hom );
            
        end


        % Check if veg_method value is Virtual
        function result = is_popup_veg_method_virtual(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_veg_method );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_veg_method) == Constants.id_veg_vir );
            
        end


        % Check if veg_vir_orientations value is Random-spread
        function result = is_popup_veg_vir_orientations_random_spread(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_veg_vir_orientations );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_veg_vir_orientations) == Constants.id_veg_vir_random_spread );
            
        end


        % Check if veg_vir_orientations value is Row-crop
        function result = is_popup_veg_vir_orientations_row_crop(obj)
            
            isOn = obj.isEnabled( obj.uiIDs.popup_veg_vir_orientations );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_veg_vir_orientations) == Constants.id_veg_vir_row_crop );
            
        end
        
        
        function funout = outputFun(obj)
            
            funout{1} = obj.inputStruct;
            
        end
        
        
        % EVENT MANAGER
        % When an element is modified in the GUI, this must be called!
        % Sync the internal copy of the interface with the GUI
        function syncFromGUI(obj, idEl)
            
            
            %% GET GLOBAL DIRECTORIES
            dir_gui_veg = Directories.getInstance.scobi_gui_veg;
            
            
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
            
              % If selected sim_mode is Snapshot
              if obj.is_popup_sim_mode_snapshot()
                      
                  obj.setElStatus(obj.uiGroups.on_popup_sim_mode, 1, 0);
                  
              % Else if selected sim_mode is Time-series
              elseif obj.is_popup_sim_mode_time_series()
    
                  obj.setElStatus([obj.uiGroups.on_popup_sim_mode], 0, 0);
                  
              else
                      
                  % TO-DO: Handle Exception?

              end
              
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
          
          
          %% SIMULATION INPUTS
          % If veg_method popup value is changed            
          if sum(intersect(idEl, obj.uiIDs.popup_veg_method)) > 0                       
                  
              % If selected veg_method is Homogenous
              if obj.is_popup_veg_method_homogenous()

                  obj.setElStatus(obj.uiGroups.on_popup_veg_method, 0, 0);                  

              % Else if selected veg_method is Virtual
              elseif obj.is_popup_veg_method_virtual()

                  obj.setElStatus([obj.uiGroups.on_popup_veg_method], 1, 0);

              else

                  % TO-DO: Handle Exception?

              end
              
              obj.updateGUI();
              
          end
          
          % If veg_vir_orientations popup value is changed            
          if sum(intersect(idEl, obj.uiIDs.popup_veg_vir_orientations)) > 0                
                  
              % If selected diel_model is Dobson
              if obj.is_popup_veg_vir_orientations_random_spread()

                  % Display a warning that the method is not
                  % implemented yet
                  waitfor(msgbox('WARNING: This method is not yet implemented!'));

                  % Change the popup value to default Dobson model
                  obj.setElVal(obj.uiIDs.popup_veg_vir_orientations, Constants.id_veg_vir_row_crop);  

              end
              
              obj.updateGUI();
              
          end
          
          
          %% RECEIVER INPUTS
          % If orientation_Rx popup value is changed            
          if sum(intersect(idEl, obj.uiIDs.popup_orientation_Rx)) > 0                
                  
              % If selected orientation_Rx is Fixed
              if obj.is_popup_orientation_Rx_fixed()

                  obj.setElStatus(obj.uiGroups.on_popup_orientation_Rx, 1, 0);                  

              % Else if selected orientation_Rx is Specular-facing
              elseif obj.is_popup_orientation_Rx_specular_facing()

                  obj.setElStatus([obj.uiGroups.on_popup_orientation_Rx], 0, 0);

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

                  % TO-DO: GUI window for GG .id_Rx_GG);            

              % Else if selected ant_pat_Rx is User-defined
              elseif obj.is_popup_ant_pat_Rx_user_defined()
  
                  last_ant_pat_res_Rx_val = str2double( obj.getElVal( obj.uiIDs.edit_ant_pat_res_Rx ) );

                  obj.setElStatus([obj.uiGroups.on_popup_ant_pat_Rx_GG], 0, 0); 

                  obj.setElStatus([obj.uiGroups.on_popup_ant_pat_Rx_user_defined], 1, 0);                  

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
          
            
          %% PUSH BUTTONS                        
            % on Load Settings
            if sum(intersect(idEl, obj.uiIDs.pb_load_inputs)) > 0
	
                filter = {'*.mat'};
                [file, path] = uigetfile(filter, 'Load GUI From Input File');
                
                
                % If file and path are not empty
                if length(path)>1 && length(file)>1
                
                    obj.inputStruct = obj.loadGUIFromInputFile( file, path );
        
                    filename = strcat( dir_gui_veg, '/', obj.lastInputFile);
                    lastInput.lastInputFileName = strcat( path, file );
                    save(filename, 'lastInput');
                    
                end
                
            end
            
            % on Save Settings
            if sum(intersect(idEl, obj.uiIDs.pb_save_inputs)) > 0
                    
                obj.saveGUIToInputFile();
                
            end
            
            % on SCoBi
            if sum(intersect(idEl, obj.uiIDs.pb_SCoBi)) > 0

                % Try to save GUI values to an input file first
                % It will perform a check for any changes to the loaded or recently 
                % saved inputs   
                savingResult = obj.saveGUIToInputFile();

                if savingResult == 1 || savingResult == 2

                    uiresume(obj.handles.panel_main);

                end
                
            end   
                     
            % on Browse Dynamic Inputs File
            if sum(intersect(idEl, obj.uiIDs.pb_dyn_inputs_file)) > 0
	
                filter = {'*.xlsx'};
                [file, path] = uigetfile(filter, 'Load Dynamic Inputs File');
                
                
                % If file and path are not empty
                if length(path)>1 && length(file)>1
                
                   filename = strcat( path, '\',file );

                   obj.setElVal(obj.uiIDs.edit_dyn_inputs_file, filename, 0);  
                    
                end
                
            end  
                     
            % on Browse Antenna Pattern Inputs File
            if sum(intersect(idEl, obj.uiIDs.pb_ant_pat_Rx_file)) > 0
	
                filter = {'*.xlsx'};
                [file, path] = uigetfile(filter, 'Load Antenna Pattern Inputs File');
                
                
                % If file and path are not empty
                if length(path)>1 && length(file)>1
                
                   filename = strcat( path, '\',file );

                   obj.setElVal(obj.uiIDs.edit_ant_pat_Rx_file, filename, 0);  
                    
                end
                
            end
                     
            % on Browse Antenna Pattern Inputs File
            if sum(intersect(idEl, obj.uiIDs.pb_veg_inputs_file)) > 0
	
                filter = {'*.xlsx'};
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
        function init( obj, handles )
            
            tic;
            
            % Call superclass's init() function
            init@SCoBiGUIManagers( obj, handles );
            
            t0 = toc;
            
            fprintf('SCoBi-Veg GUI initialization completed in %.2f seconds\n', t0);
        
        end
        
        function initPopupMenus(obj)
            obj.init_popup_sim_mode();
            obj.init_popup_gnd_cover();
            obj.init_popup_veg_method();
            obj.init_popup_veg_vir_orientation();
            obj.init_popup_pol_Tx();
            obj.init_popup_pol_Rx();
            obj.init_popup_ant_pat_Rx();
            obj.init_popup_diel_model();
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
                str = Constants.polarizations;
            end
            
            value = get(obj.handles.popup_pol_Rx,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_pol_Rx,'Value', value);
            set(obj.handles.popup_pol_Rx,'String', str);
        end

        
        % Fill the popup_pol_Tx
        function init_popup_pol_Tx(obj, str)
            
            if nargin < 2
                str = Constants.polarizations;
            end
            
            value = get(obj.handles.popup_pol_Tx,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_pol_Tx,'Value', value);
            set(obj.handles.popup_pol_Tx,'String', str);
        end

        
        % Fill the popup_veg_vir_orientations
        function init_popup_veg_vir_orientation(obj, str)
            
            if nargin < 2
                str = Constants.veg_vir_orientations;
            end
            
            value = get(obj.handles.popup_veg_vir_orientations,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_veg_vir_orientations,'Value', value);
            set(obj.handles.popup_veg_vir_orientations,'String', str);
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

        
        % Fill the popup_veg_method
        function init_popup_veg_method(obj, str)
            
            if nargin < 2
                str = Constants.veg_methods;
            end
            
            value = get(obj.handles.popup_veg_method,'Value');
            value = min(1,max(length(str), value));
            set(obj.handles.popup_veg_method,'Value', value);
            set(obj.handles.popup_veg_method,'String', str);
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
            
          i = i+1;        id.text_title = i;                pointers(i) = obj.handles.text_title;
          i = i+1;        id.text_description = i;          pointers(i) = obj.handles.text_description;
          
          i = i+1;        id.pb_load_inputs = i;            pointers(i) = obj.handles.pb_load_inputs;          
          i = i+1;        id.pb_save_inputs = i;            pointers(i) = obj.handles.pb_save_inputs;        
          i = i+1;        id.pb_SCoBi = i;                  pointers(i) = obj.handles.pb_SCoBi;       
          i = i+1;        id.pb_exit = i;                   pointers(i) = obj.handles.pb_exit;
            
            
          %% PANELS            
          i = i+1;        id.panel_sim_settings = i;        pointers(i) = obj.handles.panel_sim_settings;
          i = i+1;        id.panel_sim_inputs = i;          pointers(i) = obj.handles.panel_sim_inputs;
          i = i+1;        id.panel_Tx_inputs = i;           pointers(i) = obj.handles.panel_Tx_inputs;
          i = i+1;        id.panel_Rx_inputs = i;           pointers(i) = obj.handles.panel_Rx_inputs;
          i = i+1;        id.panel_gnd_inputs = i;          pointers(i) = obj.handles.panel_gnd_inputs;
          i = i+1;        id.panel_input_files = i;         pointers(i) = obj.handles.panel_input_files;
                        
          
          %% SIMULATION SETTINGS PANEL ELEMENTS           
          i = i+1;        id.text_sim_mode = i;             pointers(i) = obj.handles.text_sim_mode;  
          i = i+1;        id.popup_sim_mode = i;            pointers(i) = obj.handles.popup_sim_mode;
          i = i+1;        id.text_gnd_cover = i;            pointers(i) = obj.handles.text_gnd_cover;
          i = i+1;        id.popup_gnd_cover = i;           pointers(i) = obj.handles.popup_gnd_cover;
          
          i = i+1;        id.panel_preferences = i;         pointers(i) = obj.handles.panel_preferences;   % Another panel inside a panel
          i = i+1;        id.cb_write_attenuation = i;      pointers(i) = obj.handles.cb_write_attenuation;
          i = i+1;        id.cb_calc_direct_term = i;       pointers(i) = obj.handles.cb_calc_direct_term;
          i = i+1;        id.cb_calc_specular_term = i;     pointers(i) = obj.handles.cb_calc_specular_term;
          i = i+1;        id.cb_calc_diffuse_term = i;      pointers(i) = obj.handles.cb_calc_diffuse_term;
          
            
          groupIDs.sim_settings = [id.panel_sim_settings id.text_sim_mode : id.cb_calc_diffuse_term];
          
          
          %% SIMULATION INPUTS PANEL ELEMENTS           
          i = i+1;        id.text_sim_name = i;                 pointers(i) = obj.handles.text_sim_name;
          i = i+1;        id.edit_sim_name = i;                 pointers(i) = obj.handles.edit_sim_name;
          i = i+1;        id.text_campaign = i;                 pointers(i) = obj.handles.text_campaign;
          i = i+1;        id.edit_campaign = i;                 pointers(i) = obj.handles.edit_campaign;
          i = i+1;        id.text_campaign_date = i;            pointers(i) = obj.handles.text_campaign_date;
          i = i+1;        id.edit_campaign_date = i;            pointers(i) = obj.handles.edit_campaign_date;
          i = i+1;        id.text_plot = i;                     pointers(i) = obj.handles.text_plot;
          i = i+1;        id.edit_plot = i;                     pointers(i) = obj.handles.edit_plot;
          i = i+1;        id.text_veg_plant = i;                pointers(i) = obj.handles.text_veg_plant;
          i = i+1;        id.edit_veg_plant = i;                pointers(i) = obj.handles.edit_veg_plant;
          i = i+1;        id.text_veg_method = i;               pointers(i) = obj.handles.text_veg_method;
          i = i+1;        id.popup_veg_method = i;              pointers(i) = obj.handles.popup_veg_method;
          i = i+1;        id.text_veg_vir_orientations = i;     pointers(i) = obj.handles.text_veg_vir_orientations;
          i = i+1;        id.popup_veg_vir_orientations = i;    pointers(i) = obj.handles.popup_veg_vir_orientations;
          i = i+1;        id.text_Nr = i;                       pointers(i) = obj.handles.text_Nr;
          i = i+1;        id.edit_Nr = i;                       pointers(i) = obj.handles.edit_Nr;
          i = i+1;        id.text_Nfz = i;                      pointers(i) = obj.handles.text_Nfz;
          i = i+1;        id.edit_Nfz = i;                      pointers(i) = obj.handles.edit_Nfz;          

          groupIDs.sim_inputs       = [id.panel_sim_inputs id.text_sim_name : id.edit_Nfz];
          
          
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
            
          groupIDs.Tx_inputs = [id.panel_Tx_inputs id.text_f_MHz : id.popup_pol_Tx];
          
          
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
          i = i+1;        id.text_th0_Rx = i;               pointers(i) = obj.handles.text_th0_Rx;
          i = i+1;        id.edit_th0_Rx = i;               pointers(i) = obj.handles.edit_th0_Rx;
          i = i+1;        id.text_deg_th0_Rx = i;           pointers(i) = obj.handles.text_deg_th0_Rx;
          i = i+1;        id.text_ph0_Rx = i;               pointers(i) = obj.handles.text_ph0_Rx;
          i = i+1;        id.edit_ph0_Rx = i;               pointers(i) = obj.handles.edit_ph0_Rx;
          i = i+1;        id.text_deg_ph0_Rx = i;           pointers(i) = obj.handles.text_deg_ph0_Rx;
          i = i+1;        id.text_ant_pat_Rx = i;           pointers(i) = obj.handles.text_ant_pat_Rx;
          i = i+1;        id.popup_ant_pat_Rx = i;          pointers(i) = obj.handles.popup_ant_pat_Rx;
          i = i+1;        id.text_ant_pat_res_Rx = i;       pointers(i) = obj.handles.text_ant_pat_res_Rx;
          i = i+1;        id.edit_ant_pat_res_Rx = i;       pointers(i) = obj.handles.edit_ant_pat_res_Rx;
          i = i+1;        id.text_deg_ant_pat_res_Rx = i;   pointers(i) = obj.handles.text_deg_ant_pat_res_Rx;
          % Generalized-Gaussian Antenna Pattern Input
          i = i+1;        id.panel_Rx_GG = i;               pointers(i) = obj.handles.panel_Rx_GG;
          i = i+1;        id.text_hpbw_deg = i;             pointers(i) = obj.handles.text_hpbw_deg;
          i = i+1;        id.edit_hpbw_deg = i;             pointers(i) = obj.handles.edit_hpbw_deg;
          i = i+1;        id.text_deg_hpbw = i;             pointers(i) = obj.handles.text_deg_hpbw;
          i = i+1;        id.text_SLL_dB = i;               pointers(i) = obj.handles.text_SLL_dB;
          i = i+1;        id.edit_SLL_dB = i;               pointers(i) = obj.handles.edit_SLL_dB;
          i = i+1;        id.text_dB_SLL = i;               pointers(i) = obj.handles.text_dB_SLL;
          i = i+1;        id.text_XPL_dB = i;               pointers(i) = obj.handles.text_XPL_dB;
          i = i+1;        id.edit_XPL_dB = i;               pointers(i) = obj.handles.edit_XPL_dB;
          i = i+1;        id.text_dB_XPL = i;               pointers(i) = obj.handles.text_dB_XPL;
            
          groupIDs.Rx_inputs = [id.panel_Rx_inputs id.text_hr_m : id.text_deg_ant_pat_res_Rx];
          groupIDs.ant_pat_Rx_GG_inputs = [id.panel_Rx_GG : id.text_dB_XPL];
          
          
          %% GROUND INPUTS PANEL ELEMENTS           
          i = i+1;        id.text_sand_ratio = i;           pointers(i) = obj.handles.text_sand_ratio;
          i = i+1;        id.edit_sand_ratio = i;           pointers(i) = obj.handles.edit_sand_ratio;
          i = i+1;        id.text_clay_ratio = i;           pointers(i) = obj.handles.text_clay_ratio;
          i = i+1;        id.edit_clay_ratio = i;           pointers(i) = obj.handles.edit_clay_ratio;
          i = i+1;        id.text_rhob_gcm3 = i;            pointers(i) = obj.handles.text_rhob_gcm3;
          i = i+1;        id.edit_rhob_gcm3 = i;            pointers(i) = obj.handles.edit_rhob_gcm3;
          i = i+1;        id.text_gcm3_rhob = i;            pointers(i) = obj.handles.text_gcm3_rhob;         
          i = i+1;        id.text_diel_model = i;           pointers(i) = obj.handles.text_diel_model;  
          i = i+1;        id.popup_diel_model = i;          pointers(i) = obj.handles.popup_diel_model;
            
          groupIDs.gnd_inputs = [id.panel_gnd_inputs id.text_sand_ratio:id.popup_diel_model];
          
          
          %% INPUT FILES PANEL ELEMENTS
          i = i+1;        id.text_LED_dyn_inputs_file = i;  pointers(i) = obj.handles.text_LED_dyn_inputs_file;
          i = i+1;        id.text_dyn_inputs_file = i;      pointers(i) = obj.handles.text_dyn_inputs_file;
          i = i+1;        id.edit_dyn_inputs_file = i;      pointers(i) = obj.handles.edit_dyn_inputs_file;
          i = i+1;        id.pb_dyn_inputs_file = i;        pointers(i) = obj.handles.pb_dyn_inputs_file;
          i = i+1;        id.text_LED_ant_pat_Rx_file = i;  pointers(i) = obj.handles.text_LED_ant_pat_Rx_file;
          i = i+1;        id.text_ant_pat_Rx_file = i;      pointers(i) = obj.handles.text_ant_pat_Rx_file;
          i = i+1;        id.edit_ant_pat_Rx_file = i;      pointers(i) = obj.handles.edit_ant_pat_Rx_file;
          i = i+1;        id.pb_ant_pat_Rx_file = i;        pointers(i) = obj.handles.pb_ant_pat_Rx_file;
          i = i+1;        id.text_LED_veg_inputs_file = i;  pointers(i) = obj.handles.text_LED_veg_inputs_file;
          i = i+1;        id.text_veg_inputs_file = i;      pointers(i) = obj.handles.text_veg_inputs_file;
          i = i+1;        id.edit_veg_inputs_file = i;      pointers(i) = obj.handles.edit_veg_inputs_file;
          i = i+1;        id.pb_veg_inputs_file = i;        pointers(i) = obj.handles.pb_veg_inputs_file;
            

          %% GUI ENABLE/DISABLE GROUPS
          % Virtual vegetation related group
          groupIDs.veg_vir_inputs   = [id.cb_calc_diffuse_term ...
                                       id.text_veg_method : id.edit_Nfz ];
          
          % On sim_mode change
          % Because there are currently only two sim modes, it can be 
          % adjusted by only on_popup_sim_mode. If three or more items exist in 
          % the future, each should have its own enable/disable group) 
          groupIDs.on_popup_sim_mode = [id.cb_calc_diffuse_term ...
                                  id.text_veg_method ...
                                  id.popup_veg_method ...
                                  id.text_veg_vir_orientations ... 
                                  id.popup_veg_vir_orientations ...
                                  id.text_Nr ...
                                  id.edit_Nr ...
                                  id.text_Nfz ...
                                  id.edit_Nfz ];
                              
          % On gnd_cover change
          % Because there are currently only two ground covers, it can be 
          % adjusted by only on_popup_gnd_cover. If three or more items exist in 
          % the future, each should have its own enable/disable group) 
          groupIDs.on_popup_gnd_cover = [id.cb_calc_diffuse_term ...
                                   id.text_veg_plant ...
                                   id.edit_veg_plant ...
                                   id.text_veg_method ...
                                   id.popup_veg_method ...
                                   id.text_veg_vir_orientations ... 
                                   id.popup_veg_vir_orientations ...
                                   id.text_Nr ...
                                   id.edit_Nr ...
                                   id.text_Nfz ...
                                   id.edit_Nfz ...
                                   id.text_veg_inputs_file ...
                                   id.edit_veg_inputs_file ...
                                   id.pb_veg_inputs_file ];
                               
          % On veg_method change
          % Because there are currently only two vegetation methods, it can
          % be adjusted by only on_vege_method. If more items exist in the
          % future, each should have its own enable/disable group) 
          groupIDs.on_popup_veg_method = [id.text_veg_vir_orientations ... 
                                    id.popup_veg_vir_orientations ];
                                
          % On Rx_orientation change
          % Because there are currently only two Rx orientations, it can be
          % adjusted by only on_Rx_orientation. If more items exist in the
          % future, each should have its own enable/disable group)             
          groupIDs.on_popup_orientation_Rx = [id.text_th0_Rx id.edit_th0_Rx ...
                                        id.text_deg_th0_Rx id.text_ph0_Rx ...
                                        id.edit_ph0_Rx id.text_deg_ph0_Rx];
          
          % On ant_pat_Rx change
          groupIDs.on_popup_ant_pat_Rx_GG = [id.text_ant_pat_res_Rx ...
                                             id.edit_ant_pat_res_Rx ...
                                             id.text_deg_ant_pat_res_Rx ...
                                             groupIDs.ant_pat_Rx_GG_inputs ];
                                      
          groupIDs.on_popup_ant_pat_Rx_user_defined = [id.text_LED_ant_pat_Rx_file ...
                                                       id.text_ant_pat_Rx_file ...
                                                       id.edit_ant_pat_Rx_file ...
                                                       id.pb_ant_pat_Rx_file ];
                                      
          
          groupIDs.input_files = [id.panel_input_files id.text_dyn_inputs_file:id.text_veg_inputs_file];          
          groupIDs.gLED = [id.text_LED_dyn_inputs_file id.text_LED_ant_pat_Rx_file id.text_LED_veg_inputs_file];

            
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
            
            % sim_mode while gnd_cover
            if obj.is_popup_sim_mode_snapshot()

                  % If the gnd_cover is Bare-soil, then virtual
                  % vegetation related fields should not be enabled
                  if obj.is_popup_gnd_cover_bare_soil()

                      obj.setElStatus(obj.uiGroups.veg_vir_inputs, 0, 0);
                      
                      obj.setElVal( obj.uiIDs.cb_calc_diffuse_term, 0 );
                      obj.setElVal( obj.uiIDs.popup_veg_method, Constants.id_veg_hom );
                      obj.setElVal( obj.uiIDs.edit_Nr, num2str(1) );
                      obj.setElVal( obj.uiIDs.edit_Nfz, num2str(1) );

                  end

            end
              
            
            % gnd_cover while sim_mode
            if obj.is_popup_gnd_cover_vegetation()

                  % If the sim_mode is Time-series, then virtual
                  % vegetation related fields should not be enabled
                  if obj.is_popup_sim_mode_time_series()

                      obj.setElStatus(obj.uiGroups.veg_vir_inputs, 0, 0);
                      
                      obj.setElVal( obj.uiIDs.cb_calc_diffuse_term, 0 );
                      obj.setElVal( obj.uiIDs.popup_veg_method, Constants.id_veg_hom );
                      obj.setElVal( obj.uiIDs.edit_Nr, num2str(1) );
                      obj.setElVal( obj.uiIDs.edit_Nfz, num2str(1) );

                  end
                  
            end
              
            
            % gnd_cover
            if obj.is_popup_gnd_cover_bare_soil()

                  obj.setElStatus(obj.uiGroups.on_popup_gnd_cover, 0, 0);
                  
            end
              
            
            % popup_veg_method
            if obj.is_popup_veg_method_homogenous()

                  obj.setElStatus(obj.uiGroups.on_popup_veg_method, 0, 0);
                  
            end
              
            
            % popup_orientation_Rx
            if obj.is_popup_orientation_Rx_specular_facing()

                  obj.setElStatus(obj.uiGroups.on_popup_orientation_Rx, 0, 0);
                  
            end
              
            
            % popup_ant_pat_Rx
            if obj.is_popup_ant_pat_Rx_user_defined()

                  obj.setElStatus(obj.uiGroups.on_popup_ant_pat_Rx_GG, 0, 0);
                      
                  obj.setElVal( obj.uiIDs.edit_ant_pat_res_Rx, num2str(0) );

                  obj.setElStatus(obj.uiGroups.on_popup_ant_pat_Rx_user_defined, 1, 0);
                  
            elseif obj.is_popup_ant_pat_Rx_GG()

                  obj.setElStatus(obj.uiGroups.on_popup_ant_pat_Rx_GG, 1, 0);
                      
%                   obj.setElVal( obj.uiIDs.edit_ant_pat_res_Rx, num2str(1) );

                  obj.setElStatus(obj.uiGroups.on_popup_ant_pat_Rx_user_defined, 0, 0);
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
                
            % Modified LED elements
            mLED = obj.uiGroups.gLED(obj.getFlag(obj.uiGroups.gLED));
                
            % Modified text elements
            mTextEl = setdiff(idEl(textEl & obj.getFlag), mLED);

            mIdEl = setdiff(idEl(~panels & ~textEl & obj.getFlag), mLED); % id of elements that are not panels nor text elements that have been modified
                
            % For each modified panel
            for i=1:length(mPanel)
                obj.curVal{mPanel(i)} = obj.getGuiElTitle(obj.uiPointers(mPanel(i)));
            end
                
            % For each modified txt Element
            for i=1:length(mTextEl)
                obj.curVal{mTextEl(i)} = obj.getGuiElStr(obj.uiPointers(mTextEl(i)));
            end

            % For all the LEDs
            for i=1:length(mLED)
                obj.curVal{mLED(i)} = obj.getGuiElColor(obj.uiPointers(mLED(i)));
            end

            % For all the other modified elements
            for i=1:length(mIdEl)
                obj.curVal{mIdEl(i)} = obj.getGuiElVal(obj.uiPointers(mIdEl(i)));
            end
            
        end
        
        
        % Initialize specific GUI elements to this class
        function initSpecificGUI(obj)
            
            
            %% GET GLOBAL DIRECTORIES
            dir_gui_veg = Directories.getInstance.scobi_gui_veg;
            dir_input_sys = Directories.getInstance.input_sys;
            
            
            % Read the last input file, if any
            lastInputFileName = [];
            if exist( [strcat(dir_gui_veg, '\') obj.lastInputFile], 'file' )
                filename = strcat( dir_gui_veg, '\', obj.lastInputFile);
                load(filename);
                lastInputFileName = lastInput.lastInputFileName;
                % obj.lastInputFile contains the lastInputFileName 
                % string that might be saved during the last SCoBi run
            end
            
            % If there is an input file from the last execution of SCoBi
            if ~isempty( lastInputFileName ) && exist( lastInputFileName,'file' )
                
                obj.inputStruct = obj.loadGUIFromInputFile( [], lastInputFileName );
                
            % Else if the default input exists, load it
            elseif exist([strcat(dir_input_sys, '/') obj.defaultInputFileName],'file')
                
                obj.inputStruct = obj.loadGUIFromInputFile( obj.defaultInputFileName);
                
            % Else, make a warning about it and open an empty GUI
            else
                
                waitfor(msgbox('WARNING: There is no default input, you should add input files or give inputs manually!'));
                
                obj.loadEmptyGUI();
                
            end
            
        end
        
        
        % Load the state of the gui from a matlab file.
        function loadEmptyGUI( obj )

            
            %%   SIMULATION SETTINGS
            % Simulation mode: Snapshot OR Time-series
            obj.init_popup_sim_mode();
            obj.setElVal( obj.uiIDs.popup_sim_mode, Constants.id_snapshot, 0 );
            
            % Ground cover: Bare-soil OR Vegetation
            obj.init_popup_gnd_cover();
            obj.setElVal(obj.uiIDs.popup_gnd_cover, Constants.id_veg_cover, 0);
            
            % Flag to write Attenuation to Excel file
            obj.setElVal(obj.uiIDs.cb_write_attenuation, 0, 0);
            
            % Flag to calculate direct term
            obj.setElVal(obj.uiIDs.cb_calc_direct_term, 0, 0);
            
            % Flag to calculate specular term
            obj.setElVal(obj.uiIDs.cb_calc_specular_term, 0, 0);
            
            % Flag to calculate diffuse term
            obj.setElVal(obj.uiIDs.cb_calc_diffuse_term, 0, 0);
            
            
            %% SIMULATION INPUTS
            obj.setElVal(obj.uiIDs.edit_sim_name, '', 0);
            
            obj.setElVal(obj.uiIDs.edit_campaign, '', 0);
            
            obj.setElVal(obj.uiIDs.edit_campaign_date, '', 0);
            
            obj.setElVal(obj.uiIDs.edit_plot, '', 0);
            
            obj.setElVal(obj.uiIDs.edit_veg_plant, '', 0);
            
            obj.init_popup_veg_method();
            obj.setElVal(obj.uiIDs.popup_veg_method, Constants.id_veg_hom, 0);
            
            obj.init_popup_veg_vir_orientation();
            obj.setElVal(obj.uiIDs.popup_veg_vir_orientations, Constants.id_veg_vir_row_crop, 0);
            
            obj.setElVal(obj.uiIDs.edit_Nr, [], 0);
            
            obj.setElVal(obj.uiIDs.edit_Nfz, [], 0);
            
            
            %% TRANSMITTER (Tx) INPUTS
            obj.setElVal(obj.uiIDs.edit_f_MHz, [], 0);
            
            obj.setElVal(obj.uiIDs.edit_r_Tx_km, [], 0);
            
            obj.setElVal(obj.uiIDs.edit_EIRP_dB, [], 0);
            
            obj.init_popup_pol_Tx();
            obj.setElVal(obj.uiIDs.popup_pol_Tx, Constants.id_pol_R, 0);
            
            
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
            obj.setElVal(obj.uiIDs.edit_sand_ratio, [], 0);
            
            obj.setElVal(obj.uiIDs.edit_clay_ratio, [], 0);
            
            obj.setElVal(obj.uiIDs.edit_rhob_gcm3, [], 0);
            
            obj.init_popup_diel_model();
            obj.setElVal(obj.uiIDs.popup_diel_model, Constants.id_diel_dobson, 0);
            
            
            %% INPUT FILES  
            obj.setElVal(obj.uiIDs.edit_dyn_inputs_file, '', 0);     
            
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
            dir_input_dyn = Directories.getInstance.input_dyn;   
            dir_input_veg_hom = Directories.getInstance.input_veg_hom;               
            
            
            isDefaultInput = 0;
            
            % If the default input file is used
            if nargin == 2
                
                pathname = dir_input_sys;
                
                isDefaultInput = 1;
                
            end
            
            
            if ~isempty( filename )
            % Load the input file
                filename = strcat( pathname, '\',filename );
                load( filename ); % the file contains the variable inputStruct
            
            % If filename is empy, then pathname contains all    
            else
                load( pathname );
            end
            
            
%             % TO-DO: Delete these
%             inputStruct.dyn_inputs_file = 'default_input-scobi_veg-snapshot.xlsx';
%             
%             inputStruct.ant_pat_Rx_file = 'default_ant_pat.xlsx';
%             
%             inputStruct.veg_inputs_file = 'default_input-scobi_veg.xlsx';
%             
% 
%             inputStruct.hpbw_deg = 30;
%             inputStruct.SLL_dB = 25;
%             inputStruct.XPL_dB = 25;
%             save(filename, 'inputStruct');



            % If the default input is used, make the full input path 
            if isDefaultInput
                
                inputStruct.dyn_inputs_file = strcat(dir_input_dyn, '\', inputStruct.dyn_inputs_file);
                
                inputStruct.ant_pat_Rx_file = strcat(dir_input_sys, '\', inputStruct.ant_pat_Rx_file);
                
                inputStruct.veg_inputs_file = strcat(dir_input_veg_hom, '\', inputStruct.veg_inputs_file);
                
            end
            
            %%   SIMULATION SETTINGS
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
            
            % Flag to calculate direct term
            obj.setElVal(obj.uiIDs.cb_calc_direct_term, inputStruct.calc_direct_term, 0);
            
            % Flag to calculate specular term
            obj.setElVal(obj.uiIDs.cb_calc_specular_term, inputStruct.calc_specular_term, 0);
            
            % Flag to calculate diffuse term
            obj.setElVal(obj.uiIDs.cb_calc_diffuse_term, inputStruct.calc_diffuse_term, 0);
            
                    
            
            
            %% SIMULATION INPUTS
            obj.setElVal(obj.uiIDs.edit_sim_name, inputStruct.sim_name, 0);
            
            obj.setElVal(obj.uiIDs.edit_campaign, inputStruct.campaign, 0);
            
            obj.setElVal(obj.uiIDs.edit_campaign_date, inputStruct.campaign_date, 0);
            
            obj.setElVal(obj.uiIDs.edit_plot, inputStruct.plot, 0);
            
            obj.setElVal(obj.uiIDs.edit_veg_plant, inputStruct.veg_plant, 0);
            
            obj.init_popup_veg_method();
            veg_method_id = findElementIdInCell( Constants.veg_methods, inputStruct.veg_method );
            obj.setElVal(obj.uiIDs.popup_veg_method, veg_method_id, 0);
            
            if veg_method_id == Constants.id_veg_vir
                obj.init_popup_veg_vir_orientation();
                veg_vir_orientation_id = findElementIdInCell( Constants.veg_vir_orientations, inputStruct.veg_vir_orientation );
                obj.setElVal(obj.uiIDs.popup_veg_vir_orientations, veg_vir_orientation_id, 0);
            end
            
            obj.setElVal(obj.uiIDs.edit_Nr, num2str(inputStruct.Nr), 0);
            
            obj.setElVal(obj.uiIDs.edit_Nfz, num2str(inputStruct.Nfz), 0);
            
            
            %% TRANSMITTER (Tx) INPUTS
            obj.setElVal(obj.uiIDs.edit_f_MHz, num2str(inputStruct.f_MHz), 0);
            
            obj.setElVal(obj.uiIDs.edit_r_Tx_km, num2str(inputStruct.r_Tx_km), 0);
            
            obj.setElVal(obj.uiIDs.edit_EIRP_dB, num2str(inputStruct.EIRP_dB), 0);
            
            obj.init_popup_pol_Tx();
            pol_Tx_id = findElementIdInCell( Constants.polarizations, inputStruct.pol_Tx );
            obj.setElVal(obj.uiIDs.popup_pol_Tx, pol_Tx_id, 0);
            
            
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
                
            elseif ant_pat_Rx_id == Constants.id_Rx_user_defined
                
                % TO-DO : Load file
                
            end
            
            
            %% GROUND INPUTS
            obj.setElVal(obj.uiIDs.edit_sand_ratio, num2str(inputStruct.sand_ratio), 0);
            
            obj.setElVal(obj.uiIDs.edit_clay_ratio, num2str(inputStruct.clay_ratio), 0);
            
            obj.setElVal(obj.uiIDs.edit_rhob_gcm3, num2str(inputStruct.rhob_gcm3), 0);
            
            obj.init_popup_diel_model();
            diel_model_id = findElementIdInCell( Constants.diel_models, inputStruct.diel_model );
            obj.setElVal(obj.uiIDs.popup_diel_model, diel_model_id, 0);
            
            
            %% INPUT FILES  
            obj.setElVal(obj.uiIDs.edit_dyn_inputs_file, inputStruct.dyn_inputs_file, 0);     
            
            obj.setElVal(obj.uiIDs.edit_ant_pat_Rx_file, inputStruct.ant_pat_Rx_file, 0);
            
            obj.setElVal(obj.uiIDs.edit_veg_inputs_file, inputStruct.veg_inputs_file, 1);

            obj.updateGUI();
            
            %% LOAD GUI
            % Check all the dependencies
            obj.syncFromGUI(obj.uiIDs.popup_sim_mode);
            
        end
        
        % Save the state of the GUI to a matlab file
        function savingResult = saveGUIToInputFile(obj )
        % savingResult: 0 - Didn't save. Not good to go
        % savingResult: 1 - Saved. Good to go
        % savingResult: 2 - Didn't save, but good to go


        %% GET GLOBAL DIRECTORIES
        dir_gui_veg = Directories.getInstance.scobi_gui_veg;
     
            
        %% SIMULATION SETTINGS PANEL ELEMENTS
        sim_mode_id                     = obj.getElVal(obj.uiIDs.popup_sim_mode);            
        inputStruct.sim_mode            = Constants.sim_modes{ 1, sim_mode_id };

        gnd_cover_id                    = obj.getElVal(obj.uiIDs.popup_gnd_cover);            
        inputStruct.gnd_cover           = Constants.gnd_covers{ 1, gnd_cover_id };

        inputStruct.write_attenuation	= obj.getElVal(obj.uiIDs.cb_write_attenuation);

        inputStruct.calc_direct_term	= obj.getElVal(obj.uiIDs.cb_calc_direct_term);

        inputStruct.calc_specular_term	= obj.getElVal(obj.uiIDs.cb_calc_specular_term);

        inputStruct.calc_diffuse_term	= obj.getElVal(obj.uiIDs.cb_calc_diffuse_term);


        %% SIMULATION INPUTS PANEL ELEMENTS 
        inputStruct.sim_name        = obj.getElVal(obj.uiIDs.edit_sim_name);

        inputStruct.campaign        = obj.getElVal(obj.uiIDs.edit_campaign);

        inputStruct.campaign_date	= obj.getElVal(obj.uiIDs.edit_campaign_date);

        inputStruct.plot            = obj.getElVal(obj.uiIDs.edit_plot);

        inputStruct.veg_plant       = obj.getElVal(obj.uiIDs.edit_veg_plant);

        veg_method_id               = obj.getElVal(obj.uiIDs.popup_veg_method);
        inputStruct.veg_method      = Constants.veg_methods{ 1, veg_method_id };

        % If vegetation method is Virtual, then assign Virtual
        % orientation
        if veg_method_id == Constants.id_veg_vir
            veg_vir_orientation_id          = obj.getElVal(obj.uiIDs.popup_veg_vir_orientations);
            inputStruct.veg_vir_orientation	= Constants.veg_vir_orientations{ 1, veg_vir_orientation_id };
        end

        inputStruct.Nr	= str2double(obj.getElVal(obj.uiIDs.edit_Nr));

        if isnan(inputStruct.Nr)
            inputStruct.Nr = 1;
        end

        inputStruct.Nfz	= str2double(obj.getElVal(obj.uiIDs.edit_Nfz));

        if isnan(inputStruct.Nfz)
            inputStruct.Nfz = 1;
        end


        %% TRANSMITTER (Tx) INPUTS PANEL ELEMENTS  
        inputStruct.f_MHz     = str2double(obj.getElVal(obj.uiIDs.edit_f_MHz));

        inputStruct.r_Tx_km	= str2double(obj.getElVal(obj.uiIDs.edit_r_Tx_km));

        inputStruct.EIRP_dB	= str2double(obj.getElVal(obj.uiIDs.edit_EIRP_dB));

        pol_Tx_id       = obj.getElVal(obj.uiIDs.popup_pol_Tx);
        inputStruct.pol_Tx	= Constants.polarizations{ 1, pol_Tx_id };


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

        elseif ant_pat_Rx_id == Constants.id_Rx_user_defined

            % TO-DO: User0-defined filename

        end

        inputStruct.ant_pat_res_deg_Rx	= str2double(obj.getElVal(obj.uiIDs.edit_ant_pat_res_Rx));


        %% GROUND INPUTS PANEL ELEMENTS
        inputStruct.sand_ratio	= str2double(obj.getElVal(obj.uiIDs.edit_sand_ratio));

        inputStruct.clay_ratio	= str2double(obj.getElVal(obj.uiIDs.edit_clay_ratio));

        inputStruct.rhob_gcm3     = str2double(obj.getElVal(obj.uiIDs.edit_rhob_gcm3));

        diel_model_id       = obj.getElVal(obj.uiIDs.popup_diel_model);
        inputStruct.diel_model	= Constants.diel_models{ 1, diel_model_id };


        %% INPUT FILES PANEL ELEMENTS
        inputStruct.dyn_inputs_file	= obj.getElVal(obj.uiIDs.edit_dyn_inputs_file);

        inputStruct.ant_pat_Rx_file	= obj.getElVal(obj.uiIDs.edit_ant_pat_Rx_file);

        inputStruct.veg_inputs_file     = obj.getElVal(obj.uiIDs.edit_veg_inputs_file);            


        %% SAVE ALL TO FILE                
        if ~isequaln( obj.inputStruct, inputStruct )

            filter = {'*.mat'};
            [file, path] = uiputfile(filter, 'Save Input File');

            % If file and path are not empty
            if length(path)>1 && length(file)>1

                savingResult = 1;

                filename = strcat( path, file );

                save(filename, 'inputStruct');

                obj.inputStruct = inputStruct;

                filename = strcat( dir_gui_veg, '/', obj.lastInputFile);
                lastInput.lastInputFileName = strcat( path, file );
                save(filename, 'lastInput');

            else

                % No valid save file info. Do not save and not good to go
                savingResult = 0;

            end

        else
            % No change on the input file; so, do not save but good to go.
            savingResult = 2;

        end
            
        end
        
        
        function setSpecificElContent(obj)
                
            idEl = 1:length(obj.uiPointers);              % All the elements

            % Sets of panels
            panels = false(length(obj.uiPointers),1);     % init logical panels group
            panels(obj.uiGroups.gPanels) = true;           % set logical indexes of the panels
            % modified panels
            mPanel = idEl(panels & obj.setFlag);
                
            % Modified LED elements
            mLED = obj.uiGroups.gLED(obj.setFlag(obj.uiGroups.gLED));            
                
            % Sets of text elements
            textEl = false(length(obj.uiPointers),1);     % init logical text elements group
            textEl(obj.uiGroups.strEl) = true;              % set logical indexes of the text elements
            % Modified text elements
            mTextEl = setdiff(idEl(textEl & obj.setFlag), mLED);

            mIdEl = setdiff(idEl(~panels & ~textEl & obj.setFlag), mLED); % id of elements that are not panels nor text elements that have been modified
                
            % For each modified panel
            for i=1:length(mPanel)
                obj.setGuiElTitle(obj.uiPointers(mPanel(i)), obj.newVal{mPanel(i)});
            end
                
            % For all the LEDs
            for i=1:length(mLED)
                obj.setGuiElColor(obj.uiPointers(mLED(i)), obj.newVal{mLED(i)});
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
        
        
        % Test if the active file/dir paths
        % contain valid file/dir
        function updateLEDstate(obj)
            
            % Check edit_dyn_inputs_file
            if obj.isEnabled(obj.uiIDs.edit_dyn_inputs_file)
                
                filename = obj.getElVal(obj.uiIDs.edit_dyn_inputs_file);
                
                if isempty(filename)
                    
                    obj.setGUILedStatus(obj.uiIDs.text_LED_dyn_inputs_file, obj.ledKo, 0);
                    
%                     % If there is NO dyn_inputs_file, every LED should be RED
%                     for  i = 1:length(obj.uiGroups.gInINILED)
%                         obj.setGUILedStatus(obj.uiGroups.gInINILED(i), obj.ledKo, 0);
%                     end
                else
                    
                    if exist(filename,'file')
                        
                        obj.setGUILedStatus(obj.uiIDs.text_LED_dyn_inputs_file, obj.ledOk, 0);
                    else
                        obj.setGUILedStatus(obj.uiIDs.text_LED_dyn_inputs_file, obj.ledCk, 0);
                    end
                    
                end
            end
            
            % Check edit_ant_pat_Rx_file
            if obj.isEnabled(obj.uiIDs.edit_ant_pat_Rx_file)
                
                filename = obj.getElVal(obj.uiIDs.edit_ant_pat_Rx_file);
                
                if isempty(filename)
                    
                    obj.setGUILedStatus(obj.uiIDs.text_LED_ant_pat_Rx_file, obj.ledKo, 0);
                    
                else
                    
                    if exist(filename,'file')
                        
                        obj.setGUILedStatus(obj.uiIDs.text_LED_ant_pat_Rx_file, obj.ledOk, 0);
                    else
                        obj.setGUILedStatus(obj.uiIDs.text_LED_ant_pat_Rx_file, obj.ledCk, 0);
                    end
                    
                end
            end
            
            % Check edit_veg_inputs_file
            if obj.isEnabled(obj.uiIDs.edit_veg_inputs_file)
                
                filename = obj.getElVal(obj.uiIDs.edit_veg_inputs_file);
                
                if isempty(filename)
                    
                    obj.setGUILedStatus(obj.uiIDs.text_LED_veg_inputs_file, obj.ledKo, 0);
                    
                else
                    
                    if exist(filename,'file')
                        
                        obj.setGUILedStatus(obj.uiIDs.text_LED_veg_inputs_file, obj.ledOk, 0);
                    else
                        obj.setGUILedStatus(obj.uiIDs.text_LED_veg_inputs_file, obj.ledCk, 0);
                    end
                    
                end
            end
            
        end 
        
    end
    
end

