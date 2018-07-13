% TO-DO: Check comments, copyrights, etc.
classdef SCoBiGUIManager < handle
    %SCOBIGUIMANAGER This class implements handles for GUI elements, and 
    % performs the GUI-related functionalities
    
    properties(Constant)
        
        %% COLORS
        colorDisable = [0.502 0.502 0.502];   % Grey (disabled color)
        colorEnable = [0 0 0];                % Black (enabled color)
    
    end
    
    properties(GetAccess = 'private', SetAccess = 'private')
        
        handles = [];	% GUI handles
        
        
        %% UI STATES
        curState = [];      % This array [n x 1] contains the current status of abilitation of each element of the interface
        newState = [];      % This array [n x 1] contains the future status of abilitation of each element of the interface
        initialState = [];  % This array [n x 1] contains a saved status of abilitation of each element of the interface
        
        getFlag = [];       % Logical flag array [n x 1] contain 1 where a change of the value of the element is needed (from UI to current status)
        setFlag = [];       % Logical flag array [n x 1] contain 1 where a change of the value of the element is needed (from current status to UI)
       
        curVal = {};        % contains the value of each element of the interface, if one of the element has no value the cell is empty
        newVal = {};        % contains the value of the element of the interface that require to be changed
        
        
        initialized = 0;
        
        
        %% INPUT FILES        
        defaultInputFileName = 'default_input.xml';   % Name of the file containing default_settings
        lastInput = 'last_input.mat';         % Name of the file containing last settings

        
        %% LEDs    
        ledOff = 0;          % Led status off
        ledOn  = 1;          % Led status on
        ledOk  = 3;          % Led status ok
        ledKo  = 4;          % Led status ko
        ledCk  = 5;          % Led status check
        ledOp  = 6;          % Led status optional parameter
           
    end
    
    properties(GetAccess = 'public', SetAccess = 'private')
        
        %% TO ACCESS UI ELEMENTS 
        uiIDs;               % Structure containing all the UI identifiers
        uiGroups;            % Structure containing groups of elements        
        
        uiPointers;          % Array n x 2 containing at the uiIDs index the handle of the object in the GUI
    end
    
    methods
        
        % Constructor
        function obj = SCoBiGUIManager( handles )            
            
            obj.init( handles );
            
        end 
        
%         function handleGUIAction( obj, hObject, uiID )
%             
%             if strcmp( get( hObject, 'Style' ), 'pushbutton' )
%                 obj.handlePushButton( hObject, uiID );
%             elseif strcmp( get( hObject, 'Style' ), 'popupmenu' )
%                 obj.handlePopupMenu( hObject, uiID );
%                 %set( obj.handles.popup_gnd_cover, 'Enable','off') 
%             elseif strcmp( get( hObject, 'Style' ), 'edit' )
%                 obj.handleLineEdit( hObject, uiID );                
%             elseif strcmp( get( hObject, 'Style' ), 'checkbox' )
%                 obj.handleCheckBox( hObject, uiID );
%             end
%             
%         end
%         
%         function handlePushButton( obj, hObject, uiID )
%                         
%             if uiID == obj.uiIDs.pb_load_inputs
%                 
%                 [file, path] = uigetfile;
%                 
%                 loadInputsFromFile(path, file);
%                 
%             elseif uiID == obj.uiIDs.pb_save_inputs
%                  
%             elseif uiID == obj.uiIDs.pb_SCoBi
%                  
%             elseif uiID == obj.uiIDs.pb_exit
%                 obj.closeGUI;
%             end
%             
%         end

        
        % Return the status of activation of the element with id = idEl
        function isOn = isActive(obj, idEl)
            
            if ischar(obj.getElVal(idEl))
                
                if strcmp(obj.getElVal(idEl),'on')
                    isOn = true;
                else
                    isOn = false;
                end
                
            else
                isOn = logical(obj.getElVal(idEl));
            end
            
            isOn = isOn & obj.isEnabled(idEl);  % To be active an elment must be Enabled and its value = 'on'
        end

        
        % Return the status of abilitation of the element with id = idEl
        function isOn = isEnabled(obj, idEl)
            isOn = obj.newState(idEl);
        end


        % Get sim_mode
        function result = is_popup_sim_mode_snapshot(obj)
            isOn = obj.isEnabled( obj.uiIDs.popup_sim_mode );
            result = isOn && (obj.getElVal(obj.uiIDs.popup_sim_mode) == Constants.id_snapshot );
        end

        
        % Set content of the element of the interface
        function setAllElContent(obj)
            obj.setFlag(1) = false; % the figure doesn't change its status;
            
            if (sum(obj.setFlag) > 0)
%                 if (obj.setFlag(obj.uiIDs.sUPass))
%                     % Password field is the only one to be managed independently
%                     obj.setFlag(obj.uiIDs.sUPass) = false;
%                     obj.showPassword();
%                 end
                
                idEl = 1:length(obj.uiPointers);              % All the elements
                
                % Sets of panels
                panels = false(length(obj.uiPointers),1);     % init logical panels group
                panels(obj.uiGroups.gPanels) = true;           % set logical indexes of the panels
                % modified panels
                mPanel = idEl(panels & obj.setFlag);
                
                % Sets of text elements
                textEl = false(length(obj.uiPointers),1);     % init logical text elements group
                textEl(obj.uiGroups.strEl) = true;              % set logical indexes of the text elements
                % Modified LED elements
                mLED = obj.uiGroups.gLED(obj.setFlag(obj.uiGroups.gLED));
                % Modified text elements
                mTextEl = setdiff(idEl(textEl & obj.setFlag), mLED);
                
                mIdEl = setdiff(idEl(~panels & ~textEl & obj.setFlag), mLED); % id of elements that are not panels nor text elements that have been modified
                
                % For each modified panel
                for i=1:length(mPanel)
                    obj.setGuiElTitle(obj.uiPointers(mPanel(i)), obj.newVal{mPanel(i)});
                end
                
                % For each modified txt element
                for i=1:length(mTextEl)
                    obj.setGuiElStr(obj.uiPointers(mTextEl(i)), obj.newVal{mTextEl(i)});
                end
                
                % For all the LEDs
                for i=1:length(mLED)
                    obj.setGuiElColor(obj.uiPointers(mLED(i)), obj.newVal{mLED(i)});
                end
                
                % For all the other modified elements
                for i=1:length(mIdEl)
                    obj.setGuiElVal(obj.uiPointers(mIdEl(i)), obj.newVal{mIdEl(i)});
                end
                
                % Save the status in the object
                obj.curVal(obj.setFlag) = obj.newVal(obj.setFlag);
                obj.setFlag = false(size(obj.setFlag()));
            end
        end 
        

        % Set a value of an element of the interface
        function setGuiElColor(obj, hObject, color)
            set(hObject, 'ForegroundColor', color);
        end
        
        
        % Set a string of an element of the interface
        function setGuiElStr(obj, hObject, str)
            set(hObject, 'String', str);
        end
        
        
        % Get a title from an element of the interface
        function setGuiElTitle(obj, hObject, str)
            set(hObject, 'Title', str);
        end 
        

        % Set a value of an element of the interface
        function setGuiElVal(obj, hObject, value)
            set(hObject, 'Value', value);
        end
        
        
        % EVENT MANAGER
        % When an element is modified (and launch a callback function in
        % the GUI) this function must be called!
        % Sync the internal copy of the interface with the GUI
        function syncFromGUI(obj, idEl)
            
            if (nargin == 1)
                idEl = obj.uiIDs.popup_sim_mode;
            end
            
            obj.getFlag(idEl) = true;
            
            % Read all the values of the elements
            obj.getAllElContent();
            obj.setAllElContent();

            
          %% SIMULATION SETTINGS
            if sum(intersect(idEl, obj.uiIDs.popup_sim_mode)) > 0
                if obj.is_popup_sim_mode_snapshot()
                    obj.setElStatus(obj.uiIDs.popup_gnd_cover, 1, 0);
                    % Trigger the next popup menu
                    idEl = unique([idEl obj.uiIDs.popup_gnd_cover]);
                else
                    obj.setElStatus([obj.uiIDs.popup_veg_method], 1, 0);
                    % Trigger the next popup menu
                    idEl = unique([idEl obj.uiIDs.popup_veg_method]);
                end
            end
            
%             % I'm in real time, and this popup menu is active
%             if sum(intersect(idEl, obj.uiIDs.lCaptMode)) > 0
%                 if obj.isRealTime()
%                     obj.setElStatus(obj.uiGroups.Fig, 0, 0);
%                     switch obj.getElVal(obj.uiIDs.lCaptMode)
%                         case obj.idNav
%                             obj.setElStatus(obj.uiGroups.onRT_Nav, 1, 0);
%                         case obj.idRMon
%                             obj.setElStatus(obj.uiGroups.onRT_RMon, 1, 0);
%                         case obj.idMMon
%                             obj.setElStatus(obj.uiGroups.onRT_MMon, 1, 0);
%                         case obj.idRMMon
%                             obj.setElStatus(obj.uiGroups.onRT_RMMon, 1, 0);
%                     end
%                 end
%                 obj.updateGUI();
%             end
%             
%             % I'm in post processing and these popup menus are active
%             if sum(intersect(idEl, obj.uiIDs.lAlgType)) > 0
%                 obj.initProcessingType();
%                 obj.getFlag(obj.uiIDs.lProcType) = true;
%                 obj.setElStatus([obj.uiIDs.lProcType], 1, 0);
%                 % Trigger the next popup menu
%                 idEl = unique([idEl obj.uiIDs.lProcType]);
%                 obj.getAllElContent();
%             end        
%             
%             if sum(intersect(idEl, obj.uiIDs.lProcType)) > 0
%                 if obj.isKF() && (obj.getElVal(obj.uiIDs.lProcType) == obj.idCP_DD_MR)
%                     obj.setRinex();
%                 end
% 
%                 % Enable / Disable elements
%                 if obj.isPostProc()
%                     if obj.isLS()
%                         obj.setElStatus(obj.uiGroups.Fig, 0, 0);
%                         switch obj.getElVal(obj.uiIDs.lProcType)
%                             case obj.idC_SA
%                                 obj.setElStatus(obj.uiGroups.onPP_LS_C_SA, 1, 0);
%                             case obj.idC_DD
%                                 obj.setElStatus(obj.uiGroups.onPP_LS_C_DD, 1, 0);
%                             case obj.idCP_DD_L
%                                 obj.setElStatus(obj.uiGroups.onPP_LS_CP_DD_L, 1, 0);
%                             case obj.idCP_Vel
%                                 obj.setElStatus(obj.uiGroups.onPP_LS_CP_Vel, 1, 0);
%                             case obj.idC_SA_MR
%                                 obj.setElStatus(obj.uiGroups.onPP_LS_C_SA_MR, 1, 0);
%                             case obj.idCP_DD_MR
%                                 obj.setElStatus(obj.uiGroups.onPP_LS_CP_DD_MR, 1, 0);
%                         end
%                     end
%                     if obj.isKF()
%                         obj.setElStatus(obj.uiGroups.Fig, 0, 0);
%                         switch obj.getElVal(obj.uiIDs.lProcType)
%                             case obj.idC_SA
%                                 obj.setElStatus(obj.uiGroups.onPP_KF_C_SA, 1, 0);
%                             case obj.idC_DD
%                                 obj.setElStatus(obj.uiGroups.onPP_KF_C_DD, 1, 0);
%                             case obj.idCP_SA
%                                 obj.setElStatus(obj.uiGroups.onPP_KF_CP_SA, 1, 0);
%                             case obj.idCP_DD
%                                 obj.setElStatus(obj.uiGroups.onPP_KF_CP_DD, 1, 0);
%                             case obj.idCP_DD_MR                                
%                                 obj.setElStatus(obj.uiGroups.onPP_KF_CP_DD_MR, 1, 0);
%                         end
%                     end
%                     obj.getFlag(idEl) = true;
%                     obj.updateGUI();
%                 end
%                 
%                 % Verify if there's still something to enable/disable
%                 obj.onoffUIEl();
%                 % Check dependencies
%                 obj.checkUIdependencies();
%             end         
%             
%           %   INPUT/OUTPUT FILE AND FOLDERS
%           % ---------------------------------------------------------------
% 
%             % Browse for rover file
%             if sum(intersect(idEl, obj.uiIDs.bINI)) > 0
%                 obj.browseINIFile();
%             end
%             if sum(intersect(idEl, obj.uiGroups.gINI)) > 0
%                 obj.forceINIupdate();
%             end
%                 
%             % Browse output foder fo binary data
%             if sum(intersect(idEl, obj.uiIDs.bDirGoOut)) > 0
%                 obj.browseOutDir()
%             end
% 
%           %   SETTINGS - MASTER STATION
%           % ---------------------------------------------------------------
% 
%             % Toggle show password
%             if sum(intersect(idEl, obj.uiIDs.bUPass)) > 0
%                 obj.showPassword();
%             end

            
          %% PUSH BUTTONS                        
            % on Load Settings
            if sum(intersect(idEl, obj.uiIDs.pb_load_inputs)) > 0
%                 obj.loadState();
                [file, path] = uigetfile;
                
                obj.loadGUIFromInputFile( file, path )
            end
            
            % on Save Settings
            if sum(intersect(idEl, obj.uiIDs.pb_save_inputs)) > 0
%                 obj.saveState();
            end
            
            % on GO
            if sum(intersect(idEl, obj.uiIDs.pb_SCoBi)) > 0
%                 obj.go();
            end
            
            obj.onoffUIEl();
            obj.checkUIdependencies();
            
            % on Exit
            if sum(intersect(idEl, obj.uiIDs.pb_exit)) > 0
                obj.closeGUI;
            end
        end
        
    end
    
    
    % Internal initialization functions
    methods(Access = 'private')
        
        % Initialize the instance
        function init( obj, handles )
            
            tic;
            
            obj.handles = handles;
            
            % Choose default command line output for gui_goGPS_unix
            obj.handles.output = obj.handles.panel_main;
            
            set(obj.handles.panel_main,'CloseRequestFcn',@obj.closeGUI);            
              
            % Update handles structure
            guidata(obj.handles.panel_main, obj.handles);
                        
            %pixels
            set(obj.handles.panel_main, 'Units', 'pixels');
            
            %get display size
            screenSize = get(0, 'ScreenSize');
            
            %calculate the center of the display
            position = get(obj.handles.panel_main, 'Position');
            position(1) = ( screenSize(3) - position(3) );
            position(2) = ( screenSize(4) - position(4) );
            
            %center the window
            set(obj.handles.panel_main, 'Position', position);

            % Init elements ids
            obj.initGUI();
            t0 = toc;
            fprintf('SCoBi GUI initialization completed in %.2f seconds\n', t0);
        
        end
        
        % Initialize GUI and get the status of it
        function initGUI(obj)
            
            %% GET GLOBAL DIRECTORIES
            dir_input_sys = Directories.getInstance.input_sys;

            % Set value for elements ids
            obj.initUIAccess();       

            % Read interface status (initialize structures
            obj.getAllElStatus();
            
            drawnow;
            
            % Fill pop up menus
            obj.initPopupMenus(); % Popup are also modified / reloaded in importStateMatlab
            
            % Read the last input file, if any
            lastInputFileName = [];
            if exist( [dir_input_sys obj.lastInput], 'file' )
                load(obj.lastInput);
                % obj.lastInput contains the lastInputFileName 
                % string that might be saved during the last SCoBi run
            end
               
            if ~isempty( lastInputFileName ) && exist( [strcat(dir_input_sys, '/') lastInputFileName],'file' )
                obj.loadGUIFromInputFile( lastInputFileName );
            elseif exist([strcat(dir_input_sys, '/') obj.defaultInputFileName],'file')
                obj.loadGUIFromInputFile( obj.defaultInputFileName);
            else
                waitfor(msgbox('WARNING: There is no default input, you should add input files or give inputs manually!'));
            end
            
            % Read interface status as get from file
            obj.initialState = obj.curState;            
            
            % Read current ids content
            obj.getFlag = true( size(obj.curState) );
            obj.getAllElContent();
                        
            % Disable interface
            % obj.disableAll();
            obj.checkUIdependencies();

            obj.initialized = 1;
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
            
%             % Rough pre-allocation
%             pointers = zeros(177,1);	% rough estimation of the handle array size
%                                     % at the end of the function it will contain the exact size
                        
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
          
            
          groupIDs.sim_settings = [id.panel_sim_settings id.text_sim_mode:id.cb_calc_diffuse_term];
          
          
          %% SIMULATION INPUTS PANEL ELEMENTS           
          i = i+1;        id.text_sim_name = i;             pointers(i) = obj.handles.text_sim_name;
          i = i+1;        id.edit_sim_name = i;             pointers(i) = obj.handles.edit_sim_name;
          i = i+1;        id.text_campaign = i;             pointers(i) = obj.handles.text_campaign;
          i = i+1;        id.edit_campaign = i;             pointers(i) = obj.handles.edit_campaign;
          i = i+1;        id.text_campaign_date = i;        pointers(i) = obj.handles.text_campaign_date;
          i = i+1;        id.edit_campaign_date = i;        pointers(i) = obj.handles.edit_campaign_date;
          i = i+1;        id.text_plot = i;                 pointers(i) = obj.handles.text_plot;
          i = i+1;        id.edit_plot = i;                 pointers(i) = obj.handles.edit_plot;
          i = i+1;        id.text_veg_plant = i;            pointers(i) = obj.handles.text_veg_plant;
          i = i+1;        id.edit_veg_plant = i;            pointers(i) = obj.handles.edit_veg_plant;
          i = i+1;        id.text_veg_method = i;           pointers(i) = obj.handles.text_veg_method;
          i = i+1;        id.popup_veg_method = i;          pointers(i) = obj.handles.popup_veg_method;
          i = i+1;        id.text_veg_vir_orientations = i;        pointers(i) = obj.handles.text_veg_vir_orientations;
          i = i+1;        id.popup_veg_vir_orientations = i;       pointers(i) = obj.handles.popup_veg_vir_orientations;
          i = i+1;        id.text_Nr = i;                   pointers(i) = obj.handles.text_Nr;
          i = i+1;        id.edit_Nr = i;                   pointers(i) = obj.handles.edit_Nr;
          i = i+1;        id.text_Nfz = i;                  pointers(i) = obj.handles.text_Nfz;
          i = i+1;        id.edit_Nfz = i;                  pointers(i) = obj.handles.edit_Nfz;          
            
          groupIDs.sim_inputs = [id.panel_sim_inputs id.text_sim_name:id.edit_Nfz];
          
          
          %% TRANSMITTER (Tx) INPUTS PANEL ELEMENTS           
          i = i+1;        id.text_f_MHz = i;                pointers(i) = obj.handles.text_f_MHz;
          i = i+1;        id.edit_f_MHz = i;                pointers(i) = obj.handles.edit_f_MHz;
          i = i+1;        id.text_MHz = i;                  pointers(i) = obj.handles.text_MHz;
          i = i+1;        id.text_rsat_km = i;              pointers(i) = obj.handles.text_rsat_km;
          i = i+1;        id.edit_rsat_km = i;              pointers(i) = obj.handles.edit_rsat_km;
          i = i+1;        id.text_km_rsat = i;              pointers(i) = obj.handles.text_km_rsat;
          i = i+1;        id.text_EIRP_dB = i;              pointers(i) = obj.handles.text_EIRP_dB;
          i = i+1;        id.edit_EIRP_dB = i;              pointers(i) = obj.handles.edit_EIRP_dB;
          i = i+1;        id.text_dB_EIRP = i;              pointers(i) = obj.handles.text_dB_EIRP;
          i = i+1;        id.text_pol_Tx = i;               pointers(i) = obj.handles.text_pol_Tx;
          i = i+1;        id.popup_pol_Tx = i;              pointers(i) = obj.handles.popup_pol_Tx;
            
          groupIDs.Tx_inputs = [id.panel_Tx_inputs id.text_f_MHz:id.popup_pol_Tx];
          
          
          %% RECEIVER (Rx) INPUTS PANEL ELEMENTS           
          i = i+1;        id.text_hr_m = i;                 pointers(i) = obj.handles.text_hr_m;
          i = i+1;        id.edit_hr_m = i;                 pointers(i) = obj.handles.edit_hr_m;
          i = i+1;        id.text_m_hr = i;                 pointers(i) = obj.handles.text_m_hr;
          i = i+1;        id.text_G0r_dB = i;               pointers(i) = obj.handles.text_G0r_dB;
          i = i+1;        id.edit_G0r_dB = i;               pointers(i) = obj.handles.edit_G0r_dB;
          i = i+1;        id.text_dB_G0r = i;               pointers(i) = obj.handles.text_dB_G0r;
          % TO-DO: Move hpbw, SLL, and XPL to the specific Rx type input
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
            
          groupIDs.Rx_inputs = [id.panel_Rx_inputs id.text_hr_m:id.text_deg_ant_pat_res_Rx];
          
          
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
        
        % Create different groups for each type of element
        % Element that are panels
        % Element with the main content stored in:
        %  - 'Value' field
        %  - 'String' field
        function [panel, strEl, valEl] = autoElClassification(obj, uiPointers)
            panel = [];
            strEl = [];
            valEl = [];
            for i = 1 : length(uiPointers)
                
                % A panel doesn't have a Style field
                el = get(uiPointers(i));
                
                if isfield(el,'Style')                    
                    
                    if sum(strcmp(el.Style, {'edit' 'text'})) > 0
                        strEl = [strEl i];
                    elseif sum(strcmp(el.Style, {'radiobutton', 'checkbox', 'pushbuttons', 'popupmenu'})) > 0
                        valEl = [valEl i];
                    end
                else
                    % A panel doesn't have the field Style
                    % but has the field BorderType
                    if isfield(el,'BorderType')
                        panel = [panel i];
                    end
                end
            end
        end
        
        
        % Function that runs inside onoffUIEl
        % Test every logical dependence in the GUI
        % E.g. a flag that activate other fields
        function checkUIdependencies(obj)
                        
%           %   INPUT FILE TYPE
%           % --------------------------------------------------------------- 
% 
%             % Check File input dependencies
%             obj.setElStatus( [obj.uiGroups.onRin], ~obj.isBin(), 0 );
%             
%             if obj.isRinex()
%                 if obj.isStandAlone()
%                     obj.setElStatus([obj.uiGroups.RinMaster], 0, 0);
%                 end
%             end
%             obj.setElStatus([obj.uiGroups.onBin], obj.isBin() && ~obj.isRealTime(), 0);
%             
%           %   OPTIONS
%           % --------------------------------------------------------------- 
% 
%             % Reference path file
%             isOn = obj.isActive(obj.uiIDs.cRefPath);
%             obj.setElStatus([obj.uiGroups.RefPath], isOn, 0);
%             obj.setElStatus([obj.uiIDs.cConstraint], isOn, 0);
%             
%             % Plot while processing flag
%             isOn = obj.isActive(obj.uiIDs.cPlotProc);
%             obj.setElStatus([obj.uiGroups.gPlotProc], isOn, 0);
%             
%             % If the master is not available disable the flag to plot it 
%             if ~obj.isEnabled(obj.uiIDs.fRinMaster)
%                 obj.setElStatus([obj.uiIDs.cPlotMaster], 0, 0);
%             end
%             
%             % If only code is used, no ambiguities are available
%             if obj.isProcessingType(obj.idC_SA) || obj.isProcessingType(obj.idC_DD)
%                 obj.setElStatus([obj.uiIDs.cPlotAmb], 0, 0);
%             end
% 
%             % Only for the variometric approach google Earth plot is disabled
%             if obj.isLS() && obj.isProcessingType(obj.idCP_Vel)
%                 obj.setElStatus([obj.uiIDs.cGEarth], 0, 0);
%             end
%             
%             % NTRIP flag
%             isOn = obj.isActive(obj.uiIDs.cUseNTRIP);
%             obj.setElStatus([obj.uiGroups.gNTRIP], isOn, 0);
%             
%           %   INTEGER AMBIGUITY RESOLUTION
%           % --------------------------------------------------------------- 
%           
%             % LAMBDA flag
%             isOn = obj.isActive(obj.uiIDs.cLAMBDA);
%             obj.setElStatus([obj.uiGroups.gLAMBDA], isOn, 0);
% 
%             % LAMBDA version check
%             if (obj.isLambda2())
%                 obj.setElStatus([obj.uiGroups.gLAMBDA3], 0, 0);
%             end
%             if (~obj.isLambdaIls())
%                 obj.setElStatus([obj.uiGroups.gLAMBDAILS], 0, 0);
%                 if (obj.isLambda3Par())
%                     obj.setElStatus([obj.uiIDs.tP0], 1, 0);
%                     obj.setElStatus([obj.uiIDs.nP0], 1, 0);
%                     obj.setElStatus([obj.uiIDs.cP0], 1, 0);
%                 end
%             end
%             
%             if (obj.isLambda3Par())
%                 set(obj.handles.tP0,'String','Min. success rate (P0):');
%             else
%                 set(obj.handles.tP0,'String','Fixed failure rate (P0):');
%             end
% 
%             % Automatic mu flag
%             isOn = obj.isEnabled(obj.uiIDs.cMu) && obj.isActive(obj.uiIDs.cMu);
%             if (isOn)
%                 obj.setElStatus([obj.uiIDs.nMu], ~isOn, 0);
%             end
%             
%             % Default P0 flag
%             isOn = obj.isEnabled(obj.uiIDs.cP0) && obj.isActive(obj.uiIDs.cP0);
%             if (isOn)
%                 obj.setElStatus([obj.uiIDs.nP0], ~isOn, 0);
%             end
% 
%           %   SETTINGS - KALMAN FILTER - STD
%           % --------------------------------------------------------------- 
% 
%             % Error standard deviation Phase
%             isOn = obj.isActive(obj.uiIDs.bStdPhase);
%             obj.setElStatus([obj.uiIDs.nStdPhase], isOn, 0);
%             obj.setElStatus([obj.uiIDs.uStdPhase], isOn, 0);
%             
%             % DTM toggle
%             isOn = obj.isActive(obj.uiIDs.bStdDTM);
%             obj.setElStatus([obj.uiGroups.gDTM], isOn, 0);
% 
%           %   SETTINGS - KALMAN FILTER
%           % --------------------------------------------------------------- 
% 
%             % Stop Go Stop
%             if obj.isEnabled(obj.uiIDs.cStopGoStop)
%                 isOn = obj.isActive(obj.uiIDs.cStopGoStop);
%                 obj.setElStatus(obj.uiGroups.pDynModel, ~isOn, 0);
%             end
%                         
%           %   SETTINGS - PORTS
%           % ---------------------------------------------------------------
% 
%             % Ports
%             if obj.isEnabled(obj.uiIDs.lnPorts)
%                 nPorts = obj.getElVal(obj.uiIDs.lnPorts);
%                 obj.setElStatus([obj.uiGroups.lPort0], nPorts > 0, 0); nPorts = nPorts-1;
%                 obj.setElStatus([obj.uiGroups.lPort1], nPorts > 0, 0); nPorts = nPorts-1;
%                 obj.setElStatus([obj.uiGroups.lPort2], nPorts > 0, 0); nPorts = nPorts-1;
%                 obj.setElStatus([obj.uiGroups.lPort3], nPorts > 0, 0);
%             end
%             
%           %   SETTINGS - KALMAN FILTER - DYNAMIC MODEL
%           % ---------------------------------------------------------------
% 
%             obj.resetDynModel();
%             
%             % Check if static dynamic model
%             if (obj.isDynModelStatic())
%                 obj.setElStatus([obj.uiGroups.pKF_ENU], 0, 0);
%             else
%                 if (obj.isKF())
%                     obj.setElStatus([obj.uiGroups.pKF_ENU], 1, 0);
%                 end
%             end
%                 
%           %   SETTINGS - MASTER STATION
%           % ---------------------------------------------------------------
%             if obj.isEnabled(obj.uiIDs.cMPos)
%                 % If the panel is enabled
%                 isOn = obj.isActive(obj.uiIDs.cMPos);
%                 % if isOn the rest should be disabled
%                 obj.setElStatus(obj.uiIDs.lCRS, ~isOn, 0);                
%                 isXYZ = obj.getElVal(obj.uiIDs.lCRS) == obj.idXYZ;
%                 obj.setElStatus([obj.uiGroups.gMXYZ], isXYZ && ~isOn, 0);
%                 obj.setElStatus([obj.uiGroups.gMGeodetic], ~isXYZ && ~isOn, 0);                
%             else
%                 obj.setElStatus([obj.uiGroups.pMSt], 0, 0);
%             end
% 
%             
%           %   SETTINGS - MASTER SERVER
%           % ---------------------------------------------------------------
%             
%             % Password
%             obj.showPassword();
%             
%           %   MODE
%           % --------------------------------------------------------------- 
%             
%           %  % Check list boxes
%           %  isOn = obj.isRealTime();
%           %  obj.setElStatus([obj.uiIDs.lCaptMode], isOn, 0);
%           %  obj.setElStatus([obj.uiIDs.lAlgType obj.uiIDs.lProcType], ~isOn, 1);
%             
%           %   GO BUTTON AND LEDS
%           % --------------------------------------------------------------- 
%           % For each file field enabled, I have to check the existence of
%           % the folder / file to enable the go Button
%             obj.onoffUIEl();
%             obj.updateLEDstate();
%             goOk = obj.test4Go();
%             obj.setElStatus([obj.uiIDs.bSave obj.uiIDs.bGo] , goOk, 1);
        end

        % Get content of the element of the interface 
        function getAllElContent(obj)
            obj.getFlag( obj.uiIDs.panel_main) = false; % the figure doesn't change its status;
            
            if sum( obj.getFlag > 0 )                
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
                
                % Update newVal if the size of it is different from curVal
                % It may appen in the initialization process
                if sum(size(obj.newVal) == size(obj.curVal)) < 2
                    obj.newVal = obj.curVal;
                end
                obj.newVal(obj.getFlag) = obj.curVal(obj.getFlag);
                obj.getFlag = false(size(obj.getFlag()));
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
        
        
        % Get the value of an element (from the internal copy of the object)
        function value = getElVal(obj, idEl)
            value = obj.newVal{idEl};
        end
        
        
        % Get a value from an element of the interface
        function val = getGuiElColor(obj, hObject)
            val = get(hObject, 'ForegroundColor');
        end
        
        
        % Get a string from an element of the interface
        function str = getGuiElStr(obj, hObject)
            str = get(hObject, 'String');
        end
        

        % Get a title from an element of the interface
        function str = getGuiElTitle(obj, hObject)
            str = get(hObject, 'Title');
        end
        
        
        % Get a value from an element of the interface
        function val = getGuiElVal(obj, hObject)
            val = get(hObject, 'Value');
        end
        
        % Load the state of the gui from a matlab file.
        function loadGUIFromInputFile( obj, filename, pathname )
            
            
            %% GET GLOBAL DIRECTORIES
            dir_input_sys = Directories.getInstance.input_sys;
                
            
            if nargin == 2
                pathname = dir_input_sys;
            end
            
%             load(filename); % the file contains the variable state
%             obj.status = state;
            
            [sim_settings, sim_inputs, Tx_inputs, Rx_inputs, gnd_inputs, veg_inputs, dyn_inputs] = obj.loadInputsFromFile(pathname, filename);
            
            
            %%   SIMULATION SETTINGS
            obj.init_popup_sim_mode();
            obj.setElVal(obj.uiIDs.popup_sim_mode, sim_settings.sim_mode_id, 0);
            
            obj.init_popup_gnd_cover();
            obj.setElVal(obj.uiIDs.popup_gnd_cover, sim_settings.gnd_cover_id, 0);
            
            obj.setElVal(obj.uiIDs.cb_write_attenuation, sim_settings.write_attenuation, 0);
            obj.setElVal(obj.uiIDs.cb_calc_direct_term, sim_settings.calc_direct_term, 0);
            obj.setElVal(obj.uiIDs.cb_calc_specular_term, sim_settings.calc_specular_term, 0);
            obj.setElVal(obj.uiIDs.cb_calc_diffuse_term, sim_settings.calc_diffuse_term, 0);
            
            
            %% SIMULATION INPUTS
            obj.setElVal(obj.uiIDs.edit_sim_name, sim_inputs.sim_name, 0);
            obj.setElVal(obj.uiIDs.edit_campaign, sim_inputs.campaign, 0);
            obj.setElVal(obj.uiIDs.edit_campaign_date, sim_inputs.campaign_date, 0);
            obj.setElVal(obj.uiIDs.edit_plot, sim_inputs.plot, 0);
            obj.setElVal(obj.uiIDs.edit_veg_plant, sim_inputs.veg_plant, 0);
            
            obj.init_popup_veg_method();
            obj.setElVal(obj.uiIDs.popup_veg_method, sim_inputs.veg_method_id, 0);
            
            obj.init_popup_veg_vir_orientation();
            obj.setElVal(obj.uiIDs.popup_veg_vir_orientations, sim_inputs.veg_vir_orientation_id, 0);
            
            obj.setElVal(obj.uiIDs.edit_Nr, sim_inputs.Nr, 0);
            obj.setElVal(obj.uiIDs.edit_Nfz, sim_inputs.Nfz, 0);
            
            
            %% TRANSMITTER (Tx) INPUTS
            obj.setElVal(obj.uiIDs.edit_f_MHz, Tx_inputs.f_MHz, 0);
            obj.setElVal(obj.uiIDs.edit_rsat_km, Tx_inputs.rsat_km, 0);
            obj.setElVal(obj.uiIDs.edit_EIRP_dB, Tx_inputs.EIRP_dB, 0);
            
            obj.init_popup_pol_Tx();
            obj.setElVal(obj.uiIDs.popup_pol_Tx, Tx_inputs.pol_Tx_id, 0);
            
            
            %% RECEIVER (Rx) INPUTS
            obj.setElVal(obj.uiIDs.edit_hr_m, Rx_inputs.hr_m, 0);
            obj.setElVal(obj.uiIDs.edit_G0r_dB, Rx_inputs.G0r_dB, 0);
            
            obj.init_popup_pol_Rx();
            obj.setElVal(obj.uiIDs.popup_pol_Rx, Rx_inputs.pol_Rx_id, 0);
            
            obj.init_popup_orientation_Rx();
            obj.setElVal(obj.uiIDs.popup_orientation_Rx, Rx_inputs.orientation_Rx_id, 0);
            
            obj.setElVal(obj.uiIDs.edit_th0_Rx, Rx_inputs.th0_Rx_deg, 0);
            obj.setElVal(obj.uiIDs.edit_ph0_Rx, Rx_inputs.ph0_Rx_deg, 0);
            
            obj.init_popup_ant_pat_Rx();
            obj.setElVal(obj.uiIDs.popup_ant_pat_Rx, Rx_inputs.ant_pat_Rx_id, 0);
            
            obj.setElVal(obj.uiIDs.edit_ant_pat_res_Rx, Rx_inputs.ant_pat_res_deg, 0);
            
            
            %% GROUND INPUTS
            obj.setElVal(obj.uiIDs.edit_sand_ratio, gnd_inputs.sand_ratio, 0);
            obj.setElVal(obj.uiIDs.edit_clay_ratio, gnd_inputs.clay_ratio, 0);
            obj.setElVal(obj.uiIDs.edit_rhob_gcm3, gnd_inputs.rhob_gcm3, 0);
            
            obj.init_popup_diel_model();
            obj.setElVal(obj.uiIDs.popup_diel_model, gnd_inputs.diel_model, 0);
            
            
            %% INPUT FILES  
            obj.setElVal(obj.uiIDs.edit_dyn_inputs_file, dyn_inputs.filename, 0);          
            obj.setElVal(obj.uiIDs.edit_ant_pat_Rx_file, Rx_inputs.ant_pat_file, 0);
            obj.setElVal(obj.uiIDs.edit_veg_inputs_file, veg_inputs.filename, 1);
            
            
%             %   SETTINGS - KALMAN FILTER - DYNAMIC MODEL
%             % ===============================================================            
%             obj.resetDynModel();
%             obj.setElVal(obj.uiIDs.lDynModel, state.dyn_mod, 0);

            
            % Check all the dependencies
            obj.syncFromGUI(obj.uiIDs.popup_sim_mode);
        end
        
        function [sim_settings, sim_inputs, Tx_inputs, Rx_inputs, gnd_inputs, veg_inputs, dyn_inputs] = loadInputsFromFile(obj, pathname, filename)
            
            if ~isempty(pathname) && ~isempty(filename)
                
                xDoc = xmlread( strcat( pathname, '\', filename ) );

                %% SIMULATION SETTINGS
                sim_mode = getStringFromXML(xDoc, ConstantNames.set_simMode);        % 
                is_sim_mode = cellfun(@(x)isequal(x, sim_mode), Constants.sim_modes );
                [~, sim_mode_id] = find(is_sim_mode);
                sim_settings.sim_mode_id = sim_mode_id;       % 
                
                gnd_cover = getStringFromXML(xDoc, ConstantNames.set_gndCover);        	% Bare-soil OR Vegetation
                is_gnd_cover = cellfun(@(x)isequal(x, gnd_cover), Constants.gnd_covers );
                [~, gnd_cover_id] = find(is_gnd_cover);
                sim_settings.gnd_cover_id = gnd_cover_id;  
                
                sim_settings.write_attenuation = getDoubleFromXML(xDoc, ConstantNames.set_writeAttenuation);      % Flag to write Attenuation to Excel file 
                sim_settings.calc_direct_term = getDoubleFromXML(xDoc, ConstantNames.set_calcDirectTerm);  % Flag to calculate direct term  
                sim_settings.calc_specular_term = getDoubleFromXML(xDoc, ConstantNames.set_calcSpecularTerm);  % Flag to calculate specular term     
                sim_settings.calc_diffuse_term = getDoubleFromXML(xDoc, ConstantNames.set_calcDiffuseTerm);       % Flag to calculate diffuse term 


                %% SIMULATION INPUTS
                sim_inputs.sim_name = getStringFromXML(xDoc, ConstantNames.sim_simName);  % Plot of the campaign  
                sim_inputs.campaign = getStringFromXML(xDoc, ConstantNames.sim_campaign);      % The campaign of the simulation 
                sim_inputs.campaign_date = getStringFromXML(xDoc, ConstantNames.sim_campaignDate);  % Date of the campaign  
                sim_inputs.plot = getStringFromXML(xDoc, ConstantNames.sim_plot);  % Plot of the campaign  
                
                veg_method = getStringFromXML(xDoc, ConstantNames.sim_vegMethod);        % 
                is_veg_method = cellfun(@(x)isequal(x, veg_method), Constants.veg_methods );
                [~, veg_method_id] = find(is_veg_method);
                sim_inputs.veg_method_id = veg_method_id; 
                             
                veg_vir_orientation = getStringFromXML(xDoc, ConstantNames.sim_vegVirOrientation);                      % If veg_method is Virtual, then is orientation Row-crop or Random-spread?   % 
                is_veg_vir_orientation = cellfun(@(x)isequal(x, veg_vir_orientation), Constants.veg_vir_orientations );
                [~, veg_vir_orientation_id] = find(is_veg_vir_orientation);
                sim_inputs.veg_vir_orientation_id = veg_vir_orientation_id; 
                
                sim_inputs.veg_plant = getStringFromXML(xDoc, ConstantNames.sim_vegetationPlant);  % The vegetation_plant that simulations to be run 
                sim_inputs.Nr = getDoubleFromXML(xDoc, ConstantNames.sim_Nr);       % Number of Realizations  
                sim_inputs.Nfz = getDoubleFromXML(xDoc, ConstantNames.sim_Nfz);    % Number of Fresnel Zones


                %% TRANSMITTER INPUTS
                Tx_inputs.f_MHz = getDoubleFromXML(xDoc, ConstantNames.Tx_f_MHz);      % Operating frequncy (MHz)
                Tx_inputs.rsat_km = getDoubleFromXML(xDoc, ConstantNames.Tx_rsat_km);    % Transmitter radius (km) 
                % TO-DO: This is for satGeometryManual. There should be an option for satGeometry
                Tx_inputs.EIRP_dB = getDoubleFromXML(xDoc, ConstantNames.Tx_EIRP_dB);  % Equivalent Isotropic Radiated Power
                               
                pol_Tx = getStringFromXML(xDoc, ConstantNames.Tx_pol_Tx);               % Transmitter polarization
                is_pol_Tx = cellfun(@(x)isequal(x, pol_Tx), Constants.polarizations );
                [~, pol_Tx_id] = find(is_pol_Tx);
                Tx_inputs.pol_Tx_id = pol_Tx_id; 


                %% RECEIVER INPUTS
                Rx_inputs.hr_m = getDoubleFromXML(xDoc, ConstantNames.Rx_hr_m);        % Antenna Height (m)
                Rx_inputs.G0r_dB = getDoubleFromXML(xDoc, ConstantNames.Rx_G0r_dB);    % Receive Antenna Gain (dB) 
                               
                pol_Rx = getStringFromXML(xDoc, ConstantNames.Rx_pol_Rx);               % Receiver polarization
                is_pol_Rx = cellfun(@(x)isequal(x, pol_Rx), Constants.polarizations );
                [~, pol_Rx_id] = find(is_pol_Rx);
                Rx_inputs.pol_Rx_id = pol_Rx_id;
                               
                orientation_Rx = getStringFromXML(xDoc, ConstantNames.Rx_orientation_Rx);                   % Receiver antenna orientation
                is_orientation_Rx = cellfun(@(x)isequal(x, orientation_Rx), Constants.Rx_orientation );
                [~, orientation_Rx_id] = find(is_orientation_Rx);
                Rx_inputs.orientation_Rx_id = orientation_Rx_id;
                
                Rx_inputs.th0_Rx_deg = getDoubleFromXML(xDoc, ConstantNames.Rx_th0_Rx_deg);    % Receive Antenna theta (incidence) angle
                Rx_inputs.ph0_Rx_deg = getDoubleFromXML(xDoc, ConstantNames.Rx_ph0_Rx_deg);    % Receive Antenna phi (azimuth) angle
                
                ant_pat_Rx = getStringFromXML(xDoc, ConstantNames.Rx_ant_pat_Rx);                   % Receive Antenna Pattern generation method
                is_ant_pat_Rx = cellfun(@(x)isequal(x, ant_pat_Rx), Constants.Rx_ant_pats );
                [~, ant_pat_Rx_id] = find(is_ant_pat_Rx);
                Rx_inputs.ant_pat_Rx_id = ant_pat_Rx_id;
                
                Rx_inputs.ant_pat_res_deg = getDoubleFromXML(xDoc, ConstantNames.Rx_ant_pat_res_deg);    % Receive Antenna Pattern file
                Rx_inputs.ant_pat_file = getStringFromXML(xDoc, ConstantNames.rec_ant_pat_file);  % Antenna polarization


                %% GROUND INPUTS
    %             gnd_inputs.VSM_list_cm3cm3 = getDoubleArrayFromXML(xDoc, ConstantNames.dyn_VSM_list_cm3cm3);        % Theta probe   
    %             gnd_inputs.RMSH_list_cm = getDoubleArrayFromXML(xDoc, ConstantNames.dyn_RMSH_list_cm);  % Surface rms height (cm)
                gnd_inputs.sand_ratio = getDoubleFromXML(xDoc, ConstantNames.gnd_sand_ratio);  % Sand ratio of the soil texture
                gnd_inputs.clay_ratio = getDoubleFromXML(xDoc, ConstantNames.gnd_clay_ratio);  % Clay ratio of the soil texture 
                gnd_inputs.rhob_gcm3 = getDoubleFromXML(xDoc, ConstantNames.gnd_rhob_gcm3);  % Soil bulk density 
                gnd_inputs.diel_model = getDoubleFromXML(xDoc, ConstantNames.gnd_diel_model);    % Receive Antenna Pattern file



                %% VEGETATION INPUTS
                veg_inputs.filename = getStringFromXML(xDoc, ConstantNames.veg_file);  % Vegetation inputs file name


                %% DYNAMIC INPUTS
                dyn_inputs.filename = getStringFromXML(xDoc, ConstantNames.dyn_input_file);  % Vegetation inputs file name
                
            end
            
                        
        end
        
        % Enable/Disable a generic element of the interface
        function onoffGuiEl(hObject, state)
            if nargin < 2
                state = 'on';
            end
            set(hObject, 'Enable', state);
        end

        % Enable/Disable a panel element of the interface
        function onoffGuiPanel(hObject, state)
            if nargin < 2
                state = 'on';
            end
            if strcmp(state,'off')
                set(hObject, 'ForegroundColor', SCoBiGUIManager.disableCol);
            else
                set(hObject, 'ForegroundColor', SCoBiGUIManager.enableCol);
            end
        end  
        
        % Enable / Disable various elements in the interface
        % the variable newStatus will decide the future status of the
        % interface (function also known as setGuiElStatus)
        function onoffUIEl(obj)
            % Detect modified state (logical array)
            idModified = xor(obj.curState, obj.newState);
            idModified(obj.uiIDs.panel_main) = false; % the main panel doesn't change its status;
            
            if (sum(idModified) > 0)
                obj.curState = obj.newState;
                
                % ids of the modified elements
                % state id of the modified element
                listState = uint8(obj.newState)+1;
                state = {'off','on'};
                
                idEl = 1:length(obj.uiPointers);              % All the elements
                
                % Sets of panels
                panels = false(length(obj.uiPointers),1);     % init logical panels group
                panels(obj.uiGroups.gPanels) = true;           % set logical indexes of the panels
                
                % modified panels
                mPanel = idEl(panels & idModified);
                
                mIdEl = idEl(~panels & idModified);          % id of elements that are not panels that have been modified
                
                % For each modified panel
                for i=1:length(mPanel)
                    obj.onoffGuiPanel(obj.uiPointers(mPanel(i)), state{listState(mPanel(i))});
                end
                
                % For all the other modified elements
                for i=1:length(mIdEl)
                    obj.onoffGuiEl(obj.uiPointers(mIdEl(i)), state{listState(mIdEl(i))});
                end                
            end
        end
        
        % Save the stati of the gui to a matlab file
        function saveGUIToInputFile(obj,filename)
             
            %   MODE
            % ===============================================================
            state.mode              = obj.getElVal(obj.uiIDs.lProcMode);
            state.nav_mon           = obj.getElVal(obj.uiIDs.lCaptMode);
            state.kalman_ls         = obj.getElVal(obj.uiIDs.lAlgType);
            state.code_dd_sa        = obj.getElVal(obj.uiIDs.lProcType);
            
            save(filename, 'state');
        end
        
                
        % Set the value of an element
        %  - idEl           is the identifier of the object (or group of
        %                   objects (see idUI/idGroup for the list of ids)
        %  - value          could be 0/1 'on'/'off'
        %  - <autoapply>    if set obj.setAllElContent() is automatically
        %                   called after to show immediatly the modification 
        %                   in the GUI
        function setElStatus(obj, idEl, status, autoapply)
            if nargin == 3
                autoapply = true;
            end
            if ischar(status)
                if strcmp(status,'on')
                    status = true;
                else 
                    status = false;
                end
            end
            obj.newState(idEl) = logical(status);
            if autoapply
                obj.onoffUIEl();
            end
        end
        
        
        % Set the value of an element
        %  - idEl           is the identifier of the object 
        %                   (see uiIDs for the list of ids)
        %  - value          could be of any type string/boolean/number
        %  - <autoapply>    if set, obj.setAllElContent() is automatically
        %                   called after to show immediatly the modification 
        %                   in the GUI
        function setElVal(obj, idEl, value, autoapply)
           
            if nargin == 3
                autoapply = true;
            end
            
            obj.setFlag(idEl) = true;
            obj.newVal{idEl} = value;
            
            if autoapply
                obj.setAllElContent();
            end
            
        end
        
        
        % Led status
        % Set green / red status of the UI and optionally lock the UI
        % -------------------------------------------------------------------------
        function setGUILedStatus(obj, idEl, status, autoapply)
            
            if nargin == 3
                autoapply = true;
            end
            
            if (status == obj.ledOk)    % Led Ok
                if ~obj.isColor(idEl, obj.green)
                    obj.setElVal(idEl, obj.green)
                end
            elseif (status == obj.ledKo) % Led Ko
                if ~obj.isColor(idEl, obj.red)
                    obj.setElVal(idEl,obj.red)
                end
            elseif (status == obj.ledCk) % Led Check
                if ~obj.isColor(idEl, obj.yellow)
                    obj.setElVal(idEl,obj.yellow)
                end
            elseif (status == obj.ledOp) % Led Optional parameter
                if ~obj.isColor(idEl, obj.blue)
                    obj.setElVal(idEl,obj.blue)
                end
            end
            obj.setAllElContent();
            %obj.setElStatus(idEl, status > 1, autoapply)
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
    
    
    methods(Static, Access = private)
        
               
     
        % Get enabled status of a generic element
        function state = isGuiElOn(hObject)
            state = strcmp(get(hObject, 'Enable'),'on');
        end
        
        
        % Get enabled status of a panel element
        function state = isGuiPanelOn(hObject)
            state = isequal(get(hObject, 'ForegroundColor'), SCoBiGUIManager.colorEnable);
        end
        
        
        % Close box
        function closeGUI(src,evnt)
            global scobiGUI
            
            selection = questdlg('Are you sure to exit?',...
                '',...
                'Yes','No','Yes');
            
            switch selection
                case 'Yes'
                    if isempty(scobiGUI)  % Something went wrong in the initialization phase
                        delete(gcf);
                    else
                        if isfield(scobiGUI.handles, 'panel_main')
                            delete(scobiGUI.handles.panel_main);
                        else % Something went wrong in the initialization phase
                            delete(gcf);
                        end
                    end
                case 'No'
                    return
            end
        end
        
    end
    
end

