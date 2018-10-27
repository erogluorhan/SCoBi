% TO-DO: Check comments, copyrights, etc.
classdef SCoBiGUIManagers < handle
    %SCOBIGUIMANAGERS This class implements handles for GUI elements, and 
    % performs the GUI-related functionalities
    
        
    
    %% PROPERTIES
    
    %% CONSTANT PROPERTIES
    properties(Constant)
        
        % Colors
        colorDisable = [0.502 0.502 0.502];   % Grey (disabled color)
        colorEnable = [0 0 0];                % Black (enabled color)
    
    end
    
    
    %% PROTECTED PROPERTIES
    properties(GetAccess = 'protected', SetAccess = 'protected')
        
        % Simulator ID to handle simulator type specific functionalities
        simulator_id
        
        % GUI handles
        handles = [];	
        
        
        % UI STATES
        curState = [];      % This array [n x 1] contains the current status of abilitation of each element of the interface
        newState = [];      % This array [n x 1] contains the future status of abilitation of each element of the interface
        initialState = [];  % This array [n x 1] contains a saved status of abilitation of each element of the interface
        
        getFlag = [];       % Logical flag array [n x 1] contain 1 where a change of the value of the element is needed (from UI to current status)
        setFlag = [];       % Logical flag array [n x 1] contain 1 where a change of the value of the element is needed (from current status to UI)
       
        curVal = {};        % contains the value of each element of the interface, if one of the element has no value the cell is empty
        newVal = {};        % contains the value of the element of the interface that require to be changed
        
        
        initialized = 0;
           
    end
    
    
    %% PUBLIC GET, PROTECTED SET PROPERTIES
    properties(GetAccess = 'public', SetAccess = 'protected')
        
        % TO ACCESS UI ELEMENTS 
        uiIDs;               % Structure containing all the UI identifiers
        
        uiGroups;            % Structure containing groups of elements  
        
        uiPointers;          % Array n x 2 containing at the uiIDs index the handle of the object in the GUI
    end
    
    
        
    %% METHODS
    
    %% ABSTRACT METHODS
    methods( Access = 'public', Abstract = true )
        
        % Function to return values 
        outputFun(obj)
        
        % When an element is modified in the GUI, this must be called!
        % Sync the internal copy of the interface with the GUI
        syncFromGUI(obj, idEl)
        
    end
    
    methods( Access = 'protected', Abstract = true )
        
        
        % Test logical UI dependencies (a flag that activate other fields)
        checkUIdependencies(obj)
        
        % Get enable / disable status of the element of the interface
        % This function should be called only once, later in the code
        % the status is kept updated
        getAllElStatus(obj)
        
        % Get specific content in a GUI, if any
        getSpecificElContent(obj)
        
        % Initialize pop-up menus if exist in a GUI
        initPopupMenus(obj)
                
        % Initialize specific GUI elents of the sub-classes
        initSpecificGUI(obj)
        
        % Initialize UI IDs: Assign an id (integer) to each UI element
        initUIAccess(obj)
        
        % Set specific content in a GUI, if any
        setSpecificElContent(obj) 
    
    end
    
    
    %% STATIC METHODS
    methods(Static, Access = protected)
     
        % Get enabled status of a generic element
        function state = isGuiElOn(hObject)
            state = strcmp(get(hObject, 'Enable'),'on');
        end
        
        
        % Get enabled status of a panel element
        function state = isGuiPanelOn(hObject)
            state = isequal(get(hObject, 'ForegroundColor'), SCoBiGUIManagers.colorEnable);
        end
        
        
        % Close box
        function closeGUI(src,evnt)
            
            selection = questdlg('Are you sure to exit?',...
                '',...
                'Yes','No','Yes');
            
            switch selection
                case 'Yes'
                    delete(gcf);
                case 'No'
                    return
            end
        end
        
    end
    
    
    %% COMMON METHODS
    methods
        
        % Constructor
        function obj = SCoBiGUIManagers( handles, simulator_id )            
                        
            obj.init( handles, simulator_id );
            
        end 

        
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
        
        function openAboutDialog(obj)
            
           %get display size
            screenSize = get(0, 'ScreenSize');

            %calculate the center of the display
            position = [0, 0, 500, 370];
            position(1) = round( ( screenSize(3) - position(3) ) / 2 );
            position(2) = round( ( screenSize(4) - position(4) ) / 2 );

            d = dialog('Position', position, 'Name', 'About SCoBi Simulator');

            txtSCoBi = uicontrol('Parent',d,...
                                 'Style','text',...
                                 'HorizontalAlignment', 'left', ...
                                 'Position',[20 325 480 25],...
                                 'String', 'SCoBi - Signals of Opportunity Coherent Bistatic Scattering Simulator');

            txtVersion = uicontrol('Parent',d,...
                                 'Style','text',...
                                 'HorizontalAlignment', 'left', ...
                                 'Position',[20 280 480 25],...
                                 'String',strcat( 'Version: ', Constants.VERSION ));

            txtDevelopers = uicontrol('Parent',d,...
                                 'Style','text',...
                                 'HorizontalAlignment', 'left', ...
                                 'Position',[20 235 480 25],...
                                 'String', 'Developers: Orhan Eroglu, Mehmet Kurum, Dylan R. Boyd');
                             
            txtReferences = uicontrol('Parent',d,...
                                 'Style','text',...
                                 'HorizontalAlignment', 'left', ...
                                 'Position',[20 195 480 20],...
                                 'String', 'References:');

            txtReference1 = uicontrol('Parent',d,...
                                 'Style','text',...
                                 'HorizontalAlignment', 'left', ...
                                 'Position',[20 145 480 45],...
                                 'String', '[1] M. Kurum, M. Deshpande, A. T. Joseph, P. E. O’Neill, R. Lang, and O. Eroglu, “SCoBi-Veg: A generalized bistatic scattering model of reflectometry from vegetation for Signals of Opportunity applications,” IEEE Trans. Geosci. Remote Sensing, Press.');

            txtReference2 = uicontrol('Parent',d,...
                                 'Style','text',...
                                 'HorizontalAlignment', 'left', ...
                                 'Position',[20 95 480 45],...
                                 'String', '[2] M. Kurum, M. Deshpande, A. T. Joseph, P. E. O’Neill, R. H. Lang, and O. Eroglu, “Development of a coherent bistatic vegetation model for Signal of Opportunity applications at VHF/UHF-bands,” in Geoscience and Remote Sensing Symposium (IGARSS), 2017');

            txtReference3 = uicontrol('Parent',d,...
                                 'Style','text',...
                                 'HorizontalAlignment', 'left', ...
                                 'Position',[20 45 480 45],...
                                 'String', '[3] O. Eroglu, D. R. Boyd, and M. Kurum, “SCoBi: A Free, Open-Source, SoOp Coherent Bistatic Scattering Simulator Framework," IEEE Geosci. and Remote Sensing Mag., Review.');

            btnOK = uicontrol('Parent',d,...
                       'Position',[240 15 40 25],...
                       'String','OK',...
                       'Callback','delete(gcf)'); 
                   
        end
        
        function openDocumentsDialog(obj)
            
           % Create a dialog and display the URL   
            
            %get display size
            screenSize = get(0, 'ScreenSize');

            %calculate the center of the display
            position = [0, 0, 300, 150];
            position(1) = round( ( screenSize(3) - position(3) ) / 2 );
            position(2) = round( ( screenSize(4) - position(4) ) / 2 );

            d = dialog('Position', position, 'Name', 'SCoBi Documents');

            btnUserManual = uicontrol('Parent',d,...
               'Position',[50 110 200 25],...
               'String','SCoBi User''s Manual',...
               'Callback', strcat('web(''', Constants.URL_SCOBI_USER_MANUAL, ''', ''-browser'')'));


            btnDeveloperManual = uicontrol('Parent',d,...
               'Position',[50 80 200 25],...
               'String','SCoBi Developer''s Manual',...
               'Callback', strcat('web(''', Constants.URL_SCOBI_DEVELOPER_MANUAL, ''', ''-browser'')'));


            btnDesign = uicontrol('Parent',d,...
               'Position',[50 50 200 25],...
               'String','SCoBi Design (Enterprise Architect)',...
               'Callback', strcat('web(''', Constants.URL_SCOBI_DESIGN, ''', ''-browser'')'));

            btnOK = uicontrol('Parent',d,...
                       'Position',[130 20 40 25],...
                       'String','OK',...
                       'Callback','delete(gcf)'); 
                   
        end

        
        % Set content of the element of the interface
        function setAllElContent(obj)
            obj.setFlag(1) = false; % the figure doesn't change its status;
            
            if (sum(obj.setFlag) > 0)
                
                % Set specific content in a GUI, if any
                obj.setSpecificElContent();
                
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
        
        
        % Set a string of an element of the interface
        function setGuiElVisibility( obj, idEl, isVisible )
            
             set( obj.uiPointers(idEl), 'Visible', isVisible );
        end
                
    end
    
    
    %% PROTECTED METHODS
    methods(Access = 'protected')
        
        % Initialize the instance
        function init( obj, handles, simulator_id )
            
            obj.simulator_id = simulator_id; 
            
            obj.handles = handles;
            
            % Choose default command line output for gui_goGPS_unix
            obj.handles.output = obj.handles.panel_main;

            set(obj.handles.panel_main,'CloseRequestFcn',@obj.closeGUI);            

            % Update handles structure
            guidata(obj.handles.panel_main, obj.handles);

            % pixels
            set(obj.handles.panel_main, 'Units', 'pixels');
            
            %get display size
            screenSize = get(0, 'ScreenSize');
            
            %calculate the center of the display
            position = get(obj.handles.panel_main, 'Position');
            position(1) = round( ( screenSize(3) - position(3) ) / 2 );
            position(2) = round( ( screenSize(4) - position(4) ) / 2 );
            
            %center the window
            set(obj.handles.panel_main, 'Position', position);

            % Init elements ids
            obj.initGUI();
        
        end
        
        
        % Initialize GUI and get the status of it
        function initGUI(obj)

            % Set value for elements ids
            obj.initUIAccess();       

            % Read interface status (initialize structures
            obj.getAllElStatus();
            
            drawnow;
            
            % Fill pop up menus
            obj.initPopupMenus(); % Popup are also modified / reloaded in importStateMatlab
            
            % Initialize specific GUI elemnts of sub-classes, if any
            obj.initSpecificGUI();
            
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

        
        % Get content of the element of the interface 
        function getAllElContent(obj)
            obj.getFlag( obj.uiIDs.panel_main) = false; % the figure doesn't change its status;
            
            if sum( obj.getFlag > 0 ) 
                
                obj.getSpecificElContent();
                
                % Update newVal if the size of it is different from curVal
                % It may appen in the initialization process
                if sum(size(obj.newVal) == size(obj.curVal)) < 2
                    obj.newVal = obj.curVal;
                end
                
                obj.newVal(obj.getFlag) = obj.curVal(obj.getFlag);
                obj.getFlag = false(size(obj.getFlag()));
            end

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
        
        % Enable/Disable a generic element of the interface
        function onoffGuiEl(obj, hObject, state)
            if nargin < 2
                state = 'on';
            end
            set(hObject, 'Enable', state);
        end

        
        % Enable/Disable a panel element of the interface
        function onoffGuiPanel(obj, hObject, state)
            if nargin < 2
                state = 'on';
            end
            if strcmp(state,'off')
                set(hObject, 'ForegroundColor', SCoBiGUIManagers.colorDisable);
            else
                set(hObject, 'ForegroundColor', SCoBiGUIManagers.colorEnable);
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
        
                
        % Set the value of an element
        %  - idEl           is the identifier of the object (or group of
        %                   objects (see idUI/idGroup for the list of ids)
        %  - status          could be 0/1 'on'/'off'
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
        
        % Set new enable / disable status
        % Show all the new values stored in the internal state on the GUI
        % Get new values from the GUI
        function updateGUI(obj)
            obj.onoffUIEl();
            obj.getAllElContent();
            obj.checkUIdependencies();
            obj.setAllElContent();
        end
        
    end
    
end

