% TO-DO: Check comments, copyrights, etc.
classdef AntPatRxGgGUIManager < SCoBiGUIManagers
    %ANTPATRXGGGUIMANAGER This class implements handles for GUI elements, and 
    % performs the GUI-related functionalities
    
    
    properties(GetAccess = 'private', SetAccess = 'private')
               
        % Half-power beamwidth in degrees
        hpbw_deg
        
        % Side-lobe level in decibels
        SLL_dB
        
        % Cross-polarization level in decibels
        XPL_dB
           
    end
    
    
    %% METHODS
    
    %% COMMON METHODS
    methods
        
        % Constructor
        function obj = AntPatRxGgGUIManager( handles )            
            
            % Call superclass constructor
            obj = obj@SCoBiGUIManagers( handles );
            
        end
        
        
        function funout = outputFun(obj)
            
            funout{1} = obj.hpbw_deg;
            funout{2} = obj.SLL_dB;
            funout{3} = obj.XPL_dB;
            
        end
        
        
        % EVENT MANAGER
        % When an element is modified in the GUI, this must be called!
        % Sync the internal copy of the interface with the GUI
        function syncFromGUI(obj, idEl)
            
            if (nargin == 1)
                idEl = obj.uiIDs.text_description;
            end
            
            obj.getFlag(idEl) = true;
            
            % Read all the values of the elements
            obj.getAllElContent();
            obj.setAllElContent();
          
            
            %% PUSH BUTTONS                        
            % on Load Settings
            if sum(intersect(idEl, obj.uiIDs.pb_ok)) > 0
                
                obj.hpbw_deg = obj.getElVal( obj.uiIDs.edit_hpbw_deg );
                
                obj.SLL_dB = obj.getElVal( obj.uiIDs.edit_SLL_dB );
                
                obj.XPL_dB = obj.getElVal( obj.uiIDs.edit_XPL_dB );
	
                uiresume(obj.handles.panel_main);
                
            end
            
            % on SCoBi
            if sum(intersect(idEl, obj.uiIDs.pb_cancel)) > 0
                
                obj.hpbw_deg = [];
                
                obj.SLL_dB = [];
                
                obj.XPL_dB = [];
	
                uiresume(obj.handles.panel_main);
                
            end
            
            obj.onoffUIEl();
            obj.checkUIdependencies();
            
            % on Exit
            if sum(intersect(idEl, obj.uiIDs.pb_cancel)) > 0
%                 obj.closeGUI(obj);
                delete(gcf);
            end
            
        end
        
    end
    
    
    % Internal initialization functions
    methods(Access = 'protected')
        
        % Initialize the instance
        function init( obj, handles )
            
            % Call superclass's init() function
            init@SCoBiGUIManagers( obj, handles );
        
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
          
          
          %% EDIT FIELDS
          i = i+1;        id.text_hpbw_deg = i;                 pointers(i) = obj.handles.text_hpbw_deg;
          i = i+1;        id.edit_hpbw_deg = i;                 pointers(i) = obj.handles.edit_hpbw_deg;
          i = i+1;        id.text_deg_hpbw = i;                 pointers(i) = obj.handles.text_deg_hpbw;   
          
          i = i+1;        id.text_SLL_dB = i;                   pointers(i) = obj.handles.text_SLL_dB;
          i = i+1;        id.edit_SLL_dB = i;                   pointers(i) = obj.handles.edit_SLL_dB;
          i = i+1;        id.text_dB_SLL = i;                   pointers(i) = obj.handles.text_dB_SLL; 
          
          i = i+1;        id.text_XPL_dB = i;                   pointers(i) = obj.handles.text_XPL_dB;
          i = i+1;        id.text_dB_XPL = i;                   pointers(i) = obj.handles.text_dB_XPL;
          i = i+1;        id.edit_XPL_dB = i;                   pointers(i) = obj.handles.edit_XPL_dB; 
          
          
          %% PUSH BUTTONS
          i = i+1;        id.pb_ok = i;            pointers(i) = obj.handles.pb_ok;       
          i = i+1;        id.pb_cancel = i;                 pointers(i) = obj.handles.pb_cancel;
          

          groupIDs.all = [id.text_description : id.pb_cancel];
          
          
          %% Initialize object properties
          [groupIDs.gPanels, groupIDs.strEl, groupIDs.valEl] = obj.autoElClassification(pointers);
          
            
          %% SAVE
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
            
            % May not be any dependency
            
        end
        
        % Get enable / disable status of the element of the interface
        % This function should be called only once, later in the code
        % the status is kept updated
        function getAllElStatus(obj)    
                       
            panels = false( length(obj.uiPointers), 1 );
            panels( obj.uiGroups.gPanels) = true;          % logical indexes of the panels
            
            idEl = 1 : length( obj.uiPointers );
            idEl = idEl(~panels);						% id of elements that are not panels
            idEl = idEl(obj.uiIDs.text_description : end);	% The only elements to be considered starts 
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
            
            % May not be needed
            idEl = 1 : length(obj.uiPointers);              % All the elements

            % Sets of panels
            panels = false(length(obj.uiPointers),1);     % init logical panels group
            panels(obj.uiGroups.gPanels) = true;           % set logical indexes of the panels
                
            % Sets of text elements                
            textEl = false(length(obj.uiPointers),1);     % init logical text elements group
            textEl(obj.uiGroups.strEl) = true;              % set logical indexes of the text elements

            % TO-DO: Resolve edit or text field index problem
            mIdEl = setdiff(idEl(~panels & ~textEl & obj.getFlag), []); % id of elements that are not panels nor text elements that have been modified

            % For all the other modified elements
            for i=1:length(mIdEl)
                obj.curVal{mIdEl(i)} = obj.getGuiElVal(obj.uiPointers(mIdEl(i)));
            end
                        
        end
        
        
        function initPopupMenus(obj)
            
        end
        
        
        % Initialize specific GUI elements to this class
        function initSpecificGUI(obj)
            
            % May not be needed
            
        end
        
        
        function setSpecificElContent(obj)
            
            % May not be needed
            
            idEl = 1:length(obj.uiPointers);              % All the elements

            % Sets of panels
            panels = false(length(obj.uiPointers),1);     % init logical panels group
            panels(obj.uiGroups.gPanels) = true;           % set logical indexes of the panels
                
            % Sets of text elements
            textEl = false(length(obj.uiPointers),1);     % init logical text elements group
            textEl(obj.uiGroups.strEl) = true;              % set logical indexes of the text elements

            mIdEl = setdiff(idEl(~panels & ~textEl & obj.setFlag), []); % id of elements that are not panels nor text elements that have been modified

            % For all the other modified elements
            for i=1:length(mIdEl)
                obj.setGuiElVal(obj.uiPointers(mIdEl(i)), obj.newVal{mIdEl(i)});
            end
            
        end
        
        
        % Test if the active file/dir paths
        % contain valid file/dir
        function updateLEDstate(obj)
                        
            % May not be needed
            
        end 
        
    end
    
end

