
classdef gui_SCoBiMain_Manager < SCoBiGUIManagers
% class gui_SCoBiMain_Manager 
%
%   Implements handles for gui_SCoBiMain.fig and gui_SCoBiMain.m.
%   
%   simulator_id = gui_SCoBiMain;
%
%   See also SCoBiGUIManagers, gui_SCoBi_Manager.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd
%   Adapted from gui_goGPS.m class of goGPS v0.4.3 software

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0

    
    
    %% METHODS
    
    %% COMMON METHODS
    methods
        
        % Constructor
        function obj = gui_SCoBiMain_Manager( handles )            
            
            % Call superclass constructor
            obj = obj@SCoBiGUIManagers( handles, -1 );
            
        end
        
        function funout = outputFun(obj)
            
            funout{1} = obj.simulator_id;
            
        end
        
        
        % EVENT MANAGER
        % When an element is modified in the GUI, this must be called!
        % Sync the internal copy of the interface with the GUI
        function syncFromGUI(obj, idEl)
            
            if (nargin == 1)
                idEl = obj.uiIDs.pb_SCoBi_main;
            end
            
            obj.getFlag(idEl) = true;
            
            % Read all the values of the elements
            obj.getAllElContent();
            obj.setAllElContent();
            
          %% PUSH BUTTONS                        
            % on SCoBi main illustration 
            if sum(intersect(idEl, obj.uiIDs.pb_SCoBi_main)) > 0
                % Do not require any operation
            end
            
            % on the sim: Snow
            if sum(intersect(idEl, obj.uiIDs.pb_SCoBi_Snow)) > 0
                
                % Display a warning that the method is not implemented yet
                waitfor(msgbox('WARNING: This is not yet implemented! Please choose another one!'));
                  
            end
            
            % on the sim: Soils
            if sum(intersect(idEl, obj.uiIDs.pb_SCoBi_Soil)) > 0
                
                obj.simulator_id = Constants.ID_SIM_SOIL;
                
                uiresume(obj.handles.panel_main);
                  
            end
            
            % on the sim: Root-zone 
            if sum(intersect(idEl, obj.uiIDs.pb_SCoBi_Root_Zone)) > 0
                
                obj.simulator_id = Constants.ID_SIM_ROOT_ZONE;
                
                uiresume(obj.handles.panel_main);
                
            end
            
            % on the sim: Topography
            if sum(intersect(idEl, obj.uiIDs.pb_SCoBi_Topography)) > 0
                
                % Display a warning that the method is not implemented yet
                waitfor(msgbox('WARNING: This is not yet implemented! Please choose another one!'));
                
            end
            
            % on the sim: Agriculture 
            if sum(intersect(idEl, obj.uiIDs.pb_SCoBi_Agriculture)) > 0
                
                obj.simulator_id = Constants.ID_SIM_AGRICULTURE;
                
                uiresume(obj.handles.panel_main);
                
            end
            
            % on the sim: Forest 
            if sum(intersect(idEl, obj.uiIDs.pb_SCoBi_Forest)) > 0
                
                obj.simulator_id = Constants.ID_SIM_FOREST;
                
                uiresume(obj.handles.panel_main);
                
            end
            
            % on the sim: Permafrost
            if sum(intersect(idEl, obj.uiIDs.pb_SCoBi_Permafrost)) > 0
                
                % Display a warning that the method is not implemented yet
                waitfor(msgbox('WARNING: This is not yet implemented! Please choose another one!'));
                
            end
            
            % on the sim: About
            if sum(intersect(idEl, obj.uiIDs.pb_about)) > 0
                
                obj.openAboutDialog();
                  
            end
            
            % on the sim: Documents
            if sum(intersect(idEl, obj.uiIDs.pb_documents)) > 0                
                
                obj.openDocumentsDialog();
                  
            end
            
            obj.onoffUIEl();
            obj.checkUIdependencies();
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
            
            % No pop-up menus in main SCOBi window
            
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
          
          i = i+1;        id.pb_SCoBi_main = i;             pointers(i) = obj.handles.pb_SCoBi_main;          
          i = i+1;        id.pb_SCoBi_Snow = i;             pointers(i) = obj.handles.pb_SCoBi_Snow;         
          i = i+1;        id.pb_SCoBi_Soil = i;             pointers(i) = obj.handles.pb_SCoBi_Soil;          
          i = i+1;        id.pb_SCoBi_Root_Zone = i;        pointers(i) = obj.handles.pb_SCoBi_Root_Zone;    
          i = i+1;        id.pb_SCoBi_Topography = i;       pointers(i) = obj.handles.pb_SCoBi_Topography;       
          i = i+1;        id.pb_SCoBi_Agriculture = i;      pointers(i) = obj.handles.pb_SCoBi_Agriculture;       
          i = i+1;        id.pb_SCoBi_Forest = i;           pointers(i) = obj.handles.pb_SCoBi_Forest;        
          i = i+1;        id.pb_SCoBi_Permafrost = i;       pointers(i) = obj.handles.pb_SCoBi_Permafrost;    
          i = i+1;        id.pb_about = i;                  pointers(i) = obj.handles.pb_about;    
          i = i+1;        id.pb_documents = i;              pointers(i) = obj.handles.pb_documents; 
          
          groupIDs.all = id.text_title : id.pb_documents;

            
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
            idEl = idEl(obj.uiIDs.pb_SCoBi_main : end);	% The only elements to be considered starts 
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

            mIdEl = setdiff(idEl(~panels & ~textEl & obj.getFlag), []); % id of elements that are not panels nor text elements that have been modified

            % For all the other modified elements
            for i=1:length(mIdEl)
                obj.curVal{mIdEl(i)} = obj.getGuiElVal(obj.uiPointers(mIdEl(i)));
            end
            
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
        
    end
    
end

