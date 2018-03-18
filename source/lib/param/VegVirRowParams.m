classdef VegVirRowParams < handle
    %VEGVIRROWPARAMS Maintains virtual vegetation with row structures ...
    % parameters
    % It inherits VegVirParams.
    % It keeps the parameters that are specific to virtual vegetations with
    % row structures and each simulation. It can have only one instance 
    % throughout the whole simulation thanks to Singleton Pattern. Its 
    % properties should be initialized once in the simulation and then used
    % by other entities by using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        
        % The distance without vegetation between two rows (m)
        row_space
        
        % The distance between two adjacent plants within a row (m)
        col_space
        
        % Azimuth angle of field rows from local North (degrees)
        phi_row
        
        % Max scattering dist. of a plant pos between rows (m)
        plant_row_spread
        
        % Max scattering dist. of a plant pos within a row (m)
        plant_col_spread
        
        % Plugin class instance that will be handling the virtual 
        % vegetation plant generation
        plugin
        
    end
    
    methods (Access = private)
    
        function obj = VegVirRowParams
            % VEGVIRROWPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = VegVirRowParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods        
        
        function initialize(obj, row_space, col_space, phi_row, ...
                plant_row_spread, plant_col_spread, plugin )
            % INITIALIZE - Initializes all the properties
            
            obj.row_space = row_space;
            obj.col_space = col_space;
            obj.phi_row = phi_row;
            obj.plant_row_spread = plant_row_spread;
            obj.plant_col_spread = plant_col_spread;
            
            % Call the virtual vegetation plugin's initialize function
            obj.plugin = feval(plugin);
            obj.plugin.initialize;
            
            VegParams.getInstance.initializeStage( obj.plugin.stage );
                      
        end 
         
        
        function out = get.row_space(obj)
            out = obj.row_space;        
        end
        
        function out = get.col_space(obj)
            out = obj.col_space;
        end
        
        function out = get.phi_row(obj)
            out = obj.phi_row;
        end
        
        function out = get.plant_row_spread(obj)
            out = obj.plant_row_spread;
        end
        
        function out = get.plant_col_spread(obj)
            out = obj.plant_col_spread;
        end
        
        function out = get.plugin(obj)
            out = obj.plugin;        
        end
        
    end 
    
end

