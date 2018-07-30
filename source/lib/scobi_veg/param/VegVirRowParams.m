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
        row_space_m
        
        % The distance between two adjacent plants within a row (m)
        col_space_m
        
        % Azimuth angle of field rows from local North (degrees)
        phi_row_deg
        
        % Max scattering dist. of a plant pos from periodic row structure (m)
        seed_fluctuation_m
        
        % Plugin class instance that will be handling the virtual 
        % vegetation plant generation
        plugin          
        
        % Vegetation stage to separate simulations for different stages of 
        % the same vegetation 
        vegetation_stage;  
        
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
        
        function initialize(obj, row_space_m, col_space_m, phi_row_deg, ...
                seed_fluctuation_m, plugin, vegetation_stage )
            % INITIALIZE - Initializes all the properties         
        
                        
            obj.vegetation_stage = vegetation_stage;
            
            obj.row_space_m = row_space_m;
            obj.col_space_m = col_space_m;
            obj.phi_row_deg = phi_row_deg;
            obj.seed_fluctuation_m = seed_fluctuation_m;
            
            % Call the virtual vegetation plugin's initialize function
            obj.plugin = feval(plugin);
            obj.plugin.initialize( vegetation_stage );
                      
        end 
         
        
        function out = get.row_space_m(obj)
            out = obj.row_space_m;        
        end
        
        function out = get.col_space_m(obj)
            out = obj.col_space_m;
        end
        
        function out = get.phi_row_deg(obj)
            out = obj.phi_row_deg;
        end
        
        function out = get.seed_fluctuation_m(obj)
            out = obj.seed_fluctuation_m;
        end
        
        function out = get.plugin(obj)
            out = obj.plugin;        
        end
        
    end 
    
end

