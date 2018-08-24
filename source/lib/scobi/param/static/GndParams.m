classdef GndParams < handle
    %% GNDPARAMS CLASS - Maintains ground parameters
    % It keeps the parameters that are specific to the ground and each 
    % simulation. It can have only one instance throughout the whole
    % simulation thanks to Singleton Pattern. Its properties should be 
    % initialized once in the simulation and then used by other entities by 
    % using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
                
        % Gorund layer structure: Single-layered (1), or Multi-layered (2)
        gnd_structure_id
        
        % Sand content ratio of the soil texture [0,1]
        sand_ratio
        
        % Clay content ratio of the soil texture [0,1]
        clay_ratio
        
        % Soil bulk density (g/cm3)
        rhob_gcm3
        
        % Ground polarization
        polG = 'V';
        
        % Dielectric permittivity model index
        diel_model_id
        
    end
    
    
    methods (Access = private)
    
        function obj = GndParams
            % GROUNDPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = GndParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function initialize(obj, gnd_structure_id, sand_ratio, clay_ratio, rhob_gcm3, diel_model_id )
            % INITIALIZE - Initializes all the properties
                
            obj.gnd_structure_id = gnd_structure_id;
            
            obj.sand_ratio = sand_ratio;
            
            obj.clay_ratio = clay_ratio;
            
            obj.rhob_gcm3 = rhob_gcm3;  
            
            obj.diel_model_id = diel_model_id;
            
        end
        
        function out = get.gnd_structure_id(obj)
            out = obj.gnd_structure_id;
        end
        
        function out = get.diel_model_id(obj)
            out = obj.diel_model_id;
        end
        
        function out = get.sand_ratio(obj)
            out = obj.sand_ratio;
        end
        
        function out = get.clay_ratio(obj)
            out = obj.clay_ratio;
        end
        
        function out = get.rhob_gcm3(obj)
            out = obj.rhob_gcm3;
        end
        
        function out = get.polG(obj)
            out = obj.polG;
        end
    end
    
end

