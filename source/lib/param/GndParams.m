classdef GndParams < handle
    %% GNDPARAMS CLASS - Maintains ground parameters
    % It keeps the parameters that are specific to the ground and each 
    % simulation. It can have only one instance throughout the whole
    % simulation thanks to Singleton Pattern. Its properties should be 
    % initialized once in the simulation and then used by other entities by 
    % using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        % Volumetric soil moisture (cm3/cm3) [0,1] 
        % The amount of moisture content (cm3) in 1 cm3 of soil
        VSM_cm3cm3
        
        % Sand content ratio of the soil texture [0,1]
        sand_ratio
        
        % Clay content ratio of the soil texture [0,1]
        clay_ratio
        
        % Soil bulk density (g/cm3)
        rhob_gcm3
        
        % Surface rms height (cm)
        RMSH_cm
        
        % Ground polarization
        polG = 'V';
        
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
        
        function initialize(obj, VSM_cm3cm3, sand_ratio, clay_ratio, rhob_gcm3, RMSH_cm )
            % INITIALIZE - Initializes all the properties
            
            obj.VSM_cm3cm3 = VSM_cm3cm3;
            obj.sand_ratio = sand_ratio;
            obj.clay_ratio = clay_ratio;
            obj.rhob_gcm3 = rhob_gcm3;
            obj.RMSH_cm = RMSH_cm;    
        end
        
        function out = get.VSM_cm3cm3(obj)
            out = obj.VSM_cm3cm3;        
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
        
        function out = get.RMSH_cm(obj)
            out = obj.RMSH_cm;
        end
        
        function out = get.polG(obj)
            out = obj.polG;
        end
        
        function rep_VSM( obj, factor )
            obj.VSM_cm3cm3 = repmat(obj.VSM_cm3cm3, 1, factor);
        end
        
        function rep_RMSH( obj, factor )
            obj.RMSH_cm = repmat(obj.RMSH_cm, 1, factor);
        end
    end
    
end

