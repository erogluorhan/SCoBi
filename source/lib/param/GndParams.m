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
        VSM
        
        % Sand content ratio of the soil texture [0,1]
        sand
        
        % Clay content ratio of the soil texture [0,1]
        clay
        
        % Soil bulk density
        rho_b
        
        % Surface rms height (cm)
        RMSH
        
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
        
        function initialize(obj, VSM, sand, clay, rho_b, RMSH )
            % INITIALIZE - Initializes all the properties
            
            obj.VSM = VSM;
            obj.sand = sand;
            obj.clay = clay;
            obj.rho_b = rho_b;
            obj.RMSH = RMSH;    
        end
        
        function out = get.VSM(obj)
            out = obj.VSM;        
        end
        
        function out = get.sand(obj)
            out = obj.sand;
        end
        
        function out = get.clay(obj)
            out = obj.clay;
        end
        
        function out = get.rho_b(obj)
            out = obj.rho_b;
        end
        
        function out = get.RMSH(obj)
            out = obj.RMSH;
        end
        
        function out = get.polG(obj)
            out = obj.polG;
        end
        
        function rep_VSM( obj, factor )
            obj.VSM = repmat(obj.VSM, 1, factor);
        end
        
        function rep_RMSH( obj, factor )
            obj.RMSH = repmat(obj.RMSH, 1, factor);
        end
    end
    
end

