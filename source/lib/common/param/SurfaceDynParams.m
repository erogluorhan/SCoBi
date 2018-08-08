classdef SurfaceDynParams < handle
    %% SURFACEDYNPARAMS CLASS - Maintains ground parameters
    % It keeps the parameters that are specific to the dynamic surface and 
    % each simulation. It can have only one instance throughout the whole
    % simulation thanks to Singleton Pattern. Its properties should be 
    % initialized once in the simulation and then used by other entities by 
    % using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
                
        % Effective surface roughness parameter
        h
        
        % Dielectric permittivity
        eps_g
        
        
    end
    
    
    methods (Access = private)
    
        function obj = SurfaceDynParams
            % GROUNDPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = SurfaceDynParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function initialize(obj, h, eps_g )
            % INITIALIZE - Initializes all the properties
            
            obj.h = h;
                        
            obj.eps_g = eps_g;
            
        end
        
        
        function out = get.eps_g(obj)
            out = obj.eps_g;
        end
        
        
        function out = get.h(obj)
            out = obj.h;
        end
        
    end
    
end

