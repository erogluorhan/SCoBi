classdef GndMLParams < handle
    % GNDMLPARAMS Maintains specific Ground parameters for MultiLayer soil
    % simulations
    % It inherits GndParams.
    % It keeps the parameters that are specific to multi-layer ground with
    % row structures and each simulation. It can have only one instance 
    % throughout the whole simulation thanks to Singleton Pattern. Its 
    % properties should be initialized once in the simulation and then used
    % by other entities by using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        
        % Layer discretization
        delZ_m
        
        % Air layer
        zA_m
        
        % The bottom-most layer
        zB_m
        
    end
    
    methods (Access = private)
    
        function obj = GndMLParams
            % GNDMLPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = GndMLParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods        
        
        function initialize(obj, delZ_m, zA_m, zB_m )
            % INITIALIZE - Initializes all the properties
                        
            obj.delZ_m = delZ_m;
            
            obj.zA_m = zA_m;
            
            obj.zB_m = zB_m;
                      
        end 
         
        
        function out = get.delZ_m(obj)
            out = obj.delZ_m;        
        end
        
        
        function out = get.zA_m(obj)
            out = obj.zA_m;        
        end
        
        
        function out = get.zB_m(obj)
            out = obj.zB_m;        
        end
        
    end 
    
end

