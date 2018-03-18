classdef VegVirParams < handle
    %VEGVIRPARAMS Maintains virtual vegetation parameters
    % It inherits VegParams.
    % It keeps the parameters that are specific to virtual veg. and each 
    % simulation. It can have only one instance throughout the whole
    % simulation thanks to Singleton Pattern. Its properties should be 
    % initialized once in the simulation and then used by other entities by 
    % using the get() functions provided by it.
    
    
    properties (SetAccess = private, GetAccess = public)        

    end
    
    methods
        
        function initialize(obj)
                      
        end  
                
    end
    
    methods (Access = protected)
    
        function obj = VegVirParams
            % VEGVIRPARAMS - Protected constructor            
        end
    
    end
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = VegVirParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    
end

