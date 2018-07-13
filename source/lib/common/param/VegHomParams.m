classdef VegHomParams < VegParams
    %VEGHOMPARAMS Maintains homogenous vegetation parameters
    % It inherits VegParams.
    % It keeps the parameters that are specific to homogenous veg. and each 
    % simulation. It can have only one instance throughout the whole
    % simulation thanks to Singleton Pattern. Its properties should be 
    % initialized once in the simulation and then used by other entities by 
    % using the get() functions provided by it.
    
    methods (Access = private)
    
        function obj = VegHomParams
            % VEGHOMPARAMS - Private constructor            
        end
    
    end
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = VegHomParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    
end

