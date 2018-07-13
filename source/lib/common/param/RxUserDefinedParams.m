classdef RxUserDefinedParams < handle
    %% RXUSERDEFINEDPARAMS CLASS - Maintains the user-defined receiver 
    % antenna parameters
    % It keeps the parameters that are specific to the user-defined 
    % receiver antenna pattern. It can have only one instance throughout 
    % the whole simulation thanks to Singleton Pattern. 
    % Its properties should be initialized once in the simulation and then 
    % used by other entities by using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        
        % Antenna pattern file
        ant_pat_fullfilename 
    end
    
    
    methods (Access = private)
    
        function obj = RxUserDefinedParams
            % RXUSERDEFINEDPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = RxUserDefinedParams;
             end
             
             singleObj = localObj;
        end
             
    end
    
    methods
        
        function initialize(obj, ant_pat_fullfilename)
            % INITIALIZE - Initializes all the properties
            
            obj.ant_pat_fullfilename = ant_pat_fullfilename;
            
        end
        
        function out = get.ant_pat_fullfilename(obj)
            out = obj.ant_pat_fullfilename;
        end
        
    end
    
end

