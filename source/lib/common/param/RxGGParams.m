classdef RxGGParams < handle
    %% RXGGPARAMS CLASS - Maintains Generalized-Gaussian receiver antenna 
    % parameters
    % It keeps the parameters that are specific to the receiver antenna in 
    % the bistatic configuration of each simulation. It can have only one 
    % instance throughout the whole simulation thanks to Singleton Pattern. 
    % Its properties should be initialized once in the simulation and then 
    % used by other entities by using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        
        % Beamwidth (degrees)
        hpbw_deg 
        
        % Sidelobe Level (dB)
        SLL_dB 
        
        % X-pol level (dB)
        XPL_dB 
    end
    
    
    methods (Access = private)
    
        function obj = RxGGParams
            % RXGGPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = RxGGParams;
             end
             
             singleObj = localObj;
        end
             
    end
    
    methods
        
        function initialize(obj, hpbw_deg, SLL_dB, XPL_dB)
            % INITIALIZE - Initializes all the properties
            
            obj.hpbw_deg = hpbw_deg;
            obj.SLL_dB = SLL_dB;
            obj.XPL_dB = XPL_dB;
            
        end
        
        function out = get.hpbw_deg(obj)
            out = obj.hpbw_deg;
        end
        
        function out = get.SLL_dB(obj)
            out = obj.SLL_dB;
        end
        
        function out = get.XPL_dB(obj)
            out = obj.XPL_dB;
        end
        
    end
    
end

