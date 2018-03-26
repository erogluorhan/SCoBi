classdef RecParams < handle
    %% RECPARAMS CLASS - Maintains receiver antenna parameters
    % It keeps the parameters that are specific to the receiver antenna in 
    % the bistatic configuration of each simulation. It can have only one 
    % instance throughout the whole simulation thanks to Singleton Pattern. 
    % Its properties should be initialized once in the simulation and then 
    % used by other entities by using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        % Antenna Height (m)
        hr_m 
        
        % Receive Antenna Gain (dB)
        G0r_dB
        
        % Beamwidth (degrees)
        hpbw_deg 
        
        % Sidelobe Level (dB)
        SLL_dB 
        
        % X-pol level (dB)
        XPL_dB 
        
        % Antenna polarization
        polR
    end
    
    
    methods (Access = private)
    
        function obj = RecParams
            % RECPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = RecParams;
             end
             
             singleObj = localObj;
        end
             
    end
    
    methods
        
        function initialize(obj, hr_m, G0r_dB, hpbw_deg, SLL_dB, XPL_dB, polR)
            % INITIALIZE - Initializes all the properties
            
            obj.hr_m = hr_m;
            obj.G0r_dB = G0r_dB;
            obj.hpbw_deg = hpbw_deg;
            obj.SLL_dB = SLL_dB;
            obj.XPL_dB = XPL_dB;
            obj.polR = polR;        
        end
        
        function out = get.hr_m(obj)
            out = obj.hr_m;        
        end
        
        function out = get.G0r_dB(obj)
            out = obj.G0r_dB;
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
        
        function out = get.polR(obj)
            out = obj.polR;
        end
    end
    
end

