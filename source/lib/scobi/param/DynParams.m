classdef DynParams < handle
    %% DYNPARAMS CLASS - Maintains dynamic parameters
    % It keeps the dynamic parameters that changes for every simulation. 
    % It can have only one instance throughout the whole simulation thanks 
    % to Singleton Pattern. Its properties should be initialized once in 
    % the simulation and then used by other entities by using the get() 
    % functions provided by it.
    
    properties (SetAccess = private, GetAccess = public) 
        
        % Day-of-year for timestamping purposes
        DoYs
        
        % Transmitter Incidence Angle of incoming signal list - measured between ground zenith and the 
        % ground-Transmitter direction (degrees) 
        th0_Tx_list_deg
        
        % Transmitter Azimuth angle of incoming signal list (degrees)- the horizontal angle measured at the 
        % transmitter's position on the earth to Northpole, i.e. measured 
        % clockwise from the North pole IMPORTANT: It is used in Geo 
        % Staellite comm. and different from standard spherical azimuth 
        % direction that is measured from local East in the 
        % counter-clockwise direction.
        ph0_Tx_list_deg
        
        % Volumetric soil moisture (cm3/cm3) [0,1] 
        % The amount of moisture content (cm3) in 1 cm3 of soil
        VSM_list_cm3cm3
        
        % Surface rms height (cm)
        RMSH_list_cm
        
    end
    
    
    methods (Access = private)
    
        function obj = DynParams
            % DYNPARAMS - Private constructor
        end
    
    end
    
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = DynParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function initialize(obj, DoYs, th0_Tx_list_deg, ph0_Tx_list_deg, VSM_list_cm3cm3, RMSH_list_cm )
            % INITIALIZE - Initializes all the properties
            
            obj.DoYs = DoYs;
            obj.th0_Tx_list_deg = th0_Tx_list_deg;
            obj.ph0_Tx_list_deg = ph0_Tx_list_deg;
            obj.VSM_list_cm3cm3 = VSM_list_cm3cm3;
            obj.RMSH_list_cm = RMSH_list_cm;             
        end
        
        function out = get.DoYs(obj)
            out = obj.DoYs;
        end
        
        function out = get.th0_Tx_list_deg(obj)
            out = obj.th0_Tx_list_deg;
        end
        
        function out = get.ph0_Tx_list_deg(obj)
            out = obj.ph0_Tx_list_deg;
        end
        
        function out = get.RMSH_list_cm(obj)
            out = obj.RMSH_list_cm;
        end
        
        function out = get.VSM_list_cm3cm3(obj)
            out = obj.VSM_list_cm3cm3;        
        end
    end
    
end

