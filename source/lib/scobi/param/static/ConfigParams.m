classdef ConfigParams < handle
% class ConfigParams
%
%   Maintains configuration parameters. It keeps the configuration 
%   parameters that may change for every simulation. It can have only one 
%   instance throughout the whole simulation thanks to Singleton Pattern. 
%   Its properties should be initialized once in the simulation and then 
%   used by other entities by using the get() functions provided by it. 
%
%   See also initConfigParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


    
    properties (SetAccess = private, GetAccess = public) 
        
        % Day-of-year for timestamping purposes
        DoYs
        
        % Transmitter Elevation Angle list of incoming signal  - measured between ground plane and the 
        % ground-Transmitter direction (degrees) 
        el0_Tx_list_deg
        
        % Transmitter Incidence Angle list of incoming signal - measured between ground zenith and the 
        % ground-Transmitter direction (degrees) 
        th0_Tx_list_deg     % 90 - el0_Tx_list_deg
        
        % Transmitter Azimuth angle list of incoming signal (degrees)- the horizontal angle measured at the 
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
    
        function obj = ConfigParams
            % CONFIGPARAMS - Private constructor
        end
    
    end
    
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = ConfigParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function initialize(obj, DoYs, el0_Tx_list_deg, ph0_Tx_list_deg, VSM_list_cm3cm3, RMSH_list_cm )
            % INITIALIZE - Initializes all the properties
            
            obj.DoYs = DoYs;
            obj.el0_Tx_list_deg = el0_Tx_list_deg;
            % Compute incidence angles from the elevation angles
            obj.th0_Tx_list_deg = 90 - obj.el0_Tx_list_deg;
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

