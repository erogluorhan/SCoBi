classdef RxParams < handle
% class RxParams
%
%   Maintains receiver antenna parameters. It keeps the parameters that are
%   specific to the receiver antenna in the bistatic configuration of any 
%   simulation. It can have only one  instance throughout the whole 
%   simulation thanks to Singleton Pattern. Its properties should be 
%   initialized once in the simulation and then used by other entities by 
%   using the get() functions provided by it. 
%
%   See also initRxParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


    
    properties (SetAccess = private, GetAccess = public)
        % Antenna Height (m)
        hr_m 
        
        % Receive Antenna Gain (dB)
        G0r_dB
        
        % Recevier Antenna polarization
        pol_Rx
        
        % Recevier Antenna orientation index
        orientation_Rx_id
        
        % Receveier antenna observation angle
        th0_Rx_deg
        
        % Receveier antenna azimuth angle
        ph0_Rx_deg
        
        % Recevier Antenna pattern generation method index
        ant_pat_Rx_id
        
        % Antenna pattern struct that contains pattern lookup theta and phi
        % angles, and antenna pattern matrices G and g
        ant_pat_struct_Rx
        
        % Antenna pattern resolution in degrees
        ant_pat_res_deg
    end
    
    
    methods (Access = private)
    
        function obj = RxParams
            % RXPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = RxParams;
             end
             
             singleObj = localObj;
        end
             
    end
    
    methods
        function initialize(obj, hr_m, G0r_dB, pol_Rx, ant_pat_Rx_id, ...
                ant_pat_struct_Rx, ant_pat_res_deg, orientation_Rx_id, ...
                th0_Rx_deg, ph0_Rx_deg)
            % INITIALIZE - Initializes all the properties
            
            obj.hr_m = hr_m;
            obj.G0r_dB = G0r_dB;
            obj.pol_Rx = pol_Rx; 
            obj.orientation_Rx_id = orientation_Rx_id; 
            obj.th0_Rx_deg = th0_Rx_deg;
            obj.ph0_Rx_deg = ph0_Rx_deg;            
            obj.ant_pat_Rx_id = ant_pat_Rx_id;  
            obj.ant_pat_struct_Rx = ant_pat_struct_Rx;
            obj.ant_pat_res_deg = ant_pat_res_deg;
            
        end
        
        function out = get.hr_m(obj)
            out = obj.hr_m;        
        end
        
        function out = get.G0r_dB(obj)
            out = obj.G0r_dB;
        end
        
        function out = get.pol_Rx(obj)
            out = obj.pol_Rx;
        end
        
        function out = get.orientation_Rx_id(obj)
            out = obj.orientation_Rx_id;
        end
        
        function out = get.th0_Rx_deg(obj)
            out = obj.th0_Rx_deg;
        end
        
        function set_th0_Rx_deg(obj, val)
            obj.th0_Rx_deg = val;
        end
        
        function out = get.ph0_Rx_deg(obj)
            out = obj.ph0_Rx_deg;
        end
        
        function set_ph0_Rx_deg(obj, val)
            obj.ph0_Rx_deg = val;
        end
        
        function out = get.ant_pat_Rx_id(obj)
            out = obj.ant_pat_Rx_id;
        end
        
        function out = get.ant_pat_struct_Rx(obj)
            out = obj.ant_pat_struct_Rx;        
        end
        
        function out = get.ant_pat_res_deg(obj)
            out = obj.ant_pat_res_deg;
        end
        
        function set_ant_pat_res_deg(obj, val)
            obj.ant_pat_res_deg = val;
        end
        
    end
    
end

