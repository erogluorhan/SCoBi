classdef GndMLParams < handle
% class GndMLParams
%
%   Maintains multi-layer ground parameters, if any. It keeps the parameters that 
%   are specific to the multi-layer ground, if any, of any simulation. It 
%   can have only one instance throughout the entire simulation thanks to 
%   Singleton Pattern. Its properties should be initialized once in the 
%   simulation and then used by other entities by using the get() functions
%   provided by it. 
%
%   See also initGndMLParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


    
    properties (SetAccess = private, GetAccess = public)
        
        % Used when ground is multi-layered
        layer_depth_m
        
        % Calculated by using layer_depth
        num_layers
        
        % Layer discretization
        delZ_m
        
        % Layer bottom (meters)
        layer_bottom_m
        
        % Layer thickness (meters)
        layer_thickness_m
        
        % Air layer
        zA_m
        
        % The bottom-most layer
        zB_m
        
        % Layers total (meters)
        zS_m
        
        % Layer profile
        z_m
        
        % Flags to calculate several dielectric profile  fitting functions
        calc_diel_profile_fit_functions
        
        % Flag for calculating penetration depth
        % 0: Do not include
        % 1: Include
        calculate_penetration_depth
        
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
        
        function initialize(obj, layer_depth_m, ...
                delZ_m, zA_m, zB_m, calc_diel_profile_fit_functions, ...
                calculate_penetration_depth)
            % INITIALIZE - Initializes all the properties
                
            obj.layer_depth_m = layer_depth_m;
            
            obj.num_layers = length(layer_depth_m);
                        
            obj.delZ_m = delZ_m;
            
            obj.zA_m = zA_m;
            
            obj.zB_m = zB_m;
            
            obj.layer_bottom_m = [0; layer_depth_m(1 : end - 1) + diff(layer_depth_m) / 2 ] ;
            
            obj.layer_thickness_m = diff( obj.layer_bottom_m );
            
            obj.zS_m = obj.zA_m + obj.layer_bottom_m(end) + obj.zB_m ;

            obj.z_m = (0 : obj.delZ_m : obj.zS_m)';
            
            obj.calc_diel_profile_fit_functions = calc_diel_profile_fit_functions;
            
            obj.calculate_penetration_depth = calculate_penetration_depth;
                      
        end 
         
        
        function out = get.calc_diel_profile_fit_functions(obj)
            out = obj.calc_diel_profile_fit_functions;        
        end
         
        
        function out = get.delZ_m(obj)
            out = obj.delZ_m;        
        end
        
        function out = get.layer_depth_m(obj)
            out = obj.layer_depth_m;
        end
        
        function out = get.num_layers(obj)
            out = obj.num_layers;
        end
         
        
        function out = get.layer_bottom_m(obj)
            out = obj.layer_bottom_m;        
        end
         
        
        function out = get.layer_thickness_m(obj)
            out = obj.layer_thickness_m;        
        end
        
        
        function out = get.zA_m(obj)
            out = obj.zA_m;        
        end
        
        
        function out = get.zB_m(obj)
            out = obj.zB_m;        
        end
        
        
        function out = get.zS_m(obj)
            out = obj.zS_m;        
        end
        
        
        function out = get.z_m(obj)
            out = obj.z_m;        
        end
        
        function out = get.calculate_penetration_depth(obj)
            out = obj.calculate_penetration_depth;
        end
        
    end 
    
end

