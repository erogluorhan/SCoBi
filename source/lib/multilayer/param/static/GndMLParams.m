classdef GndMLParams < handle
    % GNDMLPARAMS Maintains specific Ground parameters for MultiLayer soil
    % simulations
    % It keeps the parameters that are specific to multi-layer ground with
    % row structures and each simulation. It can have only one instance 
    % throughout the whole simulation thanks to Singleton Pattern. Its 
    % properties should be initialized once in the simulation and then used
    % by other entities by using the get() functions provided by it.
    
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
                delZ_m, zA_m, zB_m, calc_diel_profile_fit_functions )
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
        
    end 
    
end

