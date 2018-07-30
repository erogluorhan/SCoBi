classdef SimSettings < handle
    %% SIMSETTINGS CLASS - Maintains simulation settings
    % It keeps the settings that are specific to each 
    % simulation. It can have only one instance throughout the whole
    % simulation thanks to Singleton Pattern. Its properties should be 
    % initialized once in the simulation and then used by other entities by 
    % using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        % Simulator to be run
        % 1: SCoBi-Veg
        % 2: SCoBi-ML
        simulator_id
        
        % Simulation mode
        % 1: Snapshot
        % 2: Time-series
        sim_mode_id
        
        % Ground cover
        % 1: Bare-soil
        % 2: Vegetation
        gnd_cover_id
        
        % Flag for whether to write attenuation to Excel file or not
        % 0: Do not write 
        % 1: Write
        write_attenuation
        
        % Flag for calculation of Direct Term
        % 0: Do not calculate
        % 1: Calculate
        calc_direct_term
        
        % Flag for calculation of Specular Term
        % 0: Do not calculate
        % 1: Calculate
        calc_specular_term
        
        % Flag for calculation of Diffuse Term
        % 0: Do not calculate
        % 1: Calculate
        calc_diffuse_term
        
        % Flag for drawing live plots during simulations
        draw_live_plots
    end
    
    
    methods (Access = private)
    
        function obj = SimSettings
            % SIMSETTINGS - Private constructor
        end
    
    end
    
    methods (Static)
        
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = SimSettings;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function initialize(obj, simulator_id, sim_mode_id, gnd_cover_id, write_attenuation, ...
                calc_direct_term, calc_specular_term, calc_diffuse_term, draw_live_plots )
            % INITIALIZE - Initializes all the properties
            
            obj.simulator_id = simulator_id;
            obj.sim_mode_id = sim_mode_id;
            obj.gnd_cover_id = gnd_cover_id;
            obj.write_attenuation = write_attenuation;
            obj.calc_direct_term = calc_direct_term;
            obj.calc_specular_term = calc_specular_term;
            obj.calc_diffuse_term = calc_diffuse_term; 
            obj.draw_live_plots = draw_live_plots; 
        end
        
        function out = get.simulator_id(obj)
            out = obj.simulator_id;        
        end
        
        function out = get.sim_mode_id(obj)
            out = obj.sim_mode_id;        
        end
        
        function out = get.gnd_cover_id(obj)
            out = obj.gnd_cover_id;        
        end
        
        function out = get.write_attenuation(obj)
            out = obj.write_attenuation;
        end        
        
        function out = get.calc_direct_term(obj)
            out = obj.calc_direct_term;
        end
        
        function out = get.calc_specular_term(obj)
            out = obj.calc_specular_term;
        end
        
        function out = get.calc_diffuse_term(obj)
            out = obj.calc_diffuse_term;
        end
        
        function out = get.draw_live_plots(obj)
            out = obj.draw_live_plots;
        end
    end
    
end

