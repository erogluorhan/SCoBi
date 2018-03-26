classdef SimSettings < handle
    %% SIMSETTINGS CLASS - Maintains simulation settings
    % It keeps the settings that are specific to each 
    % simulation. It can have only one instance throughout the whole
    % simulation thanks to Singleton Pattern. Its properties should be 
    % initialized once in the simulation and then used by other entities by 
    % using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        % Simulation mode for Incidence angle, VSM, and RMSH values
        % 0: Cross simulation
        % 1: Time series
        sim_mode
        
        % Simulator vegetation mode
        % 0: Bare soil only
        % 1: Vegtation + Bare soil
        ground_cover
        
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
        
        function initialize(obj, sim_mode, ground_cover, write_attenuation, ...
                calc_direct_term, calc_specular_term, calc_diffuse_term )
            % INITIALIZE - Initializes all the properties
            
            obj.sim_mode = sim_mode;
            obj.ground_cover = ground_cover;
            obj.write_attenuation = write_attenuation;
            obj.calc_direct_term = calc_direct_term;
            obj.calc_specular_term = calc_specular_term;
            obj.calc_diffuse_term = calc_diffuse_term; 
        end
        
        function out = get.sim_mode(obj)
            out = obj.sim_mode;        
        end
        
        function out = get.ground_cover(obj)
            out = obj.ground_cover;        
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
    end
    
end

