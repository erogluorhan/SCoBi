classdef SimParams < handle
    %% SIMPARAMS CLASS - Maintains simulation parameters
    % It keeps the parameters that are specific to the simulator and each 
    % simulation. It can have only one instance throughout the whole
    % simulation thanks to Singleton Pattern. Its properties should be 
    % initialized once in the simulation and then used by other entities by 
    % using the get() functions provided by it.
    
    properties (SetAccess = private, GetAccess = public)
        % Simulator version, e.g. "1.0"
        version = Constants.version;
        
        % Specific name of the simulation, any string
        sim_name
        
        % Campaign name, e.g. "MSU"
        campaign
        
        % Campaign date, e.g. "20171103"
        campaign_date
        
        % Plot of the campaign, e.g. "01"
        plot
        
        % Vegetation may be any name, e.g. "Corn", "Paulownia" etc.
        vegetation_plant
        
        % Number of realizations must be an integer.
        % It is needed for Monte Carlo simulations to generate the Diffuse
        % Term.
        Nr
        
        % Number of Fresnel zones
        Nfz
    end
    
    
    methods (Access = private)
    
        function obj = SimParams
            % SIMPARAMS - Private constructor
        end
    
    end
    
    methods (Static)
        
        
        function singleObj = getInstance
            % GETINSTANCE - One instance for Singleton Pattern
             persistent localObj
             
             if isempty(localObj) || ~isvalid(localObj)
                localObj = SimParams;
             end
             
             singleObj = localObj;
        end
        
    end
    
    
    methods
        
        function initialize(obj, sim_name, campaign, campaign_date, ...
                plot, vegetation_plant, Nr, Nfz )
            % INITIALIZE - Initializes all the properties
            
            obj.sim_name = sim_name;
            obj.campaign = campaign;
            obj.campaign_date = campaign_date;
            obj.plot = plot;
            obj.vegetation_plant = vegetation_plant;
            obj.Nr = Nr;  
            obj.Nfz = Nfz;
        end
        
        function updateSimName(obj, sim_name )
            obj.sim_name = sim_name;     
        end
        
        function out = get.version(obj)
            out = obj.version;        
        end
        
        function out = get.sim_name(obj)
            out = obj.sim_name;
        end
        
        function out = get.campaign(obj)
            out = obj.campaign;
        end
        
        function out = get.campaign_date(obj)
            out = obj.campaign_date;
        end
        
        function out = get.plot(obj)
            out = obj.plot;
        end
        
        function out = get.vegetation_plant(obj)
            out = obj.vegetation_plant;
        end
        
        function out = get.Nr(obj)
            out = obj.Nr;
        end
        
        function out = get.Nfz(obj)
            out = obj.Nfz;        
        end
    end
    
end

