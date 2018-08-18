classdef SimulationFolders < handle
    % SIMULATIONFOLDERS Class to keep track of all output directories
    %   This class has one attribute for each source code or input
    %   directory. Every attribute can be reached by a static getter method.
    
    
properties (SetAccess = private, GetAccess = public)


    sims_main_dir    

    
    %% ANALYSIS
    analysis


    % OUTPUT
    output
    sim

    % Input
    sim_input

    % Metadata
    metadata
    afsa
    geo
    position
    fzones
    config
    rot_lookup

    % Products
    products
    products_direct
    products_direct_field
    products_direct_power
    products_specular
    products_specular_field
    products_specular_power

    % Figures
    fig
    fig_direct
    fig_specular
    fig_specular_P
    fig_specular_P_vsTH

end
    
    
methods (Access = private)

    function obj = SimulationFolders
    end

    end


    methods (Static)

    function singleObj = getInstance
    persistent localObj

    if isempty(localObj) || ~isvalid(localObj)
       localObj = SimulationFolders;
    end

     singleObj = localObj;
    end
    end
        
    
methods
        
    function initializeStaticDirs(obj)


    %% GET GLOBAL DIRECTORIES
    main_dir = Directories.getInstance.main_dir;


    %% GET GLOBAL PARAMETERS
    % Simulation Settings
    include_in_master_sim_file = SimSettings.getInstance.include_in_master_sim_file;
    % Simulation Parameters
    sim_name = SimParams.getInstance.sim_name;


    %% INITIALIZE DIRECTORIES
    % Initialize static simulation directories
    if include_in_master_sim_file

        obj.sims_main_dir = strcat( main_dir, '\sims\master');

    else

        obj.sims_main_dir = strcat( main_dir, '\sims\temp');
    end


    %% ANALYSIS
    obj.analysis = strcat( obj.sims_main_dir, '\', 'analysis') ; 


    %% OUTPUT
    obj.output = strcat( obj.sims_main_dir , '\output');

    % TO-DO: create a unique name with timestamp
    obj.sim = strcat( obj.output, '\', sim_name);

    % Meta-Data
    obj.sim_input = strcat( obj.sim, '\input');

    % Meta-Data
    obj.metadata = strcat( obj.sim, '\metadata');

    % Average Forward Scattering Amplitude
    obj.afsa = strcat(obj.metadata, '\', 'afsa');  

    % Products
    obj.products = strcat(obj.sim, '\', 'products') ;
    obj.products_direct = strcat(obj.products, '\', 'direct') ;
    obj.products_direct_field = strcat(obj.products_direct, '\', 'field') ;
    obj.products_direct_power = strcat(obj.products_direct, '\', 'power') ;
    obj.products_specular = strcat(obj.products, '\', 'specular') ;
    obj.products_specular_field = strcat(obj.products_specular, '\', 'field') ;
    obj.products_specular_power = strcat(obj.products_specular, '\', 'power') ;

    % Figure
    obj.fig = strcat(obj.sim, '\', 'figure' ) ;
    obj.fig_direct = strcat(obj.fig, '\', 'direct') ;
    obj.fig_specular = strcat(obj.fig, '\', 'specular') ;
    obj.fig_specular_P = strcat(obj.fig_specular, '\', 'P') ;
    obj.fig_specular_P_vsTH = strcat(obj.fig_specular_P, '\', 'vs_TH') ;


    end



    function initializeDynamicDirs(obj)
    % Initialize dynamically changing simulation directories


    %% GET GLOBAL PARAMETERS
    sim_counter = ParamsManager.sim_counter;
    % Configuration Parameters
    th0_Tx_list_deg = ConfigParams.getInstance.th0_Tx_list_deg;
    th0_Tx_deg = th0_Tx_list_deg( sim_counter );
    ph0_Tx_list_deg = ConfigParams.getInstance.ph0_Tx_list_deg;
    ph0_Tx_deg = ph0_Tx_list_deg( sim_counter );


    %% Angle of incidence
    th0_Tx_deg_folder = strcat('th0_Tx_', num2str( th0_Tx_deg ), ...
                         '-ph0_Tx_', num2str( ph0_Tx_deg ) ) ;
    dir_th0_ph0_Tx_deg = strcat(obj.metadata, '\', th0_Tx_deg_folder) ;


    %% Geometry
    obj.geo = strcat(dir_th0_ph0_Tx_deg, '\', 'geo') ;
    obj.position = strcat(obj.geo, '\', 'position') ;                        
    obj.fzones = strcat(obj.geo, '\', 'fZones') ; 

    %Configuration
    obj.config = strcat(dir_th0_ph0_Tx_deg, '\', 'configuration') ;    


    %% Rotation
    obj.rot_lookup = strcat(dir_th0_ph0_Tx_deg, '\', 'rotation');


    obj.makeDynamicDirs();

    end


    function makeStaticDirs(obj)


        %% Simulations main directory
        if ~exist(obj.sims_main_dir, 'dir')
            mkdir(obj.sims_main_dir);
        end


        %% Analysis
        if ~exist(obj.analysis, 'dir')
            mkdir(obj.analysis);
        end

        %% Input folder that includes sim. report and input file used
        if ~exist(obj.sim_input, 'dir')
            mkdir(obj.sim_input)
        end

        %% Average Forward Scattering Amplitude
        if ~exist(obj.afsa, 'dir')
            mkdir(obj.afsa)
        end

        %% Products
        if ~exist(obj.products, 'dir')
            mkdir(obj.products)
        end

        if ~exist(obj.products_direct, 'dir')
            mkdir(obj.products_direct)
        end

        if ~exist(obj.products_direct_field, 'dir')
            mkdir(obj.products_direct_field)
        end

        if ~exist(obj.products_direct_power, 'dir')
            mkdir(obj.products_direct_power)
        end

        if ~exist(obj.products_specular, 'dir')
            mkdir(obj.products_specular)
        end

        if ~exist(obj.products_specular_field, 'dir')
            mkdir(obj.products_specular_field)
        end

        if ~exist(obj.products_specular_power, 'dir')
            mkdir(obj.products_specular_power)
        end

        % Figure
        if ~exist(obj.fig, 'dir')
            mkdir(obj.fig)
        end

        if ~exist(obj.fig_direct, 'dir')
            mkdir(obj.fig_direct)
        end

        if ~exist(obj.fig_specular, 'dir')
            mkdir(obj.fig_specular)
        end

        if ~exist(obj.fig_specular_P, 'dir')
            mkdir(obj.fig_specular_P)
        end

        if ~exist(obj.fig_specular_P_vsTH, 'dir')
            mkdir(obj.fig_specular_P_vsTH)
        end

    end


    function makeDynamicDirs(obj)

        %% Geometry
        if ~exist(obj.geo, 'dir')
            mkdir(obj.geo)
        end

        if ~exist(obj.position, 'dir')
            mkdir(obj.position)
        end

        if ~exist(obj.fzones, 'dir')
            mkdir(obj.fzones)
        end

        %% Configuration
        if ~exist(obj.config, 'dir')
            mkdir(obj.config)
        end

        %% Rotation
        if ~exist(obj.rot_lookup, 'dir')
            mkdir(obj.rot_lookup)
        end

    end

    function out = get.sims_main_dir(obj)
        out = obj.sims_main_dir;
    end

    function out = get.analysis(obj)
        out = obj.analysis;
    end

    function out = get.sim(obj)
        out = obj.sim;
    end

    function out = get.sim_input(obj)
        out = obj.sim_input;
    end 

    function out = get.afsa(obj)
        out = obj.afsa;
    end 

    function out = get.position(obj)
        out = obj.position;
    end 

    function out = get.fzones(obj)
        out = obj.fzones;
    end 

    function out = get.config(obj)
        out = obj.config;
    end

    function out = get.rot_lookup(obj)
        out = obj.rot_lookup;
    end

    function out = get.products_direct(obj)
        out = obj.products_direct;
    end

    function out = get.products(obj)
        out = obj.products;
    end

    function out = get.products_direct_field(obj)
        out = obj.products_direct_field;
    end

    function out = get.products_direct_power(obj)
        out = obj.products_direct_power;
    end

    function out = get.products_specular(obj)
        out = obj.products_specular;
    end

    function out = get.products_specular_field(obj)
        out = obj.products_specular_field;
    end

    function out = get.products_specular_power(obj)
        out = obj.products_specular_power;
    end

    function out = get.fig(obj)
        out = obj.fig;
    end

    function out = get.fig_direct(obj)
        out = obj.fig_direct;
    end

    function out = get.fig_specular(obj)
        out = obj.fig_specular;
    end

    function out = get.fig_specular_P(obj)
        out = obj.fig_specular_P;
    end

    function out = get.fig_specular_P_vsTH(obj)
        out = obj.fig_specular_P_vsTH;
    end
        
end
    
end

