
classdef SimulationFolders < handle
% class specularTerm 
%
%   Handles the output directories. This class has a property for each 
%   simulation output directory that is to be managed, accessed, etc.
%   . Its properties can be reached by static getter methods.
%
%   See also Directories.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0

    
    
properties (SetAccess = private, GetAccess = public)


    sims_main_dir  


    % OUTPUT
    sim

    % Input
    sim_input
    sim_input_used_files

    % Metadata
    metadata
    afsa

    % Products
    products
    products_direct
    products_direct_field
    products_direct_power
    products_specular
    products_specular_reff_coeff
    products_specular_reff_coeff_diel_profiles
    products_specular_reflectivity
    products_specular_reflectivity_diel_profiles
    products_specular_pen_dep
    products_specular_pen_dep_diel_profiles

    % Figures
    fig
    fig_direct
    fig_specular
    fig_specular_reflectivity
    fig_specular_reflectivity_vsEL

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
    sim_name = SimSettings.getInstance.sim_name;
    include_in_master_sim_file = SimSettings.getInstance.include_in_master_sim_file;
    gnd_cover_id = SimSettings.getInstance.gnd_cover_id;


    %% INITIALIZE DIRECTORIES
    % Initialize static simulation directories
    if include_in_master_sim_file

        obj.sims_main_dir = strcat( main_dir, '\sims\master');

    else

        obj.sims_main_dir = strcat( main_dir, '\sims\temp');
    end

    % TO-DO: create a unique name with timestamp
    obj.sim = strcat( obj.sims_main_dir, '\', sim_name);

    % Simulation input
    obj.sim_input = strcat( obj.sim, '\input');
    obj.sim_input_used_files = strcat( obj.sim_input, '\used_files');

    % Meta-Data
    obj.metadata = strcat( obj.sim, '\metadata');

    % Average Forward Scattering Amplitude
    if gnd_cover_id == Constants.ID_VEG_COVER
    
        obj.afsa = strcat(obj.metadata, '\', 'afsa');  
        
    end

    % Products
    obj.products = strcat(obj.sim, '\', 'products') ;
    obj.products_direct = strcat(obj.products, '\', 'direct') ;
    obj.products_direct_field = strcat(obj.products_direct, '\', 'field') ;
    obj.products_direct_power = strcat(obj.products_direct, '\', 'power') ;
    obj.products_specular = strcat(obj.products, '\', 'specular') ;
    obj.products_specular_reff_coeff = strcat(obj.products_specular, '\', 'reflection_coefficient') ;
    obj.products_specular_reflectivity = strcat(obj.products_specular, '\', 'reflectivity') ;
    obj.products_specular_pen_dep = strcat(obj.products_specular, '\', 'pen_dep') ;
    
    % Products for multiple dielectric profiles
    diel_profiles = Constants.DIEL_PROFILES;
    [~, num_diel_profiles] = size( diel_profiles );
    
    for ii = 1 : num_diel_profiles
        
        current_diel_profile = diel_profiles{1, ii};
        obj.products_specular_reff_coeff_diel_profiles{1, ii}  = strcat(obj.products_specular_reff_coeff, '\', current_diel_profile );
        obj.products_specular_reflectivity_diel_profiles{1, ii}  = strcat(obj.products_specular_reflectivity, '\', current_diel_profile );
        obj.products_specular_pen_dep_diel_profiles{1, ii}  = strcat(obj.products_specular_pen_dep, '\', current_diel_profile );
    end
        

    % Figure
    obj.fig = strcat(obj.sim, '\', 'figure' ) ;
    obj.fig_direct = strcat(obj.fig, '\', 'direct') ;
    obj.fig_specular = strcat(obj.fig, '\', 'specular') ;
    obj.fig_specular_reflectivity = strcat(obj.fig_specular, '\', 'reflectivity') ;
    obj.fig_specular_reflectivity_vsEL = strcat(obj.fig_specular_reflectivity, '\', 'vs_EL') ;


    end


    function makeStaticDirs(obj)

        
        %% GET GLOBAL PARAMETERS
        % Ground Parameters
        gnd_structure_id = GndParams.getInstance.gnd_structure_id;
        

        %% Simulations main directory
        if ~exist(obj.sims_main_dir, 'dir')
            mkdir(obj.sims_main_dir);
        end

        %% Input folder that includes sim. report and input file used
        if ~exist(obj.sim_input, 'dir')
            mkdir(obj.sim_input)
        end

        %% Input folder that contains the used input files for the sim.
        if ~exist(obj.sim_input_used_files, 'dir')
            mkdir(obj.sim_input_used_files)
        end

        %% Average Forward Scattering Amplitude
        if ~exist(obj.afsa, 'dir')
            if ~ isempty(obj.afsa)
                mkdir(obj.afsa)
            end
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

        if ~exist(obj.products_specular_reff_coeff, 'dir')
            mkdir(obj.products_specular_reff_coeff)
        end

        if ~exist(obj.products_specular_reflectivity, 'dir')
            mkdir(obj.products_specular_reflectivity)
        end
        
        if ~exist(obj.products_specular_pen_dep, 'dir')
            mkdir(obj.products_specular_pen_dep)
        end
    
        % Products for multiple dielectric profiles
        if gnd_structure_id == Constants.ID_GND_MULTI_LAYERED
            
            diel_profiles = Constants.DIEL_PROFILES;
            [~, num_diel_profiles] = size( diel_profiles );

            for ii = 1 : num_diel_profiles

                if ~exist(obj.products_specular_reff_coeff_diel_profiles{1, ii}, 'dir')

                    mkdir(obj.products_specular_reff_coeff_diel_profiles{1, ii})

                end

                if ~exist(obj.products_specular_reflectivity_diel_profiles{1, ii}, 'dir')

                    mkdir(obj.products_specular_reflectivity_diel_profiles{1, ii})

                end
                
                if ~exist(obj.products_specular_pen_dep_diel_profiles{1, ii}, 'dir')

                    mkdir(obj.products_specular_pen_dep_diel_profiles{1, ii})

                end

            end
            
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

        if ~exist(obj.fig_specular_reflectivity, 'dir')
            mkdir(obj.fig_specular_reflectivity)
        end

        if ~exist(obj.fig_specular_reflectivity_vsEL, 'dir')
            mkdir(obj.fig_specular_reflectivity_vsEL)
        end

    end

    function out = get.sims_main_dir(obj)
        out = obj.sims_main_dir;
    end

    function out = get.sim(obj)
        out = obj.sim;
    end

    function out = get.sim_input(obj)
        out = obj.sim_input;
    end 

    function out = get.sim_input_used_files(obj)
        out = obj.sim_input_used_files;
    end 

    function out = get.afsa(obj)
        out = obj.afsa;
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

    function out = get.products_specular_reff_coeff(obj)
        out = obj.products_specular_reff_coeff;
    end

    function out = get.products_specular_reff_coeff_diel_profiles(obj)
        out = obj.products_specular_reff_coeff_diel_profiles;
    end

    function out = get.products_specular_reflectivity(obj)
        out = obj.products_specular_reflectivity;
    end

    function out = get.products_specular_reflectivity_diel_profiles(obj)
        out = obj.products_specular_reflectivity_diel_profiles;
    end
    
    function out = get.products_specular_pen_dep(obj)
        out = obj.products_specular_pen_dep;
    end
    
    function out = get.products_specular_pen_dep_diel_profiles(obj)
        out = obj.products_specular_pen_dep_diel_profiles;
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

    function out = get.fig_specular_reflectivity(obj)
        out = obj.fig_specular_reflectivity;
    end

    function out = get.fig_specular_reflectivity_vsEL(obj)
        out = obj.fig_specular_reflectivity_vsEL;
    end
        
end
    
end

