classdef Constants
    %CONSTANTS Global constants
    %   Keeps the global constants of the simulation
    
    properties (Constant)
        
        %% GENERAL PURPOSE CONSTANTS (Maths, EM, etc.)
         
        km2m = 1e3;
         
        m2cm = 1e2;
         
        m2mm = 1e3;
        
        MHz2Hz = 1e6;
        
        GHz2Hz = 1e9;
        
        % Radius of earth - km -> m via multiplying by 1e3
        re = 6378 * 1e3;
        
        % Speed of light - m/s
        c = 3e8;
        
        % Air Dielectric Constant
        eps_diel_air = 1.0 - 0.0 * 1i ;
        
        version = '1.0';
        
        ant_pat_th_range_deg = 180;
        
        ant_pat_ph_range_deg = 360;
        

        %% CELL LISTS         
        % Simulators
        id_veg_agr = 1;
        id_veg_for = 2;
        id_multi_layer = 3;
        id_soil = 4;
        simulators = {'Agriculture', 'Forest', 'Root-zone', 'Soil'};
        
        % lastInputFileNames
        lastInputFileNames = {'last_input-agriculture.mat', 'last_input-forest.mat', 'last_input-root_zone.mat', 'last_input-soil.mat'};
        
        % defaultInputFileNames
        defaultInputFileNames = {'default-agriculture.mat', 'default-forest.mat', 'default-root_zone.mat', 'default-soil.mat'};
        
        % Sim-Modes
        id_snapshot = 1;
        id_time_series = 2;
        sim_modes = {'Snapshot', 'Time-series'};
        
        % Ground-cover
        id_bare_soil = 1;
        id_veg_cover = 2;  
        gnd_covers = {'Bare-soil', 'Vegetation'};
        
        % Polarizations
        id_pol_R = 1;
        id_pol_L = 2;
        id_pol_X = 3;
        id_pol_Y = 4;
        id_pol_H = 5;
        id_pol_V = 6;
        polarizations = {'R', 'L', 'X', 'Y', 'H', 'V'};
        
        % Receiver antenna orientations
        id_Rx_fixed = 1;
        id_Rx_specular_facing = 2;
        Rx_orientations = {'Fixed', 'Specular-facing'};
        
        % Receiver antenna patterns
        id_Rx_GG = 1;
        id_Rx_user_defined = 2;
        id_Rx_cos_pow_n = 3;
        Rx_ant_pats = {'Generalized-Gaussian', 'User-defined', 'Cosine to the power n'};
        
        % Dielcetric models
        id_diel_dobson = 1;
        id_diel_mironov = 2;
        id_diel_wang = 3;
        diel_models = {'Dobson', 'Mironov', 'Wang'};
        
        % Dielectric profile generation methods 
        id_diel_slab = 1;
        id_diel_logistic = 2;
        id_diel_2nd_order = 3;
        id_diel_3rd_order = 4;
        diel_profiles = {'Discrete-slab', 'Logistic', '2nd-order', '3rd-order'};
        
        % Receiver antenna orientations
        id_GUI_save = 1;
        id_GUI_save_as = 2;
        GUI_save_options = {'Save', 'Save-as'};
        
        
        %% STRUCTS
        % To determine the need for a function run in mainSCoBi
        need_for_run = struct('NO', 0, 'PARTIAL', 1, 'FULL', 2 );
        
        % For use of virtual vegetation
        particleDataStruct = struct('Type', 1, 'Kind', 2, ...
            'posX', 3, 'posY', 4, 'posZ', 5, ...
            'downAngle', 6, 'azimuthAngle', 7, ...
            'dim1', 8, 'dim2', 9, 'dim3', 10, ...
            'epsrRe', 11, 'epsrIm', 12, 'fzIndex', 13) ;

        % For use of virtual vegetation
        scattererParamsStruct = struct('posX', 1, 'posY', 2, 'posZ', 3, ...
            'downAngle', 4, 'azimuthAngle', 5, ...
            'dim1', 6, 'dim2', 7, 'dim3', 8, ...
            'epsrRe', 9, 'epsrIm', 10) ;
        
        
        %% SIMULATION FOLDERS CONTENT COUNTS
        num_afsa = 10;       % Except the Attenuation.xls
        
        factor_fscat = 8;    % 4x2 (4: # mechanisms, 2: real and imag.)
        
        factor_ant_real = 4; % 2x2 (2: Real and imag., 2: Rec and RecI)
        
        num_rot_lookup = 4;  % u_tr, u_sr, u_ts, u_tIs
        
        num_products_specular_bare = 8;  % 2x2x2 (2: Ideal or sim. ant. pat.
                                    %        2: real and imag. parts, 
                                    %        2: #e_t's ) 
                                    
        num_products_specular_veg = 16; % 2x2x2x2 (2: Ideal/ sim. ant. pat.
                                    %              2: real and imag. parts, 
                                    %              2: #e_t's 
                                    %              2: Bare and vegetation) 
                                    
                                    
        %% INPUTS
        num_gnd_single_params = 3;  % sand_ratio, clay_ratio, and rhob_gcm3
        
    end
    
end

