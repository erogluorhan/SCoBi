classdef Constants
    %CONSTANTS Global constants
    %   Keeps the global constants of the simulation
    
    properties (Constant)
        
        %% GENERAL PURPOSE CONSTANTS (Maths, EM, etc.)
         
        km2m = 1e3;
         
        m2cm = 1e2;
         
        m2mm = 1e3;
        
        deg2rad = pi / 180 ;
        
        rad2deg = 180 / pi ;
        
        MHz2Hz = 1e6;
        
        % Radius of earth - km -> m via multiplying by 1e3
        re = 6378 * 1e3;
        
        % Speed of light - m/s
        c = 3e8;
        
        version = '1.0';
        
        ant_pat_th_range_deg = 180;
        
        ant_pat_ph_range_deg = 360;
        

        %% CELL LISTS 
        % Simulators
        id_veg = 1;
        id_multi_layer = 2;
        simulators = {'SCoBi-Veg', 'SCoBi-ML'};
        
        id_bare_soil = 1;
        id_veg_cover = 2;  
        gnd_covers = {'Bare-soil', 'Vegetation'};
        
        id_veg_hom = 1;
        id_veg_vir = 2;  
        veg_methods = {'Homogenous', 'Virtual'};
        
        id_veg_vir_row_crop = 1;
        id_veg_vir_random_spread = 2;
        veg_vir_orientations = {'Row-crop', 'Random-spread'};
        
        id_pol_R = 1;
        id_pol_L = 2;
        id_pol_X = 3;
        id_pol_Y = 4;
        id_pol_H = 5;
        id_pol_V = 6;
        polarizations = {'R', 'L', 'X', 'Y', 'H', 'V'};
        
        id_Rx_fixed = 1;
        id_Rx_specular_facing = 2;
        Rx_orientations = {'Fixed', 'Specular-facing'};
        
        id_Rx_GG = 1;
        id_Rx_cos_pow_n = 2;
        id_Rx_user_defined = 3;
        Rx_ant_pats = {'Generalized-Gaussian', 'Cosine to the power n', 'User-defined'};
        
        id_mironov = 1;
        id_dobson = 2;
        diel_models = {'Mironov', 'Dobson'};
        
        
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
        
        num_ant_lookup = 1;  % only AntPat.mat
        
        factor_ant_real = 4; % 2x2 (2: Real and imag., 2: Rec and RecI)
        
        num_rot_lookup = 6;  % u_tr, u_gar, u_garI, u_sr, u_ts, u_tIs
        
        factor_rot_real = 4;    % 2x2 (2: Real and imag., 2: Rec and RecI)
        
        num_out_specular = 24;  % 2x2(2x2 + 2) (2: Bare and Veg, 2: Ideal 
                                % or modelled Rec. pattern, 2: real and 
                                % imag. parts, 2: #e_t's, 2: Power terms ) 
        
        factor_frediff_b1 = 16; % 4x2x2 (4: #incoh. mechanisms, 2:#e_t's 
                                % 2: real and imag. parts)
        
    end
    
end

