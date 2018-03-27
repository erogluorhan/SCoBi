classdef Constants
    %CONSTANTS Global constants
    %   Keeps the global constants of the simulation
    
    properties (Constant)
        
        %% GENERAL PURPOSE CONSTANTS (Maths, EM, etc.)
         
        km2m = 1e3;
         
        m2cm = 1e2;
        
        deg2rad = pi / 180 ;
        
        rad2deg = 180 / pi ;
        
        MHz2Hz = 1e6;
        
        % Radius of earth - km -> m via multiplying by 1e3
        re = 6378 * 1e3;
        
        % Speed of light - m/s
        c = 3e8;
        
        version = '1.0';
        

        %% INPUT CONSTANT STRUCTS 
        
        sim_mode = struct('SNAPSHOT', 0, 'TIME_SERIES', 1 );
        
        veg_methods = struct( 'VIRTUAL', 'vir', 'HOMOGENOUS', 'hom' );
        
        veg_vir_types = struct( 'ROW', 'row', 'RANDOM', 'rnd' );
        
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
        
        
        % SIMULATION FOLDERS CONTENTS
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

