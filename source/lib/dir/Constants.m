classdef Constants
    %CONSTANTS Global constants
    %   Keeps the global constants of the simulation
    
    properties (Constant)
        
        km2m = 1e3;
        
        deg2rad = pi / 180 ;
        
        rad2deg = 180 / pi ;
        
        MHz2Hz = 1e6;
        
        % Radius of earth - km -> m via multiplying by 1000
        re = 6378 * 1e3;
        
        % Speed of light - m/s
        c = 3e8;   
        
        version = '1.0';
        

        % INPUT CONSTANT STRUCTS 
        
        sim_mode = struct('CROSS', 0, 'TIME_SERIES', 1 );
        
        veg_methods = struct( 'VIRTUAL', 'vir', 'HOMOGENOUS', 'hom' );
        
        veg_vir_types = struct( 'ROW', 'row', 'RANDOM', 'rnd' );
        
        % To determine the need for a function run in mainSCoBi
        need_for_run = struct('NO', 0, 'PARTIAL', 1, 'FULL', 2 );
        
        particleDataStruct = struct('Type', 1, 'Kind', 2, ...
            'posX', 3, 'posY', 4, 'posZ', 5, ...
            'downAngle', 6, 'azimuthAngle', 7, ...
            'dim1', 8, 'dim2', 9, 'dim3', 10, ...
            'epsrRe', 11, 'epsrIm', 12, 'fzIndex', 13) ;

        % Below struct will be used where saved scatterer positions are used
        scattererParamsStruct = struct('posX', 1, 'posY', 2, 'posZ', 3, ...
            'downAngle', 4, 'azimuthAngle', 5, ...
            'dim1', 6, 'dim2', 7, 'dim3', 8, ...
            'epsrRe', 9, 'epsrIm', 10) ;
        
        
        % SIMULATION FOLDERS CONTENTS
        num_afsa = 10;
        factor_fscat = 8;
        num_ant_lookup = 1;
        factor_ant_real = 4;
        num_rot_lookup = 6;
        factor_rot_real = 4;
        num_out_specular = 24;
        factor_frediff_b1 = 16; % 4: #incoherent terms, 2:#t's 2: imag and real parts
        
    end
    
end

