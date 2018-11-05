classdef Constants
% class Constants
%
%   Global constant variables. This class only consists of constant (final) 
%   properties that keep the global values for any purpose.   
%
%   See also ConstantNames, Directories.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0

    
    properties (Constant)
        
        
        %% SCOBI VERSION        
        
        VERSION = '1.0.0';
        
        URL_MSSTATE = 'https://www.msstate.edu';
        
        URL_IMPRESS = 'http://impress.ece.msstate.edu/';
        
        URL_SCOBI_GITHUB = 'https://github.com/impresslab/SCoBi';
        
        URL_SCOBI_DESIGN = 'https://github.com/impresslab/SCoBi/blob/master/design/UML-EA/SCoBi-v1_0.EAP';
        
        URL_SCOBI_USER_MANUAL = 'https://github.com/impresslab/SCoBi/blob/master/docs/manuals/SCoBi-User_Manual-v1_0.pdf';
        
        URL_SCOBI_DEVELOPER_MANUAL = 'https://github.com/impresslab/SCoBi/blob/master/docs/manuals/SCoBi-Developer_Manual-v1_0.pdf';
        
        
        %% GENERAL PURPOSE CONSTANTS (Maths, EM, etc.)
         
        KM_TO_M = 1e3;
         
        M_TO_CM = 1e2;
        
        MHZ_TO_HZ = 1e6;
        
        GHZ_TO_HZ = 1e9;
        
        % Radius of earth - km -> m via multiplying by 1e3
        R_EARTH = 6378 * 1e3;
        
        % Speed of light - m/s
        LIGHTSPEED = 3e8;
        
        % Air Dielectric Constant
        EPS_DIEL_AIR = 1.0 - 0.0 * 1i ;
        
        ANT_PAT_TH_RANGE_DEG = 180;
        
        ANT_PAT_PH_RANGE_DEG = 360;
        
        RX_ANT_PAT_GG_HPBW_MAX = 60;    % degrees
        
        RX_ANT_PAT_GG_SLL_MIN = 10;     % dB
        RX_ANT_PAT_GG_SLL_MAX = 40;     % dB
        
        RX_ANT_PAT_GG_XPL_MIN = 10;     % dB
        RX_ANT_PAT_GG_XPL_MAX = 40;     % dB
        

        %% CELL LISTS         
        % Simulators
        ID_SIM_AGRICULTURE = 1;
        ID_SIM_FOREST = 2;
        ID_SIM_ROOT_ZONE = 3;
        ID_SIM_SOIL = 4;
        SIMULATORS = {'Agriculture', 'Forest', 'Root-zone', 'Soil'};
        
        % lastInputFileNames
        LAST_INPUT_FILENAMES = {'last_input-agriculture.mat', 'last_input-forest.mat', 'last_input-root_zone.mat', 'last_input-soil.mat'};
        
        % defaultInputFileNames
        DEFAULT_INPUT_FILENAMES = {'default-agriculture.mat', 'default-forest.mat', 'default-root_zone.mat', 'default-soil.mat'};
        
        % Sim-Modes
        ID_SNAPSHOT = 1;
        ID_TIME_SERIES = 2;
        SIM_MODES = {'Snapshot', 'Time-series'};
        
        % Ground-cover
        ID_BARE_SOIL = 1;
        ID_VEG_COVER = 2;  
        GND_COVERS = {'Bare-soil', 'Vegetation'};
        
        % Polarizations
        ID_POL_R = 1;
        ID_POL_L = 2;
        ID_POL_X = 3;
        ID_POL_Y = 4;
        ID_POL_H = 5;
        ID_POL_V = 6;
        POLARIZATIONS = {'R', 'L', 'X', 'Y', 'H', 'V'};
        
        % Transmitter orientations
        ID_TX_GEOSTATIONARY = 1;
        ID_TX_VARIABLE = 2;
        TX_ORIENTATIONS = {'Geo-stationary', 'Variable'};
        
        % Receiver antenna orientations
        ID_RX_FIXED = 1;
        ID_RX_SPECULAR_FACING = 2;
        RX_ORIENTATIONS = {'Fixed', 'Specular-facing'};
        
        % Receiver antenna patterns
        ID_RX_GG = 1;
        ID_RX_USER_DEFINED = 2;
        ID_RX_COS_POW_N = 3;
        RX_ANT_PATTERNS = {'Generalized-Gaussian', 'User-defined', 'Cosine to the power n'};
        
        % Dielcetric models
        ID_DIEL_DOBSON = 1;
        ID_DIEL_MIRONOV = 2;
        ID_DIEL_WANG = 3;
        DIEL_MODELS = {'Dobson', 'Mironov', 'Wang'};
        
        % Dielectric profile generation methods 
        ID_DIEL_PROFILE_2ND_ORDER = 1;
        ID_DIEL_PROFILE_3RD_ORDER = 2;
        ID_DIEL_PROFILE_LOGISTIC = 3;
        ID_DIEL_PROFILE_SLAB = 4;
        DIEL_PROFILES = {'2nd-order', '3rd-order', 'Logistic', 'Discrete-slab' };
        
        % Ground layer structure
        ID_GND_SINGLE_LAYERED = 1;
        ID_GND_MULTI_LAYERED = 2;
        GND_STRUCTURES = {'Single-layered', 'Multi-layered'};
        
        % Receiver antenna orientations
        ID_GUI_SAVE = 1;
        ID_GUI_SAVE_AS = 2;
        GUI_SAVE_OPTIONS = {'Save', 'Save-as'};
        
        
        %% STRUCTS
        % To determine the need for a function run in mainSCoBi
        NEED_TO_RUN_STRUCT = struct('NO', 0, 'PARTIAL', 1, 'FULL', 2 );
        
    end
    
end

