
function configInputFullFile = initConfigParams( inputStruct )
% function initConfigParams 
%
%   Fetches the configuration parameters' values from the input structure 
%   and initializes the ConfigParams class. It implements preprocessing of
%   data beyond intialization (such as generating the combination of 
%   Snapshot simulation mode).
%
%   configInputFullFile = initConfigParams( inputStruct )
%
%   INPUTS:
%   inputStruct: Input structure that comes from GUI
%
%   See also initAllInputParams,  initSimSettings, initTxParams, 
%   initRxParams, initGndParams, initVegParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



%% GET GLOBAL PARAMETERS
% Simulation Settings
sim_mode_id = SimSettings.getInstance.sim_mode_id;
% Ground Parameters
gnd_structure_id = GndParams.getInstance.gnd_structure_id;


% Get the configuration inputs full file path and name
configInputFullFile = inputStruct.config_inputs_file;

% Read the configuration inputs file
[num, ~, ~] = xlsread( configInputFullFile, 1 );

ind = 0;

% Timestamp - Day-of-year
DoYs = [];

% If sim_mode is Time-series, then input contains timestamps
if sim_mode_id == Constants.ID_TIME_SERIES
    
    ind = ind + 1;
    DoYs = num(:, ind);
    DoYs(any(isnan(DoYs), 2), :) = [];
    
end


orientation_Tx_id = findElementIdInCell( Constants.TX_ORIENTATIONS, inputStruct.orientation_Tx );
% If transmitter orientation is variable, get incidence and azimuth angles
% from configuration inputs file
if orientation_Tx_id == Constants.ID_TX_VARIABLE
    
    % Incidence angle
    ind = ind + 1;
    el0_Tx_list_deg = num(:, ind);    
    el0_Tx_list_deg(any(isnan(el0_Tx_list_deg), 2), :) = [];

    % Azimuth angle
    ind = ind + 1;
    ph0_Tx_list_deg = num(:, ind);    
    ph0_Tx_list_deg(any(isnan(ph0_Tx_list_deg), 2), :) = [];
    
% Else if transmitter orientation is Geo-stationary, get incidence and 
% azimuth angles from inputStruct
elseif orientation_Tx_id == Constants.ID_TX_GEOSTATIONARY
    
    % Incidence angle
    el0_Tx_list_deg = inputStruct.el0_Tx_deg;

    % Azimuth angle
    ph0_Tx_list_deg = inputStruct.ph0_Tx_deg;    
    
end

% Surface rms height (cm)
ind = ind + 1;
RMSH_list_cm = num(:, ind);       
RMSH_list_cm(any(isnan(RMSH_list_cm), 2), :) = [];

if gnd_structure_id == Constants.ID_GND_SINGLE_LAYERED
    
    num_gnd_layers = 1;
    
elseif gnd_structure_id == Constants.ID_GND_MULTI_LAYERED
   
    %% GET GLOBAL PARAMETERS
    num_gnd_layers = GndMLParams.getInstance.num_layers; 
    
end
% Volumetric soil moisture
ind = ind + 1;
VSM_list_cm3cm3 = num(:, ind : ind + (num_gnd_layers - 1) );
VSM_list_cm3cm3(any(isnan(VSM_list_cm3cm3), 2), :) = [];


%% PRE-PROCESS TIME-SERIES DATA
% If sim_mode is Time-series, then preprocess data
if sim_mode_id == Constants.ID_TIME_SERIES
    
    [num_DoY, ~]    = size( DoYs );
    [num_Th, ~]     = size( th0_Tx_list_deg );
    [num_Ph, ~]     = size( ph0_Tx_list_deg );
    [num_VSM, ~]    = size( VSM_list_cm3cm3 );
    [num_RMSH, ~]   = size( RMSH_list_cm );
    
    maxnum = max( [num_DoY, num_El, num_Ph, num_VSM, num_RMSH] );

    % The only validity for time-series input: Length of parameters are
    % equal or one. If one, replicate those.
    if (num_DoY == maxnum || num_DoY == 1) && ...
       (num_El == maxnum || num_El == 1) && ...
       (num_Ph == maxnum || num_Ph == 1) && ...
       (num_VSM == maxnum || num_VSM == 1 ) && ...
       (num_RMSH == maxnum || num_RMSH == 1)

        if num_El == 1
            el0_Tx_list_deg = repmat(el0_Tx_list_deg, 1, maxnum);
        end

        if num_DoY == 1
            DoYs = repmat(DoYs, 1, maxnum);
        end

        if num_Ph == 1
            ph0_Tx_list_deg = repmat(ph0_Tx_list_deg, 1, maxnum);
        end

        if num_VSM == 1
            VSM_list_cm3cm3 = repmat(VSM_list_cm3cm3, 1, maxnum);
        end     

        if num_RMSH == 1
            RMSH_list_cm = repmat(RMSH_list_cm, 1, maxnum);
        end 

    % input lengths are not valid for time-series, make all empty
    else
        
        DoYs = [];
        el0_Tx_list_deg = [];
        ph0_Tx_list_deg = [];
        VSM_list_cm3cm3 = [];
        RMSH_list_cm = [];
        
    end
    
elseif sim_mode_id == Constants.ID_SNAPSHOT
    
    DOYs = [];

    % Generate combinations
    combinations = allcomb(el0_Tx_list_deg, ph0_Tx_list_deg, VSM_list_cm3cm3, RMSH_list_cm);
    
    % Assign combination columns to corresponding parameters
    el0_Tx_list_deg = combinations(:, 1);
    ph0_Tx_list_deg = combinations(:, 2);
    VSM_list_cm3cm3 = combinations(:, 3:(3+num_gnd_layers-1));
    RMSH_list_cm = combinations(:, end);
    
end


% INITIALIZE RECEIVER PARAMETERS
ConfigParams.getInstance.initialize( DoYs, el0_Tx_list_deg, ph0_Tx_list_deg, VSM_list_cm3cm3, RMSH_list_cm );


% SET THE TOTAL NUMBER OF SIMULATIONS
% Simply, length of el0_Tx_list_deg can be used (But not DoYs since it may 
% be empty for Snapshot simulations)
ParamsManager.num_sims( length(el0_Tx_list_deg) );

end