
function initConfigParams( inputStruct )


%% GET GLOBAL PARAMETERS
% Simulation Settings
sim_mode_id = SimSettings.getInstance.sim_mode_id;
% Ground Parameters
num_gnd_layers = GndParams.getInstance.num_layers;


% Get the configuration inputs full file path and name
configInputFullFile = inputStruct.config_inputs_file;

% Read the configuration inputs file
[num, ~, ~] = xlsread( configInputFullFile, 1 );

ind = 0;

% Timestamp - Day-of-year
DoYs = [];

% If sim_mode is Time-series, then input contains timestamps
if sim_mode_id == Constants.id_time_series
    
    ind = ind + 1;
    DoYs = num(:, ind);
    DoYs(any(isnan(DoYs), 2), :) = [];
    
end


% TO-DO: This is for satGeometryManual. There should be an option for
% transmitterGeometry
% Incidence angle
ind = ind + 1;
th0_Tx_list_deg = num(:, ind);    
th0_Tx_list_deg(any(isnan(th0_Tx_list_deg), 2), :) = [];

% Azimuth angle
ind = ind + 1;
ph0_Tx_list_deg = num(:, ind);    
ph0_Tx_list_deg(any(isnan(ph0_Tx_list_deg), 2), :) = [];

% Surface rms height (cm)
ind = ind + 1;
RMSH_list_cm = num(:, ind);       
RMSH_list_cm(any(isnan(RMSH_list_cm), 2), :) = [];

% Volumetric soil moisture
ind = ind + 1;
VSM_list_cm3cm3 = num(:, ind : ind + (num_gnd_layers - 1) );
VSM_list_cm3cm3(any(isnan(VSM_list_cm3cm3), 2), :) = [];


%% PRE-PROCESS TIME-SERIES DATA
% If sim_mode is Time-series, then preprocess data
if sim_mode_id == Constants.id_time_series
    
    num_DoY = length( DoYs );
    num_Th = length( th0_Tx_list_deg );
    num_Ph = length( ph0_Tx_list_deg );
    num_VSM = length( VSM_list_cm3cm3 );
    num_RMSH = length( RMSH_list_cm );
    
    maxnum = max( [num_DoY, num_Th, num_Ph, num_VSM, num_RMSH] );

    % The only validity for time-series input: Length of parameters are
    % equal or one. If one, replicate those.
    if (num_DoY == maxnum || num_DoY == 1) && ...
       (num_Th == maxnum || num_Th == 1) && ...
       (num_Ph == maxnum || num_Ph == 1) && ...
       (num_VSM == maxnum || num_VSM == 1 ) && ...
       (num_RMSH == maxnum || num_RMSH == 1)

        if num_Th == 1
            th0_Tx_list_deg = repmat(th0_Tx_list_deg, 1, maxnum);
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
        th0_Tx_list_deg = [];
        ph0_Tx_list_deg = [];
        VSM_list_cm3cm3 = [];
        RMSH_list_cm = [];
        
    end
    
elseif sim_mode_id == Constants.id_snapshot
    
    DOYs = [];

    % Generate combinations
    combinations = allcomb(th0_Tx_list_deg, ph0_Tx_list_deg, VSM_list_cm3cm3, RMSH_list_cm);
    
    th0_Tx_list_deg = combinations(:, 1);
    ph0_Tx_list_deg = combinations(:, 2);
    VSM_list_cm3cm3 = combinations(:, 3);
    RMSH_list_cm = combinations(:, 4);
    
end


% Initialize Configuration Parameters
ConfigParams.getInstance.initialize( DoYs, th0_Tx_list_deg, ph0_Tx_list_deg, VSM_list_cm3cm3, RMSH_list_cm );

% Set the total number of simulations
ParamsManager.num_sims( length(th0_Tx_list_deg) );

end