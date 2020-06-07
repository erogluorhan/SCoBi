
function initGndParams( inputStruct ) 
% function initGndParams 
%
%   Fetches the ground parameters' values from the inputstructure and 
%   initializes the GndParams class
%
%   initGndParams( inputStruct ) 
%
%   INPUTS:
%   inputStruct: Input structure that comes from GUI
%
%   See also initAllInputParams, initSimSettings, initTxParams, initRxParams,
%   initConfigParams, initVegParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


% Fetch the receiver parameters from the input structure
% Get the string diel_model and convert it to an integer index
diel_model_id = findElementIdInCell( Constants.DIEL_MODELS, inputStruct.diel_model );

% Get the ground kayer structure: Single or Multi-layered
gnd_structure_id = findElementIdInCell( Constants.GND_STRUCTURES, inputStruct.gnd_structure );

% Get the configuration inputs file
configInputFullFile = inputStruct.config_inputs_file;

% Read file
[num, ~, ~] = xlsread( configInputFullFile, 2);

ind = 0;

% Ground layer depth (meters)
layer_depth_m = [];
% If number of greater than standard single ground parameter numbsers,
% ground is multi-layer
if gnd_structure_id == Constants.ID_GND_MULTI_LAYERED
    
    ind = ind + 1;
    layer_depth_m = num(:, ind);
    layer_depth_m( isnan(layer_depth_m) ) = [];
    
end

% Sand ratio of the soil texture
ind = ind + 1;
sand_ratio = num(:, ind);
sand_ratio( isnan(sand_ratio) ) = [];

% Clay ratio of the soil texture 
ind = ind + 1;
clay_ratio = num(:, ind);
clay_ratio( isnan(clay_ratio) ) = [];

% Soil bulk density
ind = ind + 1;
rhob_gcm3 = num(:, ind);
rhob_gcm3( isnan(rhob_gcm3) ) = []; 


% INITIALIZE GROUND PARAMETERS
GndParams.getInstance.initialize( gnd_structure_id, sand_ratio, clay_ratio, rhob_gcm3, diel_model_id );



%% SPECIFIC MULTI-LAYERED GROUND PARAMETERS, IF ANY
% If number of greater than standard single ground parameter numbsers,
% ground is multi-layer
if gnd_structure_id == Constants.ID_GND_MULTI_LAYERED    
    
    % Get the configuration inputs file
    calculate_penetration_depth = inputStruct.calculate_penetration_depth;
    
    % Flags to calculate several dielectric profile fit functions
    calc_diel_profile_fit_functions = inputStruct.calc_diel_profile_fit_functions;  
    
    % Initialize workspace for ground-multilayer
    initWSMultilayer();
    
    % Layer discritization
    ind = ind + 1;
    delZ_m = num(1, ind);
    delZ_m( isnan(delZ_m) ) = []; 

    % Air layer
    ind = ind + 1;
    zA_m = num(1, ind);
    zA_m( isnan(zA_m) ) = []; 

    % Bottom-most layer
    ind = ind + 1;
    zB_m = num(1, ind);
    zB_m( isnan(zB_m) ) = []; 

    % INITIALIZE GROUND MULTI-LAYER PARAMETERS
    GndMLParams.getInstance.initialize( layer_depth_m, delZ_m, zA_m, zB_m, ...
        calc_diel_profile_fit_functions, calculate_penetration_depth );
end

end