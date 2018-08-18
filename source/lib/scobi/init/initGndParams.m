
function initGndParams( inputStruct ) 

% TO-DO: Handle multi-layer ground for SCoBi-Veg
% Get the string diel_model and convert it to an integer index
diel_model_id = findElementIdInCell( Constants.diel_models, inputStruct.diel_model );

% Get the configuration inputs file
configInputFullFile = inputStruct.config_inputs_file;

% Read file
[num, ~, ~] = xlsread( configInputFullFile, 2);

[~, num_gnd_params] = size(num);

ind = 0;

% Ground layer depth (meters)
layer_depth_m = [];
% If number of greater than standard single ground parameter numbsers,
% ground is multi-layer
if num_gnd_params > Constants.num_gnd_single_params
    
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

% Initialize Ground Parameters
GndParams.getInstance.initialize( layer_depth_m, sand_ratio, clay_ratio, rhob_gcm3, diel_model_id );


%% SPECIFIC MULTI-LAYERED GROUND PARAMETERS, IF ANY
% If number of greater than standard single ground parameter numbsers,
% ground is multi-layer
if num_gnd_params > Constants.num_gnd_single_params
    
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

    % Initialize Ground Parameters
    GndMLParams.getInstance.initialize( layer_depth_m, delZ_m, zA_m, zB_m );
    
end

end