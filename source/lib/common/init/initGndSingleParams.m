
function initGndSingleParams( inputStruct ) 


% Get the string diel_model and convert it to an integer index
diel_model_id = findElementIdInCell( Constants.diel_models, inputStruct.diel_model );

% Sand ratio of the soil texture
sand_ratio = inputStruct.sand_ratio;

% Clay ratio of the soil texture 
clay_ratio = inputStruct.clay_ratio;  

% Soil bulk density   
rhob_gcm3 = inputStruct.rhob_gcm3;   

% Default ground is a single-layered infinite soil. No need for layer depth
layer_depth_m = [];

% Initialize Ground Parameters
GndParams.getInstance.initialize( layer_depth_m, sand_ratio, clay_ratio, rhob_gcm3, diel_model_id );

end