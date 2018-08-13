
function initSimParams( inputStruct )


%% GET GLOBAL PARAMETERS
% Simulation Settings
simulator_id = SimSettings.getInstance.simulator_id;
sim_mode_id = SimSettings.getInstance.sim_mode_id;
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;


% Simulation name
sim_name = inputStruct.sim_name;

% Simulation campaign
campaign = inputStruct.campaign; 

% Simulation campaign date
campaign_date = inputStruct.campaign_date; 

% Simulation campaign plot
plot = inputStruct.plot;

% Vegetation method
veg_method_id = [];
% If gnd_cover is Vegetation, then veg_method is taken from inputs
if gnd_cover_id == Constants.id_veg_cover
    
    if simulator_id == Constants.id_veg_agr ...
            || simulator_id == Constants.id_veg_for
        
        % If sim_mode is Snapshot, then veg_method depends on input
        if sim_mode_id == Constants.id_snapshot

            veg_method_id = findElementIdInCell( Constants.veg_methods, inputStruct.veg_method );

        % Else if sim_mode is Time-series, veg_method can only be Homogenous
        elseif sim_mode_id == Constants.id_time_series

            veg_method_id = Constants.id_veg_hom;

        end
        
    elseif simulator_id == Constants.id_multi_layer
        
        veg_method_id = Constants.id_veg_hom;
        
    end
    
end    

% Virtual Vegetation Orientation
veg_vir_orientation_id = [];
% If veg_method is 'Virtual' only
if veg_method_id == Constants.id_veg_vir
    veg_vir_orientation_id = findElementIdInCell( Constants.veg_vir_orientations, inputStruct.veg_vir_orientation );
end

% Vegetation plant name
vegetation_plant = [];

% If gnd_cover is Vegetation, then veg_plant is used and taken from inputs
if gnd_cover_id == Constants.id_veg_cover
    
    vegetation_plant = inputStruct.veg_plant;    
    
end

% Number of Realizations
Nr = inputStruct.Nr;

% Number of Fresnel Zones
Nfz = inputStruct.Nfz;

% Initialize Simulation Parameters
SimParams.getInstance.initialize(sim_name, campaign, campaign_date, ...
                plot, veg_method_id, veg_vir_orientation_id, vegetation_plant, Nr, Nfz );

end