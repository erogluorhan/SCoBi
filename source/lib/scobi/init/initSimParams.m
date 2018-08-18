
function initSimParams( inputStruct )


%% GET GLOBAL PARAMETERS
% Simulation Settings
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;


% Simulation name
sim_name = inputStruct.sim_name;

% Simulation campaign
campaign = inputStruct.campaign; 

% Simulation campaign date
campaign_date = inputStruct.campaign_date; 

% Simulation campaign plot
plot = inputStruct.plot;
  
% Vegetation plant name
% If gnd_cover is Vegetation, then veg_plant is used and taken from inputs
if gnd_cover_id == Constants.id_veg_cover
    
    vegetation_plant = inputStruct.veg_plant; 

elseif gnd_cover_id == Constants.id_bare_soil
    
    vegetation_plant = [];
    
end

% Number of Realizations
Nr = inputStruct.Nr;

% Number of Fresnel Zones
Nfz = inputStruct.Nfz;

% Initialize Simulation Parameters
SimParams.getInstance.initialize(sim_name, campaign, campaign_date, ...
                plot, vegetation_plant, Nr, Nfz );

end