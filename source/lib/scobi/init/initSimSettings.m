
function initSimSettings( simulator_id, inputStruct )

campaign = inputStruct.campaign;

sim_mode_id = findElementIdInCell( Constants.sim_modes, inputStruct.sim_mode );

gnd_cover_id = findElementIdInCell( Constants.gnd_covers, inputStruct.gnd_cover );

include_in_master_sim_file = inputStruct.include_in_master_sim_file;

% If gnd_cover is "Vegetation", then get write_attenuation from input
if gnd_cover_id == Constants.id_veg_cover
    
    % Flag to write Attenuation to Excel file 
    write_attenuation = inputStruct.write_attenuation;
    
elseif gnd_cover_id == Constants.id_bare_soil
    
    write_attenuation = 0;
    
end
    
draw_live_plots = 0;

% Initialize Simulation Settings
SimSettings.getInstance.initialize(campaign, simulator_id, sim_mode_id, ...
                gnd_cover_id, write_attenuation, include_in_master_sim_file, draw_live_plots );

end