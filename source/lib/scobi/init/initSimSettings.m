
function initSimSettings( simulator_id, inputStruct )

sim_mode_id = findElementIdInCell( Constants.sim_modes, inputStruct.sim_mode );

gnd_cover_id = findElementIdInCell( Constants.gnd_covers, inputStruct.gnd_cover );

if simulator_id == Constants.id_veg_agr ...
        || simulator_id == Constants.id_veg_for
    
    % Flag to write Attenuation to Excel file 
    write_attenuation = inputStruct.write_attenuation;      

    % Flag to calculate direct term  
    calc_direct_term = inputStruct.calc_direct_term;

    % Flag to calculate specular term    
    calc_specular_term = inputStruct.calc_specular_term; 

    % Flag to calculate diffuse term
    % If sim_mode is Snapshot AND gnd_cover is Vegetation, then get the flag from inputs
    if sim_mode_id == Constants.id_snapshot ...
            && gnd_cover_id == Constants.id_veg_cover

        calc_diffuse_term = inputStruct.calc_diffuse_term;

    % Else, no way for diffuse calculations
    else

        calc_diffuse_term = 0;

    end
    
    draw_live_plots = 0;

elseif simulator_id == Constants.id_multi_layer
    
    calc_direct_term = 0;
    calc_specular_term = 1;
    calc_diffuse_term = 0;
    write_attenuation = 0;
    
    draw_live_plots = 1;
    
    if sim_mode_id == Constants.id_snapshot
        draw_live_plots = 0;
    end
    
end

% Initialize Simulation Settings
SimSettings.getInstance.initialize(simulator_id, sim_mode_id, ...
                gnd_cover_id, write_attenuation, calc_direct_term, ...
                calc_specular_term, calc_diffuse_term, draw_live_plots );

end