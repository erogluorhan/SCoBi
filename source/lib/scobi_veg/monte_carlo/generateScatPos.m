%% Mehmet Kurum
% 03/16/2017

function generateScatPos(ind_realization, Nr_current)
% GENERATESCATPOS: Chooses the appropriate method to create a realization 
% of a vegetation field of the selected veg_method. 
% Defined over Fresnel Zones
% PARAMS:
%   ind_realization: Current index of the realization that is performed
%   Nr_current: Existing number of realizations in the simulation directory


%% GET GLOBAL PARAMETERS
% Simulation Parameters
veg_method_id = SimParams.getInstance.veg_method_id;
veg_vir_orientation_id = SimParams.getInstance.veg_vir_orientation_id;


%% DECISION
% Select the proper vegetation function depending on the veg_method
if veg_method_id == Constants.id_veg_hom
    
    generateHomScatPositions(ind_realization);

elseif veg_method_id == Constants.id_veg_vir
    
    if veg_vir_orientation_id == Constants.id_veg_vir_row_crop

        generateVirRowScatPositions( ind_realization, Nr_current );
    
    elseif veg_vir_orientation_id == Constants.id_veg_vir_random_spread
        
        % TO-DO: Should be implemented when Virtual Random-spread is added
        generateVirRndScatPositions();
    
    end
    
end

end