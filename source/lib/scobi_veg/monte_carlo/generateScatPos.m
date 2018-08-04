%% Mehmet Kurum
% 03/16/2017

function generateScatPos( Nr_current )
% GENERATESCATPOS: Chooses the appropriate method to create a realization 
% of a vegetation field of the selected veg_method. 
% Defined over Fresnel Zones
% PARAMS:
%   Nr_current: Existing number of realizations in the simulation directory


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nr = SimParams.getInstance.Nr;
veg_method_id = SimParams.getInstance.veg_method_id;
veg_vir_orientation_id = SimParams.getInstance.veg_vir_orientation_id;


%% REALIZATIONS
for rr = Nr_current + 1 : Nr % Number of Realization
    
    
    %% DECISION
    % Select the apropriate vegetation function depending on the veg_method
    if veg_method_id == Constants.id_veg_hom

        generateHomScatPositions(rr);

    elseif veg_method_id == Constants.id_veg_vir

        % Select the apropriate virtual vegetation orientation:Row-crop OR
        % Random-spread
        if veg_vir_orientation_id == Constants.id_veg_vir_row_crop

            generateVirRowScatPositions( rr, Nr_current );

        elseif veg_vir_orientation_id == Constants.id_veg_vir_random_spread

            % TO-DO: Should be implemented when Virtual Random-spread is added
            generateVirRndScatPositions();

        end

    end
    
end % Realization

end