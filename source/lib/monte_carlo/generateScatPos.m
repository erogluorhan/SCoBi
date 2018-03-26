%% Mehmet Kurum
% 03/16/2017

function generateScatPos(ind_realization, Nr_current)
% GENERATESCATPOS: Chooses the appropriate method to create a realization 
% of a vegetation field of the selected vegetation_method. 
% Defined over Fresnel Zones
% PARAMS:
% ind_realization: Current index of the realization that is performed
% Nr_current: Existing number of realizations in the simulation directory


%% GET GLOBAL PARAMETERS
% Simulation Parameters
vegetation_method = SimParams.getInstance.vegetation_method;
vegetation_isRow = SimParams.getInstance.vegetation_isRow;


%% DECISION
% Select the proper vegetation function depending on the vegetation_method
if strcmp( vegetation_method, Constants.veg_methods.HOMOGENOUS )
    
    generateHomScatPositions(ind_realization);

elseif strcmp( vegetation_method, Constants.veg_methods.VIRTUAL )
    
    if vegetation_isRow

        generateVirRowScatPositions( ind_realization, Nr_current );
    
    else
        
        generateVirRndScatPositions();
    
    end
    
end

end