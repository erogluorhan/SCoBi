%% Mehmet Kurum
% 03/16/2017

function generateScatPos( Nr_current )
% GENERATESCATPOS: Creates realizations of a vegetation layer homogenously 
% Defined over Fresnel Zones
% PARAMS:
%   Nr_current: Existing number of realizations in the simulation directory


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nr = SimParams.getInstance.Nr;


%% REALIZATIONS
for rr = Nr_current + 1 : Nr % Number of Realization
    
    generateHomScatPositions(rr);
    
end % Realization

end