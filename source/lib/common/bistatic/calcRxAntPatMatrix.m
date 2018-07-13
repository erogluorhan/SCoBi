% Mehmet Kurm
% March 2, 2017

function calcRxAntPatMatrix


%% GET GLOBAL PARAMETERS
% Receiver Parameters
ant_pat_Rx_id = RxParams.getInstance.ant_pat_Rx_id;


%% DECISION
% If the antenna pattern method is Generalized-Gaussian
if ant_pat_Rx_id == Constants.id_Rx_GG
    
    calcRxAntPatMatrixGG();

% Else if it is Cosine to the power n 
elseif ant_pat_Rx_id == Constants.id_Rx_cos_pow_n

    % TO-DO: Implement here when cosine-to-the-power-n antenna pattern 
    % method is added 
    
% Else if it is a user-defined (custom) antenna pattern   
elseif ant_pat_Rx_id == Constants.id_Rx_user_defined
    
    calcRxAntPatMatrixUserDefined();
    
end

end