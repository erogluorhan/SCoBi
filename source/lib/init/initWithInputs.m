% Orhan Eroglu
% Jun 2, 2017

function isInputValid = initWithInputs

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();

[isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
if isInputValid
    
    % Create but not creating simulations' directories
    SimulationFolders.getInstance.makeStaticDirs();
    
    ParamsManager.saveSimParams();
    
else
    if isTerminate
        uiwait( msgbox(terminateMsg,'Input Error', 'error', 'modal') );
        return
    end
end

end

