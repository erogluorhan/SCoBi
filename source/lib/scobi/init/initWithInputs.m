% Orhan Eroglu
% Jun 2, 2017

function [isInputValid, sim_counter_start] = initWithInputs( varin )

sim_counter_start = -1;

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();

[isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
if isInputValid
    
    % Create but not creating simulations' directories
    SimulationFolders.getInstance.makeStaticDirs();
    
    ParamsManager.copyInputFiles( varin );
    
    ParamsManager.saveInputParams();
    
    sim_counter_start = ParamsManager.determineSimCounterStart();
    
else
    if isTerminate
        uiwait( msgbox(terminateMsg,'Input Error', 'error', 'modal') );
        return
    end
end

end

