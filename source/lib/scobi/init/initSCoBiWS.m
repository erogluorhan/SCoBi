function initSCoBiWS( simulator_id )

% Add all subdirectories to the path
addpath( genpath( Directories.getInstance.input ) );

% If SCoBi-Veg (Agriculture)
if simulator_id == Constants.id_veg_agr

    initSCoBiVegWS();

% Else if SCoBi-Veg (Forest)    
elseif simulator_id == Constants.id_veg_for

    initSCoBiVegWS();
    
% Else if SCoBi-ML (MutiLayer)
elseif simulator_id == Constants.id_multi_layer

    initSCoBiMLWS();
    
else
    
    % Workspace initialization for new simulators should be done here!
    
end

end


% Initialize SCoBi-Veg workspace
function initSCoBiVegWS

addpath( genpath( Directories.getInstance.scobi_veg ) );

end


% Initialize SCoBi-ML workspace
function initSCoBiMLWS

addpath( genpath( Directories.getInstance.scobi_ml ) );

end