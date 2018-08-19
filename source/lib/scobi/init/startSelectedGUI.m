function inputStruct = startSelectedGUI( simulator_id )

inputStruct = [];

% TO-DO: Check this!
% If the OS is not UNIX OR it is MAC and Matlab version below 7.14
if (~isunix || (ismac && verLessThan('matlab', '7.14')))
    
    % If simulator is SCoBi-Veg (Agriculture)
    if simulator_id == Constants.id_veg_agr
        
        inputStruct = gui_SCoBi(simulator_id);
    
    % Else if simulator is SCoBi-Veg (Forest)
    elseif simulator_id == Constants.id_veg_for
        
        inputStruct = gui_SCoBi(simulator_id);
    
    % Else if simulator is SCoBi-ML (MultiLayer)
    elseif simulator_id == Constants.id_multi_layer
        
        inputStruct = gui_SCoBi(simulator_id);
        
    else
        
        % GUIs for new simulators should be added here.
     
    end
    
end

end