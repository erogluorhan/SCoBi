
function initGndParams( simulator_id, inputStruct ) 

% TO-DO: Handle multi-laer ground for SCoBi-Veg

% If the simulator is SCoBi-Veg (Agriculture OR Forest), set the single
% ground layer parameters (For now only! It can be multi-layer as well in
% the future)
if simulator_id == Constants.id_veg_agr ...
        || simulator_id == Constants.id_veg_for
    
    initGndSingleParams( inputStruct );
    
% Else if simulator is SCoBi-ML (Multilayer), then set the multi-layer
% ground inputs
elseif simulator_id == Constants.id_multi_layer
    
    initGndMLParams( inputStruct );
    
end

end