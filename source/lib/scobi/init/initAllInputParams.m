function initAllInputParams(simulator_id, inputStruct)

%% SIMULATION SETTINGS 
initSimSettings( simulator_id, inputStruct );


% TO-DO: Fix the bugs with the Random-spread sim run
% TO-DO: Handle the Timestamps for the Time-series sim mode
% TO-DO: Handle the Bare-soil case
%% SIMULATION INPUTS 
initSimParams( inputStruct );


%% TRANSMITTER (Tx) INPUTS
initTxParams( inputStruct );


%% RECEIVER (Rx) ANTENNA INPUTS
initRxParams( inputStruct );


% TO-DO: Diel_model added
%% GROUND INPUTS
initGndParams( inputStruct );


%% CONFIGURATION INPUTS
initConfigParams(inputStruct);


%% VEGETATION INPUTS
% It checks if ground cover is Vegetation
initVegParams( inputStruct );

end
