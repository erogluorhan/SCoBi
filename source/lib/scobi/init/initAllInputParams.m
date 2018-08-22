function varout = initAllInputParams(simulator_id, inputStruct)

varout = [];

%% SIMULATION SETTINGS 
initSimSettings( simulator_id, inputStruct );


% TO-DO: Fix the bugs with the Random-spread sim run
% TO-DO: Handle the Timestamps for the Time-series sim mode
% TO-DO: Handle the Bare-soil case

%% TRANSMITTER (Tx) INPUTS
initTxParams( inputStruct );


%% RECEIVER (Rx) ANTENNA INPUTS
varout{1,1} = initRxParams( inputStruct );


% TO-DO: Diel_model added
%% GROUND INPUTS
initGndParams( inputStruct );


%% CONFIGURATION INPUTS
varout{2,1} = initConfigParams(inputStruct);


%% VEGETATION INPUTS
% It checks if ground cover is Vegetation
varout{3,1} = initVegParams( inputStruct );

end
