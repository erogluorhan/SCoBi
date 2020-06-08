
function [varout, inputStruct] = initAllInputParams(simulator_id, inputStruct)
% function initAllInputParams 
%
%   Initializes the static parameter classes by using the user inputs.
%   Calls:
%   - initBackwardsControl to initialize inputStruct from a previous
%       version of SCoBi
%   - initSimSettings to initialize SimSettings class
%   - initTxParams to initialize TxParams class
%   - initRxParams to initialize RxParams class
%   - initGndParams to initialize GndParams class
%   - initConfigParams to initialize ConfigParams class
%   - initVegParams to initialize VegParams class
%
%   varout = initAllInputParams(simulator_id, inputStruct)
%
%   INPUTS:
%   simulator_id: Integer simulator ID that comes from GUI
%   inputStruct: Input structure that comes from GUI
%
%   See also initSimSettings, initTxParams, initRxParams, initGndParams,
%   initConfigParams, initVegParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


varout = [];

%% BACKWARDS COMPATIBILITY CONTROL
inputStruct = initBackwardsControl(inputStruct);


%% SIMULATION SETTINGS 
initSimSettings( simulator_id, inputStruct );


%% TRANSMITTER (Tx) INPUTS
initTxParams( inputStruct );


%% RECEIVER (Rx) ANTENNA INPUTS
varout{1,1} = initRxParams( inputStruct );


%% GROUND INPUTS
initGndParams( inputStruct );


%% CONFIGURATION INPUTS
varout{2,1} = initConfigParams(inputStruct);


%% VEGETATION INPUTS
% It checks if ground cover is Vegetation
varout{3,1} = initVegParams( inputStruct );

end
