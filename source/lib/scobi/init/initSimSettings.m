
function initSimSettings( simulator_id, inputStruct )
% function initSimSettings 
%
%   Fetches the simulation setting parameters' values from the input
%   structure and initializes the SimSettings class
%
%   initSimSettings( simulator_id, inputStruct )
%
%   INPUTS:
%   simulator_id: Integer simulator ID that comes from GUI
%   inputStruct: Input structure that comes from GUI
%
%   See also initAllInputParams, initTxParams, initRxParams, initGndParams,
%   initConfigParams, initVegParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


% Fetch the simulation settings parameters from the input structure
campaign = inputStruct.campaign;

sim_mode_id = findElementIdInCell( Constants.SIM_MODES, inputStruct.sim_mode );

gnd_cover_id = findElementIdInCell( Constants.GND_COVERS, inputStruct.gnd_cover );

include_in_master_sim_file = inputStruct.include_in_master_sim_file;

% If gnd_cover is "Vegetation", then get write_attenuation from input
if gnd_cover_id == Constants.ID_VEG_COVER
    
    % Flag to write Attenuation to Excel file 
    write_attenuation = inputStruct.write_attenuation;
    
elseif gnd_cover_id == Constants.ID_BARE_SOIL
    
    write_attenuation = 0;
    
end

% INITIALIZE SIMULATION SETTINGS
SimSettings.getInstance.initialize( ...
    campaign, simulator_id, sim_mode_id, ...
    gnd_cover_id, write_attenuation, include_in_master_sim_file);

end