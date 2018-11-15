
function [ant_pat_struct_Rx, ant_pat_res_deg] = initRxUserDefinedParams( inputStruct )
% function initRxUserDefinedParams 
%
%   Fetches the parameters' values for receivers with a custom antenna 
%   pattern from the inputstructure and initializes the 
%   initRxUserDefinedParams class.
%
%   [ant_pat_struct_Rx, ant_pat_res_deg] = initRxUserDefinedParams( inputStruct )
%
%   INPUTS:
%   inputStruct     : Input structure that comes from GUI
%
%   See also initRxParams, initRxGGParams, initAllInputParams, 
%   initTxParams, initSimSettings, initGndParams, initConfigParams, 
%   initVegParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


% Filename with the full path
ant_pat_Rx_file = inputStruct.ant_pat_Rx_file;    

% Initialize Generalized-Gaussian Receiver Parameters
RxUserDefinedParams.getInstance.initialize( ant_pat_Rx_file );

% Generate the antenna pattern struct
[ant_pat_struct_Rx, ant_pat_res_deg] = RxUserDefinedParams.getInstance.calcRxAntPatMatrix();

end