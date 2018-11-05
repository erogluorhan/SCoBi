
function ant_pat_struct_Rx = initRxGGParams( inputStruct, ant_pat_res_deg )
% function initRxGGParams 
%
%   Fetches the parameters' values for receivers with Generalized-Gaussian 
%   antenna pattern from the inputstructure and initializes the RxGGParams 
%   class.
%
%   ant_pat_struct_Rx = initRxGGParams( inputStruct, ant_pat_res_deg )
%
%   INPUTS:
%   inputStruct     : Input structure that comes from GUI
%   ant_pat_res_deg : Pattern resolution (Double value in degrees)
%
%   See also initRxParams, initRxUserDefinedParams, initAllInputParams, 
%   initTxParams, initSimSettings, initGndParams, initConfigParams, 
%   initVegParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



% Beamwidth (degrees)
hpbw_deg = inputStruct.hpbw_deg;      

% Sidelobe Level (dB)
SLL_dB = inputStruct.SLL_dB;      

% X-pol level (dB)
XPL_dB = inputStruct.XPL_dB;   

% Initialize Generalized-Gaussian Receiver Parameters
RxGGParams.getInstance.initialize( hpbw_deg, SLL_dB, XPL_dB );

% Generate the antenna pattern struct
ant_pat_struct_Rx = RxGGParams.getInstance.calcRxAntPatMatrix( ant_pat_res_deg );

end