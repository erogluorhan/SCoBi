
function initTxParams( inputStruct )
% function initTxParams 
%
%   Fetches the transmitter parameters' values from the input
%   structure and initializes the TxParams class
%
%   initTxParams( inputStruct )
%
%   INPUTS:
%   inputStruct: Input structure that comes from GUI
%
%   See also initAllInputParams, initSimSettings, initRxParams, initGndParams,
%   initConfigParams, initVegParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


% Fetch the transmitter parameters from the input structure
f_MHz = inputStruct.f_MHz;      % Operating frequncy (MHz)

r_Tx_m = inputStruct.r_Tx_km * Constants.KM_TO_M;    % Transmitter range from Earth"s center (km - > m) 

EIRP_dB = inputStruct.EIRP_dB;    % Equivalent Isotropic Radiated Power

pol_Tx = inputStruct.pol_Tx;  % Transmitter polarization                        


% INITIALIZE TRANSMITTER PARAMETERS
TxParams.getInstance.initialize( f_MHz, r_Tx_m, EIRP_dB, pol_Tx );

end