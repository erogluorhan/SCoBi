
function ant_pat_fullfile = initRxParams( inputStruct )
% function initRxParams 
%
%   Fetches the receiver parameters' values from the inputstructure and 
%   initializes the RxParams class. It implements interpretation of
%   inputs beyond intialization (such as deciding which antenna pattern is 
%   used and initializng that as well).
%
%   ant_pat_fullfile = initRxParams( inputStruct )
%
%   INPUTS:
%   inputStruct: Input structure that comes from GUI
%
%   See also initRxGGParams, initRxUserDefinedParams, initAllInputParams, 
%   initTxParams, initSimSettings, initGndParams, initConfigParams, 
%   initVegParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


% Fetch the receiver parameters from the input structure
ant_pat_fullfile = [];

% Antenna Height (m)
hr_m = inputStruct.hr_m;

% Receive Antenna Gain (dB) 
G0r_dB = inputStruct.G0r_dB;

% Receiver antenna polarization
pol_Rx = inputStruct.pol_Rx;

% Get the string orientation_Rx and convert it to an integer index
orientation_Rx_id = findElementIdInCell( Constants.RX_ORIENTATIONS, inputStruct.orientation_Rx );


%% DETERMINE RECEIVER ANGLES
% If receiver has a fixed orientation (observation and azimuh angles)
if orientation_Rx_id == Constants.ID_RX_FIXED

    th0_Rx_deg = inputStruct.th0_Rx_deg;   % Receiver observation (theta) angle
    ph0_Rx_deg = inputStruct.ph0_Rx_deg;   % Receiver azimuth (phi) angle

% Else if receiver always faces the specular point. Then it will always 
% have equal theta and azimuth angles with the transmitter. Ignore here!
elseif orientation_Rx_id == Constants.ID_RX_SPECULAR_FACING

    th0_Rx_deg = [];
    ph0_Rx_deg = [];
    
end


%% DETERMINE ANTENNA PATTERN
% Get the string ant_pat_Rx and convert it to an integer index    
ant_pat_Rx_id = findElementIdInCell( Constants.RX_ANT_PATTERNS, inputStruct.ant_pat_Rx );

% If receiver antenna pattern is Generalized-Gaussian
if ant_pat_Rx_id == Constants.ID_RX_GG
    
    ant_pat_res_deg_Rx = inputStruct.ant_pat_res_deg_Rx;   % Antenna pattern resolution in degrees
    
    % Set GG pattern parameters
    ant_pat_struct_Rx = initRxGGParams( inputStruct, ant_pat_res_deg_Rx );
    
elseif ant_pat_Rx_id == Constants.ID_RX_USER_DEFINED
    
    ant_pat_fullfile = inputStruct.ant_pat_Rx_file;  
    
    % Set User-defined pattern parameters
    [ant_pat_struct_Rx, ant_pat_res_deg_Rx] = initRxUserDefinedParams( inputStruct );
  
% Else if Antenna pattern is Cosine to the power n
elseif ant_pat_Rx_id == Constants.ID_RX_COS_POW_N 
    
    % Should be implemented when this pattern added
    
end


% INITIALIZE RECEIVER PARAMETERS
RxParams.getInstance.initialize( hr_m, G0r_dB, pol_Rx, ant_pat_Rx_id, ...
                               ant_pat_struct_Rx, ant_pat_res_deg_Rx, ...
                               orientation_Rx_id, th0_Rx_deg, ph0_Rx_deg );

end