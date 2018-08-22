
function ant_pat_fullfile = initRxParams( inputStruct )

ant_pat_fullfile = [];

% Antenna Height (m)
hr_m = inputStruct.hr_m;

% Receive Antenna Gain (dB) 
G0r_dB = inputStruct.G0r_dB;

% Receiver antenna polarization
pol_Rx = inputStruct.pol_Rx;

% Get the string orientation_Rx and convert it to an integer index
orientation_Rx_id = findElementIdInCell( Constants.Rx_orientations, inputStruct.orientation_Rx );


%% DETERMINE RECEIVER ANGLES
% If receiver has a fixed orientation (observation and azimuh angles)
if orientation_Rx_id == Constants.id_Rx_fixed

    th0_Rx_deg = inputStruct.th0_Rx_deg;   % Receiver observation (theta) angle
    ph0_Rx_deg = inputStruct.ph0_Rx_deg;   % Receiver azimuth (phi) angle

% Else if receiver always faces the specular point. Then it will always 
% have equal theta and azimuth angles with the transmitter. Ignore here!
elseif orientation_Rx_id == Constants.id_Rx_specular_facing

    th0_Rx_deg = [];
    ph0_Rx_deg = [];
    
end


%% DETERMINE ANTENNA PATTERN
% Get the string ant_pat_Rx and convert it to an integer index    
ant_pat_Rx_id = findElementIdInCell( Constants.Rx_ant_pats, inputStruct.ant_pat_Rx );

% If receiver antenna pattern is Generalized-Gaussian
if ant_pat_Rx_id == Constants.id_Rx_GG
    
    ant_pat_res_deg_Rx = inputStruct.ant_pat_res_deg_Rx;   % Antenna pattern resolution in degrees
    
    % Set GG pattern parameters
    ant_pat_struct_Rx = initRxGGParams( inputStruct, ant_pat_res_deg_Rx );
    
elseif ant_pat_Rx_id == Constants.id_Rx_user_defined
    
    ant_pat_fullfile = inputStruct.ant_pat_Rx_file;  
    
    % Set User-defined pattern parameters
    [ant_pat_struct_Rx, ant_pat_res_deg_Rx] = initRxUserDefinedParams( inputStruct );
  
% Else if Antenna pattern is Cosine to the power n
elseif ant_pat_Rx_id == Constants.id_Rx_cos_pow_n 
    
    % Should be implemented when this pattern added
    
end


%% INITIALIZE
% Initialize Receiver Parameters
RxParams.getInstance.initialize( hr_m, G0r_dB, pol_Rx, ant_pat_Rx_id, ...
                               ant_pat_struct_Rx, ant_pat_res_deg_Rx, ...
                               orientation_Rx_id, th0_Rx_deg, ph0_Rx_deg );

end