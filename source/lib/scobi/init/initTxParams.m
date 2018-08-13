
function initTxParams( inputStruct )

f_MHz = inputStruct.f_MHz;      % Operating frequncy (MHz)

r_Tx_m = inputStruct.r_Tx_km * Constants.km2m;    % Transmitter range from Earth"s center (km - > m) 

EIRP_dB = inputStruct.EIRP_dB;    % Equivalent Isotropic Radiated Power

pol_Tx = inputStruct.pol_Tx;  % Transmitter polarization                        

% Initialize Transmitter Parameters
TxParams.getInstance.initialize( f_MHz, r_Tx_m, EIRP_dB, pol_Tx );

end