
function ant_pat_struct_Rx = initRxGGParams( inputStruct, ant_pat_res_deg )
    
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