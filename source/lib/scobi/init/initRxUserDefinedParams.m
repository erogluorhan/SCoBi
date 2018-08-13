
function [ant_pat_struct_Rx, ant_pat_res_deg] = initRxUserDefinedParams( inputStruct )

% Filename with the full path
ant_pat_Rx_file = inputStruct.ant_pat_Rx_file;    

% Initialize Generalized-Gaussian Receiver Parameters
RxUserDefinedParams.getInstance.initialize( ant_pat_Rx_file );

% Generate the antenna pattern struct
[ant_pat_struct_Rx, ant_pat_res_deg] = RxUserDefinedParams.getInstance.calcRxAntPatMatrix();

end