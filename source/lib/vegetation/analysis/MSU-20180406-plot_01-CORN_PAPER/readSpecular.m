



%% read results
function [b_coh1veg, b_coh2veg, b0_coh1veg, b0_coh2veg, ...
    b_coh1bare, b_coh2bare, b0_coh1bare, b0_coh2bare, ...
    P_coh1veg, P_coh2veg, P0_coh1veg, P0_coh2veg, ...
    P_coh1bare, P_coh2bare, P0_coh1bare, P0_coh2bare] = readSpecular


%% GET GLOBAL DIRECTORIES
dir_out_specular_tuple = SimulationFolders.getInstance.out_specular_tuple;


%% GET GLOBAL PARAMETERS
pol_Tx = TxParams.getInstance.pol_Tx;
pol_Rx = RxParams.getInstance.pol_Rx;


%% READ SAVED OUTPUT

% 2 X 2
filename1 = strcat('Veg1', pol_Tx, pol_Rx) ;
filename2 = strcat('Veg2', pol_Tx, pol_Rx) ;
filename01 = strcat('Veg01', pol_Tx, pol_Rx) ;
filename02 = strcat('Veg02', pol_Tx, pol_Rx) ;
b_coh1veg = readComplexVar(dir_out_specular_tuple, filename1) ;
b_coh2veg = readComplexVar(dir_out_specular_tuple, filename2) ;
b0_coh1veg = readComplexVar(dir_out_specular_tuple, filename01) ;
b0_coh2veg = readComplexVar(dir_out_specular_tuple, filename02) ;

filename1 = strcat('Bare1', pol_Tx, pol_Rx) ;
filename2 = strcat('Bare2', pol_Tx, pol_Rx) ;
filename01 = strcat('Bare01', pol_Tx, pol_Rx) ;
filename02 = strcat('Bare02', pol_Tx, pol_Rx) ;
b_coh1bare = readComplexVar(dir_out_specular_tuple, filename1) ;
b_coh2bare = readComplexVar(dir_out_specular_tuple, filename2) ;
b0_coh1bare = readComplexVar(dir_out_specular_tuple, filename01) ;
b0_coh2bare = readComplexVar(dir_out_specular_tuple, filename02) ;

% 4 X 4
filename1 = strcat('P_Veg1', pol_Tx, pol_Rx) ;
filename2 = strcat('P_Veg2', pol_Tx, pol_Rx) ;
filename01 = strcat('P_Veg01', pol_Tx, pol_Rx) ;
filename02 = strcat('P_Veg02', pol_Tx, pol_Rx) ;
P_coh1veg = readVar(dir_out_specular_tuple, filename1) ;
P_coh2veg = readVar(dir_out_specular_tuple, filename2) ;
P0_coh1veg = readVar(dir_out_specular_tuple, filename01) ;
P0_coh2veg = readVar(dir_out_specular_tuple, filename02) ;

filename1 = strcat('P_Bare1', pol_Tx, pol_Rx) ;
filename2 = strcat('P_Bare2', pol_Tx, pol_Rx) ;
filename01 = strcat('P_Bare01', pol_Tx, pol_Rx) ;
filename02 = strcat('P_Bare02', pol_Tx, pol_Rx) ;
P_coh1bare = readVar(dir_out_specular_tuple, filename1) ;
P_coh2bare = readVar(dir_out_specular_tuple, filename2) ;
P0_coh1bare = readVar(dir_out_specular_tuple, filename01) ;
P0_coh2bare = readVar(dir_out_specular_tuple, filename02) ;

end

