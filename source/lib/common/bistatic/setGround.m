% Mehmet Kurum
% Feb 22, 2017

function setGround

%% GET GLOBAL DIRECTOIRES
dir_gnd = SimulationFolders.getInstance.gnd;


%% GET GLOBAL PARAMETERS
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
f_Hz = f_MHz * Constants.MHz2Hz ;
% Dynamic Parameters
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
VSM_cm3cm3 = VSM_list_cm3cm3( ParamsManager.index_VSM, : )';
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;
RMSH_cm = RMSH_list_cm( ParamsManager.index_RMSH );
% Ground Parameters
% TO-DO: Implement for different ground layers
sand_ratio = GndParams.getInstance.sand_ratio;
clay_ratio = GndParams.getInstance.clay_ratio;
rhob_gcm3 = GndParams.getInstance.rhob_gcm3;


%% CALCULATIONS
% Surface Roughness
lambda_cm = Constants.c / f_Hz * Constants.m2cm ;        % in cm
ko = 2 * pi / lambda_cm ;
h = (2 * RMSH_cm * ko) ^ 2 ;        % effective roughness parameter

% Soil Dielectric Constant
eps_g = dielg( VSM_cm3cm3, f_Hz, sand_ratio, clay_ratio, rhob_gcm3) ; % eps_g = eps_gp - j * eps_gpp
eps_g = conj(eps_g) ; % eps_g = eps_gp + i * eps_gpp
eps_g = round(eps_g * 10) / 10 ;

% Ground Parameters
% grnd_par = [h; real(eps_g); imag(eps_g)] ;
grnd_par = [h; eps_g] ;


%% SAVE
filename = 'G' ;
writeComplexVar( dir_gnd, filename, grnd_par) ;

end