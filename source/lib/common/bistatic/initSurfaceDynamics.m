% Mehmet Kurum
% Feb 22, 2017

function initSurfaceDynamics


%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Simulation Settings
simulator_id = SimSettings.getInstance.simulator_id;
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
f_Hz = f_MHz * Constants.MHz2Hz ;
% Dynamic Parameters
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
VSM_cm3cm3 = VSM_list_cm3cm3( sim_counter, : )';
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;
RMSH_cm = RMSH_list_cm( sim_counter );
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

% TO-DO: Justify the conjugate operation for SCoBi-Veg
if simulator_id == Constants.id_veg_agr ...
        || simulator_id == Constants.id_veg_for

    eps_g = conj(eps_g) ; % eps_g = eps_gp + i * eps_gpp
    
end

eps_g = round(eps_g * 10) / 10 ;


% Initialize Surface Dynamic Parameters
SurfaceDynParams.getInstance.initialize(h, eps_g);

end