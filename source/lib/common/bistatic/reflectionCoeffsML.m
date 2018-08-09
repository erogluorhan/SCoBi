function [r_g] = reflectionCoeffsML()


%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
f_Hz = f_MHz * Constants.MHz2Hz;
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
th0_Tx_deg = th0_Tx_list_deg( sim_counter );
% Dielectric Parameters
eps_diel_z2nd = DielParams.getInstance.eps_diel_z2nd;
eps_diel_z3rd = DielParams.getInstance.eps_diel_z3rd;
eps_diel_zL = DielParams.getInstance.eps_diel_zL;
eps_diel_zS = DielParams.getInstance.eps_diel_zS;
% Surface Dynamic Params
eps_g = SurfaceDynParams.getInstance.eps_g;
% Ground MultiLayer Parameters
layer_bottom_m = GndMLParams.getInstance.layer_bottom_m;
zA_m = GndMLParams.getInstance.zA_m;    % Air layer
z_m = GndMLParams.getInstance.z_m;    % Layer profile


%% CALCULATE REFLECTION COEFFICIENTS
% Wavelength
lambda_m = Constants.c / f_Hz ;


% Reflection Coefficient for Discrete Slab
zzb = (zA_m + layer_bottom_m) ;
Lzb = diff(zzb)' ; % / lambda_m ; % complex optical length in units of lambda_m
nA = sqrte(Constants.eps_diel_air) ;
nS = sqrte(eps_g) ;
nS = nS.';

na = [nA; nA; nA] ;
ns = [nS; nS; nS] ;

% input to multidiel
n = [na, ns] ;

rh = multidiel(n, Lzb, 1, th0_Tx_deg, 'te') ;
rv = multidiel(n, Lzb, 1, th0_Tx_deg, 'th') ;
r_g{Constants.id_diel_slab, 1} = [rv 0; 0 rh];


% Reflection Coefficient for Logistic Profile
[rh_L, rv_L] = calcSpecularReflectionCoeffML(lambda_m, th0_Tx_deg, z_m, eps_diel_zL) ;
r_g{Constants.id_diel_logistic, 1} = [rv_L 0; 0 rh_L] ;


% Reflection Coefficient for 2nd Order Profile
[rh_2nd, rv_2nd] = calcSpecularReflectionCoeffML(lambda_m, th0_Tx_deg, z_m, eps_diel_z2nd ) ;
r_g{Constants.id_diel_2nd_order, 1} = [rv_2nd 0; 0 rh_2nd] ;


% Reflection Coefficient for 3rd Order Profile
[rh_3rd, rv_3rd] = calcSpecularReflectionCoeffML(lambda_m, th0_Tx_deg, z_m, eps_diel_z3rd ) ;
r_g{Constants.id_diel_3rd_order, 1} = [rv_3rd 0; 0 rh_3rd] ;

end


function [rh, rv] = calcSpecularReflectionCoeffML(lambda_m, th0_Tx_deg, z_m, eps_diel_z)

Lz = diff(z_m)' / lambda_m ; % complex optical length in units of lambda_m

nAz = sqrte(Constants.eps_diel_air) ;
nmz = sqrte(eps_diel_z(2 : end, :)) ;
nSz = sqrte(eps_diel_z(end, :)) ;

% Air - % isotropic
na = [nAz; nAz; nAz] ;

% Dielectric Profile : isotropic
nm = [nmz(:, 1).'; nmz(:, 1).'; nmz(:, 1).'] ;

% Soil - isotropic
nb = [nSz(:, 1); nSz(:, 1); nSz(:, 1)] ;

%% input to multidiel
n = [na, nm, nb] ;

%% Reflection Coeffficient
rh = multidiel(n, Lz, 1, th0_Tx_deg, 'te') ;
rv = multidiel(n, Lz, 1, th0_Tx_deg, 'th') ;

end

