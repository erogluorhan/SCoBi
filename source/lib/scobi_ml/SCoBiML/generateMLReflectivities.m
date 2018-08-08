function generateMLReflectivities


%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
f_Hz = f_MHz * Constants.MHz2Hz;
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
DoYs = DynParams.getInstance.DoYs;
DoY = DoYs( sim_counter );
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
th0_Tx_deg = th0_Tx_list_deg( sim_counter );
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;   
VSM_cm3cm3 = VSM_list_cm3cm3(sim_counter,:)';
% Dielectric Parameters
eps_diel_z2nd = DielParams.getInstance.eps_diel_z2nd;
eps_diel_z3rd = DielParams.getInstance.eps_diel_z3rd;
eps_diel_zL = DielParams.getInstance.eps_diel_zL;
eps_diel_zS = DielParams.getInstance.eps_diel_zS;
% Ground Parameters
sand_ratio = GndParams.getInstance.sand_ratio;  % Texture at various depths 
clay_ratio = GndParams.getInstance.clay_ratio;  % Texture at various depths 
rhob_gcm3 = GndParams.getInstance.rhob_gcm3;    % Soil bulk density at various depths 
% Ground MultiLayer Parameters
layer_bottom_m = GndMLParams.getInstance.layer_bottom_m;
zA_m = GndMLParams.getInstance.zA_m;    % Air layer
z_m = GndMLParams.getInstance.z_m;    % Layer profile


%% CALCULATE REFLECTION COEFFICIENTS
% Wavelength
lambda_m = Constants.c / f_Hz ;

% Soil Dielectric Constant
eps_diel_soil = round(10 * dielg(VSM_cm3cm3, f_Hz, sand_ratio, clay_ratio, rhob_gcm3)) / 10 ;

zzb = (zA_m + layer_bottom_m) ;
Lzb = diff(zzb)' ; % / lambda_m ; % complex optical length in units of lambda_m
nA = sqrte(Constants.eps_diel_air) ;
nS = sqrte(eps_diel_soil) ;
nS = nS.';

na = [nA; nA; nA] ;
ns = [nS; nS; nS] ;

% input to multidiel
n = [na, ns] ;

% Reflection Coefficient for 
rh = multidiel(n, Lzb, 1, th0_Tx_deg, 'te') ;
rv = multidiel(n, Lzb, 1, th0_Tx_deg, 'th') ;
r_s = [rv 0; 0 rh] ;

[r0_coh1b, r0_coh2b] = SpecularReflection(r_s) ;


% Reflection Coefficient for Logistic Profile
[rh_L, rv_L] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z_m, eps_diel_zL) ;
r_s_L = [rv_L 0; 0 rh_L] ;

[r0_coh1b_L, r0_coh2b_L] = SpecularReflection(r_s_L) ;


% Reflection Coefficient for 2nd Order Profile
[rh_2nd, rv_2nd] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z_m, eps_diel_z2nd ) ;
r_s_2nd = [rv_2nd 0; 0 rh_2nd] ;

[r0_coh1b_2nd, r0_coh2b_2nd] = SpecularReflection(r_s_2nd) ;


% Reflection Coefficient for 3rd Order Profile
[rh_3rd, rv_3rd] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z_m, eps_diel_z3rd ) ;
r_s_3rd = [rv_3rd 0; 0 rh_3rd] ;

[r0_coh1b_3rd, r0_coh2b_3rd] = SpecularReflection(r_s_3rd) ;


if pol_Tx == 'X' && pol_Rx == 'X'
    % Reflectivity
    Rp1 = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
    Rp2 = abs(r0_coh2b(2)) .^ 2 ; % at ?/?0 = 1
    Rp1_L = abs(r0_coh1b_L(1)) .^ 2 ; % at ?/?0 = 1
    Rp2_L = abs(r0_coh2b_L(2)) .^ 2 ; % at ?/?0 = 1
    Rp1_2nd = abs(r0_coh1b_2nd(1)) .^ 2 ; % at ?/?0 = 1
    Rp2_2nd = abs(r0_coh2b_2nd(2)) .^ 2 ; % at ?/?0 = 1
    Rp1_3rd = abs(r0_coh1b_3rd(1)) .^ 2 ; % at ?/?0 = 1
    Rp2_3rd = abs(r0_coh2b_3rd(2)) .^ 2 ; % at ?/?0 = 1
else
    % Reflectivity
    Rp1 = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
    Rp2 = abs(r0_coh1b(2)) .^ 2 ; % at ?/?0 = 1
    Rp1_L = abs(r0_coh1b_L(1)) .^ 2 ; % at ?/?0 = 1
    Rp2_L = abs(r0_coh1b_L(2)) .^ 2 ; % at ?/?0 = 1
    Rp1_2nd = abs(r0_coh1b_2nd(1)) .^ 2 ; % at ?/?0 = 1
    Rp2_2nd = abs(r0_coh1b_2nd(2)) .^ 2 ; % at ?/?0 = 1
    Rp1_3rd = abs(r0_coh1b_3rd(1)) .^ 2 ; % at ?/?0 = 1
    Rp2_3rd = abs(r0_coh1b_3rd(2)) .^ 2 ; % at ?/?0 = 1
end

plotReflectivityForProfiles( DoY, Rp1, Rp2, Rp1_L, Rp2_L, Rp1_2nd, Rp2_2nd, Rp1_3rd, Rp2_3rd );

plotReflectivityVsTh( Rp1, Rp2, Rp1_L, Rp2_L, Rp1_2nd, Rp2_2nd, Rp1_3rd, Rp2_3rd );

[M1, M2] = plotAddSMpoint( DoY, VSM_cm3cm3 );

% plotMovie( M1 );

end


function [rh, rv] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z_m, eps_diel_z)

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

%% Reflection Coeffficeint
rh = multidiel(n, Lz, 1, th0_Tx_deg, 'te') ;
rv = multidiel(n, Lz, 1, th0_Tx_deg, 'th') ;

end

