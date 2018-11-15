
function [r_g] = reflectionCoeffsML
% function reflectionCoeffsML 
%
%   Calculates the equivalent reflection coefficients of the rough, 
%   multi-layered ground for the chosen ones of four dielectric profiles.
%
%   [r_g] = reflectionCoeffsML
%
%   See also specularTerm, reflectionCoeffsSingle

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.1
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%  UPDATE HISTORY  %%%%%%%%%%%%%%%%%%%%%%%%%  %
%   Version 1.0.1
%
%   November 14, 2018
%
%   Refitted the complex conjugate to make the physics-oriented multidiel 
%   function compatible with SCoBi. In engineering, j =sqrt(-1).
%   In physics, i = (-1)sqrt(-1).
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %


%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
f_Hz = f_MHz * Constants.MHZ_TO_HZ;
% Configuration Parameters
th0_Tx_list_deg = ConfigParams.getInstance.th0_Tx_list_deg;
th0_Tx_deg = th0_Tx_list_deg( sim_counter );
% Dielectric Parameters
eps_diel_z2nd = DielMLDynParams.getInstance.eps_diel_z2nd;
eps_diel_z3rd = DielMLDynParams.getInstance.eps_diel_z3rd;
eps_diel_zL = DielMLDynParams.getInstance.eps_diel_zL;
eps_diel_zS = DielMLDynParams.getInstance.eps_diel_zS;
% Ground Dynamic Params
eps_g = GndDynParams.getInstance.eps_g;
eps_g = conj(eps_g); % i --> j
% Ground MultiLayer Parameters
layer_bottom_m = GndMLParams.getInstance.layer_bottom_m;
zA_m = GndMLParams.getInstance.zA_m;    % Air layer
z_m = GndMLParams.getInstance.z_m;    % Layer profile
calc_diel_profile_fit_functions = GndMLParams.getInstance.calc_diel_profile_fit_functions;


%% CALCULATE REFLECTION COEFFICIENTS
% Wavelength
lambda_m = Constants.LIGHTSPEED / f_Hz ;

[num_diel_profiles, ~] = size(calc_diel_profile_fit_functions);

r_g = cell(num_diel_profiles, 1);

% Reflection Coefficient for Discrete Slab
if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_SLAB, 1)
    
    zzb = (zA_m + layer_bottom_m) ;
    Lzb = diff(zzb)' ; % / lambda_m ; % complex optical length in units of lambda_m
    nA = sqrte(Constants.EPS_DIEL_AIR) ;
    nS = sqrte(eps_g) ;
    nS = nS.';

    na = [nA; nA; nA] ;
    ns = [nS; nS; nS] ;

    % input to multidiel
    n = [na, ns] ;

    rh = multidiel(n, Lzb, 1, th0_Tx_deg, 'te') ;
    rv = multidiel(n, Lzb, 1, th0_Tx_deg, 'th') ;
    r_g{Constants.ID_DIEL_PROFILE_SLAB, 1} = [rv 0; 0 rh];

end

% Reflection Coefficient for Logistic Profile
if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_LOGISTIC, 1)

    [rh_L, rv_L] = calcSpecularReflectionCoeffML(lambda_m, th0_Tx_deg, z_m, eps_diel_zL) ;
    r_g{Constants.ID_DIEL_PROFILE_LOGISTIC, 1} = [rv_L 0; 0 rh_L] ;
    
end


% Reflection Coefficient for 2nd Order Profile
if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_2ND_ORDER, 1)
    
    [rh_2nd, rv_2nd] = calcSpecularReflectionCoeffML(lambda_m, th0_Tx_deg, z_m, eps_diel_z2nd ) ;
    r_g{Constants.ID_DIEL_PROFILE_2ND_ORDER, 1} = [rv_2nd 0; 0 rh_2nd] ;
    
end


% Reflection Coefficient for 3rd Order Profile
if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_3RD_ORDER, 1)
    
    [rh_3rd, rv_3rd] = calcSpecularReflectionCoeffML(lambda_m, th0_Tx_deg, z_m, eps_diel_z3rd ) ;
    r_g{Constants.ID_DIEL_PROFILE_3RD_ORDER, 1} = [rv_3rd 0; 0 rh_3rd] ;
    
end

end


function [rh, rv] = calcSpecularReflectionCoeffML(lambda_m, th0_Tx_deg, z_m, eps_diel_z)

Lz = diff(z_m)' / lambda_m ; % complex optical length in units of lambda_m

nAz = sqrte(Constants.EPS_DIEL_AIR) ;
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

