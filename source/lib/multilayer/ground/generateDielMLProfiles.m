
function generateDielMLProfiles
% function generateDielMLProfiles 
%
%   Generates the selected dielectric profiles (2nd order, 3rd order, 
%   Logistic regression, Discrete slab) if the ground structure is 
%   Multi-layered and updates the dynamic Multi-layer dielectric parameters 
%   (DielMLDynParams) with those parameter values in each simulation 
%   iteration.
%
%   - Uses the information from GndMLParams and GndDynParams
%   - Updates DielMLDynParams with the calculated profiles 
%
%   See also mainSCoBi, GndMLParams, GndDynParams, DielMLDynParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.1
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%  UPDATE HISTORY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Version 1.0.1  
%
%   November 14, 2018
%
%   Refitted the complex conjugate to make the physics-oriented multidiel 
%   function compatible with SCoBi. In engineering, j =sqrt(-1).
%   In physics, i = (-1)sqrt(-1).
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% GET GLOBAL PARAMETER
% Ground Parameters
gnd_layer_depth_m = GndMLParams.getInstance.layer_depth_m;
% Ground Dynamic Params
eps_g = GndDynParams.getInstance.eps_g;
eps_g = conj(eps_g); % i --> j
% Ground-ML Parameters
layer_bottom_m = GndMLParams.getInstance.layer_bottom_m;
layer_thickness_m = GndMLParams.getInstance.layer_thickness_m;
zA_m = GndMLParams.getInstance.zA_m;
z_m = GndMLParams.getInstance.z_m;    % Layer profile
calc_diel_profile_fit_functions = GndMLParams.getInstance.calc_diel_profile_fit_functions;


%% CALCULATIONS
% Soil Parameters
% standard variation of layer boundaries: 10%
sL_depth = 0.1 * gnd_layer_depth_m ;


%% 2ND ORDER POLYFIT
eps_diel_z2nd = [];

if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_2ND_ORDER, 1)
    
    zz = zA_m + gnd_layer_depth_m ;
    p2 = polyfit(zz, eps_g, 2) ;
    eps_diel_z2nd = polyval(p2, z_m) ;
    eps_diel_z2nd(z_m <= zA_m) = Constants.EPS_DIEL_AIR ;
    
end


%% 3RD ORDER POLYFIT
eps_diel_z3rd = [];

if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_3RD_ORDER, 1)

    zz = zA_m + gnd_layer_depth_m ;
    p3 = polyfit(zz, eps_g, 3) ;
    eps_diel_z3rd = polyval(p3, z_m) ;
    eps_diel_z3rd(z_m <= zA_m) = Constants.EPS_DIEL_AIR ;
    eps_diel_z3rd(real(eps_diel_z3rd) < 1) = 1 ;

end


%% LOGISTIC REGRESSION
eps_diel_zL = [];

if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_LOGISTIC, 1)
    
    [z_m, eps_diel_zL] ...
        = dielProfileLogistic(eps_g, layer_thickness_m, sL_depth, zA_m, z_m) ;
    
end


%% DISCRETE SLAB
eps_diel_zS = [];

if calc_diel_profile_fit_functions(Constants.ID_DIEL_PROFILE_SLAB, 1)

    eps_diel_zS = eps_diel_zL ; 
    eps_diel_zS(z_m < zA_m) = Constants.EPS_DIEL_AIR ; 

    for ii = 1 : length(layer_bottom_m)-1
        
        eps_diel_zS(z_m > (zA_m + layer_bottom_m(ii)) & z_m < (zA_m + layer_bottom_m(ii+1))) = eps_g(ii) ;
        
    end

    eps_diel_zS(z_m > (zA_m + layer_bottom_m(end))) = eps_g(end) ;
    
end


% Initialize Dielectric Parameters
DielMLDynParams.getInstance.initialize( eps_diel_z2nd, eps_diel_z3rd, eps_diel_zL, eps_diel_zS );


end


function [z_m, eps_diel_z] ...
    = dielProfileLogistic(eps_g, layer_thickness_m, sL_depth, zA_m, z_m)


%% GET GLOBAL PARAMETERS
% Ground Parameters
num_gnd_layers = GndMLParams.getInstance.num_layers;

eps1 = Constants.EPS_DIEL_AIR ;
eps_diel_z = eps1  ;
zz = zeros(1, num_gnd_layers) ;
for ii = 1 : num_gnd_layers
    
    sdel = sL_depth(ii) ;
    if ii == 1
        zz(ii) = zA_m ;
    else
        zz(ii) = zz(ii - 1) + layer_thickness_m(ii - 1) ;
    end
    eps2 = eps_g(ii) ;
    
    eps_diel_z = eps_diel_z ...
        + (eps2 - eps1) ./ (1 + exp(-2.197 * (z_m - zz(ii)) / sdel)) ;
    
    eps1 = eps2 ;
    
end

% to exclude roughness
eps_diel_z(z_m <= zA_m) = Constants.EPS_DIEL_AIR ;

end

