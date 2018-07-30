function [ eps_diel_z2, eps_diel_z3, eps_diel_zL, eps_diel_zS ] = generateDielProfiles( eps_diel_soil, z )


%% GET GLOBAL PARAMETER
% Ground Parameters
gnd_layer_depth_m = GndParams.getInstance.layer_depth_m;
% Ground-ML Parameters
zA_m = GndMLParams.getInstance.zA_m;


% Layer bottom
L_bottom = [0; gnd_layer_depth_m(1 : end - 1) + diff(gnd_layer_depth_m) / 2 ] ;

% Layer thickness
L_thickness = diff(L_bottom) ;

% Soil Parameters
% standard variation of layer boundaries: 10%
sL_depth = 0.1 * gnd_layer_depth_m ;


%% 2ND ORDER POLYFIT
zz = zA_m + gnd_layer_depth_m ;
p2 = polyfit(zz, eps_diel_soil, 2) ;
eps_diel_z2 = polyval(p2, z) ;
eps_diel_z2(z <= zA_m) = Constants.eps_diel_air ;


%% 3RD ORDER POLYFIT
zz = zA_m + gnd_layer_depth_m ;
p3 = polyfit(zz, eps_diel_soil, 3) ;
eps_diel_z3 = polyval(p3, z) ;
eps_diel_z3(z <= zA_m) = Constants.eps_diel_air ;
eps_diel_z3(real(eps_diel_z3) < 1) = 1 ;


%% SMOOTHER TRANSITION
[z, eps_diel_zL] ...
    = DielPrf0(eps_diel_soil, L_thickness, sL_depth, zA_m, z) ;


%% DISCRETE SLAB
eps_diel_zS = eps_diel_zL ; 
eps_diel_zS(z < zA_m) = Constants.eps_diel_air ; 

for ii = 1 : length(L_bottom)-1
    eps_diel_zS(z > (zA_m + L_bottom(ii)) & z < (zA_m + L_bottom(ii+1))) = eps_diel_soil(ii) ;
end

eps_diel_zS(z > (zA_m + L_bottom(end))) = eps_diel_soil(end) ;


end


function [z, eps_diel_z] ...
    = DielPrf0(eps_diel_soil, L_thickness, sL_depth, zA_m, z)


%% GET GLOBAL PARAMETERS
% Ground Parameters
num_gnd_layers = GndParams.getInstance.num_layers;

eps1 = Constants.eps_diel_air ;
eps_diel_z = eps1  ;
zz = zeros(1, num_gnd_layers) ;
for ii = 1 : num_gnd_layers
    
    sdel = sL_depth(ii) ;
    if ii == 1
        zz(ii) = zA_m ;
    else
        zz(ii) = zz(ii - 1) + L_thickness(ii - 1) ;
    end
    eps2 = eps_diel_soil(ii) ;
    
    eps_diel_z = eps_diel_z ...
        + (eps2 - eps1) ./ (1 + exp(-2.197 * (z - zz(ii)) / sdel)) ;
    
    eps1 = eps2 ;
    
end

% to exclude roughness
eps_diel_z(z <= zA_m) = Constants.eps_diel_air ;

end

