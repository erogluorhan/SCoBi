function generateDielMLProfiles


%% GET GLOBAL PARAMETER
% Ground Parameters
gnd_layer_depth_m = GndMLParams.getInstance.layer_depth_m;
% Ground Dynamic Params
eps_g = GndDynParams.getInstance.eps_g;
% Ground-ML Parameters
layer_bottom_m = GndMLParams.getInstance.layer_bottom_m;
layer_thickness_m = GndMLParams.getInstance.layer_thickness_m;
zA_m = GndMLParams.getInstance.zA_m;
z_m = GndMLParams.getInstance.z_m;    % Layer profile
calc_diel_profile_fit_functions = GndMLParams.getInstance.calc_diel_profile_fit_functions;


%% PLOT VSM data first
[fig1, fig2] = plotSMdata();


%% CALCULATIONS
% Soil Parameters
% standard variation of layer boundaries: 10%
sL_depth = 0.1 * gnd_layer_depth_m ;


%% 2ND ORDER POLYFIT
eps_diel_z2nd = [];

if calc_diel_profile_fit_functions(Constants.id_diel_2nd_order, 1)
    
    zz = zA_m + gnd_layer_depth_m ;
    p2 = polyfit(zz, eps_g, 2) ;
    eps_diel_z2nd = polyval(p2, z_m) ;
    eps_diel_z2nd(z_m <= zA_m) = Constants.eps_diel_air ;
    
end


%% 3RD ORDER POLYFIT
eps_diel_z3rd = [];

if calc_diel_profile_fit_functions(Constants.id_diel_3rd_order, 1)

    zz = zA_m + gnd_layer_depth_m ;
    p3 = polyfit(zz, eps_g, 3) ;
    eps_diel_z3rd = polyval(p3, z_m) ;
    eps_diel_z3rd(z_m <= zA_m) = Constants.eps_diel_air ;
    eps_diel_z3rd(real(eps_diel_z3rd) < 1) = 1 ;

end

%% LOGISTIC REGRESSION
eps_diel_zL = [];

if calc_diel_profile_fit_functions(Constants.id_diel_logistic, 1)
    
    [z_m, eps_diel_zL] ...
        = DielPrf0(eps_g, layer_thickness_m, sL_depth, zA_m, z_m) ;
    
end


%% DISCRETE SLAB
eps_diel_zS = [];

if calc_diel_profile_fit_functions(Constants.id_diel_slab, 1)

    eps_diel_zS = eps_diel_zL ; 
    eps_diel_zS(z_m < zA_m) = Constants.eps_diel_air ; 

    for ii = 1 : length(layer_bottom_m)-1
        
        eps_diel_zS(z_m > (zA_m + layer_bottom_m(ii)) & z_m < (zA_m + layer_bottom_m(ii+1))) = eps_g(ii) ;
        
    end

    eps_diel_zS(z_m > (zA_m + layer_bottom_m(end))) = eps_g(end) ;
    
end


%% PLOT DIELECTRIC PROFILES
plotDielProfiles( fig1, eps_diel_z2nd, eps_diel_z3rd, eps_diel_zL, eps_diel_zS, z_m );


% Initialize Dielectric Parameters
DielMLDynParams.getInstance.initialize( fig1, fig2, eps_diel_z2nd, eps_diel_z3rd, eps_diel_zL, eps_diel_zS );


end


function [z_m, eps_diel_z] ...
    = DielPrf0(eps_g, layer_thickness_m, sL_depth, zA_m, z_m)


%% GET GLOBAL PARAMETERS
% Ground Parameters
num_gnd_layers = GndMLParams.getInstance.num_layers;

eps1 = Constants.eps_diel_air ;
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
eps_diel_z(z_m <= zA_m) = Constants.eps_diel_air ;

end

