%% Mehmet Kurum
% March 15, 2017

function [S0x_m, x1_m, ax1_m, by1_m] = calcFresnelZones(ht_m, hr_m, Nfz)


%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Transmitter Parameters
f_Hz = TxParams.getInstance.f_MHz * Constants.MHz2Hz ;
% Configuration Parameters
th0_Tx_list_deg = ConfigParams.getInstance.th0_Tx_list_deg;
th0_Tx_deg = th0_Tx_list_deg( sim_counter );
el0_Tx_deg = 90 - th0_Tx_deg;


%% CALCULATIONS
lambda_m = Constants.c / f_Hz ;     % Wavelength

% Angle of Incidence of incoming signal
% TxParams.getInstance.el0_Tx_deg = 36.3287 ;        % Elevation angle
th0_Tx_rad = degtorad(th0_Tx_deg) ;

% Transmitter/Reciever ground range
Dist_m = ht_m * tan(th0_Tx_rad) ;

% Distance of specular point away from the receiver ground projection.
S0x_m = Dist_m / (1 + ht_m / hr_m) ;

% Parameters for Fresnel zone calculation
hd_m = ht_m - hr_m ;
hs_m = ht_m + hr_m ;
rd2_m = sqrt(hd_m ^ 2 + Dist_m ^ 2) ;
rs_m = sqrt(hs_m ^ 2 + Dist_m ^ 2) ;
del0_m = rs_m - rd2_m ; % the shortest distance after the direct path

tanth = hd_m / Dist_m ; %#ok<NASGU>
sinth = hd_m / rd2_m ;
costh = Dist_m / rd2_m ;
secth = 1 / costh ;

x1_m = zeros(Nfz, 1) ;
ax1_m = zeros(Nfz, 1) ;
by1_m = zeros(Nfz, 1) ;

for nn = 1 : Nfz
    
    % the path for the nth Fresnel Zone
    del_m = del0_m + nn * lambda_m / 2 ;
    
    a_m = (Dist_m * secth + del_m) / 2 ;
    b_m = sqrt(del_m ^ 2 + 2 * Dist_m * del_m * secth) / 2 ;
    c_m = (ht_m + hr_m) / 2 ;
    
    nump_m3 = c_m * (a_m ^ 2 - b_m ^ 2) * sinth * costh ;
    denp_m2 = (b_m ^ 2 * costh ^ 2 + a_m ^ 2 * sinth ^ 2) ;
    p_m = -nump_m3 / denp_m2 ;
    x1_m(nn) = Dist_m / 2 + p_m ; % distance to the center of the ellipse
    
    % semi-minor axis
    by1_m(nn) = b_m * sqrt(1 - c_m ^ 2 / denp_m2) ;
    % semi-major axis
    ax1_m(nn) = by1_m(nn) * a_m / sqrt(denp_m2) ;
    
end


end