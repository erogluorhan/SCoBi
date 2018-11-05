
function [S0x_m, x1_m, ax1_m, by1_m] = calcFresnelZones(ht_m, hr_m, Nfz)
% function calcFresnelZones 
%
%   Calculates and outputs the Fresnel zones for the specular reflection 
%   point and to use them in the calculation of the other points.  
%
%   [S0x_m, x1_m, ax1_m, by1_m] = calcFresnelZones(ht_m, hr_m, Nfz)
%
%   INPUTS:
%   ht_m: Transmitter altitude (in meters)
%   hr_m: Receiver altitude (in meters)
%   Nfz: Number of Fresnel zones
%
%   See also bistaticGeometry.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Transmitter Parameters
f_Hz = TxParams.getInstance.f_MHz * Constants.MHZ_TO_HZ ;
% Configuration Parameters
th0_Tx_list_deg = ConfigParams.getInstance.th0_Tx_list_deg;
th0_Tx_deg = th0_Tx_list_deg( sim_counter );
el0_Tx_deg = 90 - th0_Tx_deg;


%% CALCULATIONS
% Wavelength
lambda_m = Constants.LIGHTSPEED / f_Hz ;

% Angle of Incidence of incoming signal
th0_Tx_rad = deg2rad(th0_Tx_deg) ;

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

tanth = hd_m / Dist_m ;
sinth = hd_m / rd2_m ;
costh = Dist_m / rd2_m ;
secth = 1 / costh ;

x1_m = zeros(Nfz, 1) ;
ax1_m = zeros(Nfz, 1) ;
by1_m = zeros(Nfz, 1) ;

% Iterate for the decided number of Fresnel zones
% This iteration may be meaningful, when incoherent contribution is 
% considered. Nfz is always 1 for the current version. 
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