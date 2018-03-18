%% Mehmet Kurum
% March 15, 2017

function [S0x, x1, ax1, by1] = calcFresnelZones(ht, hr)


%% Get Global Parameters
f0hz = SatParams.getInstance.fMHz * 1e6 ;
EL0_deg = SatParams.getInstance.EL0_deg( ParamsManager.index_Th );

lambda = Constants.c / f0hz ;     % Wavelength

% Angle of Incidence of incoming signal
% SatParams.getInstance.EL0_deg = 36.3287 ;        % Elevation angle
th0t_deg = 90 - EL0_deg ;
th0t = degtorad(th0t_deg) ;

% Transmitter Height with respect to local planar ground plane
% ht = rd * cos(th0t) ;
% Transmitter/Reciever ground range
% Dist = rd * sin(th0t) ;
Dist = ht * tan(th0t) ;

%% Calculations


% Distance of specular point away from the receiver ground projection.
S0x = Dist / (1 + ht / hr) ;

% Parameters for Fresnel zone calculation
hd = ht - hr ;
hs = ht + hr ;
rd2 = sqrt(hd ^ 2 + Dist ^ 2) ;
rs = sqrt(hs ^ 2 + Dist ^ 2) ;
del0 = rs - rd2 ; % the shortest distance after the direct path

tanth = hd / Dist ; %#ok<NASGU>
sinth = hd / rd2 ;
costh = Dist / rd2 ;
secth = 1 / costh ;

x1 = zeros(SimParams.getInstance.Nfz, 1) ;
ax1 = zeros(SimParams.getInstance.Nfz, 1) ;
by1 = zeros(SimParams.getInstance.Nfz, 1) ;

for nn = 1 : SimParams.getInstance.Nfz
    
    % the path for the nth Fresnel Zone
    del = del0 + nn * lambda / 2 ;
    
    a = (Dist * secth + del) / 2 ;
    b = sqrt(del ^ 2 + 2 * Dist * del * secth) / 2 ;
    cc = (ht + hr) / 2 ;    % OE: Named 'cc' to prevent confusion with the other c - speed of light
    
    nump = cc * (a ^ 2 - b ^ 2) * sinth * costh ;
    denp = (b ^ 2 * costh ^ 2 + a ^ 2 * sinth ^ 2) ;
    p = -nump / denp ;
    x1(nn) = Dist / 2 + p ; % distance to the center of the ellipse
    
    % minor axis
    by1(nn) = b * sqrt(1 - cc ^ 2 / denp) ;
    % major axis
    ax1(nn) = by1(nn) * a / sqrt(denp) ;
    
end


end