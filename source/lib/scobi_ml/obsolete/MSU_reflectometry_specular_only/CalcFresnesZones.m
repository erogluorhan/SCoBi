%% Mehmet Kurum
% March 15, 2017

function [S0x, x1, ax1, by1] = CalcFresnesZones(Nfz, ht, hr)

global fMHz EL0d

% Nfz : number of fresnel zones

%% Input

f0hz = fMHz * 1e6 ;
c = 3e8 ;
lambda = c / f0hz ;     % Wavelength

% Angle of Incidence of incoming signal
% EL0d = 36.3287 ;        % Elevation angle
th0td = 90 - EL0d ;
th0t = degtorad(th0td) ;

% Transmitter Height with respect to local planar ground plane
% ht = rd * cos(th0t) ;
% Transmitter/Reciever ground range
% D = rd * sin(th0t) ;
D = ht * tan(th0t) ;

%% Calculations


% Distance of specular point away from the receiver ground projection.
S0x = D / (1 + ht / hr) ;

% Parameters for Fresnel zone calculation
hd = ht - hr ;
hs = ht + hr ;
rd2 = sqrt(hd ^ 2 + D ^ 2) ;
rs = sqrt(hs ^ 2 + D ^ 2) ;
del0 = rs - rd2 ; % the shortest distance after the direct path

tanth = hd / D ; %#ok<NASGU>
sinth = hd / rd2 ;
costh = D / rd2 ;
secth = 1 / costh ;

x1 = zeros(Nfz, 1) ;
ax1 = zeros(Nfz, 1) ;
by1 = zeros(Nfz, 1) ;

for nn = 1 : Nfz
    
    % the path for the nth Fresnel Zone
    del = del0 + nn * lambda / 2 ;
    
    a = (D * secth + del) / 2 ;
    b = sqrt(del ^ 2 + 2 * D * del * secth) / 2 ;
    c = (ht + hr) / 2 ;
    
    nump = c * (a ^ 2 - b ^ 2) * sinth * costh ;
    denp = (b ^ 2 * costh ^ 2 + a ^ 2 * sinth ^ 2) ;
    p = -nump / denp ;
    x1(nn) = D / 2 + p ; % distance to the center of the ellipse
    
    % minor axis
    by1(nn) = b * sqrt(1 - c ^ 2 / denp) ;
    % major axis
    ax1(nn) = by1(nn) * a / sqrt(denp) ;
    
end


end