% Mehmet Kurum
% Feb 25, 2017
% Mehmet Kurum
% April 6, 2017

function SetConfig3

%% Global
global hr hpbw EL0d PH0d

%% Transmitter Satellite Geometry
% geostationary satellite radius - m
rsat = 42164 * 1e3 ;   
re = 6378 * 1e3 ;        % radius of earth - m
rd = sqrt(rsat ^ 2 - (re * cos(deg2rad(EL0d))) ^ 2) - re * sin(deg2rad(EL0d)) ; % ~ Slant range

SatGeoManuel3(rd) ;
    
%% Incoming Signal
% Angle of Incidence of incoming signal
th0td = 90 - EL0d ;
th0t = degtorad(th0td) ;

% Azimuth Angle of incoming signal
ph0td = 90 - PH0d ;
ph0t = degtorad(ph0td) ;
    
%% Receiver Orientation

% Antenna Look Angle (angle of Incidence)
th0rd = th0td - 0 ;
th0r = degtorad(th0rd) ;

% Azimuth Angle rotation of receive antenna
ph0rd = ph0td + 0 ;
ph0r = degtorad(ph0rd) ;

%% Bistatic Geomtery

BiGeo3(hr, rd, th0t, ph0t, th0r, ph0r, hpbw) ;


end