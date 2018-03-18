% Mehmet Kurum
% Feb 25, 2017
% Mehmet Kurum
% April 6, 2017

function setConfig

% Get global parameters
EL0_deg = SatParams.getInstance.EL0_deg( ParamsManager.index_Th );
PH0_deg = SatParams.getInstance.PH0_deg( ParamsManager.index_Ph );
rsat = SatParams.getInstance.rsat;
hpbw_deg = RecParams.getInstance.hpbw_deg;


%% Transmitter Satellite Geometry
rd = sqrt( rsat ^ 2 - (Constants.re * cos( deg2rad( EL0_deg ) )) ^ 2) - Constants.re * sin( deg2rad( EL0_deg ) ) ; % ~ Slant range

satGeometryManuel(rd) ;
    
%% Incoming Signal
% Angle of Incidence of incoming signal
th0t_deg = 90 - EL0_deg;
th0t = degtorad(th0t_deg) ;

% Azimuth Angle (standard spherical coords) of transmitter's position
% If it was the incoming signal's azimuth, we would add 180 degrees
ph0t_deg = 90 - PH0_deg ;
ph0t = degtorad(ph0t_deg) ;
    
%% Receiver Orientation

% Antenna Look Angle (angle of Incidence)
th0r_deg = th0t_deg - 0 ;
th0r = degtorad(th0r_deg) ;

% Azimuth Angle rotation of receive antenna
ph0r_deg = ph0t_deg + 0 ;
ph0r = degtorad(ph0r_deg) ;

%% Bistatic Geomtery

bistaticGeometry(rd, th0t, ph0t, th0r, ph0r, hpbw_deg ) ;


end