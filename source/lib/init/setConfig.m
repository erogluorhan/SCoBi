% Mehmet Kurum
% Feb 25, 2017
% Mehmet Kurum
% April 6, 2017

function setConfig

%% GET GLOBAL PARAMETERS
% Satellite Parameters
EL0_deg = SatParams.getInstance.EL0_deg( ParamsManager.index_Th );
rsat_m = SatParams.getInstance.rsat_m;

% Transmitter Satellite Geometry
rd_m = sqrt( rsat_m ^ 2 - (Constants.re * cos( deg2rad( EL0_deg ) )) ^ 2) - Constants.re * sin( deg2rad( EL0_deg ) ) ; % ~ Slant range

%% SATELLITE GEOMETRY
satGeometryManuel(rd_m) ;

%% BISTATIC GEOMETRY
bistaticGeometry( rd_m ) ;


end