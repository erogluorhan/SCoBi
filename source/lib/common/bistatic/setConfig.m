% Mehmet Kurum
% Feb 25, 2017
% Mehmet Kurum
% April 6, 2017

function setConfig

%% GET GLOBAL PARAMETERS
% Transmitter Parameters
r_Tx_m = TxParams.getInstance.r_Tx_m;
% Dynamic Parameters
el0_Tx_list_deg = DynParams.getInstance.el0_Tx_list_deg;
el0_Tx_deg = el0_Tx_list_deg( ParamsManager.index_Th );

% Transmitter Geometry
rd_m = sqrt( r_Tx_m ^ 2 - (Constants.re * cos( deg2rad( el0_Tx_deg ) )) ^ 2) - Constants.re * sin( deg2rad( el0_Tx_deg ) ) ; % ~ Slant range

%% TO-DO: Work on the realistic SatGeo
%% TRANSMITTER GEOMETRY
satGeometryManuel(rd_m) ;

%% BISTATIC GEOMETRY
bistaticGeometry( rd_m ) ;


end