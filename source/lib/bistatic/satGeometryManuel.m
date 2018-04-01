% Mehmet Kurum
% April 6, 2017


function satGeometryManuel(rd_m)

%%GET GLOBAL DIRECTORIES
dir_config = SimulationFolders.getInstance.config;


%% GET GLOBAL PARAMETERS
% Satellite Parameters
EL0_deg = SatParams.getInstance.EL0_list_deg( ParamsManager.index_Th );
PH0_deg = SatParams.getInstance.PH0_list_deg( ParamsManager.index_Ph );

%% Incoming Signal
% Angle of Incidence of incoming signal
th0td = 90 - EL0_deg;
th0t = degtorad(th0td) ;

% Azimuth Angle (in standard spherical coords) of incoming signal
ph0td = 90 - PH0_deg + 180 ;  % 180 degrees added to face transmitter and receiver aligned to each other.
ph0t = degtorad(ph0td) ;


%% Antenna Rotation Matrix for Transmitter

% Rotation about z-axis (Azimuth rotation)
AntRotZt = [cos(ph0t) -sin(ph0t) 0;
    sin(ph0t) cos(ph0t)   0    ;
    0 0 1] ;

% Rotation about y-axis (Elevation rotation)
AntRotYt = [-cos(th0t) 0  sin(th0t);
    0       1     0    ;
    -sin(th0t) 0 -cos(th0t)] ;

% Alternative if ph0td = 90 - SatParams.PH0_deg
% % Rotation about y-axis (Elevation rotation)
% AntRotYt = [-cos(th0t) 0  -sin(th0t);
%     0       1     0    ;
%     sin(th0t) 0 -cos(th0t)] ;

% Rotation in both azimuth and elevation
AntRott = AntRotZt  * AntRotYt ;

%% Ground Reference Coordinate System (East-North-Up)
% located on the ground where the receiver is projected
ux = [1 0 0] ;
uy = [0 1 0] ;
uz = [0 0 1] ;

%% Transmitter Antenna Coordinate System
% transfroming antenna to reference
uxt = (AntRott * ux')' ;
uyt = (AntRott * uy')' ;
uzt = (AntRott * uz')' ;

sum(uxt .* uyt);
sum(uxt .* uzt);
sum(uzt .* uyt);

%% Image of Transmitter Antenna Coordinate System
% transfroming antenna to reference
uxtI = [uxt(1), uxt(2), -uxt(3)] ;
uytI = [uyt(1), uyt(2), -uyt(3)] ;
uztI = [uzt(1), uzt(2), -uzt(3)] ;

%% Transmfromation
% Transformation matrix for transforming a vector from the ground frame
% to transmitter system
Tgt = [uxt ; uyt ; uzt] ; % G -> T

% Transformation matrix for transforming a vector from the ground frame
% to Image transmitter system
TgtI = [uxtI ; uytI ; uztI] ; % G -> TI


%% SAVE ALL
% Tgt - Transformation G -> T
filename = 'Tgt' ;
writeVar(dir_config, filename, Tgt) ;

% TgtI - Transformation G -> TI
filename = 'TgtI' ;
writeVar(dir_config, filename, TgtI) ;

% rd_m : slant range - distance between ground and satellite
filename = 'rd_m' ;
writeVar(dir_config, filename, rd_m) ;

% PH0 : Azimuth angle from local north axis - clockwise
filename = 'PH0_deg' ;
writeVar(dir_config, filename, PH0_deg ) ;

% EL0 : elevation angle from the local horizon
filename = 'EL0_deg' ;
writeVar(dir_config, filename, EL0_deg ) ;


end