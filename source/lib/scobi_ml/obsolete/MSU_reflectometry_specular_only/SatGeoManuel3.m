% Mehmet Kurum
% April 6, 2017


function SatGeoManuel3(rd)

global EL0d PH0d FolderPath_Config

%% Incoming Signal
% Angle of Incidence of incoming signal
th0td = 90 - EL0d ;
th0t = degtorad(th0td) ;

% Azimuth Angle of incoming signal
ph0td = 90 - PH0d + 180 ;  % 180 degrees added to face transmitter and receiver aligned to each other.
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

% Alternative if ph0td = 90 - PH0d
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

%% Receiver Antenna Coordinate System
% transfroming antenna to reference
uxt = (AntRott * ux')' ;
uyt = (AntRott * uy')' ;
uzt = (AntRott * uz')' ;

sum(uxt .* uyt) ;
sum(uxt .* uzt) ;
sum(uzt .* uyt) ;

%% Image of Receiver Antenna Coordinate System
% transfroming antenna to reference
uxtI = [uxt(1), uxt(2), -uxt(3)] ;
uytI = [uyt(1), uyt(2), -uyt(3)] ;
uztI = [uzt(1), uzt(2), -uzt(3)] ;

%% Transmfromation
% Transformation matrix for transforming a vector from the ground frame
% to receiver system
Tgt = [uxt ; uyt ; uzt] ; % G -> T

% Transformation matrix for transforming a vector from the ground frame
% to Image receiver system
TgtI = [uxtI ; uytI ; uztI] ; % G -> TI

%% Saving. . . 

% Tgt - Transformation G -> T
% TgrI - Transformation G -> TI
% rd : slant range - distance between ground and satellite
% PH0 : Azimuth angle from local north axis - clockwise
% EL0 : elevation angle from the local horizon

filename = 'Tgt' ;
write_var(FolderPath_Config, filename, Tgt) ;

filename = 'TgtI' ;
write_var(FolderPath_Config, filename, TgtI) ;

filename = 'rd' ;
write_var(FolderPath_Config, filename, rd) ;

filename = 'PH0d' ;
write_var(FolderPath_Config, filename, PH0d) ;

filename = 'EL0d' ;
write_var(FolderPath_Config, filename, EL0d) ;


end