%% Mehmet Kurum
% Feb 25, 2017


function [rd, PH0_deg, EL0_deg, Tgt, TgtI] ...
    = satGeometry( LatSd, LonSd, LatBd, LonBd, LatGd, LonGd)

%% GET GLOBAL DIRECTORIES
dir_config = SimulationFolders.getInstance.config;


%% GET GLOBAL PARAMETERS
% Transmitter Parameters
r_Tx_m = TxParams.getInstance.r_Tx_m;


%% INITIALIZE REQUIRED PARAMETERS
% Earth Center
pE = [0, 0, 0] ;

% Earth Center Coordinate System
uxe = [1 0 0] ;
uye = [0 1 0] ;
uze = [0 0 1] ;



%% CALCULATIONS
k = r_Tx_m / Constants.re ;  % ~6.61072 : distance of Transmitter from the earth center in terms of re

% Position Vectors in eart centern frame
% Receiver ground position
LatG = degtorad(LatGd) ;
LonG = degtorad(LonGd) ;
pG = [cos(LatG) * cos(LonG), cos(LatG) * sin(LonG), sin(LatG)] ;

% Transmitter position
LatS = degtorad(LatSd) ;
LonS = degtorad(LonSd) ;
pT = k * [cos(LatS) * cos(LonS), cos(LatS) * sin(LonS), sin(LatS)] ;

% Boresight Position
LatB = degtorad(LatBd) ;
LonB = degtorad(LonBd) ;
pB = [cos(LatB) * cos(LonB), cos(LatB) * sin(LonB), sin(LatB)] ;


%% COORDINATE SYSTEMS
% Transmitter Transmit Antenna Coordinate System
TB = pB - pT ; % Transmitter to boresight

% zt
MagTS = vectorMagnitude(TB) ;
uzt = TB / MagTS ; % unit vector

% For antennas on a Transmitter, the reference plane is usually the
% equatorial plane. In most cases, linear-polarization
% is either horizontal, where the electric field is parallel to
% the plane of the equator, or vertical.

% yt
uyt = crossProduct(uzt, uze) ; % east or west?
Magya = vectorMagnitude(uyt) ;
uyt = uyt / Magya ; % unit vector

% xt
uxt = crossProduct(uyt, uzt) ;

% Local (Ground) Coordinate System in earth frame
% zz
MagpG = vectorMagnitude(pG) ;
uz = pG / MagpG ; % unit vector / normal to local surface

% xx
ux = crossProduct(uze, uz) ; % east
Magxp = vectorMagnitude(ux) ;
ux = ux / Magxp ; % unit vector

% yy
uy = crossProduct(uz, ux) ;  % north

% Transmitter to Ground (specular point) in  earh frame (center)
GT = pG - pT ;
rdk = vectorMagnitude(GT) ;
isn_ef = GT / rdk ;  % propagation vector (i_s^-)

rd = rdk * Constants.re ; % slant distance in m


%% TRANSFORMATIONS
% Transformation matrix for transforming a vector from the Earth center
% system to local ground system
Teg = [ux ; uy ; uz] ; % E -> G

% Transformation matrix for transforming a vector from the Earth center
% system to Transmitter antenna system
Tet = [uxt ; uyt ; uzt] ; % E -> T

% Antenna Coordinates in ground system
Tgt = Tet * Teg .' ;   % G -> T

% Image Tranmist Antenna Coordinates in ground system
TgtI = Tgt ;
TgtI(:, 3) = - TgtI(:, 3) ;  %  G -> TI


%% PROPAGATION VECTORS
% propagation vector from Transmitter to ground in ground system (reference)
isn = Teg * isn_ef' ;  % propagation vector (i_s^-)
isp = isn ;  
isp(3) = -isn(3) ; % propagation vector (i_s^+)

% propagation vector from Transmitter to ground in Transmitter transmit frame
isn_tf = Tet * isn_ef' ;   % propagation vector (i_s^-)

% propagation vector from image Transmitter to ground in image Transmitter transmit frame
isp_tIf = TgtI * isp ;% propagation vector (i_s^+)


%% ELEVATION AND AZIMUTH ANGLES
% Transmit Antenna Elevation
th0t = acos(-isn(3)) ;
th0t_deg = radtodeg(th0t) ;
EL0_deg = 90 - th0t_deg ;

% Transmit Antenna Azimuth
ph0t = atan2(-isn(2), -isn(1)) ;
ph0td = radtodeg(ph0t) ;
PH0_deg = 90 - ph0td ;

% % alternatives
% thd2 = acos(dot(-isn_ef, uz)) * 180 / pi ;
% phd2 = atan2(dot(-isn_ef, uy), dot(-isn_ef, ux)) * 180 / pi ;


%% ANGULAR POSITION OF THE GROUND IN TRANSMITTER SYSTEM
% off-axis angle of zt towards ground receiver
tht0d = acos(isn_tf(3)) * 180 / pi ;
% ground reciever orientation - azimuth
pht0d = atan2(isn_tf(2), isn_tf(1)) * 180 / pi ;

AngT2R_tf = [tht0d; pht0d] ;


%% SAVE ALL
% Tgt - Transformation G -> T
filename = 'Tgt' ;
writeVar(dir_config, filename, Tgt) ;

% TgrI - Transformation G -> TI
filename = 'TgtI' ;
writeVar(dir_config, filename, TgtI) ;

% rd : slant range - distance between ground and Transmitter
filename = 'rd' ;
writeVar(dir_config, filename, rd) ;

% PH0 : Azimuth angle from local north axis - clockwise
filename = 'PH0_deg' ;
writeVar(dir_config, filename, PH0_deg) ;

% EL0 : elevation angle from the local horizon
filename = 'EL0_deg' ;
writeVar(dir_config, filename, EL0_deg) ;


end


%% Cross-Product

function result = crossProduct(A,B)

if length(A) == 3 && length(B) == 3
    result = [A(2)*B(3)-A(3)*B(2), A(3)*B(1)-A(1)*B(3), A(1)*B(2)-A(2)*B(1)];
else
    disp('Error - one or more vector components missing!');
    result = [0, 0, 0] ;
end

end
