%% Mehmet Kurum
% Feb 25, 2017 

function bistaticGeometry( rd_m )

%% GET GLOBAL DIRECTORIES
dir_config = SimulationFolders.getInstance.config;


%% GET GLOBAL PARAMETERS
% Satellite Parameters
th0t_deg = SatParams.getInstance.th0_list_deg( ParamsManager.index_Th );
PH0_deg = SatParams.getInstance.PH0_list_deg( ParamsManager.index_Ph );
% Receiver Parameters
hpbw_deg = RecParams.getInstance.hpbw_deg;
hr_m = RecParams.getInstance.hr_m;


%% INITIALIZE REQUIRED PARAMETERS
% Ground Reference Coordinate System (East-North-Up)
% located on the ground where the receiver is projected
ux = [1 0 0] ;
uy = [0 1 0] ;
uz = [0 0 1] ;

% G : Antenna projection point on the ground
% Center of reference coordinate system
pG_m = [0; 0; 0] ;


%% CALCULATIONS
% Antenna Parameters/Orientation - Receiver
% 3dB Beamwidths
% Bthrd = 40 ;
Bthr_rad = degtorad(hpbw_deg) ;
% Bphrd = 40 ;
Bphr_rad = degtorad(hpbw_deg) ;


% Incoming Signal
% Angle of Incidence of incoming signal
th0t_rad = degtorad(th0t_deg) ;

% Azimuth Angle (standard spherical coords) of transmitter's position
% If it was the incoming signal's azimuth, we would add 180 degrees
ph0t_deg = 90 - PH0_deg ;
ph0t_rad = degtorad(ph0t_deg) ;


% Receiver Orientation
% Antenna Look Angle (angle of Incidence)
th0r_deg = th0t_deg - 0 ;
th0r_rad = degtorad(th0r_deg) ;

% Azimuth Angle rotation of receive antenna
ph0r_deg = ph0t_deg + 0 ;
ph0r_rad = degtorad(ph0r_deg) ;


% Antenna Rotation Matrix for reciever
% Rotation about z-axis (Azimuth rotation)
AntRotZr = [cos(ph0r_rad) -sin(ph0r_rad) 0;
    sin(ph0r_rad) cos(ph0r_rad)   0    ;
    0 0 1] ;

% Rotation about y-axis (Elevation rotation)
AntRotYr = [-cos(th0r_rad) 0  sin(th0r_rad);
    0       1     0    ;
    -sin(th0r_rad) 0 -cos(th0r_rad)] ;

% Rotation in both azimuth and elevation
AntRotr = AntRotZr  * AntRotYr ;


% Azimuth direction of transmitted signal.
AntRotZt = [cos(ph0t_rad) -sin(ph0t_rad) 0;
    sin(ph0t_rad) cos(ph0t_rad)   0    ;
    0 0 1] ;


% Transmitter position
% T : Transmitter Antenna position
pT_m = rd_m * [sin(th0t_rad) * cos(ph0t_rad); sin(th0t_rad) * sin(ph0t_rad); cos(th0t_rad)] ;
ht_m = pT_m(3) ;
% TI : Transmitter Image Antenna position
pTI_m = [pT_m(1); pT_m(2); -ht_m] ;

% A : Antenna projection point on the ground
% % pHt = [pT_m(1), pT_m(2), 0] ;

% Specular point and 1st Fresnel zone ellipse
[S0x_m, x1_m, ax1_m, by1_m] = calcFresnelZones(ht_m, hr_m) ;

pS_m = [S0x_m; 0; 0] ;
pS2_m = AntRotZt * pS_m  ; % specular point location

% Center of first Fresnel ellipse
pSc_m = [x1_m(1); 0; 0] ;
pSc2_m = AntRotZt * pSc_m  ;

ellipse_s_m = [ax1_m, by1_m] ; % specular point Fresnel zone
ellipse_s_centers_m = AntRotZt * [x1_m'; zeros(1,10); zeros(1,10)] ;

%% Reciever footprint - ellipse

% major axis
dar_m = hr_m * (tan(th0r_rad + Bthr_rad / 2) - tan(th0r_rad - Bthr_rad / 2)) ;

% angle of incidnce at the center of the ellipse
thcr_rad = atan(tan(th0r_rad + Bthr_rad / 2) - dar_m / hr_m / 2) ;
thcr_deg = radtodeg(thcr_rad) ; %#ok<NASGU>

% minor axis
dbr_m = 2 * hr_m * sec(thcr_rad) * tan(Bphr_rad / 2) ;

ellipse_r_m = [dar_m; dbr_m] ; % receiver footprint

%% Receiver position and pointing locations

% R : Receiver Antenna position
pR_m = [0; 0; hr_m] ;
% pR2 = AntRotZr * pR_m ;

% RI : Image Receiver Antenna position
pRI_m = [0; 0; -hr_m] ;

% G : Antenna projection point on the ground
% Center of reference coordinate system
pG2_m = AntRotZr * pG_m ;

% B : Boresight point
pBr_m = [hr_m * tan(th0r_rad); 0; 0] ;
pBr2_m = AntRotZr * pBr_m ;

% Ellipse Center - C
pCr_m = [hr_m * tan(thcr_rad); 0; 0] ;
pCr2_m = AntRotZr * pCr_m ;

%% all points
% pT_m: Transmitter, pS2_m: specular point, pR2: receiver, pG2_m: ground (reference), pBr2_m:boresight,
% pCr2_m: center of footprint, pSc2_m: center of fresnel zone
AllPoints_m = [pT_m, pTI_m, pS2_m, pR_m, pRI_m, pG2_m, pBr2_m, pCr2_m, pSc2_m] ;

%% Satellite to Receiver
RT_m = pR_m - pT_m ;
magRT_m = vectorMagnitude(RT_m) ;
idn = RT_m / magRT_m ;  % propagation vector (i_d^-)

%% Satellite to Specular point
ST_m = pS2_m - pT_m ;
magST_m = vectorMagnitude(ST_m) ;
isn = ST_m / magST_m ;  % propagation vector (i_s^-)

%% Specular point to Reciever
RS_m = pR_m - pS2_m ;
magRS_m = vectorMagnitude(RS_m) ;
osp = RS_m / magRS_m ;  % propagation vector (o_s^+)

%% Specular point to Reciever
RIS_m = pRI_m - pS2_m ;
magRIS_m = vectorMagnitude(RIS_m) ;
osn = RIS_m / magRIS_m ;  % propagation vector (o_s^-)

%% Receiver Antenna Coordinate System
uxr = (AntRotr * ux')' ;
uyr = (AntRotr * uy')' ;
uzr = (AntRotr * uz')' ;

%% Image of Receiver Antenna Coordinate System
uxrI = [uxr(1), uxr(2), -uxr(3)] ;
uyrI = [uyr(1), uyr(2), -uyr(3)] ;
uzrI = [uzr(1), uzr(2), -uzr(3)] ;

%% Specular Point Coordinate System
uxs = (AntRotZt * ux')' ; 
uys = (AntRotZt * uy')' ;
uzs = (AntRotZt * uz')' ;

%% Transformations

% Transformation matrix for transforming a vector from the ground frame
% to local (specular) ground system
Tgs = [uxs ; uys ; uzs] ; % G -> S

% Transformation matrix for transforming a vector from the ground frame
% to receiver system
Tgr = [uxr ; uyr ; uzr] ; % G -> R

% Transformation matrix for transforming a vector from the ground frame
% to Image receiver system
TgrI = [uxrI ; uyrI ; uzrI] ; % G -> RI

% % Receive antenna Coordinates in local (specular) ground system
% Trs = Tgs * Tgr .' ;   % R -> S

%% The incidence angle on the receiver in receiver coordinates

% propagation vector from satellite to reciever in receiver antenna system
idn_rf = Tgr * idn ;

% FROM SAT TO REC
% off-axis angle of zr towards satellite
th0 = acos(-idn_rf(3)) * 180 / pi ;
%  orientation - azimuth
ph0 = atan2(-idn_rf(2), -idn_rf(1)) * 180 / pi ;

AngT2R_rf = [th0; convertAngleTo360Range( ph0 )]  ;

% propagation vector from specular point to receiver in receiver antenna system
osp_rf = Tgr * osp ;

% FROM SPEC TO REC
% off-axis angle of zr towards specular point
th0 = acos(-osp_rf(3)) * 180 / pi ;
%  orientation - azimuth
ph0 = atan2(-osp_rf(2), -osp_rf(1)) * 180 / pi ;

AngS2R_rf = [th0; convertAngleTo360Range( ph0 ) ]  ;


%% The incidence angle in spacular frame
% propagation vector from satellite to ground in local (specular) ground system
isn_sf = Tgs * isn ;

% Sat-to-SP
% off-axis angle of zs towards satellite
th0 = acos(-isn_sf(3)) * 180 / pi ;
%  orientation - azimuth
ph0 = atan2(-isn_sf(2), -isn_sf(1)) * 180 / pi ;

AngT2S_sf = [th0; convertAngleTo360Range( ph0 ) ]  ;


%% SAVE ALL
% idn -  propagation vector (i_d^-)
filename = 'idn' ;
writeVar(dir_config, filename, idn) ;

% isn - propagation vector (i_s^-)
filename = 'isn' ;
writeVar(dir_config, filename, isn) ;

% osp - propagation vector (o_s^+)
filename = 'osp' ;
writeVar(dir_config, filename, osp) ;

% osn - propagation vector (o_s^-)
filename = 'osn' ;
writeVar(dir_config, filename, osn) ;

% Tgs - Transformation G -> S
filename = 'Tgs' ;
writeVar(dir_config, filename, Tgs) ;

% Tgr - Transformation G -> R
filename = 'Tgr' ;
writeVar(dir_config, filename, Tgr) ;

% TgrI - Transformation G -> RI
filename = 'TgrI' ;
writeVar(dir_config, filename, TgrI) ;

% AntRotZr - Receiver Rotation about z-axis (Azimuth rotation)
filename = 'AntRotZr' ;
writeVar(dir_config, filename, AntRotZr) ;

% AntRotYr - Receiver Rotation about y-axis (Elevation rotation)
filename = 'AntRotYr' ;
writeVar(dir_config, filename, AntRotYr) ;

% AntRotr  - Receiver Rotation in both azimuth and elevation
filename = 'AntRotr' ;
writeVar(dir_config, filename, AntRotr) ;

% AntRotZt -  Rotation matrix that describes azimuth direction of transmitted signal.
filename = 'AntRotZt' ;
writeVar(dir_config, filename, AntRotZt) ;

% ellipse_s_m - specular point Fresnel zone [major and minor axes]
filename = 'ellipse_s_m' ;
writeVar(dir_config, filename, ellipse_s_m) ;

filename = 'ellipse_s_centers_m' ;
writeVar(dir_config, filename, ellipse_s_centers_m) ;

% ellipse_r_m - receiver footprint ellipse [major and minor axes]
filename = 'ellipse_r_m' ;
writeVar(dir_config, filename, ellipse_r_m) ;

% AllPoints_m - pT_m, pS2_m, pR_m, pG2_m, pBr2_m, pCr2_m, pSc2_m in ground (refrence) frame (G)
filename = 'AllPoints_m' ;
writeVar(dir_config, filename, AllPoints_m) ;

% AngT2R_rf - Incidence angle (T -> R) in receiver frame (R)
filename = 'AngT2R_rf' ;
writeVar(dir_config, filename, AngT2R_rf) ;

% AngS2R_rf - Incidence angle (S -> R) in receiver frame (R)
filename = 'AngS2R_rf' ;
writeVar(dir_config, filename, AngS2R_rf) ;

% AngT2S_sf - Incidence angle (T -> S) in specular frame (S)
filename = 'AngT2S_sf' ;
writeVar(dir_config, filename, AngT2S_sf) ;


    
end
