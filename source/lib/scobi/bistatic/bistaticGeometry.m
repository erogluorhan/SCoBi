%% Mehmet Kurum
% Feb 25, 2017 

function [rd_m, idn, isn, osp, osn, Tgs, Tgr, TgrI, AntRotZ_Rx, ...
          AntRotY_Rx, AntRot_Rx, AntRotZ_Tx, ellipse_FP_Rx_m, ...
          AllPoints_m, AngT2R_rf, AngS2R_rf, AngT2S_sf] = bistaticGeometry


%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Transmitter Parameters
r_Tx_m = TxParams.getInstance.r_Tx_m;
% Configuration Parameters
th0_Tx_list_deg = ConfigParams.getInstance.th0_Tx_list_deg;
th0_Tx_deg = th0_Tx_list_deg( sim_counter );   % Currrent theta ( Incidence Angle )
th0_Tx_rad = deg2rad(th0_Tx_deg) ;
el0_Tx_deg = 90 - th0_Tx_deg;
el0_Tx_rad = deg2rad(el0_Tx_deg) ;
ph0_Tx_list_deg = ConfigParams.getInstance.ph0_Tx_list_deg;
ph0_Tx_deg = ph0_Tx_list_deg( ParamsManager.sim_counter );   % Current phi (Azimuth Angle (standard spherical coords) of transmitter's position)
ph0_Tx_deg = 90 - ph0_Tx_deg ;                            % If it was the incoming signal's azimuth, we would add 180 degrees
ph0_Tx_rad = deg2rad(ph0_Tx_deg) ;
% Receiver Parameters
hr_m = RxParams.getInstance.hr_m;                   % Receiver altitude
ant_pat_Rx_id = RxParams.getInstance.ant_pat_Rx_id; % Receiver antenna pattern generation method (Look at Constants.Rx_ant_pats)
th0_Rx_deg = RxParams.getInstance.th0_Rx_deg;       % Antenna Looking Angle (angle of Incidence)
th0_Rx_rad = deg2rad(th0_Rx_deg) ;
ph0_Rx_deg = RxParams.getInstance.ph0_Rx_deg;       % Azimuth Angle of Receiver position
ph0_Rx_deg = 90 - ph0_Rx_deg;
ph0_Rx_rad = deg2rad(ph0_Rx_deg) ;
% If recevier antenna pattern is Generalized-Gaussian, its footprint
% ellipse will be calculated below
if ant_pat_Rx_id == Constants.id_Rx_GG
    hpbw_deg = RxGGParams.getInstance.hpbw_deg;
end


%% INITIALIZE REQUIRED PARAMETERS
% Ground Reference Coordinate System (East-North-Up)
% located on the ground where the receiver is projected
ux = [1 0 0] ;
uy = [0 1 0] ;
uz = [0 0 1] ;

% G : Antenna projection point on the ground
% Center of reference coordinate system
pos_Gnd_m = [0; 0; 0] ;


%% CALCULATIONS
% Antenna Parameters/Orientation - Receiver

% If recevier antenna pattern is Generalized-Gaussian, its footprint
% ellipse will be calculated below
if ant_pat_Rx_id == Constants.id_Rx_GG
    % 3dB Beamwidths
    % Bthrd = 40 ;
    Bthr_rad = deg2rad(hpbw_deg) ;
    % Bphrd = 40 ;
    Bphr_rad = deg2rad(hpbw_deg) ;
end


% Antenna Rotation Matrix for reciever
% Rotation about z-axis (Azimuth rotation)
AntRotZ_Rx = [ cos(ph0_Rx_rad), -sin(ph0_Rx_rad), 0 ;
             sin(ph0_Rx_rad),  cos(ph0_Rx_rad), 0 ;
             0,                0,               1 ];

% Rotation about y-axis (Elevation rotation)
AntRotY_Rx = [ -cos(th0_Rx_rad), 0,  sin(th0_Rx_rad) ;
              0,               1,  0 ;
             -sin(th0_Rx_rad), 0, -cos(th0_Rx_rad) ];

% Rotation in both azimuth and elevation
AntRot_Rx = AntRotZ_Rx  * AntRotY_Rx ;


% Azimuth direction of transmitted signal.
AntRotZ_Tx = [ cos(ph0_Tx_rad), -sin(ph0_Tx_rad), 0 ;
             sin(ph0_Tx_rad),  cos(ph0_Tx_rad), 0 ;
             0,                0,               1 ];


% Transmitter position
% Slant range (Approximated)
rd_m = sqrt( r_Tx_m ^ 2 - (Constants.re * cos( el0_Tx_rad ) ) ^ 2) - Constants.re * sin( el0_Tx_rad ) ;
% T : Transmitter Antenna position
pos_Tx_m = rd_m * [sin(th0_Tx_rad) * cos(ph0_Tx_rad); sin(th0_Tx_rad) * sin(ph0_Tx_rad); cos(th0_Tx_rad)] ;
ht_m = pos_Tx_m(3) ;
% TI : Transmitter Image Antenna position
pos_TxI_m = [pos_Tx_m(1); pos_Tx_m(2); -ht_m] ;

% A : Antenna projection point on the ground
% % pHt = [pos_Tx_m(1), pos_Tx_m(2), 0] ;

% Specular point and 1st Fresnel zone ellipse
Nfz = 1;  % number of fresnel zones is 1 for specular 
[S0x_m, x1_m, ~, ~] = calcFresnelZones(ht_m, hr_m, Nfz) ;

pos_SP_m = [S0x_m; 0; 0] ;
pos_SP_m = AntRotZ_Tx * pos_SP_m  ; % specular point location

% Center of first Fresnel ellipse
pos_FZ_m = [x1_m(1); 0; 0] ;
pos_FZ_m = AntRotZ_Tx * pos_FZ_m  ;


%% RECEIVER FOOTPRINT - ELLIPSE
% If recevier antenna pattern is Generalized-Gaussian, its footprint
% ellipse will be calculated below
ellipse_FP_Rx_m = [];

if ant_pat_Rx_id == Constants.id_Rx_GG
    
    % major axis
    dar_m = hr_m * (tan(th0_Rx_rad + Bthr_rad / 2) - tan(th0_Rx_rad - Bthr_rad / 2)) ;

    % angle of incidnce at the center of the ellipse
    thcr_rad = atan(tan(th0_Rx_rad + Bthr_rad / 2) - dar_m / hr_m / 2) ;
    thcr_deg = rad2deg(thcr_rad) ; %#ok<NASGU>

    % minor axis
    dbr_m = 2 * hr_m * sec(thcr_rad) * tan(Bphr_rad / 2) ;

    ellipse_FP_Rx_m = [dar_m; dbr_m] ; % receiver footprint
    
end


%% RECEIVER POSITION AND POINTING LOCATIONS
% R : Receiver Antenna position
pos_Rx_m = [0; 0; hr_m] ;
% pos_Rx_m = AntRotZ_Rx * pos_Rx_m ;

% RI : Image Receiver Antenna position
pos_RxI_m = [0; 0; -hr_m] ;

% G : Antenna projection point on the ground
% Center of reference coordinate system
pos_Gnd_m = AntRotZ_Rx * pos_Gnd_m ;

% B : Boresight point
pos_B_Rx_m = [hr_m * tan(th0_Rx_rad); 0; 0] ;
pos_B_Rx_m = AntRotZ_Rx * pos_B_Rx_m ;

% If recevier antenna pattern is Generalized-Gaussian, its footprint
% ellipse will be calculated below
if ant_pat_Rx_id == Constants.id_Rx_GG
    % Ellipse Center - C
    pos_FP_Rx_m = [hr_m * tan(thcr_rad); 0; 0] ;
    pos_FP_Rx_m = AntRotZ_Rx * pos_FP_Rx_m ;
else
    pos_FP_Rx_m = [];
end


%% ALL POINTS
% pos_Tx_m: Transmitter, 
% pos_SP_m: specular point, 
% pos_Rx_m: Receiver, 
% pos_Gnd_m: Ground (reference), 
% pos_B_Rx_m: Receiver Boresight,
% pos_FP_Rx_m: Center of Receiver's FootPrint, 
% pos_FZ_m: Center of Fresnel Zone
AllPoints_m = [pos_Tx_m, pos_TxI_m, pos_SP_m, pos_Rx_m, pos_RxI_m, pos_Gnd_m, pos_B_Rx_m, pos_FP_Rx_m, pos_FZ_m] ;


%% PROPAGATION VECTORS
% Transmitter to Receiver
RT_m = pos_Rx_m - pos_Tx_m ;
magRT_m = vectorMagnitude(RT_m) ;
idn = RT_m / magRT_m ;  % propagation vector (i_d^-)

% Transmitter to Specular point
ST_m = pos_SP_m - pos_Tx_m ;
magST_m = vectorMagnitude(ST_m) ;
isn = ST_m / magST_m ;  % propagation vector (i_s^-)

% Specular point to Reciever
RS_m = pos_Rx_m - pos_SP_m ;
magRS_m = vectorMagnitude(RS_m) ;
osp = RS_m / magRS_m ;  % propagation vector (o_s^+)

% Specular point to Reciever
RIS_m = pos_RxI_m - pos_SP_m ;
magRIS_m = vectorMagnitude(RIS_m) ;
osn = RIS_m / magRIS_m ;  % propagation vector (o_s^-)


%% COORDINATE SYSTEMS
% Receiver Antenna Coordinate System
uxr = (AntRot_Rx * ux')' ;
uyr = (AntRot_Rx * uy')' ;
uzr = (AntRot_Rx * uz')' ;

% Image of Receiver Antenna Coordinate System
uxrI = [uxr(1), uxr(2), -uxr(3)] ;
uyrI = [uyr(1), uyr(2), -uyr(3)] ;
uzrI = [uzr(1), uzr(2), -uzr(3)] ;

% Specular Point Coordinate System
uxs = (AntRotZ_Tx * ux')' ; 
uys = (AntRotZ_Tx * uy')' ;
uzs = (AntRotZ_Tx * uz')' ;


%% TRANSFORMATIONS
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
% propagation vector from Transmitter to reciever in receiver antenna system
idn_rf = Tgr * idn ;


%% ANGLES
% Transmitter to Receiver
% off-axis angle of zr towards Transmitter
th0 = acos(-idn_rf(3)) * 180 / pi ;
%  orientation - azimuth
ph0 = atan2(-idn_rf(2), -idn_rf(1)) * 180 / pi ;

AngT2R_rf = [th0; convertAngleTo360Range( ph0 )]  ;

% propagation vector from specular point to receiver in receiver antenna system
osp_rf = Tgr * osp ;

% Specular Point to Receiver
% off-axis angle of zr towards specular point
th0 = acos(-osp_rf(3)) * 180 / pi ;
%  orientation - azimuth
ph0 = atan2(-osp_rf(2), -osp_rf(1)) * 180 / pi ;

AngS2R_rf = [th0; convertAngleTo360Range( ph0 ) ]  ;


%% The incidence angle in spacular frame
% propagation vector from Transmitter to ground in local (specular) ground system
isn_sf = Tgs * isn ;

% Transmitter to Specular Point
% off-axis angle of zs towards Transmitter
th0 = acos(-isn_sf(3)) * 180 / pi ;
%  orientation - azimuth
ph0 = atan2(-isn_sf(2), -isn_sf(1)) * 180 / pi ;

AngT2S_sf = [th0; convertAngleTo360Range( ph0 ) ]  ;
    
end
