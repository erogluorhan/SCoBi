%% Mehmet Kurum
% Feb 25, 2017 

function BiGeo3(hr, rd, th0t, ph0t, th0r, ph0r, Bthrd)

%% Global
global FolderPath_Config 

%% Antenna Parameters/Orientation - Receiver
% 3dB Beamwidths
% Bthrd = 40 ;
Bthr = degtorad(Bthrd) ;
% Bphrd = 40 ;
Bphr = degtorad(Bthrd) ;


%% Antenna Rotation Matrix for reciever

% Rotation about z-axis (Azimuth rotation)
AntRotZr = [cos(ph0r) -sin(ph0r) 0;
    sin(ph0r) cos(ph0r)   0    ;
    0 0 1] ;

% Rotation about y-axis (Elevation rotation)
AntRotYr = [-cos(th0r) 0  sin(th0r);
    0       1     0    ;
    -sin(th0r) 0 -cos(th0r)] ;

% Rotation in both azimuth and elevation
AntRotr = AntRotZr  * AntRotYr ;

%% Azimuth direction of transmitted signal.

AntRotZt = [cos(ph0t) -sin(ph0t) 0;
    sin(ph0t) cos(ph0t)   0    ;
    0 0 1] ;

%% Transmitter position

% T : Transmitter Antenna position
pT = rd * [sin(th0t) * cos(ph0t); sin(th0t) * sin(ph0t); cos(th0t)] ;
ht = pT(3) ;
% TI : Transmitter Image Antenna position
pTI = [pT(1); pT(2); -ht] ;

% A : Antenna projection point on the ground
% % pHt = [pT(1), pT(2), 0] ;

%% Specular point and 1st Fresnel zone ellipse
Nfz = 1 ; % local
[S0x, x1, ax1, by1] = CalcFresnesZones(Nfz, ht, hr) ;

pS = [S0x; 0; 0] ;
pS2 = AntRotZt * pS  ; % specular point location

% Center of first Fresnel ellipse
pSc = [x1(1); 0; 0] ;
pSc2 = AntRotZt * pSc  ;

ellipse_s = [ax1, by1] ; % specular point Fresnel zone

%% Reciever footprint - ellipse

% major axis
dar = hr * (tan(th0r + Bthr / 2) - tan(th0r - Bthr / 2)) ;

% angle of incidnce at the center of the ellipse
thcr = atan(tan(th0r + Bthr / 2) - dar / hr / 2) ;
thcrd = radtodeg(thcr) ; %#ok<NASGU>

% minor axis
dbr = 2 * hr * sec(thcr) * tan(Bphr / 2) ;

ellipse_r = [dar; dbr] ; % receiver footprint

%% Receiver position and pointing locations

% R : Receiver Antenna position
pR = [0; 0; hr] ;
% pR2 = AntRotZr * pR ;

% RI : Image Receiver Antenna position
pRI = [0; 0; -hr] ;

% G : Antenna projection point on the ground
% Center of reference coordinate system
pG = [0; 0; 0] ;
pG2 = AntRotZr * pG ;

% B : Boresight point
pBr = [hr * tan(th0r); 0; 0] ;
pBr2 = AntRotZr * pBr ;

% Ellipse Center - C
pCr = [hr * tan(thcr); 0; 0] ;
pCr2 = AntRotZr * pCr ;

%% all points
% pT: Transmitter, pS2: specular point, pR2: receiver, pG2: ground (reference), pBr2:boresight,
% pCr2: center of footprint, pSc2: center of fresnel zone
AllPoints = [pT, pTI, pS2, pR, pRI, pG2, pBr2, pCr2, pSc2] ;

%% Satellite to Receiver
RT = pR - pT ;
magRT = vectorMag(RT) ;
idn = RT / magRT ;  % propagation vector (i_d^-)

%% Satellite to Specular point
ST = pS2 - pT ;
magST = vectorMag(ST) ;
isn = ST / magST ;  % propagation vector (i_s^-)

%% Specular point to Reciever
RS = pR - pS2 ;
magRS = vectorMag(RS) ;
osp = RS / magRS ;  % propagation vector (o_s^+)

%% Specular point to Reciever
RIS = pRI - pS2 ;
magRIS = vectorMag(RIS) ;
osn = RIS / magRIS ;  % propagation vector (o_s^-)

%% Ground Reference Coordinate System (East-North-Up)
% located on the ground where the receiver is projected
ux = [1 0 0] ;
uy = [0 1 0] ;
uz = [0 0 1] ;

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

AngT2R_rf = [th0; ph0]  ;

% propagation vector from specular point to receiver in receiver antenna system
osp_rf = Tgr * osp ;

% FROM SPEC TO REC
% off-axis angle of zr towards specular point
th0 = acos(-osp_rf(3)) * 180 / pi ;
%  orientation - azimuth
ph0 = atan2(-osp_rf(2), -osp_rf(1)) * 180 / pi ;

AngS2R_rf = [th0; ph0]  ;


%% The incidence angle in spacular frame

% propagation vector from satellite to ground in local (specular) ground system
isn_sf = Tgs * isn ;

% FROM SAT TO SPEC
% off-axis angle of zs towards satellite
th0 = acos(-isn_sf(3)) * 180 / pi ;
%  orientation - azimuth
ph0 = atan2(-isn_sf(2), -isn_sf(1)) * 180 / pi ;

AngT2S_sf = [th0; ph0]  ;


%% Saving. . . 

% AntRotZr - Receiver Rotation about z-axis (Azimuth rotation)
% AntRotYr - Receiver Rotation about y-axis (Elevation rotation)
% AntRotr  - Receiver Rotation in both azimuth and elevation
% AntRotZt -  Rotation matrix that describes azimuth direction of transmitted signal.
% ellipse_s - specular point Fresnel zone [major and minor axes]
% ellipse_r - receiver footprint ellipse [major and minor axes]
% AllPoints - pT, pS2, pR, pG2, pBr2, pCr2, pSc2 in ground (refrence) frame (G)
% AngT2R_rf - Incidence angle (T -> R) in receiver frame (R)
% AngS2R_rf - Incidence angle (S -> R) in receiver frame (R)
% AngT2S_sf - Incidence angle (T -> S) in specular frame (S)
% Tgs - Transformation G -> S
% Tgr - Transformation G -> R
% TgrI - Transformation G -> RI
% idn -  propagation vector (i_d^-)
% isn - propagation vector (i_s^-)
% osp - propagation vector (o_s^+)
% osn - propagation vector (o_s^-)

filename = 'idn' ;
write_var(FolderPath_Config, filename, idn) ;

filename = 'isn' ;
write_var(FolderPath_Config, filename, isn) ;

filename = 'osp' ;
write_var(FolderPath_Config, filename, osp) ;

filename = 'osn' ;
write_var(FolderPath_Config, filename, osn) ;

filename = 'Tgs' ;
write_var(FolderPath_Config, filename, Tgs) ;

filename = 'Tgr' ;
write_var(FolderPath_Config, filename, Tgr) ;

filename = 'TgrI' ;
write_var(FolderPath_Config, filename, TgrI) ;

filename = 'AntRotZr' ;
write_var(FolderPath_Config, filename, AntRotZr) ;

filename = 'AntRotYr' ;
write_var(FolderPath_Config, filename, AntRotYr) ;

filename = 'AntRotr' ;
write_var(FolderPath_Config, filename, AntRotr) ;

filename = 'AntRotZt' ;
write_var(FolderPath_Config, filename, AntRotZt) ;

filename = 'ellipse_s' ;
write_var(FolderPath_Config, filename, ellipse_s) ;

filename = 'ellipse_r' ;
write_var(FolderPath_Config, filename, ellipse_r) ;

filename = 'AllPoints' ;
write_var(FolderPath_Config, filename, AllPoints) ;

filename = 'AngT2R_rf' ;
write_var(FolderPath_Config, filename, AngT2R_rf) ;

filename = 'AngS2R_rf' ;
write_var(FolderPath_Config, filename, AngS2R_rf) ;

filename = 'AngT2S_sf' ;
write_var(FolderPath_Config, filename, AngT2S_sf) ;


    
end


%% Function that calculates magnitude of the given vector

function abs = vectorMag(vec)

[m,n]=size(vec);
if (m~=1)&&(n~=1)  % or unit colomn or unit row
    abs = 0;
    disp('Error - vector is not of proper dimensions');
else
    abs = sqrt(sum (vec.^2));
end;

end
