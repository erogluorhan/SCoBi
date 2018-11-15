
function updateRotMatDynParams
% function updateRotMatDynParams 
%
%   Calculates the antenna polarization rotation matrices and updates the
%   RotMatDynParams class with those values in each simulation iteration.  
%
%   See also mainSCoBi, RotMatDynParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


%% GET GLOBAL PARAMETERS
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Ground Parameters
polG = GndParams.getInstance.polG;
% Bistatic Parameters
idn = BistaticDynParams.getInstance.idn;  % propagation vector (i_d^-)
isn = BistaticDynParams.getInstance.isn;  % propagation vector (i_s^-)
osp = BistaticDynParams.getInstance.osp;  % propagation vector (o_s^+)
% Transformations
Tgs = BistaticDynParams.getInstance.Tgs;
Tgr = BistaticDynParams.getInstance.Tgr;
Tgt = BistaticDynParams.getInstance.Tgt;


%% CALCULATIONS
%% TRANSMITTER to RECEIVER
% u_t_r(i_d^-)
[ut1, ut2, ur1, ur2] = tanUnitVectors(Tgt, Tgr, idn, pol_Tx, pol_Rx ) ;

% Polarization Basis Dot Products
u11 = dot(ut1, conj(ur1)) ; u12 = dot(ut2, conj(ur1)) ;
u21 = dot(ut1, conj(ur2)) ; u22 = dot(ut2, conj(ur2)) ;

% 2 X 2
u_tr = [u11, u12; u21, u22] ;
% % % 4 X 4
% % U_tr = calcMuller(u_tr) ;


%% TRANSMITTER to SPECULAR POINT
% u_t_s(i_s^-)
[ut1, ut2, uvsi, uhsi] = tanUnitVectors(Tgt, Tgs, isn, pol_Tx, polG) ;

% Polarization Basis Dot Products
u11 = dot(ut1, conj(uvsi)) ; u12 = dot(ut2, conj(uvsi)) ;
u21 = dot(ut1, conj(uhsi)) ; u22 = dot(ut2, conj(uhsi)) ;

% 2 X 2
u_ts = [u11, u12; u21, u22] ;
% % % 4 X 4
% % U_ts = calcMuller(u_ts) ;


% 2 X 2
u_tIs = u_ts ; %% added april 29, 2017
% % % 4 X 4
% % U_tIs = U_ts ; %% added May 1, 2017


%% SPECULAR POINT to RECEIVER
% u_p_r(o_s^+)
[uvso, uhso, ur1, ur2] = tanUnitVectors(Tgs, Tgr, osp, polG, pol_Rx) ;

% Polarization Basis Dot Products
u11 = dot(uvso, conj(ur1)) ; u12 = dot(uhso, conj(ur1)) ;
u21 = dot(uvso, conj(ur2)) ; u22 = dot(uhso, conj(ur2)) ;

% 2 X 2
u_sr = [u11, u12; u21, u22] ;
% % % 4 X 4
% % U_sr = calcMuller(u_sr) ;


%% UPDATE ROTATION MATRICES DYNAMIC PARAMETERS
RotMatDynParams.getInstance.update( u_sr, u_tr, u_ts, u_tIs );


end



% Calculates unit vectors on the tangential plane
function [u1p1, u1p2, u2p1, u2p2] = tanUnitVectors(Tg1, Tg2, k, pol1, pol2)

% k is a unit vector (given in reference frame) from Frame 1 to Frame 2
% ux1 = Tg1(1, :) ; uy1 = Tg1(2, :) ; uz1 = Tg1(3, :) ;
% ux2 = Tg2(1, :) ; uy2 = Tg2(2, :) ; uz2 = Tg2(3, :) ;

%% Angles to frames # 1 and # 2
k_1f = Tg1 * k ; % unit propagation vector in frame #1
k_2f = Tg2 * k ; % unit propagation vector in frame #2

% off-axis angle of z1 towards frame 2
th1 = acos(k_1f(3)) ;
% orientation - azimuth
ph1 = atan2(k_1f(2), k_1f(1)) ;

% off-axis angle of z2 towards frame 1
th2 = acos(-k_2f(3)) ;
% orientation - azimuth
ph2 = atan2(-k_2f(2), -k_2f(1)) ;

%% unit vectors in phi direction [HPOL]
AR = 0 ; AT = 0 ; AP = 1 ;
[ux1p, uy1p, uz1p] = sph2Car2(th1, ph1, AR, AT, AP) ;
[ux2p, uy2p, uz2p] = sph2Car2(th2, ph2, AR, AT, AP) ;

% phi vector of frame 1
u1ph = Tg1' * [ux1p; uy1p; uz1p] ;  % in reference frame
% phi vector of frame 2
u2ph = Tg2' * [ux2p; uy2p; uz2p] ;  % in reference frame

%% unit vectors in theta direction [VPOL]
AR = 0 ; AT = 1 ; AP = 0 ;
[ux1t, uy1t, uz1t] = sph2Car2(th1, ph1, AR, AT, AP) ;
[ux2t, uy2t, uz2t] = sph2Car2(th2, ph2, AR, AT, AP) ;

% theta vector of frame 1
u1th = Tg1' * [ux1t; uy1t; uz1t] ;  % in reference frame
% theta vector of frame 2
u2th = Tg2' * [ux2t; uy2t; uz2t] ;  % in reference frame


%% unit vectors in Ludwig basis (reference: YPOL)
AR = 0 ; AT = sin(ph1) ; AP = cos(ph1) ;
[ux1Y, uy1Y, uz1Y] = sph2Car2(th1, ph1, AR, AT, AP) ;
AR = 0 ; AT = sin(ph2) ; AP = cos(ph2) ;
[ux2Y, uy2Y, uz2Y] = sph2Car2(th2, ph2, AR, AT, AP) ;

% reference vector (Y-axis) of frame 1
u1Y = Tg1' * [ux1Y; uy1Y; uz1Y] ;  % in reference frame
% reference vector (Y-axis) of frame 2
u2Y = Tg2' * [ux2Y; uy2Y; uz2Y] ;  % in reference frame

%% unit vectors in Ludwig basis (cross: XPOL)
AR = 0 ; AT = cos(ph1) ; AP = -sin(ph1) ;
[ux1X, uy1X, uz1X] = sph2Car2(th1, ph1, AR, AT, AP) ;
AR = 0 ; AT = cos(ph2) ; AP = -sin(ph2) ;
[ux2X, uy2X, uz2X] = sph2Car2(th2, ph2, AR, AT, AP) ;

% cross vector (X-axis) of frame 1
u1X = Tg1' * [ux1X; uy1X; uz1X] ;  % in reference frame
% cross vector (X-axis) of frame 2
u2X = Tg2' * [ux2X; uy2X; uz2X] ;  % in reference frame

%% unit vectors in circular polarization (RPOL and LPOL)

% frame 1
u1R = 1/ sqrt(2) * (u1X - 1i * u1Y) ;   % port 1    RHCP
u1L = 1/ sqrt(2) * (u1X + 1i * u1Y) ;   % port 2    LHCP

% frame 2
u2R = 1 / sqrt(2) * (u2X - 1i * u2Y) ;   % port 1    RHCP
u2L = 1 / sqrt(2) * (u2X + 1i * u2Y) ;   % port 2    LHCP


% Determine w.r.t the polarization
if pol1 == Constants.POLARIZATIONS{Constants.ID_POL_H}
    u1p1 = u1ph ;
    u1p2 = u1th ;
elseif pol1 == Constants.POLARIZATIONS{Constants.ID_POL_V}
    u1p1 = u1th ;
    u1p2 = u1ph ;
elseif pol1 == Constants.POLARIZATIONS{Constants.ID_POL_X}
    u1p1 = u1X ;
    u1p2 = u1Y ;
elseif pol1 == Constants.POLARIZATIONS{Constants.ID_POL_Y}
    u1p1 = u1Y ;
    u1p2 = u1X ;
elseif pol1 == Constants.POLARIZATIONS{Constants.ID_POL_R}
    u1p1 = u1R ;
    u1p2 = u1L ;
elseif pol1 == Constants.POLARIZATIONS{Constants.ID_POL_L}
    u1p1 = u1L ;
    u1p2 = u1R ;
    
end

% Determine w.r.t the polarization
if pol2 == Constants.POLARIZATIONS{Constants.ID_POL_H}
    u2p1 = u2ph ;
    u2p2 = u2th ;
elseif pol2 == Constants.POLARIZATIONS{Constants.ID_POL_V}
    u2p1 = u2th ;
    u2p2 = u2ph ;
elseif pol2 == Constants.POLARIZATIONS{Constants.ID_POL_X}
    u2p1 = u2X ;
    u2p2 = u2Y ;
elseif pol2 == Constants.POLARIZATIONS{Constants.ID_POL_Y}
    u2p1 = u2Y ;
    u2p2 = u2X ;
elseif pol2 == Constants.POLARIZATIONS{Constants.ID_POL_R}
    u2p1 = u2R ;
    u2p2 = u2L ;
elseif pol2 == Constants.POLARIZATIONS{Constants.ID_POL_L}
    u2p1 = u2L ;
    u2p2 = u2R ;
end

% % elseif pol2 == 'R' % Circular pol rotation has swaped to match trasnmit and receive directions
% %     u2p1 = u2L ;
% %     u2p2 = u2R ;
% % elseif pol2 == 'L'
% %     u2p1 = u2R ;
% %     u2p2 = u2L ;
% % end

end


% Conversion from Spherical to Cartesian coordinate system
function [AX, AY, AZ] = sph2Car2(theta, phi, AR, AT, AP)

st = sin(theta) ;
ct = cos(theta) ;
cp = cos(phi) ;
sp = sin(phi) ;

AX = AR .* st .* cp + AT .* ct .* cp - AP .* sp ;
AY = AR .* st .* sp + AT .* ct .* sp + AP .* cp ;
AZ = AR .* ct - AT .* st ;

end