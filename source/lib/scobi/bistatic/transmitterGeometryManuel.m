
function [Tgt, TgtI] = transmitterGeometryManuel
% function transmitterGeometryManuel 
%
%   Calculates the transmitter geometry-related  transformation matrices.  
%
%   [Tgt, TgtI] = transmitterGeometryManuel
%
%   See also updateBistaticDynParams, TxParams.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Configuration Parameters
th0_Tx_list_deg = ConfigParams.getInstance.th0_Tx_list_deg;
th0_Tx_deg = th0_Tx_list_deg( sim_counter );
th0_Tx_rad = deg2rad(th0_Tx_deg) ;
ph0_Tx_list_deg = ConfigParams.getInstance.ph0_Tx_list_deg;
ph0_Tx_deg = ph0_Tx_list_deg( sim_counter );

%% Incoming Signal
% Azimuth Angle (in standard spherical coords) of incoming signal
ph0in_Tx_deg = 90 - ph0_Tx_deg + 180 ;  % 180 degrees added to face transmitter and receiver aligned to each other.
ph0in_Tx_rad = deg2rad(ph0in_Tx_deg) ;


%% TRANSMITTER ANTENNA ROTATION MATRICES
% Rotation about z-axis (Azimuth rotation)
AntRotZ_Tx = [ cos(ph0in_Tx_rad), -sin(ph0in_Tx_rad), 0 ;
               sin(ph0in_Tx_rad),  cos(ph0in_Tx_rad), 0 ;
               0,                  0,                 1 ];

% Rotation about y-axis (Elevation rotation)
AntRotY_Tx = [ -cos(th0_Tx_rad), 0,  sin(th0_Tx_rad);
                0,               1,  0 ;
               -sin(th0_Tx_rad), 0, -cos(th0_Tx_rad) ];

% Alternative calculation if ph0in_Tx_deg = 90 - TxParams.ph0_Tx_deg
% Rotation about y-axis (Elevation rotation)
% AntRotY_Tx = [ -cos(th0_Tx_rad), 0, -sin(th0_Tx_rad) ;
%                 0,               1,  0 ;
%                 sin(th0_Tx_rad), 0, -cos(th0_Tx_rad) ];

% Rotation in both azimuth and elevation
AntRot_Tx = AntRotZ_Tx  * AntRotY_Tx ;


%% Ground Reference Coordinate System (East-North-Up)
% located on the ground where the receiver is projected
ux = [1 0 0] ;
uy = [0 1 0] ;
uz = [0 0 1] ;


%% Transmitter Antenna Coordinate System
% transfroming antenna to reference
uxt = (AntRot_Tx * ux')' ;
uyt = (AntRot_Tx * uy')' ;
uzt = (AntRot_Tx * uz')' ;


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


end