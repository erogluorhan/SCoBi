
function F = eDisc(TIN, PIN, TS, PS, TH, PH, FHZ, T, A, B, EPSDC)
% function eDisc 
%
%   The bistatic scattering amplitude from a lossy dielectric disc is
%   calculated. The phase center is at the bottom of the disc. Both
%   thin and thick routines are used. 
%
%   F = eDisc(TIN, PIN, TS, PS, TH, PH, FHZ, T, A, B, EPSDC)
%
%   INPUTS:
%   TIN,PIN = INCIDENT ANGLES (RAD)
%   TS,PS = SCATTERED ANGLES (RAD)
%   TH,PH = ROTATION ANGLES (RAD)
%   FHZ = FRQUENCY (HZ)
%   RAD = RADIUS OF CYLINDER (M)
%   L = LENGTH OF CYLINDER (M)
%   EPS = RELATIVE DIELECTRIC CONSTANT
%   F = BISTATIC SCATTERING AMPLITUDES
%   F(1) = FHH, F(2) = FVH, F(3) = FHV, F(4) = FVV
%
%   See also calcPropagation, eCylinder.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



F1 = rotDisc(TIN, PIN, TS, PS, TH, PH, FHZ, T, A, B, EPSDC) ;

F = [F1(1, 1) F1(2, 1); F1(3, 1) F1(4, 1)] ;
end


function FF = rotDisc(TIN1, PIN1, TS1, PS1, TH1, PH1, FHZ, T, A, B, EPSDC)

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%  Calculates the scattering amplitudes with respect to laboratory
%  coordinates by transforming from prime coordinates.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% Preallocation

FF = zeros(4, 1) ;

%%

DELTH = abs(TIN1 - TH1) ;

if DELTH < 0.001 % Avoids end-on scattering
    TH2 = TH1 + 0.01 ;
else
    TH2 = TH1 ;
end


%%         Define angles and trigonometric quantities

PSB = PH1 - PS1 ;

PIB = PH1 - PIN1 ;
STI = sin(TIN1) ;
CTH = cos(TH2) ;
CTI = cos(TIN1) ;
STH = sin(TH2) ;
STS = sin(TS1) ;
CTS = cos(TS1) ;
SPIB = sin(PIB) ;
SPSB = sin(PSB) ;
CPIB = cos(PIB) ;
CPSB = cos(PSB) ;

XI = STI * CTH * CPIB - CTI * STH ;
XS = STS * CTH * CPSB - CTS * STH ;
YI = -STI * SPIB ;
YS = -STS * SPSB ;
ZI = STI * STH * CPIB + CTI * CTH ;
ZS = STS * STH * CPSB + CTS * CTH ;

%%         Relation to prime angles

if TH1 == 0 && PH1 == 0
    
    PIP = PIN1 ;
    PSP = PS1 ;
    TIP = TIN1 ;
    TSP = TS1 ;
    
else
    
    TIP = acos(ZI) ;
    TSP = acos(ZS) ;
    PIP = atan2(YI, XI) ;
    
    if XS == 0 && YS == 0
        
        PSP = atan(-SPSB / CPSB) ;
        
    else
        
        PSP = atan2(YS, XS) ;
        
    end
    
end


if TIP >= 0.0 && TIP < 0.01; TIP = 0.01 ;   end
if TIP > 3.1316 && TIP <= pi; TIP = 3.1316 ; end
if TSP >= 0.0 && TSP < 0.01; TSP = 0.01 ;      end
if TSP > 3.1316 && TSP <= pi; TSP = 3.1316 ; end

%% thin_eDisc
FP = thin_eDisc(TIP, PIP, TSP, PSP, FHZ, T, A, B, EPSDC) ;

%%       Calculate dot products

SPSP = sin(PSP) ;
CPSP = cos(PSP) ;
STSP = sin(TSP) ;
CTSP = cos(TSP) ;
STIP = sin(TIP) ;
CTIP = cos(TIP) ;
SPIP = sin(PIP) ;
CPIP = cos(PIP) ;
%
HIXP = CTH * SPIB ;
HIYP = CPIB ;
HIZP = STH * SPIB ;
VIXP = -CTI * CTH * CPIB - STI * STH ;
VIYP = CTI * SPIB ;
VIZP = -CTI * STH * CPIB + STI * CTH ;
HSXP = -CTH * SPSB ;
HSYP = -CPSB ;
HSZP = -STH * SPSB ;
VSXP = -CTS * CTH * CPSB - STS * STH  ;
VSYP = CTS * SPSB ;
VSZP = -CTS * STH * CPSB + STS * CTH ;
%
HHS = SPSP * HSXP - CPSP * HSYP ;
HVS = -CTSP * CPSP * HSXP - CTSP * SPSP * HSYP + STSP * HSZP ;
VHS = SPSP * VSXP - CPSP * VSYP ;
VVS = -CTSP * CPSP * VSXP - CTSP * SPSP * VSYP + STSP * VSZP ;
%
HHI = -SPIP * HIXP + CPIP * HIYP ;
VHI = -SPIP * VIXP + CPIP * VIYP ;
HVI = -CTIP * CPIP * HIXP - CTIP * SPIP * HIYP + STIP * HIZP ;
VVI = -CTIP * CPIP * VIXP - CTIP * SPIP * VIYP + STIP * VIZP ;

%%        Transformation
%         Fhh=F(1)      Fhv=F(2)     Fvh=F(3)     Fvv=F(4)

FF(1, 1) = FP(1, 1) * HHS * HHI + FP(2, 1) * HHS * HVI ...
    + FP(3, 1) * HVS * HHI + FP(4, 1) * HVS * HVI ;

FF(2, 1) = FP(1, 1) * HHS * VHI + FP(2, 1) * HHS * VVI ...
    + FP(3, 1) * HVS * VHI + FP(4, 1) * HVS * VVI ;

FF(3, 1) = FP(1, 1) * VHS * HHI + FP(2, 1) * VHS * HVI ...
    + FP(3, 1) * VVS * HHI + FP(4, 1) * VVS * HVI ;

FF(4, 1) = FP(1, 1) * VHS * VHI + FP(2, 1) * VHS * VVI ...
    + FP(3, 1) * VVS * VHI + FP(4, 1) * VVS * VVI ;


%
end

function FP1 = thin_eDisc(TIP1, PIP1, TSP1, PSP1, FHZ, T, A, B, EPSDC)

% ***********************************************************************
% Calculates scattering amplitudes of a thin elliptic dielectric disk
%  with respect to prime coordinates. The major radius is along the x axis.
% ************************************************************************

%% Global Variables

C = 3e8 ;
AK0 = 2 * pi * FHZ / C ;

%% Preallocation

FP1 = zeros(4, 1) ;

%% Define constants

DCNST1 = (AK0 ^ 2) * (EPSDC - 1.0) / (4.0 * pi) ;
VOL = pi * (A * B) * T ;
D3 = (EPSDC - 1.0) / EPSDC ;

%%   DEFINE ANGLE

PHPD = PIP1 - PSP1 ;

%%   Define trigonometric quantities

CTIP = cos(TIP1) ;
CTSP = cos(TSP1) ;
CPIP = cos(PIP1) ;
CPSP = cos(PSP1) ;
CPHPD = cos(PHPD) ;
STIP = sin(TIP1) ;
STSP = sin(TSP1) ;
SPIP = sin(PIP1) ;
SPSP = sin(PSP1) ;
SPHPD = sin(PHPD) ;

%%   Calculate dot products

HSHI = -CPHPD ;
HSVI = CTIP * SPHPD ;
VSHI = CTSP * SPHPD ;
VSVI = CPHPD * CTIP * CTSP + STIP * STSP ;

%%   Calculate the bessel function which results from the integration
%      (making use of BESC1 subroutine)

ALPHX = -STIP * CPIP - STSP * CPSP ;
ALPHY = -STIP * SPIP - STSP * SPSP ;
A1 = (ALPHX * A) ^ 2 + (ALPHY * B) ^ 2 ;
ARG = AK0 * sqrt(A1) ;
% BES = BESC1(ARG) ;
if ARG == 0
    ARG = 1e-5 ;
    BES = besselj(1, ARG) ./ ARG ;
else
    BES = besselj(1, ARG) ./ ARG ;
end
DCNST = DCNST1 * VOL * 2.0 * BES ;

%%    CALCULATE THE SCATTERING AMPLITUDES
%
FP1(1, 1) = DCNST * HSHI ;
FP1(2, 1) = DCNST * HSVI ;
FP1(3, 1) = DCNST * VSHI ;
FP1(4, 1) = DCNST * (VSVI - D3 * STIP * STSP) ;

%
end