
function F = CYLINDER(TIN, PIN, TS, PS, TH, PH, FHZ, RAD, L, EPS)

%% ****************************************************************
%
%     The bistatic scattering amplitude from a lossy dielectric cylinder
%     is calculated. The phase center is at the bottom of the cylinder.
%     Both thin and thick routines are used.
%
%       TIN,PIN = INCIDENT ANGLES (RAD)
%       TS,PS = SCATTERED ANGLES (RAD)
%       TH,PH = ROTATION ANGLES (RAD)
%       FHZ = FRQUENCY (HZ)
%       RAD = RADIUS OF CYLINDER (M)
%       L = LENGTH OF CYLINDER (M)
%       EPS = RELATIVE DIELECTRIC CONSTANT
%       F = BISTATIC SCATTERING AMPLITUDES
%       F(1) = FHH, F(2) = FVH, F(3) = FHV, F(4) = FVV

%%
F = ROTATE2(TIN, PIN, TS, PS, TH, PH, RAD, L, FHZ, EPS) ;

return ;

end


function FF = ROTATE2(TIN, PIN, TS, PS, TH, PH, RAD, L, FHZ, EPS)

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%  Calculates the scattering amplitudes with respect to laboratory
%  coordinates by transforming from prime coordinates.
%      End of cylinder is at the origin
%  1= thick cylinder     4=thin cylinder
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%

%% Define constants
C = 2.997925E08 ;

%% Input data
AK0 = 2 * pi * FHZ / C ;
EPSR = EPS ;

%% Avoids end-on scattering
DELTH = abs(TIN - TH) ;

if DELTH < 0.001
    TH = TH + 0.01 ;
end

%% Initialize

FF = zeros(4, 1) ;
FP = zeros(4, 1) ;

%% Define angles and trigonometric quantities

PSB = PH - PS ;
PIB = PH - PIN ;

STI = sin(TIN) ;
CTH = cos(TH) ;
CTI = cos(TIN) ;
STH = sin(TH) ;
STS = sin(TS) ;
CTS = cos(TS) ;
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

%%  Relation to prime angles

if TH == 0.0 && PH == 0
    
    PIP = PIN ;
    PSP = PS ;
    TIP = TIN ;
    TSP = TS ;
    
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

if TIP >= 0.0 && TIP < 0.01; TIP = 0.01 ;  end
if TIP > 3.1316 && TIP <= pi; TIP = 3.1316 ; end
if TSP >= 0.0 && TSP < 0.01; TSP = 0.01 ;  end
if TSP > 3.1316 && TSP <= pi; TSP = 3.1316 ; end

SPSP = sin(PSP) ;
CPSP = cos(PSP) ;
STSP = sin(TSP) ;
CTSP = cos(TSP) ;
STIP = sin(TIP) ;
CTIP = cos(TIP) ;
SPIP = sin(PIP) ;
CPIP = cos(PIP) ;

% X0 = AK0 * RAD * STIP ;
EMSI = EPSR - CTIP ^ 2 ;
X1 = AK0 * RAD * sqrt(EMSI) ;
% X2R = AK0 * RAD * STSP ;
% X0ABS = abs(X0) ;
X1ABS = abs(X1) ;
% X2RABS = abs(X2R) ;

if X1ABS > 0.05
    FP  = DIECYN(TIP, PIP, TSP, PSP, RAD, L, FHZ, EPS) ;
elseif X1ABS <= 0.05
    FP = THINCYL(TIP, PIP, TSP, PSP, RAD, L, FHZ, EPS) ;
end

%%   Calculate dot products

HIXP = CTH * SPIB ;
HIYP = CPIB ;
HIZP = STH * SPIB ;
VIXP = -CTI * CTH * CPIB - STI * STH ;
VIYP = CTI * SPIB ;
VIZP = -CTI * STH * CPIB + STI * CTH ;
HSXP = -CTH * SPSB ;
HSYP = -CPSB ;
HSZP = -STH * SPSB ;
VSXP = -CTS * CTH * CPSB - STS * STH ;
VSYP = CTS * SPSB ;
VSZP = -CTS * STH * CPSB + STS * CTH ;

HHS = SPSP * HSXP - CPSP * HSYP ;
HVS = -CTSP * CPSP * HSXP - CTSP * SPSP * HSYP + STSP * HSZP ;
VHS = SPSP * VSXP - CPSP * VSYP ;
VVS = -CTSP * CPSP * VSXP - CTSP * SPSP * VSYP + STSP * VSZP ;

HHI = -SPIP * HIXP + CPIP * HIYP ;
VHI = -SPIP * VIXP + CPIP * VIYP ;
HVI = -CTIP * CPIP * HIXP - CTIP * SPIP * HIYP + STIP * HIZP ;
VVI = -CTIP * CPIP * VIXP - CTIP * SPIP * VIYP + STIP * VIZP ;

%% Transformation
%  Fhh = F(1)   Fhv = F(2)   Fvh = F(3)   Fvv = F(4)

FF(1, 1) = FP(1, 1) * HHS * HHI + FP(2, 1) * HHS * HVI ...
    + FP(3, 1) * HVS * HHI + FP(4, 1) * HVS * HVI ;


FF(2, 1) = FP(1, 1) * HHS * VHI + FP(2, 1) * HHS * VVI ...
    + FP(3, 1) * HVS * VHI + FP(4, 1) * HVS * VVI ;

FF(3, 1) = FP(1, 1) * VHS * HHI + FP(2, 1) * VHS * HVI ...
    + FP(3, 1) * VVS * HHI + FP(4, 1) * VVS * HVI ;

FF(4, 1) = FP(1, 1) * VHS * VHI + FP(2, 1) * VHS * VVI ...
    + FP(3, 1) * VVS * VHI + FP(4, 1) * VVS * VVI ;

return ;

end


function FPP = THINCYL(TIP, PIP, TSP, PSP, RAD, L, FHZ, EPS)

%% *****************************************************************
%
%      Calculates scattering amplitudes for thin cylinder wrt prime co.
%            End of cylinder at the origin
%
% *******************************************************************


%% Define constants
C = 2.997925E08 ;
EJ = 0.0 + 1i * 1.0 ;

%% Input data
AK0 = 2 * pi * FHZ / C ;
EPSR = EPS ;                                     

%%  Define constants
DCNST1 = AK0 ^ 2 * (EPSR - 1.0) / (4.0 * pi) ;
VOL = pi * RAD ^ 2 * L ;
D1 = (EPSR - 1.0) / (EPSR + 1.0) ;
D2 = 2.0 / (EPSR + 1.0) ;

%%  Initialize FPP
FPP = zeros(4, 1) ;

%%
PHPD = PIP - PSP ;

CTIP = cos(TIP) ;
STIP = sin(TIP) ;
CTSP = cos(TSP) ;
STSP = sin(TSP) ;
CPHPD = cos(PHPD) ;
SPHPD = sin(PHPD) ;

%%  Calculate dot products
HSHI = -CPHPD ;
HSVI = CTIP * SPHPD ;
VSHI = CTSP * SPHPD ;
VSVI = CPHPD * CTIP * CTSP + STIP * STSP ;

%%  Calculate the sinc function 
ALPH = -CTIP - CTSP ;
ARG = EJ * AK0 * ALPH * L ;
ABARG = abs(ARG) ;

if ABARG < 1.0E-04
    SCTPTN = 1.0 + ARG / 2.0 ;
else
    SCTPTN = (exp(ARG) - 1.0) / ARG ;
end

DCNST = DCNST1 * VOL * SCTPTN ;

%% Calculate scattering amplitudes
%
%   FPP(1)=FPPHH                FPP(2)=FPPHV
%   FPP(3)=FPPVH                FPP(4)=FPPVV
%
FPP(1, 1) = DCNST * D2 * HSHI ;
FPP(2, 1) = DCNST * D2 * HSVI ;
FPP(3, 1) = DCNST * D2 * VSHI ;
FPP(4, 1) = DCNST * (D2 * VSVI + D1 * STIP * STSP) ;

return ;

end



function FP = DIECYN(TIP, PIP, TSP, PSP, RAD, L, FHZ, EPS)
%%
%    CALCULATES SCATTERING AMPLITUDES W.R.T PRIME COORDINATES
%             End of cylinder is at the origin
%

%% Preallocation

FP = zeros(4, 1) ;
JN = zeros(25, 3) ;
DJN = zeros(25, 3) ;

HN = zeros(25, 1) ;
DHN = zeros(25, 1) ;

%%  Input data

C = 2.997925E08 ;
AK0 = 2 * pi * FHZ / C ;
EPSR = EPS ;
NMAX = 12 ;
X0 = AK0 * RAD * sin(TIP) ;
EMSI = EPSR - (cos(TIP)) ^ 2 ;
X1 = AK0 * RAD * sqrt(EMSI) ;
X2R = AK0 * RAD * sin(TSP) ;
X2 = X2R + 1i * 0.0 ;

JN(1 : 12, :) = [besselj(0 : 11, X0) ; besselj(0 : 11, X1) ; besselj(0 : 11, X2)] .' ;
DJN(1 : 12, :) = [besseljd(0 : 11, X0) ; besseljd(0 : 11, X1) ; besseljd(0 : 11, X2)] .' ;

HN(1 : 12, 1) = besselh(0 : 11, X0) .' ;
DHN(1 : 12, 1) = besselhd(0 : 11, X0) .' ;


for LJ = 1 : NMAX
    
    LJ1 = LJ - 1 ;
    
    [DH, DV] = DIEINF(TIP, LJ, RAD, FHZ, EPS, JN, DJN, HN, DHN) ;
    EIP = I6I7P(X1, X2, LJ, RAD) ;
    EI = INTBES(X1, X2, LJ, RAD, EIP, JN, DJN) ;
    
    I3 = EI(LJ, 1) ;
    I4 = EI(LJ, 2) ;
    I5 = EI(LJ, 3) ;
    I6 = EI(LJ, 4) ;
    I7 = EI(LJ, 5) ;
    
    FNP = SCATN(TIP, PIP, TSP, PSP, DH, DV, L, LJ1, FHZ, EPS, I3, I4, I5, I6, I7) ;
    
    if LJ == 1
        
        FPHH = FNP(1, 1) ;
        FPHV = FNP(2, 1) ;
        FPVH = FNP(3, 1) ;
        FPVV = FNP(4, 1) ;
        TESTHH = FPHH ;
        TESTHV = FPHV ;
        TESTVH = FPVH ;
        TESTVV = FPVV ;
        
    else
        
        FPHH = FPHH + FNP(1, 1) ;
        FPHV = FPHV + FNP(2, 1) ;
        FPVH = FPVH + FNP(3, 1) ;
        FPVV = FPVV + FNP(4, 1) ;
      
        DV(1, 1) = (-1) ^ LJ1 * DV(1, 1) ;
        DV(2, 1) = (-1) ^ LJ1 * DV(2, 1) ;
        DV(3, 1) = (-1) ^ LJ1 * DV(3, 1) ;
        DV(4, 1) = -(-1) ^ LJ1 * DV(4, 1) ;
        DV(5, 1)= -(-1) ^ LJ1 * DV(5, 1) ;
        DH(1, 1)= -(-1) ^ LJ1 * DH(1, 1) ;
        DH(2, 1)= -(-1) ^ LJ1 * DH(2, 1) ;
        DH(3, 1)= -(-1) ^ LJ1 * DH(3, 1) ;
        DH(4, 1)= (-1) ^ LJ1 * DH(4, 1) ;
        DH(5, 1)= (-1) ^ LJ1 * DH(5, 1) ;
        
        %         I3 = I3 ;
        I4 = -EI(LJ, 3) ;
        I5 = -EI(LJ, 2) ;
        I6 = -EI(LJ, 5) ;
        I7 = -EI(LJ, 4) ;
        LJ1 = -LJ1 ;
        
        FNP = SCATN(TIP, PIP, TSP, PSP, DH, DV, L, LJ1, FHZ, EPS, I3, I4, I5, I6, I7) ;
        
        FPHH = FPHH + FNP(1, 1) ;
        FPHV = FPHV + FNP(2, 1) ;
        FPVH = FPVH + FNP(3, 1) ;
        FPVV = FPVV + FNP(4, 1) ;
       
        % Test of convergence for hh pol'n
        CFHHR = real(FPHH) ;
        CFHHI = imag(FPHH) ;
        TEHHR = real(TESTHH) ;
        TEHHI = imag(TESTHH) ;
        
        if CFHHR == 0.0 || CFHHI == 0
            KONHH = 1 ;
        else
            DIFHHR = (TEHHR - CFHHR) / CFHHR ;
            DIFHHI = (TEHHI - CFHHI) / CFHHI ;
            
            if abs(DIFHHR) < 0.01 && abs(DIFHHI) < 0.01
                KONHH = 1 ;
            else
                KONHH = 0 ;
            end
            
        end
        
        % Test of convergence for hv pol'n
        CFHVR = real(FPHV) ;
        CFHVI = imag(FPHV) ;
        TEHVR = real(TESTHV) ;
        TEHVI = imag(TESTHV) ;
        
        if CFHVR == 0.0 || CFHVI  == 0.0
            KONHV = 1 ; 
        else
            DIFHVR = (TEHVR - CFHVR) / CFHVR ;
            DIFHVI = (TEHVI - CFHVI) / CFHVI ;
            
            if abs(DIFHVR) < 0.01 && abs(DIFHVI) < 0.01
                KONHV = 1  ;
            else
                KONHV = 0 ;
            end
            
        end
      
        % Test of convergence for vh pol'n
        CFVHR = real(FPVH) ;
        CFVHI = imag(FPVH) ;
        TEVHR = real(TESTVH) ;
        TEVHI = imag(TESTVH) ;
        
        if CFVHR == 0.0 || CFVHI  == 0.0
            KONVH = 1 ; 
        else
            DIFVHR = (TEVHR - CFVHR) / CFVHR ;
            DIFVHI = (TEVHI - CFVHI) / CFVHI ;
            
            if abs(DIFVHR)< 0.01 && abs(DIFVHI)< 0.01
                KONVH = 1 ;
            else
                KONVH = 0 ;
            end
            
        end
    
        % Test of convergence for vv pol'n
        CFVVR = real(FPVV) ;
        CFVVI = imag(FPVV) ;
        TEVVR = real(TESTVV) ;
        TEVVI = imag(TESTVV) ;
        
        if CFVVR == 0.0 || CFVVI  == 0.0
            KONVV = 1 ; 
        else
            DIFVVR = (TEVVR - CFVVR) / CFVVR ;
            DIFVVI = (TEVVI - CFVVI) / CFVVI ;
            
            if abs(DIFVVR) < 0.01 && abs(DIFVVI)< 0.01
                KONVV = 1 ;
            else
                KONVV = 0 ;
            end
            
        end
       
        KON = KONHH + KONHV + KONVH + KONVV ;
        
        if KON == 4
            break ; 
        end
        
        TESTHH = FPHH ;
        TESTHV = FPHV ;
        TESTVH = FPVH ;
        TESTVV = FPVV ;
        
    end
    
end

FP(1, 1) = FPHH ;
FP(2, 1) = FPHV ;
FP(3, 1) = FPVH ;
FP(4, 1) = FPVV ;

return ;

end


function [DH, DV] = DIEINF(TIP, LJ, RAD, FHZ, EPS, JN, DJN, HN, DHN)
%%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% Calculates the internal field coefficients for a vertical infinitely
% long dielectric cylinder.
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%% Preallocation
DH = zeros(5, 1) ;
DV = zeros(5, 1) ;

ANV = zeros(25, 1) ;
ANH = zeros(25, 1) ;
BNV = zeros(25, 1) ;
BNH = zeros(25, 1) ;
CNV = zeros(25, 1) ;
CNH = zeros(25, 1) ;
CCB = zeros(25, 1) ;

%% Define constants
C = 2.997925E08 ;
EJ = 0.0 + 1i * 1.0 ;

%% Input data
AK0 = 2 * pi * FHZ / C ;
EPSR = EPS ;

LJ1 = LJ - 1 ;
CIMP = 120 * pi ;

%%
STIP = sin(TIP) ;
CTIP = cos(TIP) ;
EMSI = EPSR - CTIP ^ 2 ;
EJN = (-EJ) ^ LJ1 ;
X0 = AK0 * RAD * STIP ;
S0 = 1 / STIP ;
S1 = EPSR / sqrt(EMSI) ;
R1 = 1 / sqrt(EMSI) ;
QN = -LJ1 * CTIP * (R1 ^ 2 - S0 ^ 2) / (AK0 * RAD) ;
QN2 = QN * QN ;
% QNB = QN * X0 ;
JD = JN(LJ, 1) * DJN(LJ, 2) ;
DJ = DJN(LJ, 1) * JN(LJ, 2) ;
HDJ = HN(LJ, 1) * DJN(LJ, 2) ;
DHJ = DHN(LJ, 1) * JN(LJ, 2) ;
VN = S1 * JD - S0 * DJ ;
PN = (R1 * HDJ - S0 * DHJ) * 1.0E-09 ;
SN = (S1 * HDJ - S0 * DHJ) * 1.0E-09 ;
MN = R1 * JD - S0 * DJ ;
MNN = MN * SN * 1.0E+09 ;
PNN = PN * SN ;
VP = VN * PN * 1.0E+09 ;
QHJ = QN * HN(LJ, 1) * JN(LJ, 2) * 1.0E-09 ;
QHJ2 = QHJ * QHJ ;
J22 = JN(LJ, 2) * JN(LJ, 2) ;
QJHJ = QN2 * JN(LJ, 1) * HN(LJ, 1) * J22 ;
DENOM = PNN - QHJ2 ;
SX0 = 2 * S0 / (pi * X0) ;

%%
%     Calculates coefficients excatly for argument X0 GT 0.08 ; uses
%     small argument approximations for Hankel function for X0 LT 0.08
%
%     Exact calculation :

CNV(LJ, 1) = -(VP - QJHJ) * 1.0E-05 / (DENOM * 1.0E+13) ;
CNH(LJ, 1) = -(MNN - QJHJ) * 1.0E-05 / (DENOM * 1.0E+13) ;
CCB(LJ, 1) = (SX0 * QN * J22) * 1.0E-05 / (DENOM * 1.0E+13) ;

ANV(LJ, 1) = EJN * STIP * (JN(LJ, 1) + CNV(LJ, 1) * HN(LJ, 1)) / JN(LJ, 2) ;
BNV(LJ, 1) = EJN * STIP * (HN(LJ, 1) * CCB(LJ, 1)) / (JN(LJ, 2) * CIMP) ;
ANH(LJ, 1) = EJN * STIP * (CCB(LJ, 1) * HN(LJ, 1)) / JN(LJ, 2) ;
BNH(LJ, 1) = -EJN * STIP * (JN(LJ, 1) + CNH(LJ, 1) * HN(LJ, 1)) / (JN(LJ, 2) * CIMP) ;

DV(1, 1) = ANV(LJ, 1) ;
DV(2, 1) = -EJ * R1 * CTIP * ANV(LJ, 1) ;
DV(3, 1) = -LJ1 * R1 ^ 2 * BNV(LJ, 1) * CIMP / AK0 ;
DV(4, 1) = LJ1 * R1 ^ 2 * CTIP * ANV(LJ, 1) / AK0 ;
DV(5, 1) = -EJ * R1 * BNV(LJ, 1) * CIMP ;

DH(1, 1) = ANH(LJ, 1) ;
DH(2, 1) = -EJ * R1 * CTIP * ANH(LJ, 1) ;
DH(3, 1) = -LJ1 * R1 ^ 2 * BNH(LJ, 1) * CIMP / AK0 ;
DH(4, 1) = LJ1 * R1 ^ 2 * CTIP * ANH(LJ, 1) / AK0 ;
DH(5, 1) = -EJ * R1 * BNH(LJ, 1) * CIMP ;


end


function EIP = I6I7P(XX1, XX2, LJ, RAD)

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% This subroutine evaluates I6 & I7 for N greater than zero
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%% Preallocation
EIP = zeros(25, 2) ;

%%
N = LJ - 1 ;
X1 = XX1 * XX1 / 4.0 ;
X2 = XX2 * XX2 / 4.0 ;
K = 0 ;
TERM1 =  SUM1(N, K, X2) ;
SI6 = TERM1 ;
B = 1.0 + 1i * 0.0 ;

for K = 1 : 100
    
    B = -B * X1 / K / (N + K) ;
    SNK6 = SUM1(N, K, X2) ;
    TERMK = B * SNK6 ;
    SI6 = SI6 + TERMK ;
    TEST = abs(TERMK / SI6) ;
    if TEST < 1.E-12,
        break;  
    end
    
end

I6 = SI6 * XX2 * RAD / 4.0 ;

if N == 0.0
    
    EIP(LJ, 1) = I6  ;
    if N == 0,
        I7 = -EIP(1, 1) ;
    end
    
    EIP(LJ, 2) = I7 ;
    
    return ;
    
else
    
    for J = 1 : N
        I6 = I6 * XX1 * XX2 / 4.0 / J / (J + 1) ;
    end
    
end

%  <<   I7  <<

BB = 1.0 + 1i * 0.0 ;
TERM2 =  SUM2(N, 0, X2) ;
SI7 = TERM2 ;

for KK = 1 : 100
    
    BB = -BB * X1 / KK / (N + KK) ;
    SNK7 = SUM2(N, KK, X2) ;
    TRMK = BB * SNK7 ;
    SI7 = SI7 + TRMK ;
    TES = abs(TRMK / SI7) ;
    
    if TES < 1.E-12,
        break; 
    end
    
end

%  74   *** I7* SUM HAS NOT CONVERGED *****
I7 = SI7 * XX1 * RAD / 4.0  ;

if N == 1
    
    EIP(LJ, 1) = I6  ;
    if N == 0,
        I7 = -EIP(1, 1) ;
    end
    
    EIP(LJ, 2) = I7 ;
    
    return ;
    
else
    
    for JJ = 2 : N
        I7 = I7 * XX1 * XX2 / 4.0 / JJ /(JJ - 1) ;
    end
    
end 

end


function SUM = SUM1(N, K, X2)

TERMR = 1.0 / (N + K + 1) ;
TERM = TERMR + 1i * 0.0 ;
SUM = TERM ;

for  L = 1 : 100
    
    TERM = -TERM * X2 * (N + K + L) / L /(N + K + L + 1) / (N + L + 1) ;
    SUM = SUM + TERM ;
    TEST = abs(TERM / SUM) ;
    if TEST < 1.E-12
        return ; 
    end
    
end

%  75    FORMAT(44X,'**********I6 (L)* SUM HAS NOT CONVERGED *****')

end


function SUM = SUM2(N, K, X2)

TERMR = 1.0 / (N + K) ;
TERM = TERMR + 1i * 0.0 ;
SUM = TERM ;

for L = 1 : 100
    
    TERM = -TERM * X2 * (N + K + L - 1) / L /(N + K + L) / (N + L - 1) ;
    SUM = SUM + TERM ;
    TEST = abs(TERM / SUM) ;
    if TEST < 1.E-12
        return ;
    end
    
end

%    76    FORMAT(44X,'**********I7 (L)* SUM HAS NOT CONVERGED *****')

end


function EI = INTBES(XX1, XX2, LJ, RAD, EIP, JN, DJN)

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% Integration of product of Bessel functions for N.GE.0.
%
%     FOR I4&I5 USES RECURRENCE FORMULAS
%     FOR I6&I7 USES SERIES EXPANSIONS
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% Preallocation
EI = zeros(20, 5) ;

%%
LAM1 = XX1 / RAD ;
LAM2 = XX2 / RAD ;

LJ1 = LJ - 1.0 ;
JP = LJ + 1.0 ;

if LJ >= 2
    
    I3 = RAD * (LAM1 * JN(LJ, 3) * DJN(LJ, 2) - LAM2 * JN(LJ, 2) * DJN(LJ, 3)) / (LAM2 ^ 2 - LAM1 ^ 2) ;
    I6 = EIP(LJ, 1) ;
    I7 = EIP(LJ, 2) ;
    I4 = (LJ1 * I6) / LAM1 - I3P(JP, RAD, JN, DJN, LAM1, LAM2) ;
    I5 = I3M(LJ1, RAD, JN, DJN, LAM1, LAM2) - (LJ1 * I7) / LAM1 ;
    
    EI(LJ, 1) = I3 ;
    EI(LJ, 2) = I4 ;
    EI(LJ, 3) = I5 ;
    EI(LJ, 4) = I6 ;
    EI(LJ, 5) = I7 ;
    
    return ;
end

I3 = RAD * (LAM1 * JN(1, 3) * DJN(1, 2) - LAM2 * JN(1, 2) * DJN(1, 3)) ...
    / (LAM2 ^ 2 - LAM1 ^ 2) ;

I6 = EIP(LJ, 1) ;
I7 = EIP(LJ, 2) ;
I4 = (LJ1 * I6) / LAM1 - I3P(JP, RAD, JN, DJN, LAM1, LAM2) ;
I3D = RAD * (LAM1 * JN(2, 3) * DJN(2, 2) - LAM2 * JN(2, 2) * DJN(2, 3)) / (LAM2 ^ 2 - LAM1 ^ 2) ;
I5 = +I3D - (LJ1 * I7) / LAM1 ;

EI(LJ, 1) = I3 ;
EI(LJ, 2) = I4 ;
EI(LJ, 3) = I5 ;
EI(LJ, 4) = I6 ;
EI(LJ, 5) = I7 ;


end


function I3Px = I3P(JP, RAD, JN, DJN, LAM1, LAM2)
%%

LJ = JP ;

I3Px = RAD * (LAM1 * JN(LJ, 3) * DJN(LJ, 2) - LAM2 * JN(LJ, 2) * DJN(LJ, 3)) ...
    / (LAM2 ^ 2 - LAM1 ^ 2) ;

end


function I3Mx = I3M(LJ1, RAD, JN, DJN, LAM1, LAM2)
%% 

LJ = LJ1 ;

I3Mx = RAD * (LAM1 * JN(LJ, 3) * DJN(LJ, 2) - LAM2 * JN(LJ, 2) * DJN(LJ, 3)) ...
    / (LAM2 ^ 2 - LAM1 ^ 2) ;

end


function FNP = SCATN(TIP, PIP, TSP, PSP, DH, DV, L, LJ1, FHZ, EPS, I3, I4, I5, I6, I7) 
%%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% Calculates scattering amplitudes for each N, N=0,1,..
%      End of cylinder is at the origin
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% Preallocation

FNP = zeros(4, 1) ;

%%  Define constants
C = 2.997925E08 ;
EJ = 0.0 + 1i * 1.0 ;

%%  Input data
AK0 = 2 * pi * FHZ / C ;
EPSR = EPS ;

% STIP = sin(TIP) ;
CTIP = cos(TIP) ;
% SPIP = sin(PIP) ;
% CPIP = cos(PIP) ;
SPSP = sin(PSP) ;
CPSP = cos(PSP) ;
STSP = sin(TSP) ;
CTSP = cos(TSP) ;

%%     Calculation of dot products
HSPXP = SPSP ;
HSPYP = -CPSP ;
HSPZP = 0.0 ;
VSPXP = -CTSP * CPSP ;
VSPYP = -CTSP * SPSP ;
VSPZP = STSP ;

HA0H = HSPZP * DH(1, 1) ;
HA0V = HSPZP * DV(1, 1) ;
VA0H = VSPZP * DH(1, 1) ;
VA0V = VSPZP * DV(1, 1) ;

HA1H = HSPXP * DH(2, 1) + HSPYP * DH(5, 1) ;
HA1V = HSPXP * DV(2, 1) + HSPYP * DV(5, 1) ;
VA1H = VSPXP * DH(2, 1) + VSPYP * DH(5, 1) ;
VA1V = VSPXP * DV(2, 1) + VSPYP * DV(5, 1) ;

HA2H = HSPXP * DH(3, 1) + HSPYP * DH(4, 1) ;
HA2V = HSPXP * DV(3, 1) + HSPYP * DV(4, 1) ;
VA2H = VSPXP * DH(3, 1) + VSPYP * DH(4, 1) ;
VA2V = VSPXP * DV(3, 1) + VSPYP * DV(4, 1) ;

HB1H = -HSPXP * DH(5, 1) + HSPYP * DH(2, 1) ;
HB1V = -HSPXP * DV(5, 1) + HSPYP * DV(2, 1) ;
VB1H = -VSPXP * DH(5, 1) + VSPYP * DH(2, 1) ;
VB1V = -VSPXP * DV(5, 1) + VSPYP * DV(2, 1) ;

HB2H = -HSPXP * DH(4, 1) + HSPYP * DH(3, 1) ;
HB2V = -HSPXP * DV(4, 1) + HSPYP * DV(3, 1) ;
VB2H = -VSPXP * DH(4, 1) + VSPYP * DH(3, 1) ;
VB2V = -VSPXP * DV(4, 1) + VSPYP * DV(3, 1) ;

A = AK0 ^ 2 * (EPSR - 1) / (4 * pi) ;
ALPHA = CTIP + CTSP ;
PRMT = AK0 * ALPHA * L ;
ARG = EJ * PRMT ;

if PRMT == 0.0
    I1 = L + 1i * 0.0 ;
else
    I1 = (1.0 - exp(-ARG)) * L / ARG ;  %I1 replaced by conj(I1) Lang 9/99
end

EJN = (-EJ) ^ LJ1 ;
COEFF = pi * A * I1 * EJN ;
PPP = PIP - PSP ;
EPPPM = exp(-EJ * LJ1 * PPP) ;
EPSP = exp(EJ * PSP) ;
EMPSP = exp(-EJ * PSP) ;

I3HA0H = I3 * HA0H * 2 ;
I46AHH = -EJ * EPSP * (I4 * HA1H + I6 * HA2H) ;
I57AHH = EJ * EMPSP * (I5 * HA1H + I7 * HA2H) ;
I46BHH = -EPSP * (I4 * HB1H + I6 * HB2H) ;
I57BHH = -EMPSP * (I5 * HB1H + I7 * HB2H) ;

FNP(1, 1) = COEFF * EPPPM * (I3HA0H + I46AHH + I57AHH + I46BHH + I57BHH) ;

I3HA0V = I3 * HA0V * 2 ;
I46AHV = -EJ * EPSP * (I4 * HA1V + I6 * HA2V) ;
I57AHV = EJ * EMPSP * (I5 * HA1V + I7 * HA2V) ;
I46BHV = -EPSP * (I4 * HB1V + I6 * HB2V) ;
I57BHV = -EMPSP * (I5 * HB1V + I7 * HB2V) ;

FNP(2, 1) = COEFF * EPPPM * (I3HA0V + I46AHV + I57AHV + I46BHV + I57BHV) ;

I3VA0H = I3 * VA0H * 2 ;
I46AVH = -EJ * EPSP * (I4 * VA1H + I6 * VA2H) ;
I57AVH = EJ * EMPSP * (I5 * VA1H + I7 * VA2H) ;
I46BVH = -EPSP * (I4 * VB1H + I6 * VB2H) ;
I57BVH = -EMPSP *(I5 * VB1H + I7 * VB2H) ;

FNP(2, 1) = COEFF * EPPPM * (I3VA0H + I46AVH + I57AVH + I46BVH + I57BVH) ;

I3VA0V = I3 * VA0V * 2 ;
I46AVV = -EJ * EPSP * (I4 * VA1V + I6 * VA2V) ;
I57AVV = EJ * EMPSP * (I5 * VA1V + I7 * VA2V) ;
I46BVV = -EPSP * (I4 * VB1V + I6 * VB2V) ;
I57BVV = -EMPSP * (I5 * VB1V + I7 * VB2V) ;

FNP(4, 1) = COEFF * EPPPM * (I3VA0V + I46AVV + I57AVV + I46BVV + I57BVV) ;

end


function hd = besselhd(n, z)

% Derivative of Henkel of first kind

hn = besselh(n, z) ;
hnp = besselh(n + 1, z) ;

hd = n .* hn / z - hnp ;

end

function jd = besseljd(n, z)

% Derivative of Bessel of first kind

jn = besselj(n, z) ;
jnp = besselj(n + 1, z) ;

jd = n .* jn / z - jnp ;

end