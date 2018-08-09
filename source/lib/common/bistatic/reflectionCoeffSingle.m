
function [RGHI, RGVI, RGHS, RGVS] = reflectionCoeffSingle(TIN, TS, EPSG, h)

% ************************************************************************
%
% Calculates the equivalent reflection coeff. of the rough ground from
% the Fresnel reflection coeff. of an avg. flat surface.
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% RGQIF : Fresnel reflection coeff. for a plane wave with TIN for Q poln
% RGPSF : Fresnel reflection coeff. for a plane wave with TS for P poln
% RGQI : Reflection coeff of rough surface with TIN for Q poln
% RGPS : Reflection coeff of rough surface with TS for P poln
%
% ***********************************************************************
%
%         Calculate reflection coefficients:
%         RGHI,RGVI : Fresnel ref. coeff. from gnd for incidence angle TIN
%         RGHS,RGVS : Fresnel ref. coeff. from gnd for incidence angle TS
%

STI = sin(TIN) ;
CTI = cos(TIN) ;
STS = sin(TS) ;
CTS = cos(TS) ;
%
ESQI = sqrt(EPSG - STI ^ 2) ;
ESQS = sqrt(EPSG - STS ^ 2) ;
%
RGHIF = (CTI - ESQI) ./ (CTI + ESQI) ;
% RGVIF = (EPSG * CTI - ESQI) ./ (EPSG * CTI + ESQI) ;
RGVIF = (-EPSG * CTI + ESQI) ./ (EPSG * CTI + ESQI) ; % - sign chagned 09/24/17
RGHSF = (CTS - ESQS) ./ (CTS + ESQS) ;
% RGVSF = (EPSG * CTS - ESQS) ./ (EPSG * CTS + ESQS) ;
RGVSF = (-EPSG * CTS + ESQS) ./ (EPSG * CTS + ESQS) ; % - sign chagned 09/24/17
%
QZSGMI2 = h * CTI ^ 2 / 2.0 ;
QZSGMS2 = h * CTS ^ 2 / 2.0 ;

% QZSGMI = 2. * AK0 * CTI * STDEV ;
% QZSGMS = 2. * AK0 * CTS * STDEV ;
% QZSGMI2 = QZSGMI ^ 2 / 2.0 ;
% QZSGMS2 = QZSGMS ^ 2 / 2.0 ;
%
RGHI = RGHIF * exp(-QZSGMI2) ;
RGVI = RGVIF * exp(-QZSGMI2) ;
RGHS = RGHSF * exp(-QZSGMS2) ;
RGVS = RGVSF * exp(-QZSGMS2) ;

end
