
function [RGHI, RGVI, RGHS, RGVS] = reflectionCoeffSingle(TIN, TS, EPSG, h)
% function reflectionCoeffSingle 
%
%   Calculates the equivalent reflection coeff. of the rough ground from 
%   the Fresnel reflection coeff. of an avg. flat surface. Ground is
%   considered as surace-only. Generates the following outputs:
%   
%   RGHI,RGVI : Fresnel ref. coeff. from gnd for incidence angle TIN
%   RGHS,RGVS : Fresnel ref. coeff. from gnd for scattering angle TS
%
%   [RGHI, RGVI, RGHS, RGVS] = reflectionCoeffSingle(TIN, TS, EPSG, h)
%
%   INPUTS:
%   TIN:    incidence angle
%   TS:     scattering angle
%   EPSG:   Dielectric constant
%   h:      Effective roughness parameter ((2 * RMSH_cm * k0) ^ 2)
%
%   See also specularTerm.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


STI = sin(TIN) ;
CTI = cos(TIN) ;
STS = sin(TS) ;
CTS = cos(TS) ;


ESQI = sqrt(EPSG - STI ^ 2) ;
ESQS = sqrt(EPSG - STS ^ 2) ;


RGHIF = (CTI - ESQI) ./ (CTI + ESQI) ;
% RGVIF = (EPSG * CTI - ESQI) ./ (EPSG * CTI + ESQI) ;
RGVIF = (-EPSG * CTI + ESQI) ./ (EPSG * CTI + ESQI) ;
RGHSF = (CTS - ESQS) ./ (CTS + ESQS) ;
% RGVSF = (EPSG * CTS - ESQS) ./ (EPSG * CTS + ESQS) ;
RGVSF = (-EPSG * CTS + ESQS) ./ (EPSG * CTS + ESQS) ;


QZSGMI2 = h * CTI ^ 2 / 2.0 ;
QZSGMS2 = h * CTS ^ 2 / 2.0 ;

% QZSGMI = 2. * AK0 * CTI * STDEV ;
% QZSGMS = 2. * AK0 * CTS * STDEV ;
% QZSGMI2 = QZSGMI ^ 2 / 2.0 ;
% QZSGMS2 = QZSGMS ^ 2 / 2.0 ;


RGHI = RGHIF * exp(-QZSGMI2) ;
RGVI = RGVIF * exp(-QZSGMI2) ;
RGHS = RGHSF * exp(-QZSGMS2) ;
RGVS = RGVSF * exp(-QZSGMS2) ;

end
