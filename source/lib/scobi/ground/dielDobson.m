
function diel = dielDobson(f_Hz, VSM, sand_ratio, clay_ratio, rhob_gcm3)
% function dielDobson 
%
%   Computes the complex dielectric constant of soil as a function of 
%   frequency, soil moisture, sand/clay fractions, and surface temperature.  
%   This particular implementation combines Ann Hsu's code with our 
%   existing Dobson model.  This approach starts with fixing RHOS at 2.66 
%   and computes for EPSS and RHOB. Low frequency correction is made for 
%   both the real and imaginary parts of the dielectric constant.  
%
%   diel = dielDobson(f_Hz, VSM, sand_ratio, clay_ratio, rhob_gcm3)
%
%   INPUTS:
%   f_Hz:       Frequency (Hertz)
%   VSM:        Volumetric soil moisture (cm3/cm3) [0,1]
%   rhob_gcm3:  Soil bulk density (g cm-3)
%   sand_ratio: Mass fraction of sand content in soil
%   clay_ratio: Mass fraction of clay content in soil
%
%   See also updateGndDynParams, dielMironov, dielWang.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd
%   Adapted from Steven Chan, 03/2011

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


% Cnvert freq. to GHz
f_GHz = f_Hz / Constants.GHZ_TO_HZ;

% Set physical constants and bounds
alpha   = 0.65;       % optimized coefficient
epso    = 8.854e-12;  % permittivity of free space (F/m)
epswinf = 4.9;        % high-frequency limit of free water dielectric constant
ts = 25 ;              % Soil temperature (deg C)

% Compute dielectric constant of soil solid particles
rhos = 2.66;
epss = (1.01 + 0.44 * rhos) ^ 2 - 0.062;
% por = 0.505 - 0.142*sand_ratio - 0.037*clay_ratio;
% fv = 1 - por;
% rhob_gcm3 = fv*rhos;

% Compute optimized coefficient values
betar  =  1.27480 - 0.519 * sand_ratio - 0.152 * clay_ratio ;
betai  =  1.33797 - 0.603 * sand_ratio - 0.166 * clay_ratio ;

% Compute empirical expression for effective conductivity
sigmae = -1.645 + 1.939 * rhob_gcm3 - 2.25622 * sand_ratio + 1.594 * clay_ratio ;

if (f_GHz < 1.4)
   sigmae = 0.0467 + 0.2204 * rhob_gcm3 - 0.4111 * sand_ratio + 0.6614 * clay_ratio ;
end

% Compute dielectric constant of free water (expressions for
% relaxation time and static dielectric constant obtained from
% Ulaby et al., Vol. 3)
omtau  =  f_GHz * (1.1109e-1 - 3.824e-3 * ts + 6.938e-5 * ts .^ 2 - 5.096e-7 * ts .^ 3) ;
epswo  =  88.045 - 0.4147 * ts + 6.295e-4 * ts .^ 2 + 1.075e-5 * ts .^ 3 ;
fac    = (epswo - epswinf) ./ (1.0 + omtau .^ 2) ;
epsfwr =  epswinf + fac ;
epsfwi = (omtau .* fac) + ( sigmae .* (rhos - rhob_gcm3) ./ ...\
                                      (2.0 * pi * f_GHz * 1.0e9 * epso * rhos * VSM) );

% Compute dielectric constant of soil
er_r   = ( 1.0 + (epss ^ alpha - 1.0) * rhob_gcm3 / rhos + ...\
           (VSM .^ betar) .* (epsfwr .^ alpha) - VSM ) .^ (1.0 / alpha) ;
er_i   = ( (VSM .^ betai) .* (epsfwi .^ alpha) ) .^ (1.0 / alpha) ;

% Correction factor (Peplinski et al., 1995):
if (f_GHz < 1.4)
    er_r = 1.15 * er_r - 0.68 ;
end

% Combine the dielectric constant (complex number)
diel = er_r + 1i * er_i ;

end
