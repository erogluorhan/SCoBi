
function diel = dielWang(VSM, sand_ratio, clay_ratio, rhob_gcm3)
% function dielWang 
%
%   The 1980 Wang & Schmugge soil dielectric model.  This is a straight
%   adaption of CMEM's codes.
%
%   diel = dielWang(VSM, sand_ratio, clay_ratio, rhob_gcm3)
%
%   INPUTS:
%   VSM:        Volumetric soil moisture (cm3/cm3) 
%   sand_ratio: Mass fraction of sand content in soil
%   clay_ratio: Mass fraction of clay content in soil[0,1]
%   rhob_gcm3:  Soil bulk density (g cm-3)
%
%   See also updateGndDynParams, dielMironov, dielDobson.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd
%   Adapted from Steven Chan, 03/2011

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


p = 1 - rhob_gcm3 / 2.65 ;  % Porosity = 1 - bulk density / particle density.  See Wikipedia.
ei = 3.2 + 0.1 * 1i ;       % CMEM
er = 5.5 + 0.2 * 1i ;       % CMEM
ea = 1.0 + 0.0 * 1i ;       % Common sense
ew = 79.5 + 6.63 * 1i ;     % Wang & Schmugge's paper
wp = 0.06774 - 0.064 * sand_ratio + 0.478 * clay_ratio ;
gamma = -0.57 * wp + 0.481 ;
wt = 0.49 * wp + 0.165 ;


% Add conductivity loss
alpha = 100 * wp ;
alpha(alpha > 26) = 26 ;

% Imaginary part
LF = alpha .* VSM .^ 2 ;


% Real part
DC = zeros(size(VSM)) ;

idx = VSM <= wt ;
ex = ei + (ew - ei) * (VSM ./ wt) .* gamma ;
eps = VSM .* ex + (p - VSM) * ea + (1 - p) * er ;

DC(idx,1) = eps(idx) ;


idx = VSM > wt ;
ex = ei + (ew - ei) .* gamma ;
eps = wt .* ex + (VSM - wt) * ew + (p - VSM) * ea + (1 - p) * er ;

DC(idx,1) = eps(idx) ;


% Combine the dielectric constant (complex number)
diel = DC + 1i * LF ;

end
