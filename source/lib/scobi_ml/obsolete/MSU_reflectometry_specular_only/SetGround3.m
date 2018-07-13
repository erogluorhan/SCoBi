% Mehmet Kurum
% Feb 22, 2017

function SetGround3

%% Global
global fMHz tp
global FolderPath_gnd

%% Surface Roughness
fhz = fMHz * 1e6 ;
c = 3e8 ;
lambda = c / fhz * 1e2 ;        % in cm
ko = 2 * pi / lambda ;
sig = 1 ;                       % in cm
h = (2 * sig * ko) ^ 2 ;        % effective roughness parameter

%% Soil Texture
sand = 0.80; % percentage of sand
clay = 0.07; % percentage of clay
rho_b = 1.25 ; % Soil Bulk Density

%% Soil Dielectric Constant
eps_g = dielg(tp / 100, fhz, sand, clay, rho_b) ; % eps_g = eps_gp - j * eps_gpp
eps_g = conj(eps_g) ; % eps_g = eps_gp + i * eps_gpp
eps_g = round(eps_g * 10) / 10 ;

%% Ground Parameters
grnd_par = [h, real(eps_g), imag(eps_g)] ;

%% Saving...

filename = 'G' ;
write_var(FolderPath_gnd, filename, grnd_par) ;

end


function epsoil = dielg(amv, fhz, sand, clay, bulkrho)


% c
% c     amv is fractional volumetric soil moisture
% c     computes dielectric constant of the soil
% c      (see ulaby's book, vol iii, page 2103)
% c Note: epsoil is chosen in stead of epsg because epsg is used in ground.h
% c
%       implicit none
%       include 'ground.h'
%       include 'flag.h'
% c

%       complex *8 epss,epfw,epsoil,tmp,ei,term1,term2,epalfa
%       real fghz,fhz,pi,amv,alfa,beta,rhoss,sigma,fac1,fac2,fac3

% c      real bulkrho	!changed bulkrho defined already
% c
fghz = fhz / 1.0e+9 ;
alfa = 0.65 ;
beta = 1.09 - 0.11 * sand + 0.18 * clay ;     %! from Ulaby's book
% c     beta=1.1
% c
% c     beta varies from 1.0 for sandy soil to 1.17 to silty clay
% c           beta=1.09-0.11*s+0.18*c, where
% c     sand/100 and clay/100 are fractions of soil by weight.
% c
rhoss = 2.65 ;
ei = 1i ;
epss = 4.7 ;
term1 = bulkrho * (epss^alfa - 1.0) / rhoss ;
tmp = 1.0 - ei * fghz / 18.64 ;

% c                             epfw=4.9+74.1/tmp
sigma = -1.645 + 1.939 * bulkrho - 0.02013 * 50.0 + 0.01594 * 50.0 ;
fac1 = 2.0 * pi * fghz * 1.0e+9 * 8.854e-12 ;
fac2 = sigma / fac1 ;
fac3 = (rhoss - bulkrho) ./ (rhoss * amv) ;
epfw = 4.9 + 74.1 /tmp + ei * fac2 * fac3 ;

term2 = (amv .^ beta) .* (epfw .^ alfa) - amv ;
epalfa = 1.0 + term1 + term2 ;
epsoil = epalfa .^ (1.0 / alfa) ;
epsoil = conj(epsoil) ;


end