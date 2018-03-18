% Mehmet Kurum
% Feb 22, 2017

function setGround

% Get Global Parameters
fMHz = SatParams.getInstance.fMHz;
VSM = GndParams.getInstance.VSM( ParamsManager.index_VSM );
RMSH = GndParams.getInstance.RMSH( ParamsManager.index_RMSH );
sand = GndParams.getInstance.sand;
clay = GndParams.getInstance.clay;
rho_b = GndParams.getInstance.rho_b;


%% Surface Roughness
fhz = fMHz * 1e6 ;
lambda = Constants.c / fhz * 1e2 ;        % in cm
ko = 2 * pi / lambda ;
h = (2 * RMSH * ko) ^ 2 ;        % effective roughness parameter

%% Soil Dielectric Constant
eps_g = dielg( VSM, fhz, sand, clay, rho_b) ; % eps_g = eps_gp - j * eps_gpp
eps_g = conj(eps_g) ; % eps_g = eps_gp + i * eps_gpp
eps_g = round(eps_g * 10) / 10 ;

%% Ground Parameters
grnd_par = [h, real(eps_g), imag(eps_g)] ;

%% Saving...

filename = 'G' ;
writeVar(SimulationFolders.getInstance.gnd, filename, grnd_par) ;

end


function epsoil = dielg(vsm, fhz, sand, clay, bulkrho)


% c
% c     vsm is fractional volumetric soil moisture
% c     computes dielectric constant of the soil
% c      (see ulaby's book, vol iii, page 2103)
% c Note: epsoil is chosen in stead of epsg because epsg is used in ground.h
% c
%       implicit none
%       include 'ground.h'
%       include 'flag.h'
% c

%       complex *8 epss,epfw,epsoil,tmp,ei,term1,term2,epalfa
%       real fghz,fhz,pi,vsm,alfa,beta,rhoss,sigma,fac1,fac2,fac3

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
fac3 = (rhoss - bulkrho) ./ (rhoss * vsm) ;
epfw = 4.9 + 74.1 /tmp + ei * fac2 * fac3 ;

term2 = (vsm .^ beta) .* (epfw .^ alfa) - vsm ;
epalfa = 1.0 + term1 + term2 ;
epsoil = epalfa .^ (1.0 / alfa) ;
epsoil = conj(epsoil) ;


end