% Mehmet Kurum
% Feb 22, 2017

function setGround

% Get Global Parameters
f_MHz = SatParams.getInstance.f_MHz;
VSM_cm3cm3 = GndParams.getInstance.VSM_cm3cm3( ParamsManager.index_VSM );
RMSH_cm = GndParams.getInstance.RMSH_cm( ParamsManager.index_RMSH );
sand_ratio = GndParams.getInstance.sand_ratio;
clay_ratio = GndParams.getInstance.clay_ratio;
rhob_gcm3 = GndParams.getInstance.rhob_gcm3;


%% Surface Roughness
f_Hz = f_MHz * Constants.MHz2Hz ;
lambda_cm = Constants.c / f_Hz * Constants.m2cm ;        % in cm
ko = 2 * pi / lambda_cm ;
h = (2 * RMSH_cm * ko) ^ 2 ;        % effective roughness parameter

%% Soil Dielectric Constant
eps_g = dielg( VSM_cm3cm3, f_Hz, sand_ratio, clay_ratio, rhob_gcm3) ; % eps_g = eps_gp - j * eps_gpp
eps_g = conj(eps_g) ; % eps_g = eps_gp + i * eps_gpp
eps_g = round(eps_g * 10) / 10 ;

%% Ground Parameters
grnd_par = [h, real(eps_g), imag(eps_g)] ;

%% Saving...

filename = 'G' ;
writeVar(SimulationFolders.getInstance.gnd, filename, grnd_par) ;

end


function epsoil = dielg(vsm, f_Hz, sand_ratio, clay_ratio, bulkrho)


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
%       real fghz,f_Hz,pi,vsm,alfa,beta,rhoss,sigma,fac1,fac2,fac3

% c      real bulkrho	!changed bulkrho defined already
% c

% TO-Do: Dr. Kurum: Check if sand and clay ratio OR percentage?
fghz = f_Hz / 1.0e+9 ;
alfa = 0.65 ;
beta = 1.09 - 0.11 * sand_ratio + 0.18 * clay_ratio ;     %! from Ulaby's book
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