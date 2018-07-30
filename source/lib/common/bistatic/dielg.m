
function epsoil = dielg(amv, fhz, sand, clay, bulkrho)


% c
% c     amv is fractional volumetric soil moisture
% c     computes dielectric constant of the soil
% c      (see ulaby's book, vol iii, page 2103)
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
epfw = 4.9 + 74.1 /tmp + ei * fac2 .* fac3 ;

term2 = (amv .^ beta) .* (epfw .^ alfa) - amv ;
epalfa = 1.0 + term1 + term2 ;
epsoil = epalfa .^ (1.0 / alfa) ;
epsoil = conj(epsoil) ;


end