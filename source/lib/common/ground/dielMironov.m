%
%  dielMironov  Independent implementation of Mironov soil dielectric model --
%           straight from his 2009 TGRS paper.
%
%           Mironov's MSBDM is based on GRMDM, which is in turn based on
%           RMDM.  The complex dielectric constant of a soil/water mixture
%           is given in Eqns 11-13 in the paper.
%
%  USAGE
%
%           diel = dielMironov(f_Hz, VSM, clay_ratio)
%
%           f_Hz:  Frequency in Hertz
%           VSM:    Soil moisture in cm3/cm3
%           clay_ratio:    Clay fraction (0.0-1.0)
%
%  Adapted from Steven Chan, 03/2011
%
function diel = dielMironov(f_Hz, VSM, clay_ratio)
                                                                              % Section IV
C = clay_ratio * 100 ;                                                                                  % Section VI
diel = zeros(size(VSM)) ;

% Mironov's regression expressions based on Curtis, Dobson, and Hallikainen datasets
nd = 1.634 - 0.539e-2 * C + 0.2748e-4 * C .^ 2 ;                                                % Eqn 17
kd = 0.03952 - 0.04038e-2 * C ;                                                               % Eqn 18
mvt = 0.02863 + 0.30673e-2 * C ;                                                              % Eqn 19
eps0b = 79.8 - 85.4e-2 * C + 32.7e-4 * C .^ 2 ;                                                 % Eqn 20
taub = 1.062e-11 + 3.450e-12 * 1e-2 * C ;                                                     % Eqn 21
sigb = 0.3112 + 0.467e-2 * C ;                                                                % Eqn 22
sigu = 0.3631 + 1.217e-2 * C ;                                                                % Eqn 23
eps0u = 100 ;                                                                                 % Eqn 24
tauu = 8.5e-12 ;                                                                              % Eqn 25

% Debye relaxation equations for water as a function of frequency
eps0 = 8.854e-12 ;                                                                            % Vacuum permittivity
epsinf = 4.9 ;                                                                                % Section IV
epsb_real = epsinf + ( (eps0b - epsinf) ./ (1 + (2 * pi * f_Hz * taub) .^ 2) ) ;                           % Eqn 16
epsb_imag = (eps0b - epsinf) ./ (1 + (2 * pi * f_Hz * taub) .^ 2) ...
    .* (2 * pi * f_Hz * taub) + sigb ./ (2 * pi * eps0 * f_Hz) ; % Eqn 16
epsu_real = epsinf + ( (eps0u - epsinf) ./ (1 + (2 * pi * f_Hz * tauu) .^ 2) );                           % Eqn 16
epsu_imag = (eps0u - epsinf) ./ (1 + (2 * pi * f_Hz * tauu) .^2 ) ...
    .* (2 * pi * f_Hz * tauu) + sigu ./ (2 * pi * eps0 * f_Hz) ; % Eqn 16

% Refractive indices
nb = 1/sqrt(2) * sqrt( sqrt(epsb_real.^2 + epsb_imag.^2) + epsb_real );                      % Eqn 14
kb = 1/sqrt(2) * sqrt( sqrt(epsb_real.^2 + epsb_imag.^2) - epsb_real );                      % Eqn 14
nu = 1/sqrt(2) * sqrt( sqrt(epsu_real.^2 + epsu_imag.^2) + epsu_real );                      % Eqn 14
ku = 1/sqrt(2) * sqrt( sqrt(epsu_real.^2 + epsu_imag.^2) - epsu_real );                      % Eqn 14

% n(*) are refractive indices, k(*) are normalized attenuation coefficients
%   m: moist soil
%   d: dry soil
%   b: bound soil water (BSW)
%   u: unbound (free) soil water (FSW)
idx = VSM <= mvt;
  nm = nd + (nb - 1) .* VSM;                                                                  % Eqn 12
  km = kd + kb .* VSM;                                                                        % Eqn 13
  er_r_real = nm.^2 - km.^2;                                                                 % Eqn 11
  er_r_imag = 2 * nm .* km;                                                                  % Eqn 11
  tmp = er_r_real + 1i * er_r_imag;
diel(idx) = tmp(idx);

idx = VSM > mvt;
  nm = nd + (nb - 1) .* mvt + (nu - 1) .* (VSM - mvt);                                        % Eqn 12
  km = kd + kb .* mvt + ku .* (VSM - mvt);                                                    % Eqn 13
  er_r_real = nm.^2 - km.^2;                                                                 % Eqn 11
  er_r_imag = 2 * nm .* km;                                                                  % Eqn 11
  tmp = er_r_real + 1i * er_r_imag;
diel(idx) = tmp(idx);

end
