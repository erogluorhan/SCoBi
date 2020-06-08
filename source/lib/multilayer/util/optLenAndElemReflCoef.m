%% optLenAndElemReflCoef.m - reflection response of isotropic or birefringent multilayer structure
% function optLenAndElemReflCoef
%
%   Creates the elementary reflection coefficients and optical lengths used
%   by the multilayer.m function to calculate either a multilayer reflection
%   coefficient or penetration depth. Additionally, the coefficient ci is
%   returned.
%
%   This code is originally taken from Orfanidis' multidiel.m function with
%   modificiation made for polarization-dependent optical lengths.
%
%   See also specularTerm, reflectionCoeffsSingle

%   Copyright © 2017-2020 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Notes from multidiel.m
%
%          na | n1 | n2 | ... | nM | nb
% left medium | L1 | L2 | ... | LM | right medium 
%   interface 1    2    3     M   M+1
%
% Inputs
%
% n      = isotropic 1x(M+2), uniaxial 2x(M+2), or biaxial 3x(M+2), matrix of refractive indices
% L      = vector of optical lengths of layers, in units of lambda_0
% lambda = vector of free-space wavelengths at which to evaluate the reflection response
% theta  = incidence angle from left medium (in degrees)
% pol    = for 'tm' or 'te', parallel or perpendicular, p or s, polarizations
%
%
% notes: M is the number of layers (M >= 0)
%        n = [na, n1, n2, ..., nM, nb]        = 1x(M+2) row vector of isotropic indices
%
%            [ na1  n11  n12  ...  n1M  nb1 ]   3x(M+2) matrix of birefringent indices, 
%        n = [ na2  n21  n22  ...  n2M  nb2 ] = if 2x(M+2), it is extended to 3x(M+2)
%            [ na3  n31  n32  ...  n3M  nb3 ]   by repeating the top row
%
%        optical lengths are in units of a reference free-space wavelength lambda_0:
%        for i=1,2,...,M,  L(i) = n(1,i) * l(i), for TM, 
%                          L(i) = n(2,i) * l(i), for TE,
%        TM and TE L(i) are the same in isotropic case. If M=0, use L=[].
%
%        lambda is in units of lambda_0, that is, lambda/lambda_0 = f_0/f
%
%        reflectance = |Gamma|^2, input impedance = Z = (1+Gamma)./(1-Gamma)
%
%        delta(i) = 2*pi*[n(1,i) * l(i) * sqrt(1 - (Na*sin(theta))^2 ./ n(3,i).^2))]/lambda, for TM
%        delta(i) = 2*pi*[n(2,i) * l(i) * sqrt(1 - (Na*sin(theta))^2 ./ n(2,i).^2))]/lambda, for TE
%
%        if n(3,i)=n(3,i+1)=Na, then will get NaN's at theta=90 because of 0/0, (see also FRESNEL)
%
%        it uses SQRTE, which is a modified version of SQRT approriate for evanescent waves
%
%        see also MULTIDIEL1, MULTIDIEL2

% Sophocles J. Orfanidis - 1999-2008 - www.ece.rutgers.edu/~orfanidi/ewa
%%

function [optL, r, c] = optLenAndElemReflCoef(n, L, lambda, theta, pol)

if nargin == 0, help multidiel ; return ; end
if nargin <= 4, pol = 'te' ; end
if nargin == 3, theta = 0 ; end

if size(n, 2) == 1, n = n' ; end                          % in case n is entered as column 

K = size(n, 1) ;                                          % birefringence dimension
M = size(n, 2) - 2 ;                                      % number of layers

if K == 1, n = [n; n; n] ; end                            % isotropic case
if K == 2, n = [n(1, :); n] ; end                         % uniaxial case

if M == 0, L = []; end                                    % single interface, no slabs

theta = theta * pi / 180;

if strcmp(pol, 'te')
    Nsin2 = (n(2,1) * sin(theta)) ^ 2 ;                   % (Na*sin(tha))^2              
    %c = conj(sqrt(conj(1 - Nsin2 ./ n(2,:).^2)));      % old version
    c = sqrte(1 - Nsin2 ./ n(2, :) .^ 2) ;                % coefficient ci, or cos(th(i)) in isotropic case
    nT = n(2, :) .* c ;                                   % transverse refractive indices
    r = n2r(nT) ;                                         % r(i) = (nT(i-1)-nT(i)) / (nT(i-1)+nT(i))
else
    Nsin2 = (n(1, 1) * n(3, 1) * sin(theta)) ^ 2 / (n(3, 1) ^ 2 * cos(theta) ^ 2 + n(1, 1) ^ 2 * sin(theta) ^ 2) ;
    %c = conj(sqrt(conj(1 - Nsin2 ./ n(3,:).^2)));
    c = sqrte(1 - Nsin2 ./ n(3, :) .^ 2) ;
    nTinv = c ./ n(1, :) ;                                % nTinv(i) = 1/nT(i) to avoid NaNs
    r = -n2r(nTinv) ;                                     % minus sign because n2r(n) = -n2r(1./n)
end

if M > 0
    if strcmp(pol, 'te') % modifed 10/20/09 by MK
    L = L .* n(2, 2 : M + 1) .* c(2 : M + 1) ;                                  % polarization-dependent optical lengths
    else
    L = L .* n(1, 2 : M + 1) .* c(2 : M + 1);                                  % polarization-dependent optical lengths        
    end
end


optL = L;
end

%% n2r.m - refractive indices to reflection coefficients of M-layer structure
%
% Usage: r = n2r(n)
%
% n = refractive indices = [na,n(1),...,n(M),nb] 
% r = reflection coefficients = [r(1),...,r(M+1)]
%
% notes: there are M layers, M+1 interfaces, and M+2 media

% Sophocles J. Orfanidis - 1999-2008 - www.ece.rutgers.edu/~orfanidi/ewa
%%
function r = n2r(n)

if nargin == 0, help n2r ; return ; end

r = -diff(n) ./ (2 * n(1 : end - 1) + diff(n)) ;

end
