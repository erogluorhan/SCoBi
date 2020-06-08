%% multidiel.m - reflection response of isotropic or birefringent multilayer structure
% function multidiel
%
%   Creates the multilayer reflection coefficient.
%
%   This code is originally taken from Orfanidis' multidiel.m function with
%   modificiation made for versatile use of elementary reflection 
%   coefficient calculation.
%
%   See also specularTerm, reflectionCoeffsSingle, optLenAndElemReflCoef

%   Copyright © 2017-2020 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          na | n1 | n2 | ... | nM | nb
% left medium | L1 | L2 | ... | LM | right medium 
%   interface 1    2    3     M   M+1
%
% Usage: [Gamma,Z] = multidiel(n,L,lambda,theta,pol)
%        [Gamma,Z] = multidiel(n,L,lambda,theta)       (equivalent to pol='te')
%        [Gamma,Z] = multidiel(n,L,lambda)             (equivalent to theta=0)
%
% n      = isotropic 1x(M+2), uniaxial 2x(M+2), or biaxial 3x(M+2), matrix of refractive indices
% L      = vector of optical lengths of layers, in units of lambda_0
% lambda = vector of free-space wavelengths at which to evaluate the reflection response
% theta  = incidence angle from left medium (in degrees)
% pol    = for 'tm' or 'te', parallel or perpendicular, p or s, polarizations
%
% Gamma = reflection response at interface-1 into left medium evaluated at lambda 
% Z     = transverse wave impedance at interface-1 in units of eta_a (left medium)
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
function [Gamma, Z] = multidiel(n, L, lambda, theta, pol)

[L, r, ~] = optLenAndElemReflCoef(n, L, lambda, theta, pol);

% Number of layers
M = size(n, 2) - 2; 

Gamma = r(M + 1) * ones(1, length(lambda)) ;             % initialize Gamma at right-most interface

for i = M : -1 : 1                                       % forward layer recursion 
    delta = 2 * pi * L(i) ./ lambda ;                    % phase thickness in i-th layer
    z = exp(-2 * 1i * delta) ;                          
    Gamma = (r(i) + Gamma .* z) ./ (1 + r(i) * Gamma .* z) ;
end

Z = (1 + Gamma) ./ (1 - Gamma) ;

end