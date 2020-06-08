function [depth_ind] = pendep(Gamma, n, L, lambda, theta, pol)
% function pendep
%
%   Calculates the penetration depth Pt/Po = e^-1 from a multilayer 
%   dielectric structure. 
%   Created by Mehmet Kurum - June 2, 2019
%
%   See also specularTerm, reflectionCoeffsSingle

%   Copyright © 2017-2020 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.3

% Initialize elementary reflection coefficients and optical lengths
[L, r, c] = optLenAndElemReflCoef(n, L, lambda, theta, pol);

% Number of layers
M = size(n, 2) - 2; 

% convert from degrees to radians
theta = theta * pi / 180;

%% begin forward layer
TT = zeros(M + 1, 1) ; % transmissivity at the right side of the boundary
for ii = 1 : M % backward layer recursion
    
    Matching = 1 / (1 + r(ii)) * [1 r(ii); r(ii) 1] ;
    if ii == 1 % first boundary
        Prog = [1 0; 0 1] ;
        MP = [1 0; 0 1] ;
    else
        delta = 2 * pi * L(ii) ./ lambda ;  % phase thickness in i-th layer
        Prog = [exp(1i * delta) 0; 0 exp(-1i * delta)] ;
    end
       
    MP = MP * Prog * Matching ;  % cascading layers
    invMP = inv(MP) ;    
    temp = invMP(1, 1) + invMP(1, 2) * Gamma ;
    if strcmp(pol, 'te')
        cc = real(c(ii)) / cos(theta) ;
    else
        cc =  cos(theta) / real(c(ii)) ;        
    end
    nn = real(n(1, ii+1)) / real(n(1, 1)) ;
    TT(ii) = abs(temp) ^ 2 *  cc * nn ; % transmissivity
    
    if TT(ii) < exp(-1)  % penetration depth
        depth_ind = ii ;
        return
    end
    
end



end
