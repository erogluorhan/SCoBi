
function UU = calcMuller(uu)
% function calcMuller 
%
%   Calculates and outputs the 4 x 4 Mueller matrix of a 2 x 2 input matrix.  
%
%   UU = calcMuller(uu)
%
%   INPUTS:
%   uu: 2 x 2 double matrix

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


% Get the elements of the 2 x 2 input matrix
u11 = uu(1, 1) ; u12 = uu(1, 2) ;
u21 = uu(2, 1) ; u22 = uu(2, 2) ;

% Calculate the elements of 4 x 4 Mueller matrix 
U11 = abs(u11) .^ 2 ;
U12 = abs(u12) .^ 2 ;
U13 = real(u11 .* conj(u12)) ;
U14 = -imag(u11 .* conj(u12)) ;

U21 = abs(u21) .^ 2 ;
U22 = abs(u22) .^ 2 ;
U23 = real(u21 .* conj(u22)) ;
U24 = -imag(u21 .* conj(u22)) ;

U31 = 2 * real(u11 .* conj(u21)) ; 
U32 = 2 * real(u12 .* conj(u22)) ; 
U33 = real(u11 .* conj(u22) + u12 .* conj(u21)) ;
U34 = -imag(u11 .* conj(u22) - u12 .* conj(u21)) ;

U41 = 2 * imag(u11 .* conj(u21)) ; 
U42 = 2 * imag(u12 .* conj(u22)) ; 
U43 = imag(u11 .* conj(u22) + u12 .* conj(u21)) ;
U44 = real(u11 .* conj(u22) - u12 .* conj(u21)) ;

% Combine elements into the output Mueller matrix
UU = [U11, U12, U13, U14; ...
    U21, U22, U23, U24; ...
    U31, U32, U33, U34; ...
    U41, U42, U43, U44] ;


end

