
function [result] = convertDecibelToNatural( valDB )
% function convertDecibelToNatural 
%
%   Converts any magnitude in decibels to natural value.  
%
%   [result] = convertDecibelToNatural( valDB )
%
%   INPUTS:
%   valDB: Double value that represents a magnitude in decibels

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0

% Calculate the natural value of the input given in decibels
result = power( 10, valDB / 10 );

end

