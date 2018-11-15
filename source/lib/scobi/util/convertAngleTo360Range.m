%% Orhan Eroglu
% 13 July, 2017

%% Function that converts any angle to the [0,360] range
function newAngle = convertAngleTo360Range( angle )
% function convertAngleTo360Range 
%
%   Converts any angular value (in degrees) into [0, 360) range.  
%
%   newAngle = convertAngleTo360Range( angle )
%
%   INPUTS:
%   angle: Any negative or positive double (in degrees)

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0

% Caclulate the factor that drops the angle into [0, 360)
factor = floor( angle / 360 );

% Calculate the equivalent angular value in [0, 360)
newAngle = angle - factor * 360;


end