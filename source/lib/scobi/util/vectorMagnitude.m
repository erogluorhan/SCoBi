

function abs = vectorMagnitude(vec)
% function vectorMagnitude 
%
%   Calculates magnitude of the given vector  
%
%   abs = vectorMagnitude(vec)
%
%   INPUTS:
%   vec: A row or column vector

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



% Get the size of vec
[m,n] = size(vec);

% Check if vec is of the proper dimensions or unit column or unit row
if (m ~= 1)&&(n ~= 1)  
    abs = 0;
    disp('Error - vector is not of proper dimensions');
else
    % Calculate the magnitude
    abs = sqrt(sum (vec.^2));
end

end