function [ element_id ] = findElementIdInCell( cellArr, element )
% function findElementIdInCell 
%
%   Returns the index of a value in a cell array, if exists. If not, 
%   returns -1  
%
%   [ element_id ] = findElementIdInCell( cellArr, element )
%
%   INPUTS:
%   cellArr: Cell array
%   element: Value to be found in the cell array

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


is_element = cellfun(@(x)isequal(x, element), cellArr );

[~, element_id] = find(is_element);

end

