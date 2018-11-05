
function writeVarIncremental(pathname, filename, index, var)
% function writeVarIncremental 
%
%   Writes the given number values to an existing files in an incremental
%   fashion by using the writeVar function. In other words, it uses the 
%   main writer function "writeVar" and an index to add the given "var" i
%   nto the "index" of the given "filename"   
%
%   writeVarIncremental(pathname, filename, index, var)
%
%   INPUTS:
%   pathname  : The folder path to store the files
%   filename  : The filename that will hold the numbers.
%   index     : The index that "var" will be put in
%   var       : Number values that will be stored.
%
%   See also writeVar, readVar

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



% First read the existing variable, if any    
currentVar = readVar(pathname, filename);
    
% If no current variable, write var as initial
if isnan( currentVar )

    [N, M] = size(var);

% Else if current variable exists, append var to the end
else
    % First read the current variable's size
    [N, M] = size(currentVar);
    
    if index > M
        % Enlarge the current var w.r.t var's end index
        M = index; 
    end
    
    currentVar(:, index) = var;
    var = currentVar;

end

% Create the folder if not exists, set the full filename
if ~isempty(pathname)

    if exist(pathname, 'dir') ~= 7
        mkdir( pathname )
    end
    
    filename = strcat( pathname, '\', filename, '.dat' ) ;

end

% Open the file to write
fid = fopen(filename, 'w') ;


fprintf(fid, '%6.4f\n', ndims(var)) ;
fprintf(fid, '%6.4f\n', size(var)) ;


for nn = 1 : N

    for mm = 1 : M

        fprintf(fid, '%6.16f  ', var(nn, mm)) ;

    end

    fprintf(fid, '\n') ;

end

fclose(fid) ;
    
end
    