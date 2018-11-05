
function writeVar(pathname, filename, var)
% function writeVar 
%
%   Writes the given number values to file. 
%
%   writeVar(pathname, filename, var)
%
%   INPUTS:
%   pathname: The folder path to store the files
%   filename: The filename that will hold the numbers.
%   var:Number values that will be stored.
%
%   File Structure to be written: The file has the same sized array of 
%   values with the size of "var" input). The first line of the file shows 
%   the number of dimensinons of "var", following lines (enough number of 
%   them depending on the number of dimensions) contain the length of each 
%   dimension, then the  following lines have the values corresponding to
%   the values given in "var". For example:
%   
%   var = [ 19, 25, 16, 7;
%           33, 18, 0,  4;
%           40, 28, 37, 0 ]
%
%   example.dat
%   2.0
%   3.0
%   4.0
%   19 25 16 7
%   33 18 0  4
%   40 28 37 0
%
%   See also readVar, writeVarIncremental.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


    % Check if the given pathname variable is not empty
    if ~isempty(pathname)
        
        % If there is no such folder, create it
        if exist(pathname, 'dir') ~= 7, mkdir( pathname ), end
        
        % Consrtuct the full file name
        filename = strcat( pathname, '\', filename, '.dat' ) ;
        
    end
    
    % Open the file for writing
    fid = fopen(filename, 'w') ;
    
    % Get the size of var
    [N, M] = size(var) ;
    
    % Write the number of dimensions of var into the first line
    fprintf(fid, '%6.4f\n', ndims(var)) ;
    
    % Write the dimensions of var into the next lines
    fprintf(fid, '%6.4f\n', size(var)) ;
    
    % Write the values of var into the file
    for nn = 1 : N
        
        for mm = 1 : M
    
            fprintf(fid, '%6.16f  ', var(nn, mm)) ;
        
        end
        
        fprintf(fid, '\n') ;
        
    end
    
    % Close the file
    fclose(fid) ;
    
end
    