
function writeComplexVarIncremental(pathname, filename, index, cmplxvar)
% function writeComplexVarIncremental 
%
%   Writes the given complex number values to two separate files in an 
%   incremental fashion by using writeComplexVar function. In other words, 
%   it uses the main writer function. "writeComplexVar" and an index to add
%   the given "cmplxvar" into the "index" of the given "filename"
%
%   writeComplexVar(pathname, filename, cmplxvar)
%
%   INPUTS:
%   pathname  : The folder path to store the files
%   filename  : The filename that will be combined with "_r" and "_i" 
%   additions to represent the real and imaginary parts of the complex
%   numbers.
%   index     ; The index that "cmplxvar" will be put in
%   cmplxvar  : Complex number values that will be stored as real and
%   imaginary parts.
%
%   See also writeComplexVar, readComplexVar, readComplexVarIncremental

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0

% Write the real part in an incremental fashion
filename_real = strcat(filename, '_r') ;
writeVarIncremental(pathname, filename_real, index, real(cmplxvar)) ;

% Write the imaginary part in an incremental fashion
filename_imag = strcat(filename, '_i') ;
writeVarIncremental(pathname, filename_imag, index, imag(cmplxvar)) ;

end