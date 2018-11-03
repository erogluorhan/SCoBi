
function writeComplexVar(pathname, filename, cmplxvar)
% function writeComplexVar 
%
%   Writes the given complex number values to two separate files to 
%   store the real and imaginary parts of the numbers separately, by using 
%   writeVar function. In fact, the given filename is combined with "_r" 
%   and "_i" additions, respectively. These files contain real-number 
%   values for real and imaginary parts of the complex numbers. 
%
%   writeComplexVar(pathname, filename, cmplxvar)
%
%   INPUTS:
%   pathname: The folder path to store the files
%   filename: The filename that will be combined with "_r" and "_i" 
%   additions to represent the real and imaginary parts of the complex
%   numbers.
%   cmplxvar: Complex number values that will be stored as real and
%   imaginary parts.
%
%   File Structure to be written: The real ("_r") and Imaginary ("_i") part
%   files have the same sized array of values (based on the size of 
%   cmplxvar input). The first line of the file shows the number of 
%   dimensinons of cmplxvar, following lines (enugh number of them 
%   depending on the number of dimensions) contain the length of each 
%   dimension, then the  following lines have the values corresponding to
%   the values given in cmplxvar. For example:
%   
%   cmplxvar = [ 19 + 1i, 25 + 2i, 16 + 1i, 7;
%                33 + 3i, 18 + 1i, 0,       4 + 4i;
%                40,      28 + 2i. 37 + 3i, 0 ]
%
%   example_r.dat
%   2.0
%   3.0
%   4.0
%   19 25 16 7
%   33 18 0  4
%   40 28 37 0
%   
%   example_i.dat
%   2.0
%   3.0
%   4.0
%   1 2 1 0
%   3 1 0 4
%   0 2 3 0
%
%   See also writeVar, readComplexVar, writeComplexVarIncremental, 
%   readComplexVarIncremental.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0


% Write the real part of cmplxvar
filename_real = strcat(filename, '_r') ;
writeVar(pathname, filename_real, real(cmplxvar)) ;

% Write the imaginary part of cmplxvar
filename_imag = strcat(filename, '_i') ;
writeVar(pathname, filename_imag, imag(cmplxvar)) ;

end