function cpmlxvar = readComplexVar(pathname, filename)

filename_real = strcat(filename, '_r') ;
var_real = readVar(pathname, filename_real) ;

filename_imag = strcat(filename, '_i') ;
var_imag = readVar(pathname, filename_imag) ;


cpmlxvar = var_real + 1i * var_imag ;

end