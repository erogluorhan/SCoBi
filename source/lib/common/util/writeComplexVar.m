

function writeComplexVar(pathname, filename, cmplxvar)


filename_real = strcat(filename, '_r') ;
writeVar(pathname, filename_real, real(cmplxvar)) ;

filename_imag = strcat(filename, '_i') ;
writeVar(pathname, filename_imag, imag(cmplxvar)) ;

end