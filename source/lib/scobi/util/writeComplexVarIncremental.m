

function writeComplexVarIncremental(pathname, filename, index, cmplxvar)


filename_real = strcat(filename, '_r') ;
writeVarIncremental(pathname, filename_real, index, real(cmplxvar)) ;

filename_imag = strcat(filename, '_i') ;
writeVarIncremental(pathname, filename_imag, index, imag(cmplxvar)) ;

end