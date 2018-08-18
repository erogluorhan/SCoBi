

function writeComplexVarIncremental(pathname, filename, start_index, end_index, cmplxvar)


filename_real = strcat(filename, '_r') ;
writeVarIncremental(pathname, filename_real, start_index, end_index, real(cmplxvar)) ;

filename_imag = strcat(filename, '_i') ;
writeVarIncremental(pathname, filename_imag, start_index, end_index, imag(cmplxvar)) ;

end