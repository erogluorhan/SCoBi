
function writeVar(pathname, filename, var)


    if ~isempty(pathname)
        
        if exist(pathname, 'dir') ~= 7, mkdir( pathname ), end
        filename = strcat( pathname, '\', filename, '.dat' ) ;
        
    end
    
    fid = fopen(filename, 'w') ;
    
    [N, M] = size(var) ;
    
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
    