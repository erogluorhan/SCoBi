function readExistingVegParams

% TO-DO: Test the running cases of this function

pathname = SimulationFolders.getInstance.veg;


dim_layers = readVar(pathname, 'dim_layers') ;

scat_cal_veg = readVar(pathname, 'scat_cal_veg') ;

TYPKND = readVar(pathname, 'TYPKND') ;

dsty = readVar(pathname, 'dsty') ;

dim1 = readVar(pathname, 'dim1') ;

dim2 = readVar(pathname, 'dim2') ;

dim3 = readVar(pathname, 'dim3') ;

epsr = readComplexVar(pathname, 'epsr') ;

parm1 = readVar(pathname, 'parm1') ;

parm2 = readVar(pathname, 'parm2') ;

if ( ~isnan( sum(dim_layers) ) && ~isnan( sum(sum(sum(scat_cal_veg))) ) && ...
     ~isnan( sum(sum(TYPKND)) ) && ~isnan( sum(sum(sum(dsty))) ) && ...
     ~isnan( sum(sum(sum(dim1))) ) && ~isnan( sum(sum(sum(dim2))) ) && ~isnan( sum(sum(sum(dim3))) ) && ...
     ~isnan( sum(sum(sum(epsr))) ) && ~isnan( sum(sum(sum(parm1))) ) && ~isnan( sum(sum(sum(parm2))) ) )

    VegParams.getInstance.initialize2( dim_layers, scat_cal_veg, TYPKND, ...
        dsty, dim1, dim2, dim3, epsr, parm1, parm2)


else
    disp('Error reading existing vegetation parameters!')    
end

end

