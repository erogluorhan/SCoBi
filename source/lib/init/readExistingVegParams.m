function readExistingVegParams

% TO-DO: Test the running cases of this function

%% GET GLOBAL DIRECTORIES
dir_veg = SimulationFolders.getInstance.veg;


dim_layers_m = readVar(dir_veg, ConstantNames.veg_hom_dimLayers_m) ;

scat_cal_veg = readVar(dir_veg, ConstantNames.veg_hom_scatCalVeg) ;

TYPKND = readVar(dir_veg, ConstantNames.veg_hom_TYPKND) ;

dsty = readVar(dir_veg, ConstantNames.veg_hom_dsty) ;

dim1_m = readVar(dir_veg, ConstantNames.veg_hom_dim1_m) ;

dim2_m = readVar(dir_veg, ConstantNames.veg_hom_dim2_m) ;

dim3_m = readVar(dir_veg, ConstantNames.veg_hom_dim3_m) ;

epsr = readComplexVar(dir_veg, ConstantNames.veg_hom_epsr) ;

parm1_deg = readVar(dir_veg, ConstantNames.veg_hom_parm1_deg) ;

parm2_deg = readVar(dir_veg, ConstantNames.veg_hom_parm2_deg) ;

if ( ~isnan( sum(dim_layers_m) ) && ~isnan( sum(sum(sum(scat_cal_veg))) ) && ...
     ~isnan( sum(sum(TYPKND)) ) && ~isnan( sum(sum(sum(dsty))) ) && ...
     ~isnan( sum(sum(sum(dim1_m))) ) && ~isnan( sum(sum(sum(dim2_m))) ) && ~isnan( sum(sum(sum(dim3_m))) ) && ...
     ~isnan( sum(sum(sum(epsr))) ) && ~isnan( sum(sum(sum(parm1_deg))) ) && ~isnan( sum(sum(sum(parm2_deg))) ) )

    VegParams.getInstance.initialize2( dim_layers_m, scat_cal_veg, TYPKND, ...
        dsty, dim1_m, dim2_m, dim3_m, epsr, parm1_deg, parm2_deg)


else
    disp('Error reading existing vegetation parameters!')    
end

end

