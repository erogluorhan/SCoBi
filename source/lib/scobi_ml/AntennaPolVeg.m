function [r0_coh1b, r0_coh2b] = AntennaPolVeg(r_s)


%addpath('C:\Users\kurum\Desktop\Purdue\Model\MSU_MultiLayerCPmodel\MSU_reflectometry_specular_only')

% add reflectometry file path that ignores diffuse terms
% fp = strcat(Directories.getInstance.scobi_ml, '\MSU_reflectometry_specular_only') ;
% addpath(fp) ;


% inputmaster2 ;
% 
% BiE2ES_main_specular_only ;
mainSCoBi;

[r0_coh1b, r0_coh2b] = SpecularReflection(r_s) ;

end