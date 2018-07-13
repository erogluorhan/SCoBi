

function [doy_SM, SM5, SM10, SM20, SM40] = Read_Ground_Data


clc ;

%% Directory
myDirect = pwd ;
dirSM = strcat(Directories.getInstance.scobi_ml, '\Ground Data\SM') ;

%% Date - Soil Moisture
fname = dir(fullfile(dirSM, 'DOY*'))' ;
fileName = strcat(dirSM, '\', fname.name) ;
doy_SM = load(fileName) ;

%% Soil Moisture - Corn
fname = dir(fullfile(dirSM, 'SM5*'))' ;
fileName = strcat(dirSM, '\', fname.name) ;
SM5 = load(fileName) ;
fname = dir(fullfile(dirSM, 'SM10*'))' ;
fileName = strcat(dirSM, '\', fname.name) ;
SM10 = load(fileName) ;
fname = dir(fullfile(dirSM, 'SM20*'))' ;
fileName = strcat(dirSM, '\', fname.name) ;
SM20 = load(fileName) ;
fname = dir(fullfile(dirSM, 'SM40*'))' ;
fileName = strcat(dirSM, '\', fname.name) ;
SM40 = load(fileName) ;

end