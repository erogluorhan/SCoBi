%% Mehmet Kurum
% May 02, 2017
% Modified November 10, 2017

function avgDiffuseTerm

%% Get Global Parameters
% Simulation Parameters
Nr = SimParams.getInstance.Nr ;
% Ground Parameters
VSM_cm3cm3 = GndParams.getInstance.VSM_cm3cm3( ParamsManager.index_VSM );
RMSH_cm = GndParams.getInstance.RMSH_cm( ParamsManager.index_RMSH );


%% Reading
% per kind  : P4 = 4 X Nfz X NkindMax X Ntype X Nlayer
% per type  : P3 = 4 X Nfz X Ntype X Nlayer
% per layer : P2 = 4 X Nfz X Nlayer
% medium    : P1 = 4 X Nfz


%% read Ki
filename1 = strcat('Ki') ;
Ki = readComplexVar(SimulationFolders.getInstance.freqdiff, filename1) ;
KKi = abs(Ki) ^ 2 / 4 / pi ;

%% Vegetation Information
% filename = 'TYPKND' ;
% TYPKND = readVar(FolderPath_Veg, filename) ;

% filename = 'scatcalveg' ;
% scatcalveg = readVar(FolderPath_Veg, filename) ;

%% Layer parameters
% TYPKND %#ok<NOPRT>
% [Nlayer, Ntype1] = size(TYPKND) ; %
% NkindMax = max(max(TYPKND)) ; 
% sTYPKND = sum(TYPKND) ;
% Ntype = length(sTYPKND(sTYPKND ~= 0)) ; % L, B, T

for rr = 1 : Nr
    
    disp(strcat('Realization: ', num2str(rr)))
    
    %% read output
    % P1
    pathname = strcat(SimulationFolders.getInstance.freqdiff_P1, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;
    filename1 = strcat('P1_inc1_t1', '_R', num2str(rr)) ;
    P1_inc1_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P1_inc2_t1', '_R', num2str(rr)) ;
    P1_inc2_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P1_inc3_t1', '_R', num2str(rr)) ;
    P1_inc3_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P1_inc4_t1', '_R', num2str(rr)) ;
    P1_inc4_t1 = readVar(pathname, filename1) ;
    
    filename1 = strcat('P1_inc1_t2', '_R', num2str(rr)) ;
    P1_inc1_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P1_inc2_t2', '_R', num2str(rr)) ;
    P1_inc2_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P1_inc3_t2', '_R', num2str(rr)) ;
    P1_inc3_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P1_inc4_t2', '_R', num2str(rr)) ;
    P1_inc4_t2 = readVar(pathname, filename1) ;
    
    % P2
    pathname = strcat(SimulationFolders.getInstance.freqdiff_P2, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;
    filename1 = strcat('P2_inc1_t1', '_R', num2str(rr)) ;
    P2_inc1_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P2_inc2_t1', '_R', num2str(rr)) ;
    P2_inc2_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P2_inc3_t1', '_R', num2str(rr)) ;
    P2_inc3_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P2_inc4_t1', '_R', num2str(rr)) ;
    P2_inc4_t1 = readVar(pathname, filename1) ;
    
    filename1 = strcat('P2_inc1_t2', '_R', num2str(rr)) ;
    P2_inc1_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P2_inc2_t2', '_R', num2str(rr)) ;
    P2_inc2_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P2_inc3_t2', '_R', num2str(rr)) ;
    P2_inc3_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P2_inc4_t2', '_R', num2str(rr)) ;
    P2_inc4_t2 = readVar(pathname, filename1) ;
    
    % P3
    pathname = strcat(SimulationFolders.getInstance.freqdiff_P3, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;
    filename1 = strcat('P3_inc1_t1', '_R', num2str(rr)) ;
    P3_inc1_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P3_inc2_t1', '_R', num2str(rr)) ;
    P3_inc2_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P3_inc3_t1', '_R', num2str(rr)) ;
    P3_inc3_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P3_inc4_t1', '_R', num2str(rr)) ;
    P3_inc4_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P3_inc1_t2', '_R', num2str(rr)) ;
    P3_inc1_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P3_inc2_t2', '_R', num2str(rr)) ;
    P3_inc2_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P3_inc3_t2', '_R', num2str(rr)) ;
    P3_inc3_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P3_inc4_t2', '_R', num2str(rr)) ;
    P3_inc4_t2 = readVar(pathname, filename1) ;
    
    % P4
    pathname = strcat(SimulationFolders.getInstance.freqdiff_P4, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;
    filename1 = strcat('P4_inc1_t1', '_R', num2str(rr)) ;
    P4_inc1_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P4_inc2_t1', '_R', num2str(rr)) ;
    P4_inc2_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P4_inc3_t1', '_R', num2str(rr)) ;
    P4_inc3_t1 = readVar(pathname, filename1) ;
    filename1 = strcat('P4_inc4_t1', '_R', num2str(rr)) ;
    P4_inc4_t1 = readVar(pathname, filename1) ;
    
    filename1 = strcat('P4_inc1_t2', '_R', num2str(rr)) ;
    P4_inc1_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P4_inc2_t2', '_R', num2str(rr)) ;
    P4_inc2_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P4_inc3_t2', '_R', num2str(rr)) ;
    P4_inc3_t2 = readVar(pathname, filename1) ;
    filename1 = strcat('P4_inc4_t2', '_R', num2str(rr)) ;
    P4_inc4_t2 = readVar(pathname, filename1) ;
    
    if rr == 1 % initilizing (Intensities)
        
        PP1_inc1_t1 = zeros(size(P1_inc1_t1)) ;
        PP1_inc2_t1 = PP1_inc1_t1 ;
        PP1_inc3_t1 = PP1_inc1_t1 ;
        PP1_inc4_t1 = PP1_inc1_t1 ;
        PP1_inc_t1 = PP1_inc1_t1 ;
        PP1_inc1_t2 = zeros(size(P1_inc1_t2)) ;
        PP1_inc2_t2 = PP1_inc1_t2 ;
        PP1_inc3_t2 = PP1_inc1_t2 ;
        PP1_inc4_t2 = PP1_inc1_t2 ;
        PP1_inc_t2 = PP1_inc1_t2 ;
        
        PP2_inc1_t1 = zeros(size(P2_inc1_t1)) ;
        PP2_inc2_t1 = PP2_inc1_t1 ;
        PP2_inc3_t1 = PP2_inc1_t1 ;
        PP2_inc4_t1 = PP2_inc1_t1 ;
        PP2_inc_t1 = PP2_inc1_t1 ;
        PP2_inc1_t2 = zeros(size(P2_inc1_t2)) ;
        PP2_inc2_t2 = PP2_inc1_t2 ;
        PP2_inc3_t2 = PP2_inc1_t2 ;
        PP2_inc4_t2 = PP2_inc1_t2 ;
        PP2_inc_t2 = PP2_inc1_t2 ;
        
        PP3_inc1_t1 = zeros(size(P3_inc1_t1)) ;
        PP3_inc2_t1 = PP3_inc1_t1 ;
        PP3_inc3_t1 = PP3_inc1_t1 ;
        PP3_inc4_t1 = PP3_inc1_t1 ;
        PP3_inc_t1 = PP3_inc1_t1 ;
        PP3_inc1_t2 = zeros(size(P3_inc1_t2)) ;
        PP3_inc2_t2 = PP3_inc1_t2 ;
        PP3_inc3_t2 = PP3_inc1_t2 ;
        PP3_inc4_t2 = PP3_inc1_t2 ;
        PP3_inc_t2 = PP3_inc1_t2 ;
        
        PP4_inc1_t1 = zeros(size(P4_inc1_t1)) ;
        PP4_inc2_t1 = PP4_inc1_t1 ;
        PP4_inc3_t1 = PP4_inc1_t1 ;
        PP4_inc4_t1 = PP4_inc1_t1 ;
        PP4_inc_t1 = PP4_inc1_t1 ;
        PP4_inc1_t2 = zeros(size(P4_inc1_t2)) ;
        PP4_inc2_t2 = PP4_inc1_t2 ;
        PP4_inc3_t2 = PP4_inc1_t2 ;
        PP4_inc4_t2 = PP4_inc1_t2 ;
        PP4_inc_t2 = PP4_inc1_t2 ;
    end
    
    %% summing intensities over realizations
    % medium
    PP1_inc1_t1 = PP1_inc1_t1 + P1_inc1_t1 ;
    PP1_inc2_t1 = PP1_inc2_t1 + P1_inc2_t1 ;
    PP1_inc3_t1 = PP1_inc3_t1 + P1_inc3_t1 ;
    PP1_inc4_t1 = PP1_inc4_t1 + P1_inc4_t1 ;
    P1_inc_t1 = P1_inc1_t1 + P1_inc2_t1 + P1_inc3_t1 + P1_inc4_t1 ;
    PP1_inc_t1 = PP1_inc_t1 + P1_inc_t1 ;
    
    PP1_inc1_t2 = PP1_inc1_t2 + P1_inc1_t2 ;
    PP1_inc2_t2 = PP1_inc2_t2 + P1_inc2_t2 ;
    PP1_inc3_t2 = PP1_inc3_t2 + P1_inc3_t2 ;
    PP1_inc4_t2 = PP1_inc4_t2 + P1_inc4_t2 ;
    P1_inc_t2 = P1_inc1_t2 + P1_inc2_t2 + P1_inc3_t2 + P1_inc4_t2 ;
    PP1_inc_t2 = PP1_inc_t2 + P1_inc_t2 ;
    
    % layer
    PP2_inc1_t1 = PP2_inc1_t1 + P2_inc1_t1 ;
    PP2_inc2_t1 = PP2_inc2_t1 + P2_inc2_t1 ;
    PP2_inc3_t1 = PP2_inc3_t1 + P2_inc3_t1 ;
    PP2_inc4_t1 = PP2_inc4_t1 + P2_inc4_t1 ;
    P2_inc_t1 = P2_inc1_t1 + P2_inc2_t1 + P2_inc3_t1 + P2_inc4_t1 ;
    PP2_inc_t1 = PP2_inc_t1 + P2_inc_t1 ;
    
    PP2_inc1_t2 = PP2_inc1_t2 + P2_inc1_t2 ;
    PP2_inc2_t2 = PP2_inc2_t2 + P2_inc2_t2 ;
    PP2_inc3_t2 = PP2_inc3_t2 + P2_inc3_t2 ;
    PP2_inc4_t2 = PP2_inc4_t2 + P2_inc4_t2 ;
    P2_inc_t2 = P2_inc1_t2 + P2_inc2_t2 + P2_inc3_t2 + P2_inc4_t2 ;
    PP2_inc_t2 = PP2_inc_t2 + P2_inc_t2 ;
    
    % type
    PP3_inc1_t1 = PP3_inc1_t1 + P3_inc1_t1 ;
    PP3_inc2_t1 = PP3_inc2_t1 + P3_inc2_t1 ;
    PP3_inc3_t1 = PP3_inc3_t1 + P3_inc3_t1 ;
    PP3_inc4_t1 = PP3_inc4_t1 + P3_inc4_t1 ;
    P3_inc_t1 = P3_inc1_t1 + P3_inc2_t1 + P3_inc3_t1 + P3_inc4_t1 ;
    PP3_inc_t1 = PP3_inc_t1 + P3_inc_t1 ;
    
    PP3_inc1_t2 = PP3_inc1_t2 + P3_inc1_t2 ;
    PP3_inc2_t2 = PP3_inc2_t2 + P3_inc2_t2 ;
    PP3_inc3_t2 = PP3_inc3_t2 + P3_inc3_t2 ;
    PP3_inc4_t2 = PP3_inc4_t2 + P3_inc4_t2 ;
    P3_inc_t2 = P3_inc1_t2 + P3_inc2_t2 + P3_inc3_t2 + P3_inc4_t2 ;
    PP3_inc_t2 = PP3_inc_t2 + P3_inc_t2 ;
    
    % kind
    PP4_inc1_t1 = PP4_inc1_t1 + P4_inc1_t1 ;
    PP4_inc2_t1 = PP4_inc2_t1 + P4_inc2_t1 ;
    PP4_inc3_t1 = PP4_inc3_t1 + P4_inc3_t1 ;
    PP4_inc4_t1 = PP4_inc4_t1 + P4_inc4_t1 ;
    P4_inc_t1 = P4_inc1_t1 + P4_inc2_t1 + P4_inc3_t1 + P4_inc4_t1 ;
    PP4_inc_t1 = PP4_inc_t1 + P4_inc_t1 ;
    
    PP4_inc1_t2 = PP4_inc1_t2 + P4_inc1_t2 ;
    PP4_inc2_t2 = PP4_inc2_t2 + P4_inc2_t2 ;
    PP4_inc3_t2 = PP4_inc3_t2 + P4_inc3_t2 ;
    PP4_inc4_t2 = PP4_inc4_t2 + P4_inc4_t2 ;
    P4_inc_t2 = P4_inc1_t2 + P4_inc2_t2 + P4_inc3_t2 + P4_inc4_t2 ;
    PP4_inc_t2 = PP4_inc_t2 + P4_inc_t2 ;
    
end  % Nr


%% Avergaging received power and converting dB ...
% medium
PP1_inc1_t1_dB = 10 * log10(PP1_inc1_t1(1 : 2, :, :) / Nr * KKi) ;
PP1_inc2_t1_dB = 10 * log10(PP1_inc2_t1(1 : 2, :, :) / Nr * KKi) ;
PP1_inc3_t1_dB = 10 * log10(PP1_inc3_t1(1 : 2, :, :) / Nr * KKi) ;
PP1_inc4_t1_dB = 10 * log10(PP1_inc4_t1(1 : 2, :, :) / Nr * KKi) ;
PP1_inc_t1_dB = 10 * log10(PP1_inc_t1(1 : 2, :, :) / Nr * KKi) ;
PP1_inc1_t2_dB = 10 * log10(PP1_inc1_t2(1 : 2, :, :) / Nr * KKi) ;
PP1_inc2_t2_dB = 10 * log10(PP1_inc2_t2(1 : 2, :, :) / Nr * KKi) ;
PP1_inc3_t2_dB = 10 * log10(PP1_inc3_t2(1 : 2, :, :) / Nr * KKi) ;
PP1_inc4_t2_dB = 10 * log10(PP1_inc4_t2(1 : 2, :, :) / Nr * KKi) ;
PP1_inc_t2_dB = 10 * log10(PP1_inc_t2(1 : 2, :, :) / Nr * KKi) ;

% layer
PP2_inc1_t1_dB = 10 * log10(PP2_inc1_t1(1 : 2, :, :, :) / Nr * KKi) ;
PP2_inc2_t1_dB = 10 * log10(PP2_inc2_t1(1 : 2, :, :, :) / Nr * KKi) ;
PP2_inc3_t1_dB = 10 * log10(PP2_inc3_t1(1 : 2, :, :, :) / Nr * KKi) ;
PP2_inc4_t1_dB = 10 * log10(PP2_inc4_t1(1 : 2, :, :, :) / Nr * KKi) ;
PP2_inc_t1_dB = 10 * log10(PP2_inc_t1(1 : 2, :, :, :) / Nr * KKi) ;
PP2_inc1_t2_dB = 10 * log10(PP2_inc1_t2(1 : 2, :, :, :) / Nr * KKi) ;
PP2_inc2_t2_dB = 10 * log10(PP2_inc2_t2(1 : 2, :, :, :) / Nr * KKi) ;
PP2_inc3_t2_dB = 10 * log10(PP2_inc3_t2(1 : 2, :, :, :) / Nr * KKi) ;
PP2_inc4_t2_dB = 10 * log10(PP2_inc4_t2(1 : 2, :, :, :) / Nr * KKi) ;
PP2_inc_t2_dB = 10 * log10(PP2_inc_t2(1 : 2, :, :, :) / Nr * KKi) ;

% type
PP3_inc1_t1_dB = 10 * log10(PP3_inc1_t1(1 : 2, :, :, :, :) / Nr * KKi) ;
PP3_inc2_t1_dB = 10 * log10(PP3_inc2_t1(1 : 2, :, :, :, :) / Nr * KKi) ;
PP3_inc3_t1_dB = 10 * log10(PP3_inc3_t1(1 : 2, :, :, :, :) / Nr * KKi) ;
PP3_inc4_t1_dB = 10 * log10(PP3_inc4_t1(1 : 2, :, :, :, :) / Nr * KKi) ;
PP3_inc_t1_dB = 10 * log10(PP3_inc_t1(1 : 2, :, :, :, :) / Nr * KKi) ;
PP3_inc1_t2_dB = 10 * log10(PP3_inc1_t2(1 : 2, :, :, :, :) / Nr * KKi) ;
PP3_inc2_t2_dB = 10 * log10(PP3_inc2_t2(1 : 2, :, :, :, :) / Nr * KKi) ;
PP3_inc3_t2_dB = 10 * log10(PP3_inc3_t2(1 : 2, :, :, :, :) / Nr * KKi) ;
PP3_inc4_t2_dB = 10 * log10(PP3_inc4_t2(1 : 2, :, :, :, :) / Nr * KKi) ;
PP3_inc_t2_dB = 10 * log10(PP3_inc_t2(1 : 2, :, :, :, :) / Nr * KKi) ;

% kind
PP4_inc1_t1_dB = 10 * log10(PP4_inc1_t1(1 : 2, :, :, :, :, :) / Nr * KKi) ;
PP4_inc2_t1_dB = 10 * log10(PP4_inc2_t1(1 : 2, :, :, :, :, :) / Nr * KKi) ;
PP4_inc3_t1_dB = 10 * log10(PP4_inc3_t1(1 : 2, :, :, :, :, :) / Nr * KKi) ;
PP4_inc4_t1_dB = 10 * log10(PP4_inc4_t1(1 : 2, :, :, :, :, :) / Nr * KKi) ;
PP4_inc_t1_dB = 10 * log10(PP4_inc_t1(1 : 2, :, :, :, :, :) / Nr * KKi) ;
PP4_inc1_t2_dB = 10 * log10(PP4_inc1_t2(1 : 2, :, :, :, :, :) / Nr * KKi) ;
PP4_inc2_t2_dB = 10 * log10(PP4_inc2_t2(1 : 2, :, :, :, :, :) / Nr * KKi) ;
PP4_inc3_t2_dB = 10 * log10(PP4_inc3_t2(1 : 2, :, :, :, :, :) / Nr * KKi) ;
PP4_inc4_t2_dB = 10 * log10(PP4_inc4_t2(1 : 2, :, :, :, :, :) / Nr * KKi) ;
PP4_inc_t2_dB = 10 * log10(PP4_inc_t2(1 : 2, :, :, :, :, :) / Nr * KKi) ;

%% Fresnel ellipses
filenamex = 'ellipse_s_m' ;
ellipse_s_m = readVar(SimulationFolders.getInstance.config, filenamex) ;
area_s = pi * ellipse_s_m(:, 1) .* ellipse_s_m(:, 2) ;
% P1_areas_s = repmat(area_s, 1, 2, length(VSM_cm3cm3)) ;
P1_areas_s = repmat(area_s, 1, 2 ) ;
% TO-DO: Check the following: permutation converted to transpose due to the
% removal of multiple SM values from the output 
% P1_areas_s = permute(P1_areas_s, [2, 3, 1]) ;
P1_areas_s = P1_areas_s' ;
P1_areas_dB = 10 * log10(P1_areas_s) ;

%%
KKi_dB = 10 * log10(KKi) ;

%% NBRCS - medium
NBRCS1_t1_dB = PP1_inc1_t1_dB - KKi_dB - P1_areas_dB ;
NBRCS2_t1_dB = PP1_inc2_t1_dB - KKi_dB - P1_areas_dB ;
NBRCS3_t1_dB = PP1_inc3_t1_dB - KKi_dB - P1_areas_dB ;
NBRCS4_t1_dB = PP1_inc4_t1_dB - KKi_dB - P1_areas_dB ;
NBRCS_t1_dB = PP1_inc_t1_dB - KKi_dB - P1_areas_dB ;
NBRCS1_t2_dB = PP1_inc1_t2_dB - KKi_dB - P1_areas_dB ;
NBRCS2_t2_dB = PP1_inc2_t2_dB - KKi_dB - P1_areas_dB ;
NBRCS3_t2_dB = PP1_inc3_t2_dB - KKi_dB - P1_areas_dB ;
NBRCS4_t2_dB = PP1_inc4_t2_dB - KKi_dB - P1_areas_dB ;
NBRCS_t2_dB = PP1_inc_t2_dB - KKi_dB - P1_areas_dB ;

%% Saving
% P1
pathname = strcat(SimulationFolders.getInstance.out_diffuse_P1, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;
filename1 = strcat('PP1_inc1_t1_dB') ;
writeVar(pathname, filename1, (PP1_inc1_t1_dB))
filename1 = strcat('PP1_inc2_t1_dB') ;
writeVar(pathname, filename1, (PP1_inc2_t1_dB))
filename1 = strcat('PP1_inc3_t1_dB') ;
writeVar(pathname, filename1, (PP1_inc3_t1_dB))
filename1 = strcat('PP1_inc4_t1_dB') ;
writeVar(pathname, filename1, (PP1_inc4_t1_dB))
filename1 = strcat('PP1_inc_t1_dB') ;
writeVar(pathname, filename1, (PP1_inc_t1_dB))

filename1 = strcat('PP1_inc1_t2_dB') ;
writeVar(pathname, filename1, (PP1_inc1_t2_dB))
filename1 = strcat('PP1_inc2_t2_dB') ;
writeVar(pathname, filename1, (PP1_inc2_t2_dB))
filename1 = strcat('PP1_inc3_t2_dB') ;
writeVar(pathname, filename1, (PP1_inc3_t2_dB))
filename1 = strcat('PP1_inc4_t2_dB') ;
writeVar(pathname, filename1, (PP1_inc4_t2_dB))
filename1 = strcat('PP1_inc_t2_dB') ;
writeVar(pathname, filename1, (PP1_inc_t2_dB))

% P2
pathname = strcat(SimulationFolders.getInstance.out_diffuse_P2, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;
filename1 = strcat('PP2_inc1_t1_dB') ;
writeVar(pathname, filename1, (PP2_inc1_t1_dB))
filename1 = strcat('PP2_inc2_t1_dB') ;
writeVar(pathname, filename1, (PP2_inc2_t1_dB))
filename1 = strcat('PP2_inc3_t1_dB') ;
writeVar(pathname, filename1, (PP2_inc3_t1_dB))
filename1 = strcat('PP2_inc4_t1_dB') ;
writeVar(pathname, filename1, (PP2_inc4_t1_dB))
filename1 = strcat('PP2_inc_t1_dB') ;
writeVar(pathname, filename1, (PP2_inc_t1_dB))

filename1 = strcat('PP2_inc1_t2_dB') ;
writeVar(pathname, filename1, (PP2_inc1_t2_dB))
filename1 = strcat('PP2_inc2_t2_dB') ;
writeVar(pathname, filename1, (PP2_inc2_t2_dB))
filename1 = strcat('PP2_inc3_t2_dB') ;
writeVar(pathname, filename1, (PP2_inc3_t2_dB))
filename1 = strcat('PP2_inc4_t2_dB') ;
writeVar(pathname, filename1, (PP2_inc4_t2_dB))
filename1 = strcat('PP2_inc_t2_dB') ;
writeVar(pathname, filename1, (PP2_inc_t2_dB))

% P3
pathname = strcat(SimulationFolders.getInstance.out_diffuse_P3, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;
filename1 = strcat('PP3_inc1_t1_dB') ;
writeVar(pathname, filename1, (PP3_inc1_t1_dB))
filename1 = strcat('PP3_inc2_t1_dB') ;
writeVar(pathname, filename1, (PP3_inc2_t1_dB))
filename1 = strcat('PP3_inc3_t1_dB') ;
writeVar(pathname, filename1, (PP3_inc3_t1_dB))
filename1 = strcat('PP3_inc4_t1_dB') ;
writeVar(pathname, filename1, (PP3_inc4_t1_dB))
filename1 = strcat('PP3_inc_t1_dB') ;
writeVar(pathname, filename1, (PP3_inc_t1_dB))

filename1 = strcat('PP3_inc1_t2_dB') ;
writeVar(pathname, filename1, (PP3_inc1_t2_dB))
filename1 = strcat('PP3_inc2_t2_dB') ;
writeVar(pathname, filename1, (PP3_inc2_t2_dB))
filename1 = strcat('PP3_inc3_t2_dB') ;
writeVar(pathname, filename1, (PP3_inc3_t2_dB))
filename1 = strcat('PP3_inc4_t2_dB') ;
writeVar(pathname, filename1, (PP3_inc4_t2_dB))
filename1 = strcat('PP3_inc_t2_dB') ;
writeVar(pathname, filename1, (PP3_inc_t2_dB))

% P4
pathname = strcat(SimulationFolders.getInstance.out_diffuse_P4, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;
filename1 = strcat('PP4_inc1_t1_dB') ;
writeVar(pathname, filename1, (PP4_inc1_t1_dB))
filename1 = strcat('PP4_inc2_t1_dB') ;
writeVar(pathname, filename1, (PP4_inc2_t1_dB))
filename1 = strcat('PP4_inc3_t1_dB') ;
writeVar(pathname, filename1, (PP4_inc3_t1_dB))
filename1 = strcat('PP4_inc4_t1_dB') ;
writeVar(pathname, filename1, (PP4_inc4_t1_dB))
filename1 = strcat('PP4_inc_t1_dB') ;
writeVar(pathname, filename1, (PP4_inc_t1_dB))

filename1 = strcat('PP4_inc1_t2_dB') ;
writeVar(pathname, filename1, (PP4_inc1_t2_dB))
filename1 = strcat('PP4_inc2_t2_dB') ;
writeVar(pathname, filename1, (PP4_inc2_t2_dB))
filename1 = strcat('PP4_inc3_t2_dB') ;
writeVar(pathname, filename1, (PP4_inc3_t2_dB))
filename1 = strcat('PP4_inc4_t2_dB') ;
writeVar(pathname, filename1, (PP4_inc4_t2_dB))
filename1 = strcat('PP4_inc_t2_dB') ;
writeVar(pathname, filename1, (PP4_inc_t2_dB))

% NBRCS
pathname = strcat(SimulationFolders.getInstance.out_diffuse_NBRCS, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;
filename1 = strcat('NBRCS1_t1_dB') ;
writeVar(pathname, filename1, (NBRCS1_t1_dB))
filename1 = strcat('NBRCS2_t1_dB') ;
writeVar(pathname, filename1, (NBRCS2_t1_dB))
filename1 = strcat('NBRCS3_t1_dB') ;
writeVar(pathname, filename1, (NBRCS3_t1_dB))
filename1 = strcat('NBRCS4_t1_dB') ;
writeVar(pathname, filename1, (NBRCS4_t1_dB))
filename1 = strcat('NBRCS_t1_dB') ;
writeVar(pathname, filename1, (NBRCS_t1_dB))

filename1 = strcat('NBRCS1_t2_dB') ;
writeVar(pathname, filename1, (NBRCS1_t2_dB))
filename1 = strcat('NBRCS2_t2_dB') ;
writeVar(pathname, filename1, (NBRCS2_t2_dB))
filename1 = strcat('NBRCS3_t2_dB') ;
writeVar(pathname, filename1, (NBRCS3_t2_dB))
filename1 = strcat('NBRCS4_t2_dB') ;
writeVar(pathname, filename1, (NBRCS4_t2_dB))
filename1 = strcat('NBRCS_t2_dB') ;
writeVar(pathname, filename1, (NBRCS_t2_dB))

end
