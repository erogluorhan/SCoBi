% Mehmet Kurum
% Feb 22, 2017


function SetNlayerVeg3

%% Global 
global FolderPath_Veg

%% Configuration

% Receive Antenna    
%           ^
%           |
%           h
%           |
%           |
%--------------------------------------------------------- z = 0
%                                           I = 1
%---------------------------------------------
%                                           I = 2
% --------------------------------------------
%                  
%......................................................... 
%                         
%                                           
%                   
% -----------------------------------------------
%                                           I = NLayer               
%=====GROUND=============================================== z = -d


%% Vegetation Paramters
D = [0.50 0.30]' ; % D = 80 cm
Nlayer = length(D) ;
TYPKND = [1 0 0 0; 1 1 0 0] ; % L, S
sTYPKND = sum(TYPKND) ;
Ntype = length(sTYPKND(sTYPKND ~= 0)) ; % L, S

Nkind = max(max(TYPKND)) ;

dsty = zeros(Nkind, Ntype, Nlayer) ;
dim1 = dsty ;
dim2 = dsty ;
dim3 = dsty ;
epsr = dsty ;
parm1 = dsty ;
parm2 = dsty ;

% to set if we calculate scattering
scatcalveg = zeros(Nkind, Ntype, Nlayer) ; 

%% Scatter Objects
[L1, T1] = Scat_Obj ;

%% Layer 1
% type 1 (Leaf), kind 1
dsty(1, 1, 1) = L1.DENSITY ; % in m^-3       
dim1(1, 1, 1) = L1.SEMI_MAJORX ; dim2(1, 1, 1) = L1.SEMI_MINORX ; 
dim3(1, 1, 1) = L1.THICKNESS ; epsr(1, 1, 1) = L1.EPSILON ; 
parm1(1, 1, 1) = L1.PARM1 ; parm2(1, 1, 1) = L1.PARM2 ;  
LTK{1, 1, 1} = 'L1' ;
scatcalveg(1, 1, 1) = 0 ; % NO scattering contribution

%% Layer 2
% type 1 (Leaf), kind 1
dsty(1, 1, 2) = L1.DENSITY ; % in m^-3       
dim1(1, 1, 2) = L1.SEMI_MAJORX ; dim2(1, 1, 2) = L1.SEMI_MINORX ; 
dim3(1, 1, 2) = L1.THICKNESS ; epsr(1, 1, 2) = L1.EPSILON ; 
parm1(1, 1, 2) = L1.PARM1 ; parm2(1, 1, 2) = L1.PARM2 ;  
LTK{1, 1, 2} = 'L1' ;
scatcalveg(1, 1, 2) = 0 ; % NO scattering contribution

% type 2 (Stem), kind 1
dsty(1, 2, 2) = T1.DENSITY / D(2) ; % in 1 m^-3
dim1(1, 2, 2) = T1.RADIUS ; dim2(1, 2, 2) = T1.RADIUS ; 
dim3(1, 2, 2) = T1.LENGTH ; epsr(1, 2, 2) = T1.EPSILON ;    
parm1(1, 2, 2) = T1.PARM1 ; parm2(1, 2, 2) = T1.PARM2 ;  
LTK{1, 2, 2} = 'T1' ;
scatcalveg(1, 2, 2) = 1 ; % YES scattering contribution

%% Save vegetation type kinds
save LTK LTK

%% Saving...

filename = 'scatcalveg' ;
write_var(FolderPath_Veg, filename, scatcalveg) ;

filename = 'D' ;
write_var(FolderPath_Veg, filename, D) ;

filename = 'TYPKND' ;
write_var(FolderPath_Veg, filename, TYPKND) ;

filename = 'DSTY' ;
write_var(FolderPath_Veg, filename, dsty) ;

filename = 'DIM1' ;
write_var(FolderPath_Veg, filename, dim1) ;

filename = 'DIM2' ;
write_var(FolderPath_Veg, filename, dim2) ;

filename = 'DIM3' ;
write_var(FolderPath_Veg, filename, dim3) ;

filename = 'ESPR' ;
write_cplxvar(FolderPath_Veg, filename, epsr) ;

filename = 'PARM1' ;
write_var(FolderPath_Veg, filename, parm1) ;

filename = 'PARM2' ;
write_var(FolderPath_Veg, filename, parm2) ;



end


%% Create Object Kinds
function [L1, T1] = Scat_Obj

dnsty = 2050 ;                     % density in m^-3
dim1 = 2.9 * 1e-2 ;                 % radius in m
dim2 = dim1 ;
dim3 = 0.2 * 1e-3 ;                % thichness in m
epsr = 23 + 1i * 9 ;             % dielectric Constant
prb = [0 45] ;                          % uniforma angle dustribution
ObjTyp = 'L1' ;
L1 = Gen_Object(dnsty, prb, dim1, dim2, dim3, epsr, ObjTyp) ;

dnsty = 680 ;                     % density in m^-2
dim1 = 0.127 * 1e-2 ;                   % radius in m
dim2 = dim1 ;                   % radius in m
dim3 = 26.7 * 1e-2 ;                   % length in m
epsr = 15 + 1i * 5 ;             % dielectric Constant
prb = [0 25] ;                          % uniforma angle dustribution
ObjTyp = 'T1' ;
T1 = Gen_Object(dnsty, prb, dim1, dim2, dim3, epsr, ObjTyp) ;

end


%% Generate objects for types (disk or cylider)
function Object = Gen_Object(dnsty, prob, dim1, dim2, dim3, epsr, ObjTyp)

% generate the object
if ObjTyp(1, 1) == 'L'
   
    Object = struct('PARM1', prob(1, 1), 'PARM2', prob(1, 2),...
        'SEMI_MAJORX', dim1, 'SEMI_MINORX', dim2, 'THICKNESS', dim3, ...
        'DENSITY', dnsty,'EPSILON', epsr, 'OBJ_CAT', ObjTyp) ;   
else
    
    Object = struct('PARM1', prob(1, 1), 'PARM2', prob(1, 2),...
        'RADIUS', dim1, 'LENGTH', dim3, ...
        'DENSITY', dnsty,'EPSILON', epsr, 'OBJ_CAT', ObjTyp) ;
end

end