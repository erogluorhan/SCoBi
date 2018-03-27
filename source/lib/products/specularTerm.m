%% Mehmet Kurum
% 03/12/2017
% modified - 11/09/2017
% modified - 11/15/2017

% Specular (Coherent) Term
function specularTerm

%% GET GLOBAL DIRECTORIES
dir_config = SimulationFolders.getInstance.config;
dir_ant_lookup = SimulationFolders.getInstance.ant_lookup ;
dir_rot_lookup = SimulationFolders.getInstance.rot_lookup ;
dir_out_specular = SimulationFolders.getInstance.out_specular;


%% GET GLOBAL PARAMETERS
% Receiver Parameters
G0r_dB = RecParams.getInstance.G0r_dB;
G0r = convertDecibelToNatural( G0r_dB );
polR = RecParams.getInstance.polR;
% Satellite Parameters
EIRP_dB = SatParams.getInstance.EIRP_dB;
EIRP = convertDecibelToNatural( EIRP_dB );
f_MHz = SatParams.getInstance.f_MHz;
g_t = SatParams.getInstance.g_t ; % ideal
e_t1 = SatParams.getInstance.e_t1 ; 
e_t2 = SatParams.getInstance.e_t2 ;
polT = SatParams.getInstance.polT;
% Ground Parameters
VSM_cm3cm3 = GndParams.getInstance.VSM_cm3cm3( ParamsManager.index_VSM );
RMSH_cm = GndParams.getInstance.RMSH_cm( ParamsManager.index_RMSH );


%% READ OR LOAD META-DATA
% Load Receiver Antenna Pattern
load([dir_ant_lookup '\AntPat.mat'], 'G', 'g', 'th', 'ph')

% Read All Positions
% AllPoints_m = [pT_m, pTI_m, pS2_m, pR_m, pRI_m, pG2_m, pBr2_m, pCr2_m, pSc2_m] ;
filenamex = 'AllPoints_m' ;
AllPoints_m = readVar(dir_config, filenamex) ;
pT_m = AllPoints_m(:, 1) ;        % Transmitter position
pS2_m = AllPoints_m(:, 3) ;       % Specular point
pR_m = AllPoints_m(:, 4) ;        % Receiver position

% Read SP-to-Rec Rotation Angle
filename = 'AngS2R_rf' ;
AngS2R_rf = readVar(dir_config, filename) ;
% Read Sat-to-SP Rotation Angle
filename = 'AngT2S_sf' ;
AngT2S_sf = readVar(dir_config, filename) ;

% Transmitter-Receiver Rotation Matrix
load([dir_rot_lookup '\u_ts.mat'], 'u_ts')
load([dir_rot_lookup '\u_sr.mat'], 'u_sr')


%% INITIALIZE REQUIRED PARAMETERS
% Transmitter Pol State
E_t1 = [1; 0; 0; 0] ;
E_t2 = [0; 0; 0; 1] ;
Q = [1 0 0 0; 0 0 0 1; 0 1 1 0; 0 -1i -1i 0] ;
E_t1 = Q * E_t1 ;
E_t2 = Q * E_t2 ;
% Ideal Receiver Antenna pattern
g_r0 = [1 0 ; 0 1] ; % ideal


%% CALCULATIONS
% Slant Range
ST = pS2_m - pT_m ;         % Transmitter to Specular
r_st = vectorMagnitude(ST) ;  % slant range

RS = pR_m - pS2_m ;         % Specular to Receiver
r_sr = vectorMagnitude(RS) ;  % slant range

f_Hz = f_MHz * Constants.MHz2Hz ;
lambda_m = Constants.c / f_Hz ;     % Wavelength
k0 = 2 * pi * f_Hz / Constants.c ;    % Wave number

% Factor Kc
K = 1i * sqrt(EIRP) * sqrt(G0r) * lambda_m / (4 * pi) ;

Kc = K * exp(1i * k0 * (r_st + r_sr)) / (r_st + r_sr) ; 

thrd = AngS2R_rf(1, 1) ;
phrd = AngS2R_rf(2, 1) ;

thd = round(2 * radtodeg(th)) / 2 ; % rounding operation is due to accuracy concerns
phd = round(2 * radtodeg(ph)) / 2 ;

% Receiver Antenna values in the transmitter directions
ind_th = thd == round(2 * thrd(1, 1)) / 2 ; % round is to make it a multiple of 0.5
ind_ph = phd == round(2 * phrd(1, 1)) / 2 ;

g11 = cell2mat(g(1,1)) ; g12 = cell2mat(g(1,2)) ;
g21 = cell2mat(g(2,1)) ; g22 = cell2mat(g(2,2)) ;

g_r = [g11(ind_th & ind_ph), g12(ind_th & ind_ph); ...
    g21(ind_th & ind_ph), g22(ind_th & ind_ph)] ;

% 4 X 4
G_r = calc_Muller(g_r) ;
% 4 X 4
G_r0 = calc_Muller(g_r0) ;
% 4 X 4
G_t = calc_Muller(g_t) ;

% Transmitter-Receiver Rotation Matrix
% 4 X 4
U_ts = calc_Muller(u_ts) ;
% 4 X 4
U_sr = calc_Muller(u_sr) ;

% Specular Reflection Matrix
thsd = AngT2S_sf(1, 1) ;

[R_sv, R_sb, r_sv, r_sb] = CalcSRM(thsd) ; % r_sv for vegetation, r_sb for bare soil


%% SPECULAR TERM
% P = 4 x 1
% b = 2 x 1

b_coh1v = g_r * u_sr * r_sv * u_ts * g_t * e_t1 ;   % field 
b_coh2v = g_r * u_sr * r_sv * u_ts * g_t * e_t2 ;   % field
b0_coh1v = g_r0 * u_sr * r_sv * u_ts * g_t * e_t1 ;  % field
b0_coh2v = g_r0 * u_sr * r_sv * u_ts * g_t * e_t2 ;  % field

b_coh1b = g_r * u_sr * r_sb * u_ts * g_t * e_t1 ;   % field
b_coh2b = g_r * u_sr * r_sb * u_ts * g_t * e_t2 ;   % field
b0_coh1b = g_r0 * u_sr * r_sb * u_ts * g_t * e_t1 ;  % field
b0_coh2b = g_r0 * u_sr * r_sb * u_ts * g_t * e_t2 ;  % field

P_coh1v = G_r * U_sr * R_sv * U_ts * G_t * E_t1 ;   % POWER
P_coh2v = G_r * U_sr * R_sv * U_ts * G_t * E_t2 ;   % POWER
P0_coh1v = G_r0 * U_sr * R_sv * U_ts * G_t * E_t1 ;  % POWER
P0_coh2v = G_r0 * U_sr * R_sv * U_ts * G_t * E_t2 ;  % POWER

P_coh1b = G_r * U_sr * R_sb * U_ts * G_t * E_t1 ;   % POWER
P_coh2b = G_r * U_sr * R_sb * U_ts * G_t * E_t2 ;   % POWER
P0_coh1b = G_r0 * U_sr * R_sb * U_ts * G_t * E_t1 ;  % POWER
P0_coh2b = G_r0 * U_sr * R_sb * U_ts * G_t * E_t2 ;  % POWER


%% SAVE OUTPUTS
% Kc Factor
filename3 = strcat('Kc') ;
writeComplexVar(dir_out_specular, filename3, Kc)

% 2 X 2
pathname = strcat(dir_out_specular, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) ) ;

filename1 = strcat('Veg1', polT, polR) ;
filename2 = strcat('Veg2', polT, polR) ;
filename01 = strcat('Veg01', polT, polR) ;
filename02 = strcat('Veg02', polT, polR) ;
writeComplexVar(pathname, filename1, b_coh1v)
writeComplexVar(pathname, filename2, b_coh2v)
writeComplexVar(pathname, filename01, b0_coh1v)
writeComplexVar(pathname, filename02, b0_coh2v)

filename1 = strcat('Bare1', polT, polR) ;
filename2 = strcat('Bare2', polT, polR) ;
filename01 = strcat('Bare01', polT, polR) ;
filename02 = strcat('Bare02', polT, polR) ;
writeComplexVar(pathname, filename1, b_coh1b)
writeComplexVar(pathname, filename2, b_coh2b)
writeComplexVar(pathname, filename01, b0_coh1b)
writeComplexVar(pathname, filename02, b0_coh2b)

% 4 X 4
filename1 = strcat('P_Veg1', polT, polR) ;
filename2 = strcat('P_Veg2', polT, polR) ;
filename01 = strcat('P_Veg01', polT, polR) ;
filename02 = strcat('P_Veg02', polT, polR) ;
writeVar(pathname, filename1, P_coh1v)
writeVar(pathname, filename2, P_coh2v)
writeVar(pathname, filename01, P0_coh1v)
writeVar(pathname, filename02, P0_coh2v)

filename1 = strcat('P_Bare1', polT, polR) ;
filename2 = strcat('P_Bare2', polT, polR) ;
filename01 = strcat('P_Bare01', polT, polR) ;
filename02 = strcat('P_Bare02', polT, polR) ;
writeVar(pathname, filename1, P_coh1b)
writeVar(pathname, filename2, P_coh2b)
writeVar(pathname, filename01, P0_coh1b)
writeVar(pathname, filename02, P0_coh2b)


end



%% Calculate Specular Reflection Matrix(SRM)
function  [R_sv, R_sb, r_sv, r_sb] = CalcSRM(thsd)

%% GET GLOBAL DIRECTORIES
dir_afsa = SimulationFolders.getInstance.afsa;
dir_gnd = SimulationFolders.getInstance.gnd;


%% GET GLOBAL PARAMETERS
% Vegetation Parameters
dim_layers_m = VegParams.getInstance.dim_layers_m;
num_layers = VegParams.getInstance.num_layers;


%% READ META-DATA
% Incremental Propagation Constant and look-up angles
filename = 'dKz' ;
dKz = readComplexVar(dir_afsa, filename) ;
filename = 'ANGDEG' ;
ANGDEG = readVar(dir_afsa, filename) ;
% Ground Data
filename = 'G' ;
grnd_par = readVar(dir_gnd, filename) ;
% Layer Thickness
% filename = 'D' ;
% D = readVar(dir_veg, filename) ;
% Nlayer = length(D) ;


%% CALCULATIONS
dKz_s = squeeze(dKz(:, ANGDEG == round(thsd), :)) ;

%% Transmissivity Matrix
ArgH = 0 ;
ArgV = 0 ;

for ii = 1 : num_layers
    
    ArgH = ArgH + dKz_s(1, ii) * dim_layers_m(ii, 1) ;
    ArgV = ArgV + dKz_s(2, ii) * dim_layers_m(ii, 1) ;
    
end

% vegetation transmissivity
t_sv = [exp(+1i * ArgV) 0; 0 exp(+1i * ArgH)] ;
% unity transmissivity for bare soil
t_sb = [1 0; 0 1] ;


%% GROUND REFLECTION MATRIX
h = grnd_par(1, 1) ;

epsg = grnd_par(1, 2) + 1i * grnd_par(1, 3) ;

ths = degtorad(thsd) ;
[RGHIF, RGVIF, ~, ~] = reflectionCoeff(ths, ths, epsg, h) ;

r_g = [RGVIF 0; 0 RGHIF] ;


%% SPECULAR REFLECTION MATRIX
% 2 X 2
r_sv = t_sv * r_g * t_sv ;
r_sb = t_sb * r_g * t_sb ;
% 4 X 4
R_sv = calc_Muller(r_sv) ;
R_sb = calc_Muller(r_sb) ;


end