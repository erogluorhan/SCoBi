%% Mehmet Kurum
% 03/12/2017
% modified - 11/09/2017
% modified - 11/15/2017

% Specular (Coherent) Term
function specularTerm

% Get global parameters
G0r = convertDecibelToNatural( RecParams.getInstance.G0r_dB );
EIRP = convertDecibelToNatural( SatParams.getInstance.EIRP_dB );
fMHz = SatParams.getInstance.fMHz;
VSM = GndParams.getInstance.VSM( ParamsManager.index_VSM );
RMSH = GndParams.getInstance.RMSH( ParamsManager.index_RMSH );

% TO-DO: Get global directories here 


%% Poistions
% pT: Transmitter, pS2: specular point, pR2: receiver, pG2: ground (reference), pBr2:boresight,
% pCr2: center of footprint, pSc2: center of fresnel zone
% AllPoints = [pT, pTI, pS2, pR, pRI, pG2, pBr2, pCr2, pSc2] ;
filenamex = 'AllPoints' ;
AllPoints = readVar(SimulationFolders.getInstance.config, filenamex) ;

pT = AllPoints(:, 1) ;        % Transmitter
pS2 = AllPoints(:, 3) ;       % Specular point
pR = AllPoints(:, 4) ;        % Receiver

%% Slant Range
ST = pS2 - pT ;         % Transmitter to Specular
r_st = vectorMagnitude(ST) ;  % slant range

RS = pR - pS2 ;         % Specular to Receiver
r_sr = vectorMagnitude(RS) ;  % slant range

f0hz = fMHz * Constants.MHz2Hz ;
lambda = Constants.c / f0hz ;     % Wavelength
k0 = 2 * pi * f0hz / Constants.c ;    % Wave number


%% Factor Kc

K = 1i * sqrt(EIRP) * sqrt(G0r) * lambda / (4 * pi) ;

Kc = K * exp(1i * k0 * (r_st + r_sr)) / (r_st + r_sr) ; 

%% Transmitter Pol State
e_t1 = SatParams.getInstance.e_t1 ; 
e_t2 = SatParams.getInstance.e_t2 ;
E_t1 = [1; 0; 0; 0] ;
E_t2 = [0; 0; 0; 1] ;
Q = [1 0 0 0; 0 0 0 1; 0 1 1 0; 0 -1i -1i 0] ;
E_t1 = Q * E_t1 ;
E_t2 = Q * E_t2 ;
%% Ideal Transmitter Antenna pattern
g_t = SatParams.getInstance.g_t ; % ideal
%% Ideal Receiver Antenna pattern
g_r0 = [1 0 ; 0 1] ; % ideal
%% Real Receiver Antenna Pattern
pathname = SimulationFolders.getInstance.ant_lookup ;
load([pathname '\AntPat.mat'], 'G', 'g', 'th', 'ph')

filename = 'AngS2R_rf' ;
AngS2R_rf = readVar(SimulationFolders.getInstance.config, filename) ;

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
%% Transmitter-Receiver Rotation Matrix
pathname = SimulationFolders.getInstance.rot_lookup ;
load([pathname '\u_ts.mat'], 'u_ts')
load([pathname '\u_sr.mat'], 'u_sr')
% 4 X 4
U_ts = calc_Muller(u_ts) ;
% 4 X 4
U_sr = calc_Muller(u_sr) ;

%% Scattering Specular Matrix
filename = 'AngT2S_sf' ;
AngT2S_sf = readVar(SimulationFolders.getInstance.config, filename) ;
thsd = AngT2S_sf(1, 1) ;

% r_sv for vegetation, r_sb for bare soil
[R_sv, R_sb, r_sv, r_sb] = CalcSRM(thsd) ;

%% Specular Term
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

%% save output
pathname = SimulationFolders.getInstance.out_specular;

filename3 = strcat('Kc') ;
writeComplexVar(pathname, filename3, Kc)

% 2 X 2
pathname = strcat(SimulationFolders.getInstance.out_specular, '\VSM_', num2str( VSM ), '-RMSH_', num2str( RMSH ) ) ;

filename1 = strcat('Veg1', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename2 = strcat('Veg2', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename01 = strcat('Veg01', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename02 = strcat('Veg02', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
writeComplexVar(pathname, filename1, b_coh1v)
writeComplexVar(pathname, filename2, b_coh2v)
writeComplexVar(pathname, filename01, b0_coh1v)
writeComplexVar(pathname, filename02, b0_coh2v)

filename1 = strcat('Bare1', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename2 = strcat('Bare2', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename01 = strcat('Bare01', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename02 = strcat('Bare02', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
writeComplexVar(pathname, filename1, b_coh1b)
writeComplexVar(pathname, filename2, b_coh2b)
writeComplexVar(pathname, filename01, b0_coh1b)
writeComplexVar(pathname, filename02, b0_coh2b)

% 4 X 4
filename1 = strcat('P_Veg1', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename2 = strcat('P_Veg2', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename01 = strcat('P_Veg01', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename02 = strcat('P_Veg02', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
writeVar(pathname, filename1, P_coh1v)
writeVar(pathname, filename2, P_coh2v)
writeVar(pathname, filename01, P0_coh1v)
writeVar(pathname, filename02, P0_coh2v)

filename1 = strcat('P_Bare1', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename2 = strcat('P_Bare2', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename01 = strcat('P_Bare01', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
filename02 = strcat('P_Bare02', SatParams.getInstance.polT, RecParams.getInstance.polR) ;
writeVar(pathname, filename1, P_coh1b)
writeVar(pathname, filename2, P_coh2b)
writeVar(pathname, filename01, P0_coh1b)
writeVar(pathname, filename02, P0_coh2b)
    


end

%% Calculate Specular Reflection Matrix(SRM)

function  [R_sv, R_sb, r_sv, r_sb] = CalcSRM(thsd)


%% Reading Incremental Propagation Constant
filename = 'dKz' ;
dKz = readComplexVar(SimulationFolders.getInstance.afsa, filename) ;
filename = 'ANGDEG' ;
ANGDEG = readVar(SimulationFolders.getInstance.afsa, filename) ;

dKz_s = squeeze(dKz(:, ANGDEG == round(thsd), :)) ;

%% Layer Thickness
% filename = 'D' ;
% D = readVar(FolderPath_Veg, filename) ;
% Nlayer = length(D) ;

%% Transmissivity Matrix
ArgH = 0 ;
ArgV = 0 ;

for ii = 1 : VegParams.getInstance.num_layers
    
    ArgH = ArgH + dKz_s(1, ii) * VegParams.getInstance.dim_layers(ii, 1) ;
    ArgV = ArgV + dKz_s(2, ii) * VegParams.getInstance.dim_layers(ii, 1) ;
    
end

% vegetation trasmissivity
t_sv = [exp(+1i * ArgV) 0; 0 exp(+1i * ArgH)] ;
% unity transmissivity for bare soil
t_sb = [1 0; 0 1] ;

%% Ground Reflection Matrix
filename = 'G' ;
grnd_par = readVar(SimulationFolders.getInstance.gnd, filename) ;
h = grnd_par(1, 1) ;

epsg = grnd_par(1, 2) + 1i * grnd_par(1, 3) ;

ths = degtorad(thsd) ;
[RGHIF, RGVIF, ~, ~] = reflectionCoeff(ths, ths, epsg, h) ;

r_g = [RGVIF 0; 0 RGHIF] ;

%% Specular Reflection Matrix
% 2 X 2
r_sv = t_sv * r_g * t_sv ;
r_sb = t_sb * r_g * t_sb ;
% 4 X 4
R_sv = calc_Muller(r_sv) ;
R_sb = calc_Muller(r_sb) ;

end