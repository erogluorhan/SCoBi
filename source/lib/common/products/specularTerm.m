%% Mehmet Kurum
% 03/12/2017
% modified - 11/09/2017
% modified - 11/15/2017

% Specular (Coherent) Term
function specularTerm

%% GET GLOBAL DIRECTORIES
dir_config = SimulationFolders.getInstance.config;
dir_rot_lookup = SimulationFolders.getInstance.rot_lookup ;
dir_out_specular = SimulationFolders.getInstance.out_specular;
dir_out_specular_tuple = SimulationFolders.getInstance.out_specular_tuple;


%% GET GLOBAL PARAMETERS
% Receiver Parameters
G0r_dB = RxParams.getInstance.G0r_dB;
G0r = convertDecibelToNatural( G0r_dB );
pol_Rx = RxParams.getInstance.pol_Rx;
ant_pat_Rx_id = RxParams.getInstance.ant_pat_Rx_id;
ant_pat_res_deg = RxParams.getInstance.ant_pat_res_deg;
ant_pat_struct_Rx = RxParams.getInstance.ant_pat_struct_Rx;
% Transmitter Parameters
EIRP_dB = TxParams.getInstance.EIRP_dB;
EIRP = convertDecibelToNatural( EIRP_dB );
f_MHz = TxParams.getInstance.f_MHz;
g_t = TxParams.getInstance.g_t ; % ideal
e_t1 = TxParams.getInstance.e_t1 ; 
e_t2 = TxParams.getInstance.e_t2 ;
pol_Tx = TxParams.getInstance.pol_Tx;
% Bistatic Parameters
AllPoints_m = BistaticParams.getInstance.AllPoints_m;
AngS2R_rf = BistaticParams.getInstance.AngS2R_rf; % SP->Rx Rotation Angle
AngT2S_sf = BistaticParams.getInstance.AngT2S_sf; % Tx->SP Rotation Angle


%% READ OR LOAD META-DATA
% Transmitter-Receiver Rotation Matrix
load([dir_rot_lookup '\u_ts.mat'], 'u_ts')
load([dir_rot_lookup '\u_sr.mat'], 'u_sr')


%% INITIALIZE REQUIRED PARAMETERS
% Receiver Antenna Pattern and Look-up Angles (th and ph)
g = ant_pat_struct_Rx.g;
th = ant_pat_struct_Rx.th;
ph = ant_pat_struct_Rx.ph;

% Transmitter Pol State
E_t1 = [1; 0; 0; 0] ;
E_t2 = [0; 0; 0; 1] ;
Q = [1 0 0 0; 0 0 0 1; 0 1 1 0; 0 -1i -1i 0] ;
E_t1 = Q * E_t1 ;
E_t2 = Q * E_t2 ;
% Ideal Receiver Antenna pattern
g_r0 = [1 0 ; 0 1] ; % ideal


%% CALCULATIONS
% AllPoints_m = [pos_Tx_m, pos_TxI_m, pos_SP_m, pos_Rx_m, pos_RxI_m, pos_Gnd_m, pos_B_Rx_m, pos_FP_Rx_m, pos_FZ_m] ;
pos_Tx_m = AllPoints_m(:, 1) ;        % Transmitter position
pos_SP_m = AllPoints_m(:, 3) ;       % Specular point
pos_Rx_m = AllPoints_m(:, 4) ;        % Receiver position

% Slant Range
ST = pos_SP_m - pos_Tx_m ;         % Transmitter to Specular
r_st = vectorMagnitude(ST) ;  % slant range

RS = pos_Rx_m - pos_SP_m ;         % Specular to Receiver
r_sr = vectorMagnitude(RS) ;  % slant range

f_Hz = f_MHz * Constants.MHz2Hz ;
lambda_m = Constants.c / f_Hz ;     % Wavelength
k0 = 2 * pi * f_Hz / Constants.c ;    % Wave number

% Factor Kc
K = 1i * sqrt(EIRP) * sqrt(G0r) * lambda_m / (4 * pi) ;

Kc = K * exp(1i * k0 * (r_st + r_sr)) / (r_st + r_sr) ; 

thrd = AngS2R_rf(1, 1) ;
phrd = AngS2R_rf(2, 1) ;


if ant_pat_Rx_id == Constants.id_Rx_user_defined
    
    [~, Nth] = size(th);

    % Calculate the antenna pattern resolution in degrees
    ant_pat_res_deg = Constants.ant_pat_th_range_deg / (Nth - 1);

end

ant_pat_res_factor = 1 / ant_pat_res_deg;

thd = round( ant_pat_res_factor * radtodeg(th)) / ant_pat_res_factor ; % rounding operation is due to accuracy concerns
phd = round( ant_pat_res_factor * radtodeg(ph)) / ant_pat_res_factor ; % to make the angles multiples of ant_pat_res_deg

% Receiver Antenna values in the transmitter directions
ind_th = thd == round( ant_pat_res_factor * thrd(1, 1)) / ant_pat_res_factor ; % round is to make it a multiple of ant_pat_res_deg
ind_ph = phd == round( ant_pat_res_factor * phrd(1, 1)) / ant_pat_res_factor ;

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
filename1 = strcat('Veg1', pol_Tx, pol_Rx) ;
filename2 = strcat('Veg2', pol_Tx, pol_Rx) ;
filename01 = strcat('Veg01', pol_Tx, pol_Rx) ;
filename02 = strcat('Veg02', pol_Tx, pol_Rx) ;
writeComplexVar(dir_out_specular_tuple, filename1, b_coh1v)
writeComplexVar(dir_out_specular_tuple, filename2, b_coh2v)
writeComplexVar(dir_out_specular_tuple, filename01, b0_coh1v)
writeComplexVar(dir_out_specular_tuple, filename02, b0_coh2v)

filename1 = strcat('Bare1', pol_Tx, pol_Rx) ;
filename2 = strcat('Bare2', pol_Tx, pol_Rx) ;
filename01 = strcat('Bare01', pol_Tx, pol_Rx) ;
filename02 = strcat('Bare02', pol_Tx, pol_Rx) ;
writeComplexVar(dir_out_specular_tuple, filename1, b_coh1b)
writeComplexVar(dir_out_specular_tuple, filename2, b_coh2b)
writeComplexVar(dir_out_specular_tuple, filename01, b0_coh1b)
writeComplexVar(dir_out_specular_tuple, filename02, b0_coh2b)

% 4 X 4
filename1 = strcat('P_Veg1', pol_Tx, pol_Rx) ;
filename2 = strcat('P_Veg2', pol_Tx, pol_Rx) ;
filename01 = strcat('P_Veg01', pol_Tx, pol_Rx) ;
filename02 = strcat('P_Veg02', pol_Tx, pol_Rx) ;
writeVar(dir_out_specular_tuple, filename1, P_coh1v)
writeVar(dir_out_specular_tuple, filename2, P_coh2v)
writeVar(dir_out_specular_tuple, filename01, P0_coh1v)
writeVar(dir_out_specular_tuple, filename02, P0_coh2v)

filename1 = strcat('P_Bare1', pol_Tx, pol_Rx) ;
filename2 = strcat('P_Bare2', pol_Tx, pol_Rx) ;
filename01 = strcat('P_Bare01', pol_Tx, pol_Rx) ;
filename02 = strcat('P_Bare02', pol_Tx, pol_Rx) ;
writeVar(dir_out_specular_tuple, filename1, P_coh1b)
writeVar(dir_out_specular_tuple, filename2, P_coh2b)
writeVar(dir_out_specular_tuple, filename01, P0_coh1b)
writeVar(dir_out_specular_tuple, filename02, P0_coh2b)


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
% Surface Dynamic Paramaters
h = SurfaceDynParams.getInstance.h;   % Effective roughness parameters
eps_g = SurfaceDynParams.getInstance.eps_g;   % Dielectric permittivity


%% READ META-DATA
% Incremental Propagation Constant and look-up angles
filename = 'dKz' ;
dKz = readComplexVar(dir_afsa, filename) ;
filename = 'ANGDEG' ;
ANGDEG = readVar(dir_afsa, filename) ;
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
ths = degtorad(thsd) ;
[RGHIF, RGVIF, ~, ~] = reflectionCoeff(ths, ths, eps_g, h) ;

r_g = [RGVIF 0; 0 RGHIF] ;


%% SPECULAR REFLECTION MATRIX
% 2 X 2
r_sv = t_sv * r_g * t_sv ;
r_sb = t_sb * r_g * t_sb ;
% 4 X 4
R_sv = calc_Muller(r_sv) ;
R_sb = calc_Muller(r_sb) ;


end