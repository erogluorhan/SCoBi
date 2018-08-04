%% Mehmet Kurum
% 03/12/2017
% modified - 11/09/2017

% Direct Term
function directTerm

%% GET GLOBAL DIRECTORIES
dir_config = SimulationFolders.getInstance.config;
dir_rot_lookup = SimulationFolders.getInstance.rot_lookup;
dir_out_direct = SimulationFolders.getInstance.out_direct;


%% GET GLOBAL PARAMETERS
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
EIRP_dB = TxParams.getInstance.EIRP_dB;
EIRP = convertDecibelToNatural( EIRP_dB );
g_t = TxParams.getInstance.g_t;
e_t1 = TxParams.getInstance.e_t1;
e_t2 = TxParams.getInstance.e_t2;
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
G0r_dB = RxParams.getInstance.G0r_dB;
G0r = convertDecibelToNatural( G0r_dB );
ant_pat_struct_Rx = RxParams.getInstance.ant_pat_struct_Rx;


%% READ OR LOAD META-DATA
% pos_FP_Rx_m: center of footprint, pos_FZ_m: center of fresnel zone
% AllPoints_m = [pos_Tx_m, pos_TxI_m, pos_SP_m, pos_Rx_m, pos_RxI_m, pos_Gnd_m, pos_B_Rx_m, pos_FP_Rx_m, pos_FZ_m] ;
filenamex = 'AllPoints_m' ;
AllPoints_m = readVar( dir_config, filenamex );
pos_Tx_m = AllPoints_m(:, 1) ;        % Transmitter
pos_Rx_m = AllPoints_m(:, 4) ;         % Receiver

filename = 'AngT2R_rf' ;
AngT2R_rf = readVar( dir_config, filename) ;

% Transmitter-Receiver Rotation Matrix
load([dir_rot_lookup '\u_tr.mat'], 'u_tr')


% Receiver Antenna Pattern and Look-up Angles (th and ph)
g = ant_pat_struct_Rx.g;
th = ant_pat_struct_Rx.th;
ph = ant_pat_struct_Rx.ph;


%% CALCULATIONS

% Slant range
RT_m = pos_Rx_m - pos_Tx_m ;          % Transmitter to Receiver
rd_m = vectorMagnitude(RT_m) ;    % slant range

f_Hz = f_MHz * Constants.MHz2Hz ;
lambda_m = Constants.c / f_Hz ;     % Wavelength
k0 = 2 * pi * f_Hz / Constants.c ;    % Wave number


% Factor Kd
K_m = 1i * sqrt(EIRP) * sqrt(G0r) * lambda_m / (4 * pi) ;

Kd = K_m * exp(1i * k0 * rd_m) / rd_m ;

thrd = AngT2R_rf(1, 1) ;
phrd = AngT2R_rf(2, 1) ;

thd = round(2 * th * Constants.rad2deg) / 2 ; % rounding operation is due to accuracy concerns
phd = round(2 * ph * Constants.rad2deg) / 2 ;

% Receiver Antenna values in the transmitter directions
ind_th = thd == round(2 * thrd(1, 1)) / 2 ; % round is to make it a multiple of 0.5
ind_ph = phd == round(2 * phrd(1, 1)) / 2 ;
g11 = cell2mat(g(1,1)) ; g12 = cell2mat(g(1,2)) ;
g21 = cell2mat(g(2,1)) ; g22 = cell2mat(g(2,2)) ;
g_r = [g11(ind_th & ind_ph), g12(ind_th & ind_ph); ...
    g21(ind_th & ind_ph), g22(ind_th & ind_ph)] ;


%% DIRECT TERM
b_d1 = g_r * u_tr * g_t * e_t1 ;     % field
b_d2 = g_r * u_tr * g_t * e_t2 ;     % field

% P_d = abs(b_d) .^ 2 ;                  % Power
% P_d_dB = 10 * log10(P_d) ;            % Power in dB


%% SAVE
filename1 = strcat('Dir1', pol_Tx, pol_Rx) ;
writeComplexVar(dir_out_direct, filename1, b_d1)

filename2 = strcat('Dir2', pol_Tx, pol_Rx) ;
writeComplexVar(dir_out_direct, filename2, b_d2)

filename3 = strcat('Kd') ;
writeComplexVar(dir_out_direct, filename3, Kd)

end