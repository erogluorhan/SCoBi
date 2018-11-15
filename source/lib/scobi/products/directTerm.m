
function directTerm
% function directTerm 
%
%   Calculates the direct (line-of-sight) received field and power between 
%   the transmitter and the receiver. Stores the calculated values into simulation output 
%   folders in an incremental fashion as the simulation iterations continue.
%
%   See also mainSCoBi, specularTerm.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



%% GET GLOBAL DIRECTORIES
dir_products = SimulationFolders.getInstance.products;
dir_products_direct = SimulationFolders.getInstance.products_direct;
dir_products_direct_field = SimulationFolders.getInstance.products_direct_field;
dir_products_direct_power = SimulationFolders.getInstance.products_direct_power;


%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Simulation Settings
sim_mode_id = SimSettings.getInstance.sim_mode_id;
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
EIRP_dB = TxParams.getInstance.EIRP_dB;
EIRP = convertDecibelToNatural( EIRP_dB );
g_t = TxParams.getInstance.g_t;
e_t1 = TxParams.getInstance.e_t1;
e_t2 = TxParams.getInstance.e_t2;
% Receiver Parameters
G0r_dB = RxParams.getInstance.G0r_dB;
G0r = convertDecibelToNatural( G0r_dB );
ant_pat_struct_Rx = RxParams.getInstance.ant_pat_struct_Rx;
ant_pat_res_deg = RxParams.getInstance.ant_pat_res_deg;
% Configuration Parameters
DoYs = ConfigParams.getInstance.DoYs;
% Bistatic Dynamic Parameters
AllPoints_m = BistaticDynParams.getInstance.AllPoints_m;
AngT2R_rf = BistaticDynParams.getInstance.AngT2R_rf;
% Rotation Dynamic Parameters
u_tr = RotMatDynParams.getInstance.u_tr;


%% INITIALIZE REQUIRED PARAMETERS
% AllPoints_m = [pos_Tx_m, pos_TxI_m, pos_SP_m, pos_Rx_m, pos_RxI_m, pos_Gnd_m, pos_B_Rx_m, pos_FZ_m] ;
pos_Tx_m = AllPoints_m(:, 1) ;        % Transmitter
pos_Rx_m = AllPoints_m(:, 4) ;         % Receiver

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
% Slant range
RT_m = pos_Rx_m - pos_Tx_m ;          % Transmitter to Receiver
rd_m = vectorMagnitude(RT_m) ;    % slant range

f_Hz = f_MHz * Constants.MHZ_TO_HZ ;
lambda_m = Constants.LIGHTSPEED / f_Hz ;     % Wavelength
k0 = 2 * pi * f_Hz / Constants.LIGHTSPEED ;    % Wave number


% Factor Kd
K_m = 1i * sqrt(EIRP) * sqrt(G0r) * lambda_m / (4 * pi) ;

Kd = K_m * exp(1i * k0 * rd_m) / rd_m ;

thrd = AngT2R_rf(1, 1) ;
phrd = AngT2R_rf(2, 1) ;

ant_pat_res_factor = 1 / ant_pat_res_deg;

thd = round( ant_pat_res_factor * rad2deg(th)) / ant_pat_res_factor ; % rounding operation is due to accuracy concerns
phd = round( ant_pat_res_factor * rad2deg(ph)) / ant_pat_res_factor ; % to make the angles multiples of ant_pat_res_deg

% Receiver Antenna values in the transmitter direction
ind_th = thd == round( ant_pat_res_factor * thrd(1, 1)) / ant_pat_res_factor ; % round is to make it a multiple of ant_pat_res_deg
ind_ph = phd == round( ant_pat_res_factor * phrd(1, 1)) / ant_pat_res_factor ;

g11 = cell2mat(g(1,1)) ; g12 = cell2mat(g(1,2)) ;
g21 = cell2mat(g(2,1)) ; g22 = cell2mat(g(2,2)) ;

g_r = [g11(ind_th & ind_ph), g12(ind_th & ind_ph); ...
    g21(ind_th & ind_ph), g22(ind_th & ind_ph)] ;

% 4 X 4
G_r = calcMuller(g_r) ;
% 4 X 4
G_r0 = calcMuller(g_r0) ;
% 4 X 4
G_t = calcMuller(g_t) ;

% Transmitter-Receiver Rotation  Matrix
% 4 X 4
U_tr = calcMuller(u_tr) ;


%% DIRECT TERM
% Field
b_d1 = g_r * u_tr * g_t * e_t1 ;
b_d2 = g_r * u_tr * g_t * e_t2 ;
b0_d1 = g_r0 * u_tr * g_t * e_t1 ;
b0_d2 = g_r0 * u_tr * g_t * e_t2 ;

% Power
P_d1 = G_r * U_tr * G_t * E_t1 ;
P_d2 = G_r * U_tr * G_t * E_t2 ;
P0_d1 = G_r0 * U_tr * G_t * E_t1 ;
P0_d2 = G_r0 * U_tr * G_t * E_t2 ;


%% SAVE
% If sim_mode is Time-series, write DoYs
if sim_mode_id == Constants.ID_TIME_SERIES
    
    DoY = DoYs( sim_counter );
    filename = 'DoYs';
    writeVarIncremental( dir_products, filename, sim_counter, DoY )
    
end

% Kd factor
filename = strcat('Kd') ;
writeComplexVar(dir_products_direct, filename, Kd);

filename1 = 'Dir1';
filename2 = 'Dir2';
filename01 = 'Dir01';
filename02 = 'Dir02';

% Field: 2 X 1
writeComplexVarIncremental(dir_products_direct_field, filename1, sim_counter, b_d1);
writeComplexVarIncremental(dir_products_direct_field, filename2, sim_counter, b_d2);
writeComplexVarIncremental(dir_products_direct_field, filename01, sim_counter, b0_d1);
writeComplexVarIncremental(dir_products_direct_field, filename02, sim_counter, b0_d2);

% Power: 4 X 1
writeVarIncremental(dir_products_direct_power, filename1, sim_counter, P_d1);
writeVarIncremental(dir_products_direct_power, filename2, sim_counter, P_d2);
writeVarIncremental(dir_products_direct_power, filename01, sim_counter, P0_d1);
writeVarIncremental(dir_products_direct_power, filename02, sim_counter, P0_d2);


end