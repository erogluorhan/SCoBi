%% Mehmet Kurum
% 03/12/2017
% modified - 11/09/2017

% Direct Term
function directTerm

% Get global parameters
f_MHz = SatParams.getInstance.f_MHz
EIRP_dB = SatParams.getInstance.EIRP_dB;
g_t = SatParams.getInstance.g_t;
e_t1 = SatParams.getInstance.e_t1;
e_t2 = SatParams.getInstance.e_t2;
polT = SatParams.getInstance.polT;
polR = RecParams.getInstance.polR;
G0r = convertDecibelToNatural( RecParams.getInstance.G0r_dB );
EIRP = convertDecibelToNatural( EIRP_dB );

% TO-DO: Get global directories here


%%
% pT_m: Transmitter, pS2_m: specular point, pR2: receiver, pG2_m: ground (reference), pBr2_m:boresight,
% pCr2_m: center of footprint, pSc2_m: center of fresnel zone
% AllPoints_m = [pT_m, pTI_m, pS2_m, pR_m, pRI_m, pG2_m, pBr2_m, pCr2_m, pSc2_m] ;
filenamex = 'AllPoints_m' ;
AllPoints_m = readVar(SimulationFolders.getInstance.config, filenamex) ;

pT_m = AllPoints_m(:, 1) ;        % Transmitter
pR_m = AllPoints_m(:, 4) ;         % Receiver

%% Slant range

RT_m = pR_m - pT_m ;          % Satellite to Receiver
rd_m = vectorMagnitude(RT_m) ;    % slant range

f_Hz = f_MHz * Constants.MHz2Hz ;
lambda_m = Constants.c / f_Hz ;     % Wavelength
k0 = 2 * pi * f_Hz / Constants.c ;    % Wave number


%% Factor Kd

K_m = 1i * sqrt(EIRP) * sqrt(G0r) * lambda_m / (4 * pi) ;

Kd = K_m * exp(1i * k0 * rd_m) / rd_m ;

%% Receiver Antenna Pattern
pathname = SimulationFolders.getInstance.ant_lookup ;
load([pathname '\AntPat.mat'], 'g', 'th', 'ph')

filename = 'AngT2R_rf' ;
AngT2R_rf = readVar(SimulationFolders.getInstance.config, filename) ;

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


%% Transmitter-Receiver Rotation Matrix

pathname = SimulationFolders.getInstance.rot_lookup ;
load([pathname '\u_tr.mat'], 'u_tr')


%% Direct Term

b_d1 = g_r * u_tr * g_t * e_t1 ;     % field
b_d2 = g_r * u_tr * g_t * e_t2 ;     % field

% P_d = abs(b_d) .^ 2 ;                  % Power
% P_d_dB = 10 * log10(P_d) ;            % Power in dB

% save output
pathname = SimulationFolders.getInstance.out_direct;

filename1 = strcat('Dir1', polT, polR) ;
filename2 = strcat('Dir2', polT, polR) ;
filename3 = strcat('Kd') ;


writeComplexVar(pathname, filename1, b_d1)
writeComplexVar(pathname, filename2, b_d2)
writeComplexVar(pathname, filename3, Kd)



end