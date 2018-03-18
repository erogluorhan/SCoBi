%% Mehmet Kurum
% 03/12/2017
% modified - 11/09/2017

% Direct Term
function directTerm

% Get global parameters
fMHz = SatParams.getInstance.fMHz
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
% pT: Transmitter, pS2: specular point, pR2: receiver, pG2: ground (reference), pBr2:boresight,
% pCr2: center of footprint, pSc2: center of fresnel zone
% AllPoints = [pT, pTI, pS2, pR, pRI, pG2, pBr2, pCr2, pSc2] ;
filenamex = 'AllPoints' ;
AllPoints = readVar(SimulationFolders.getInstance.config, filenamex) ;

pT = AllPoints(:, 1) ;        % Transmitter
pR = AllPoints(:, 4) ;         % Receiver

%% Slant range

RT = pR - pT ;          % Satellite to Receiver
rd = vectorMagnitude(RT) ;    % slant range

f0hz = fMHz * Constants.MHz2Hz ;
lambda = Constants.c / f0hz ;     % Wavelength
k0 = 2 * pi * f0hz / Constants.c ;    % Wave number


%% Factor Kd

K = 1i * sqrt(EIRP) * sqrt(G0r) * lambda / (4 * pi) ;

Kd = K * exp(1i * k0 * rd) / rd ;

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