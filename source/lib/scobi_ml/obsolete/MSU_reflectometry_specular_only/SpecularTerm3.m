%% Mehmet Kurum
% 03/12/2017

% Specular (Coherent) Term
function SpecularTerm3

%% Global
global fMHz EIRP G0r polT polR 
global FolderPath_Config FolderPath_Ant FolderPath_rot FolderPath_out

%% Poistions
% pT: Transmitter, pS2: specular point, pR2: receiver, pG2: ground (reference), pBr2:boresight,
% pCr2: center of footprint, pSc2: center of fresnel zone
% AllPoints = [pT, pTI, pS2, pR, pRI, pG2, pBr2, pCr2, pSc2] ;
filenamex = 'AllPoints' ;
AllPoints = read_var(FolderPath_Config, filenamex) ;

pT = AllPoints(:, 1) ;        % Transmitter
pS2 = AllPoints(:, 3) ;       % Specular point
pR = AllPoints(:, 4) ;        % Receiver

%% Free space Loss
ST = pS2 - pT ;         % Transmitter to Specular
r_st = vectorMag(ST) ;  % slant range

RS = pR - pS2 ;         % Specular to Receiver
r_sr = vectorMag(RS) ;  % slant range

f0hz = fMHz * 1e6 ;
c = 3e8 ;
lambda = c / f0hz ;     % Wavelength

Ls = (4 * pi * r_st / lambda) ^ 2 ;  % free space loss

%% Factor Ks
Ks = 1i * sqrt(EIRP) * sqrt(G0r) / (r_sr * sqrt(Ls)) ; %#ok<NASGU>

%% Transmitter Pol State
e_t1 = [1; 0] ; 
e_t2 = [0; 1] ;
E_t1 = [1; 0; 0; 0] ;
E_t2 = [0; 0; 0; 1] ;
Q = [1 0 0 0; 0 0 0 1; 0 1 1 0; 0 -1i -1i 0] ;
E_t1 = Q * E_t1 ;
E_t2 = Q * E_t2 ;
%% Ideal Transmitter Antenna pattern
g_t = [1 0 ; 0 1] ; % ideal
%% Ideal Receiver Antenna pattern
g_r0 = [1 0 ; 0 1] ; % ideal
%% Real Receiver Antenna Pattern
pathname = strcat(FolderPath_Ant, '\', 'LOOK-UP') ;
load([pathname '\AntPat.mat'], 'G', 'g', 'th', 'ph')

filename = 'AngS2R_rf' ;
AngS2R_rf = read_var(FolderPath_Config, filename) ;

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
pathname = strcat(FolderPath_rot, '\', 'LOOK-UP') ;
load([pathname '\u_ts.mat'], 'u_ts')
load([pathname '\u_sr.mat'], 'u_sr')
% 4 X 4
U_ts = calc_Muller(u_ts) ;
% 4 X 4
U_sr = calc_Muller(u_sr) ;

%% Scattering Specular Matrix
filename = 'AngT2S_sf' ;
AngT2S_sf = read_var(FolderPath_Config, filename) ;
thsd = AngT2S_sf(1, 1) ;

% s_sv for vegetation, s_sb for bare soil
[sig_sv, sig_sb, s_sv, s_sb] = CalcSSM(thsd) ;

%% Specular Term
Nsm = length(s_sb) ; % Number of SM values
b_coh1v = cell(Nsm, 1) ; P_coh1v = cell(Nsm, 1) ;
b_coh2v = b_coh1v ; P_coh2v = P_coh1v ;
b0_coh1v = b_coh1v ; P0_coh1v = P_coh1v ;
b0_coh2v = b_coh1v ; P0_coh2v = P_coh1v ;

b_coh1b = cell(Nsm, 1) ; P_coh1b = cell(Nsm, 1) ;
b_coh2b = b_coh1b ; P_coh2b = P_coh1b ;
b0_coh1b = b_coh1b ; P0_coh1b = P_coh1b ;
b0_coh2b = b_coh1b ; P0_coh2b = P_coh1b ;

for ii = 1 : Nsm    % over soil moisture

    b_coh1v{ii} = g_r * u_sr * s_sv{ii} * u_ts * g_t * e_t1 ;   % field
    b_coh2v{ii} = g_r * u_sr * s_sv{ii} * u_ts * g_t * e_t2 ;   % field
    b0_coh1v{ii} = g_r0 * u_sr * s_sv{ii} * u_ts * g_t * e_t1 ;  % field
    b0_coh2v{ii} = g_r0 * u_sr * s_sv{ii} * u_ts * g_t * e_t2 ;  % field
    
    b_coh1b{ii} = g_r * u_sr * s_sb{ii} * u_ts * g_t * e_t1 ;   % field
    b_coh2b{ii} = g_r * u_sr * s_sb{ii} * u_ts * g_t * e_t2 ;   % field
    b0_coh1b{ii} = g_r0 * u_sr * s_sb{ii} * u_ts * g_t * e_t1 ;  % field
    b0_coh2b{ii} = g_r0 * u_sr * s_sb{ii} * u_ts * g_t * e_t2 ;  % field
    
    P_coh1v{ii} = G_r * U_sr * sig_sv{ii} * U_ts * G_t * E_t1 ;   % POWER
    P_coh2v{ii} = G_r * U_sr * sig_sv{ii} * U_ts * G_t * E_t2 ;   % POWER
    P0_coh1v{ii} = G_r0 * U_sr * sig_sv{ii} * U_ts * G_t * E_t1 ;  % POWER
    P0_coh2v{ii} = G_r0 * U_sr * sig_sv{ii} * U_ts * G_t * E_t2 ;  % POWER
    
    P_coh1b{ii} = G_r * U_sr * sig_sb{ii} * U_ts * G_t * E_t1 ;   % POWER
    P_coh2b{ii} = G_r * U_sr * sig_sb{ii} * U_ts * G_t * E_t2 ;   % POWER
    P0_coh1b{ii} = G_r0 * U_sr * sig_sb{ii} * U_ts * G_t * E_t1 ;  % POWER
    P0_coh2b{ii} = G_r0 * U_sr * sig_sb{ii} * U_ts * G_t * E_t2 ;  % POWER

end

%% save output
pathname = strcat(FolderPath_out, '\SPECULAR\') ;

% 2 X 2
filename1 = strcat('Veg1', polT, polR) ;
filename2 = strcat('Veg2', polT, polR) ;
filename01 = strcat('Veg01', polT, polR) ;
filename02 = strcat('Veg02', polT, polR) ;
write_cplxvar(pathname, filename1, cell2mat(b_coh1v))
write_cplxvar(pathname, filename2, cell2mat(b_coh2v))
write_cplxvar(pathname, filename01, cell2mat(b0_coh1v))
write_cplxvar(pathname, filename02, cell2mat(b0_coh2v))

filename1 = strcat('Bare1', polT, polR) ;
filename2 = strcat('Bare2', polT, polR) ;
filename01 = strcat('Bare01', polT, polR) ;
filename02 = strcat('Bare02', polT, polR) ;
write_cplxvar(pathname, filename1, cell2mat(b_coh1b))
write_cplxvar(pathname, filename2, cell2mat(b_coh2b))
write_cplxvar(pathname, filename01, cell2mat(b0_coh1b))
write_cplxvar(pathname, filename02, cell2mat(b0_coh2b))

% 4 X 4
filename1 = strcat('P_Veg1', polT, polR) ;
filename2 = strcat('P_Veg2', polT, polR) ;
filename01 = strcat('P_Veg01', polT, polR) ;
filename02 = strcat('P_Veg02', polT, polR) ;
write_var(pathname, filename1, cell2mat(P_coh1v))
write_var(pathname, filename2, cell2mat(P_coh2v))
write_var(pathname, filename01, cell2mat(P0_coh1v))
write_var(pathname, filename02, cell2mat(P0_coh2v))

filename1 = strcat('P_Bare1', polT, polR) ;
filename2 = strcat('P_Bare2', polT, polR) ;
filename01 = strcat('P_Bare01', polT, polR) ;
filename02 = strcat('P_Bare02', polT, polR) ;
write_var(pathname, filename1, cell2mat(P_coh1b))
write_var(pathname, filename2, cell2mat(P_coh2b))
write_var(pathname, filename01, cell2mat(P0_coh1b))
write_var(pathname, filename02, cell2mat(P0_coh2b))

end

%% Calculate Scattering Specular Matrix(SSM)

function  [sig_sv, sig_sb, s_sv, s_sb] = CalcSSM(thsd)

%% Global
global FolderPath_Veg FolderPath_afsa FolderPath_gnd

%% Reading Incremental Propagation Constant
filename = 'dKz' ;
dKz = read_cplxvar(FolderPath_afsa, filename) ;
filename = 'ANGDEG' ;
ANGDEG = read_var(FolderPath_afsa, filename) ;

dKz_s = squeeze(dKz(:, ANGDEG == round(thsd), :)) ;

%% Layer Thickness
filename = 'D' ;
D = read_var(FolderPath_Veg, filename) ;
Nlayer = length(D) ;

%% Transmissivity Matrix
ArgH = 0 ;
ArgV = 0 ;

for ii = 1 : Nlayer
    
    ArgH = ArgH + dKz_s(1, ii) * D(ii, 1) ;
    ArgV = ArgV + dKz_s(2, ii) * D(ii, 1) ;
    
end

% vegetation trasmissivity
t_sv = [exp(+1i * ArgV) 0; 0 exp(+1i * ArgH)] ;
% unity transmissivity for bare soil
t_sb = [1 0; 0 1] ;

%% Reflection Matrix
filename = 'G' ;
grnd_par = read_var(FolderPath_gnd, filename) ;
h = grnd_par(1, 1) ;

Nsm = (length(grnd_par) - 1) / 2 ;  % Number of SM values
epsg = zeros(Nsm, 1) ;

for ii = 1 : Nsm
    epsg(ii, 1) = grnd_par(1, ii + 1) + 1i * grnd_par(1, ii + Nsm + 1) ;
end

ths = degtorad(thsd) ;
[RGHIF, RGVIF, ~, ~] = refcoeff(ths, ths, epsg, h) ;

r_s = cell(Nsm, 1) ;
for ii = 1 : Nsm
    r_s{ii} = [RGVIF(ii) 0; 0 RGHIF(ii)] ;
end

%% Scattering Specular Matrix
s_sv = cell(Nsm, 1) ;
s_sb = cell(Nsm, 1) ;
sig_sv = cell(Nsm, 1) ;
sig_sb = cell(Nsm, 1) ;

for ii = 1 : Nsm
    
    % 2 X 2
    s_sv{ii} = t_sv * r_s{ii} * t_sv ;
    s_sb{ii} = t_sb * r_s{ii} * t_sb ;
    % 4 X 4
    sig_sv{ii} = 4 * pi * calc_Muller(s_sv{ii}) ;
    sig_sb{ii} = 4 * pi * calc_Muller(s_sb{ii}) ;
    
end


end

%% Function that calculates magnitude of the given vector

function abs = vectorMag(vec)

[m,n]=size(vec);
if (m~=1)&&(n~=1)  % or unit colomn or unit row
    abs = 0;
    disp('Error - vector is not of proper dimensions');
else
    abs = sqrt(sum (vec.^2));
end;

end

function UU = calc_Muller(uu)

u11 = uu(1, 1) ; u12 = uu(1, 2) ;
u21 = uu(2, 1) ; u22 = uu(2, 2) ;

U11 = abs(u11) .^ 2 ;
U12 = abs(u12) .^ 2 ;
U13 = real(u11 .* conj(u12)) ;
U14 = -imag(u11 .* conj(u12)) ;

U21 = abs(u21) .^ 2 ;
U22 = abs(u22) .^ 2 ;
U23 = real(u21 .* conj(u22)) ;
U24 = -imag(u21 .* conj(u22)) ;

U31 = 2 * real(u11 .* conj(u21)) ; 
U32 = 2 * real(u12 .* conj(u22)) ; 
U33 = real(u11 .* conj(u22) + u12 .* conj(u21)) ;
U34 = -imag(u11 .* conj(u22) - u12 .* conj(u21)) ;

U41 = 2 * imag(u11 .* conj(u21)) ; 
U42 = 2 * imag(u12 .* conj(u22)) ; 
U43 = imag(u11 .* conj(u22) + u12 .* conj(u21)) ;
U44 = real(u11 .* conj(u22) - u12 .* conj(u21)) ;

UU = [U11, U12, U13, U14; ...
    U21, U22, U23, U24; ...
    U31, U32, U33, U34; ...
    U41, U42, U43, U44] ;

end