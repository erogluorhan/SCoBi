%% Mehmet Kurum
% 09/24/2017

function [r0_coh1b, r0_coh2b] = SpecularReflection(r_s)

%% GLOBAL DIRECTORIES
dir_config = SimulationFolders.getInstance.config;
dir_ant_lookup = SimulationFolders.getInstance.ant_lookup;
dir_rot_lookup = SimulationFolders.getInstance.rot_lookup;
dir_afsa = SimulationFolders.getInstance.afsa;

%% GET GLOBAL PARAMETERS
dim_layers_m = VegParams.getInstance.dim_layers_m;
num_layers = VegParams.getInstance.num_layers;

%% Transmitter Pol State
e_t1 = [1; 0] ; 
e_t2 = [0; 1] ;

%% Ideal Transmitter Antenna pattern
g_t = [1 0 ; 0 1] ; % ideal
%% Ideal Receiver Antenna pattern
g_r0 = [1 0 ; 0 1] ; % ideal
%% Real Receiver Antenna Pattern
load([dir_ant_lookup '\AntPat.mat'], 'G', 'g', 'th', 'ph')

filename = 'AngS2R_rf' ;
AngS2R_rf = readVar(dir_config, filename) ;

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

%% Transmitter-Receiver Rotation Matrix
load([dir_rot_lookup '\u_ts.mat'], 'u_ts')
load([dir_rot_lookup '\u_sr.mat'], 'u_sr')


%% Scattering Specular Matrix
filename = 'AngT2S_sf' ;
AngT2S_sf = readVar(dir_config, filename) ;
thsd = AngT2S_sf(1, 1) ;


%% Reading Incremental Propagation Constant
filename = 'dKz' ;
dKz = readComplexVar(dir_afsa, filename) ;
filename = 'ANGDEG' ;
ANGDEG = readVar(dir_afsa, filename) ;

dKz_s = squeeze(dKz(:, ANGDEG == round(thsd), :)) ;

%% Transmissivity Matrix
ArgH = 0 ;
ArgV = 0 ;

for ii = 1 : num_layers
    
    ArgH = ArgH + dKz_s(1, ii) * dim_layers_m(ii, 1) ;
    ArgV = ArgV + dKz_s(2, ii) * dim_layers_m(ii, 1) ;
    
end

% vegetation trasmissivity
t_sv = [exp(+1i * ArgV) 0; 0 exp(+1i * ArgH)] ;
% unity transmissivity for bare soil
t_sb = [1 0; 0 1] ;
% r_s = [1 0; 0 1] ;


%% Specular Term
r_coh1v = g_r * u_sr * t_sv * r_s * t_sv * u_ts * g_t * e_t1 ;   % field
r_coh2v = g_r * u_sr * t_sv * r_s * t_sv * u_ts * g_t * e_t2 ;   % field

r0_coh1v = g_r0 * u_sr * t_sv * r_s * t_sv * u_ts * g_t * e_t1 ;   % field ;
r0_coh2v = g_r0 * u_sr * t_sv * r_s * t_sv * u_ts * g_t * e_t2 ;   % field ;

r_coh1b = g_r * u_sr * t_sb * r_s * t_sb * u_ts * g_t * e_t1 ;   % field
r_coh2b = g_r * u_sr * t_sb * r_s * t_sb * u_ts * g_t * e_t2 ;   % field

r0_coh1b = g_r0 * u_sr * t_sb * r_s * t_sb * u_ts * g_t * e_t1 ;   % field ;
r0_coh2b = g_r0 * u_sr * t_sb * r_s * t_sb * u_ts * g_t * e_t2 ;   % field ;

end

