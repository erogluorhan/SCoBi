%% Mehmet Kurum
% 09/24/2017

function [r0_coh1b, r0_coh2b] = SpecularReflection(r_s)


%% GLOBAL DIRECTORIES
dir_rot_lookup = SimulationFolders.getInstance.rot_lookup;
dir_afsa = SimulationFolders.getInstance.afsa;


%% GET GLOBAL PARAMETERS
% Simulation Settings
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
% Vegetation Parameters
if gnd_cover_id == Constants.id_veg_cover
    dim_layers_m = VegParams.getInstance.dim_layers_m;
    num_layers = VegParams.getInstance.num_layers;
end
% Transmitter Parameters
g_t = TxParams.getInstance.g_t ;    % Ideal Transmitter Antenna pattern
e_t1 = TxParams.getInstance.e_t1 ;  % Transmitter Pol State 
e_t2 = TxParams.getInstance.e_t2 ;
% Receiver Parameters
ant_pat_Rx_id = RxParams.getInstance.ant_pat_Rx_id;
ant_pat_res_deg = RxParams.getInstance.ant_pat_res_deg;
ant_pat_struct_Rx = RxParams.getInstance.ant_pat_struct_Rx;
% Bistatic Parameters
AngS2R_rf = BistaticParams.getInstance.AngS2R_rf; % SP->Rx Rotation Angle
AngT2S_sf = BistaticParams.getInstance.AngT2S_sf; % Tx->SP Rotation Angle


%% INITIALIZE REQUIRED PARAMETERS
% Ideal Receiver Antenna pattern
g_r0 = [1 0 ; 0 1] ; % ideal


%% READ OR LOAD META-DATA
%% Load Transmitter-Receiver Rotation Matrix
load([dir_rot_lookup '\u_ts.mat'], 'u_ts')
load([dir_rot_lookup '\u_sr.mat'], 'u_sr')
 


% Read Incremental Propagation Constant
filename = 'dKz' ;
dKz = readComplexVar(dir_afsa, filename) ;
filename = 'ANGDEG' ;
ANGDEG = readVar(dir_afsa, filename) ;


% Receiver Antenna Pattern and Look-up Angles (th and ph)
g = ant_pat_struct_Rx.g;
th = ant_pat_struct_Rx.th;
ph = ant_pat_struct_Rx.ph;


%% CALCULATIONS
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


thsd = AngT2S_sf(1, 1) ;

% Read Incremental Propagation Constant
dKz_s = squeeze(dKz(:, ANGDEG == round(thsd), :)) ;

%% Transmissivity Matrix
ArgH = 0 ;
ArgV = 0 ;

if gnd_cover_id == Constants.id_veg_cover

    for ii = 1 : num_layers

            ArgH = ArgH + dKz_s(1, ii) * dim_layers_m(ii, 1) ;
            ArgV = ArgV + dKz_s(2, ii) * dim_layers_m(ii, 1) ;

    end

end

% vegetation trasmissivity
t_sv = [exp(+1i * ArgV) 0; 0 exp(+1i * ArgH)] ;
% unity transmissivity for bare soil
t_sb = [1 0; 0 1] ;
% r_s = [1 0; 0 1] ;

% TO-DO: Solve for vegetated covers

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

