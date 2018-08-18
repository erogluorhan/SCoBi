%% Mehmet Kurum
% 03/12/2017
% modified - 11/09/2017
% modified - 11/15/2017

% Specular (Coherent) Term
function specularTerm


%% GET GLOBAL DIRECTORIES
dir_rot_lookup = SimulationFolders.getInstance.rot_lookup ;
dir_products = SimulationFolders.getInstance.products;
dir_products_specular = SimulationFolders.getInstance.products_specular;
dir_products_specular_field = SimulationFolders.getInstance.products_specular_field;
dir_products_specular_power = SimulationFolders.getInstance.products_specular_power;


%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Simulation Settings
sim_mode_id = SimSettings.getInstance.sim_mode_id;
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
% Receiver Parameters
G0r_dB = RxParams.getInstance.G0r_dB;
G0r = convertDecibelToNatural( G0r_dB );
ant_pat_res_deg = RxParams.getInstance.ant_pat_res_deg;
ant_pat_struct_Rx = RxParams.getInstance.ant_pat_struct_Rx;
% Transmitter Parameters
EIRP_dB = TxParams.getInstance.EIRP_dB;
EIRP = convertDecibelToNatural( EIRP_dB );
f_MHz = TxParams.getInstance.f_MHz;
g_t = TxParams.getInstance.g_t ; % ideal
e_t1 = TxParams.getInstance.e_t1 ; 
e_t2 = TxParams.getInstance.e_t2 ;
% Bistatic Parameters
AllPoints_m = BistaticDynParams.getInstance.AllPoints_m;
AngS2R_rf = BistaticDynParams.getInstance.AngS2R_rf; % SP->Rx Rotation Angle
% Configuration Parameters
DoYs = ConfigParams.getInstance.DoYs;
% Ground Parameters
num_gnd_layers = GndParams.getInstance.num_layers;


%% READ OR LOAD META-DATA
% Transmitter-Receiver Rotation Matrix
load([dir_rot_lookup '\u_ts.mat'], 'u_ts')
load([dir_rot_lookup '\u_sr.mat'], 'u_sr')


%% INITIALIZE REQUIRED PARAMETERS
% AllPoints_m = [pos_Tx_m, pos_TxI_m, pos_SP_m, pos_Rx_m, pos_RxI_m, pos_Gnd_m, pos_B_Rx_m, pos_FP_Rx_m, pos_FZ_m] ;
pos_Tx_m = AllPoints_m(:, 1) ;        % Transmitter position
pos_SP_m = AllPoints_m(:, 3) ;       % Specular point
pos_Rx_m = AllPoints_m(:, 4) ;        % Receiver position

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

ant_pat_res_factor = 1 / ant_pat_res_deg;

thd = round( ant_pat_res_factor * radtodeg(th)) / ant_pat_res_factor ; % rounding operation is due to accuracy concerns
phd = round( ant_pat_res_factor * radtodeg(ph)) / ant_pat_res_factor ; % to make the angles multiples of ant_pat_res_deg

% Receiver Antenna values in the specular direction
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

[R_sv, R_sb, r_sv, r_sb] = CalcSRM() ; % r_sv for vegetation, r_sb for bare soil


%% SPECULAR TERM
% For single ground layer
num_diel_profiles = 1;

% For multiplw ground layers
if num_gnd_layers > 1
    
    % Number of dielectric profiles depends on whether simulator is SCoBi-Veg
    % or SCoBi-ML. SCoBi-ML has several diel profiles.
    [~, num_diel_profiles] = size( Constants.diel_profiles );
    
end

% P = 4 x 1
% b = 2 x 1

for ii = 1 : num_diel_profiles
    
    % Field
    % Vegetation
    b_coh1v(:,ii) = g_r * u_sr * r_sv{ii,1} * u_ts * g_t * e_t1 ;
    b_coh2v(:,ii) = g_r * u_sr * r_sv{ii,1} * u_ts * g_t * e_t2 ;
    b0_coh1v(:,ii) = g_r0 * u_sr * r_sv{ii,1} * u_ts * g_t * e_t1 ;
    b0_coh2v(:,ii) = g_r0 * u_sr * r_sv{ii,1} * u_ts * g_t * e_t2 ;

    % Bare-soil
    b_coh1b(:,ii) = g_r * u_sr * r_sb{ii,1} * u_ts * g_t * e_t1 ;
    b_coh2b(:,ii) = g_r * u_sr * r_sb{ii,1} * u_ts * g_t * e_t2 ;
    b0_coh1b(:,ii) = g_r0 * u_sr * r_sb{ii,1} * u_ts * g_t * e_t1 ;
    b0_coh2b(:,ii) = g_r0 * u_sr * r_sb{ii,1} * u_ts * g_t * e_t2 ;

    % Power
    % Vegetation
    P_coh1v(:,ii) = G_r * U_sr * R_sv{ii,1} * U_ts * G_t * E_t1 ;
    P_coh2v(:,ii) = G_r * U_sr * R_sv{ii,1} * U_ts * G_t * E_t2 ;
    P0_coh1v(:,ii) = G_r0 * U_sr * R_sv{ii,1} * U_ts * G_t * E_t1 ;
    P0_coh2v(:,ii) = G_r0 * U_sr * R_sv{ii,1} * U_ts * G_t * E_t2 ;

    % Bare-soil
    P_coh1b(:,ii) = G_r * U_sr * R_sb{ii,1} * U_ts * G_t * E_t1 ;
    P_coh2b(:,ii) = G_r * U_sr * R_sb{ii,1} * U_ts * G_t * E_t2 ;
    P0_coh1b(:,ii) = G_r0 * U_sr * R_sb{ii,1} * U_ts * G_t * E_t1 ;
    P0_coh2b(:,ii) = G_r0 * U_sr * R_sb{ii,1} * U_ts * G_t * E_t2 ;

end


%% SAVE OUTPUTS
% Start and end indices of DoY variable to be appended to the
% cumulative (incremental) direct term are the same, because number of
% dielectric profiles is not effective on direct term
start_index = sim_counter;
end_index = sim_counter;
if sim_mode_id == Constants.id_time_series
    
    DoY = DoYs( sim_counter );    
    filename = 'DoYs';
    writeVarIncremental( dir_products, filename, start_index, end_index, DoY )
    
end

% Kc Factor
filename3 = 'Kc';
writeComplexVar(dir_products_specular, filename3, Kc)

% Calculate start and end indices of each variable to be appended to the
% cumulative (incremental) specular term
start_index = 1 + (sim_counter-1) * num_diel_profiles;
end_index = sim_counter * num_diel_profiles;

% Field: 2 X 1
if gnd_cover_id == Constants.id_veg_cover
    
    filename1 = 'Veg1';
    filename2 = 'Veg2';
    filename01 = 'Veg01';
    filename02 = 'Veg02';
    writeComplexVarIncremental(dir_products_specular_field, filename1, start_index, end_index, b_coh1v)
    writeComplexVarIncremental(dir_products_specular_field, filename2, start_index, end_index, b_coh2v)
    writeComplexVarIncremental(dir_products_specular_field, filename01, start_index, end_index, b0_coh1v)
    writeComplexVarIncremental(dir_products_specular_field, filename02, start_index, end_index, b0_coh2v)
    
end

filename1 = 'Bare1';
filename2 = 'Bare2';
filename01 = 'Bare01';
filename02 = 'Bare02';
writeComplexVarIncremental(dir_products_specular_field, filename1, start_index, end_index, b_coh1b)
writeComplexVarIncremental(dir_products_specular_field, filename2, start_index, end_index, b_coh2b)
writeComplexVarIncremental(dir_products_specular_field, filename01, start_index, end_index, b0_coh1b)
writeComplexVarIncremental(dir_products_specular_field, filename02, start_index, end_index, b0_coh2b)

% Power: 4 X 1
if gnd_cover_id == Constants.id_veg_cover
    
    filename1 = 'Veg1';
    filename2 = 'Veg2';
    filename01 = 'Veg01';
    filename02 = 'Veg02';
    writeVarIncremental(dir_products_specular_power, filename1, start_index, end_index, P_coh1v);
    writeVarIncremental(dir_products_specular_power, filename2, start_index, end_index, P_coh2v);
    writeVarIncremental(dir_products_specular_power, filename01, start_index, end_index, P0_coh1v);
    writeVarIncremental(dir_products_specular_power, filename02, start_index, end_index, P0_coh2v);
    
end

filename1 = 'Bare1';
filename2 = 'Bare2';
filename01 = 'Bare01';
filename02 = 'Bare02';
writeVarIncremental(dir_products_specular_power, filename1, start_index, end_index, P_coh1b);
writeVarIncremental(dir_products_specular_power, filename2, start_index, end_index, P_coh2b);
writeVarIncremental(dir_products_specular_power, filename01, start_index, end_index, P0_coh1b);
writeVarIncremental(dir_products_specular_power, filename02, start_index, end_index, P0_coh2b);


end



%% Calculate Specular Reflection Matrix(SRM)
function  [R_sv, R_sb, r_sv, r_sb] = CalcSRM()

%% GET GLOBAL DIRECTORIES
dir_afsa = SimulationFolders.getInstance.afsa;


%% GET GLOBAL PARAMETERS
% Vegetation Parameters
dim_layers_m = VegParams.getInstance.dim_layers_m;
num_veg_layers = VegParams.getInstance.num_layers;
% Ground Parameters
num_gnd_layers = GndParams.getInstance.num_layers;
% Ground Dynamic Paramaters
h = GndDynParams.getInstance.h;   % Effective roughness parameters
eps_g = GndDynParams.getInstance.eps_g;   % Dielectric permittivity
% Bistatic Parameters
AngT2S_sf = BistaticDynParams.getInstance.AngT2S_sf; % Tx->SP Rotation Angle


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
% Specular Reflection Matrix
thsd = AngT2S_sf(1, 1) ;

dKz_s = squeeze(dKz(:, ANGDEG == round(thsd), :)) ;

%% Transmissivity Matrix
ArgH = 0 ;
ArgV = 0 ;

for ii = 1 : num_veg_layers
    
    ArgH = ArgH + dKz_s(1, ii) * dim_layers_m(ii, 1) ;
    ArgV = ArgV + dKz_s(2, ii) * dim_layers_m(ii, 1) ;
    
end

% vegetation transmissivity
t_sv = [exp(+1i * ArgV) 0; 0 exp(+1i * ArgH)] ;
% unity transmissivity for bare soil
t_sb = [1 0; 0 1] ;


%% GROUND REFLECTION MATRIX
% If single-layered ground
if num_gnd_layers == 1
    
    ths = degtorad(thsd) ;
    [RGHIF, RGVIF, ~, ~] = reflectionCoeffSingle(ths, ths, eps_g, h) ;

    r_g = [RGVIF 0; 0 RGHIF] ;


    %% SPECULAR REFLECTION MATRIX
    % 2 X 2
    r_sv{1,1} = t_sv * r_g * t_sv ;
    r_sb{1,1} = t_sb * r_g * t_sb ;
    % 4 X 4
    R_sv{1,1} = calc_Muller(r_sv{1,1}) ;
    R_sb{1,1} = calc_Muller(r_sb{1,1}) ;

% Else if multi-layered ground
else
    
    [r_g_cell] = reflectionCoeffsML();
    
    num_diel_profiles = length(r_g_cell);
    
    for ii = 1 : num_diel_profiles
        
        r_sv{ii,1} = t_sv * r_g_cell{ii,1} * t_sv ;
        r_sb{ii,1} = t_sb * r_g_cell{ii,1} * t_sb ;
        % 4 X 4
        R_sv{ii,1} = calc_Muller(r_sv{ii,1});
        R_sb{ii,1} = calc_Muller(r_sb{ii,1});
        
    end

end


end