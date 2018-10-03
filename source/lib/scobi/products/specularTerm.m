%% Mehmet Kurum
% 03/12/2017
% modified - 11/09/2017
% modified - 11/15/2017

% Specular (Coherent) Term
function specularTerm


%% GET GLOBAL DIRECTORIES
dir_products = SimulationFolders.getInstance.products;
dir_products_specular = SimulationFolders.getInstance.products_specular;
dir_products_specular_reff_coeff = SimulationFolders.getInstance.products_specular_reff_coeff;
dir_products_specular_reflectivity = SimulationFolders.getInstance.products_specular_reflectivity;
dir_products_specular_reff_coeff_diel_profiles = SimulationFolders.getInstance.products_specular_reff_coeff_diel_profiles;
dir_products_specular_reflectivity_diel_profiles = SimulationFolders.getInstance.products_specular_reflectivity_diel_profiles;


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
% Configuration Parameters
DoYs = ConfigParams.getInstance.DoYs;
% Ground Parameters
gnd_structure_id = GndParams.getInstance.gnd_structure_id;
% Ground Multi-layer Parameters
calc_diel_profile_fit_functions = [];
if gnd_structure_id == Constants.id_gnd_multi_layered
    calc_diel_profile_fit_functions = GndMLParams.getInstance.calc_diel_profile_fit_functions;
end
% Bistatic Dynamic Parameters
AllPoints_m = BistaticDynParams.getInstance.AllPoints_m;
AngS2R_rf = BistaticDynParams.getInstance.AngS2R_rf; % SP->Rx Rotation Angle
% Rotation Matrices Dynamic Parameters
u_ts = RotMatDynParams.getInstance.u_ts;
u_sr = RotMatDynParams.getInstance.u_sr;


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

% For multiple ground layers
if gnd_structure_id == Constants.id_gnd_multi_layered
    
    % Number of dielectric profiles depends on number of ground layers.
    [~, num_diel_profiles] = size( Constants.diel_profiles );
    
end

% P = 4 x 1
% b = 2 x 1

if gnd_structure_id == Constants.id_gnd_single_layered

    % Field
    % Vegetation
    b_coh1v = g_r * u_sr * r_sv{1,1} * u_ts * g_t * e_t1 ;
    b_coh2v = g_r * u_sr * r_sv{1,1} * u_ts * g_t * e_t2 ;
    b0_coh1v = g_r0 * u_sr * r_sv{1,1} * u_ts * g_t * e_t1 ;
    b0_coh2v = g_r0 * u_sr * r_sv{1,1} * u_ts * g_t * e_t2 ;

    % Bare-soil
    b_coh1b = g_r * u_sr * r_sb{1,1} * u_ts * g_t * e_t1 ;
    b_coh2b = g_r * u_sr * r_sb{1,1} * u_ts * g_t * e_t2 ;
    b0_coh1b = g_r0 * u_sr * r_sb{1,1} * u_ts * g_t * e_t1 ;
    b0_coh2b = g_r0 * u_sr * r_sb{1,1} * u_ts * g_t * e_t2 ;

    % Power
    % Vegetation
    P_coh1v = G_r * U_sr * R_sv{1,1} * U_ts * G_t * E_t1 ;
    P_coh2v = G_r * U_sr * R_sv{1,1} * U_ts * G_t * E_t2 ;
    P0_coh1v = G_r0 * U_sr * R_sv{1,1} * U_ts * G_t * E_t1 ;
    P0_coh2v = G_r0 * U_sr * R_sv{1,1} * U_ts * G_t * E_t2 ;

    % Bare-soil
    P_coh1b = G_r * U_sr * R_sb{1,1} * U_ts * G_t * E_t1 ;
    P_coh2b = G_r * U_sr * R_sb{1,1} * U_ts * G_t * E_t2 ;
    P0_coh1b = G_r0 * U_sr * R_sb{1,1} * U_ts * G_t * E_t1 ;
    P0_coh2b = G_r0 * U_sr * R_sb{1,1} * U_ts * G_t * E_t2 ;
    
elseif gnd_structure_id == Constants.id_gnd_multi_layered 

    for ii = 1 : num_diel_profiles
    
         if calc_diel_profile_fit_functions(ii, 1)
            
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

    end

end


%% SAVE OUTPUTS
% Day-of-year for Time-series simulations
if sim_mode_id == Constants.id_time_series
    
    DoY = DoYs( sim_counter );    
    filename = 'DoYs';
    writeVarIncremental( dir_products, filename, sim_counter, DoY )
    
end

% Kc Factor
filename = 'Kc';
writeComplexVar(dir_products_specular, filename, Kc)


% Determine the output fodlers w.r.t number of ground layers
if num_diel_profiles == 1
    
    pathname_reff_coeff{1,1} = dir_products_specular_reff_coeff;
    pathname_reflectivity{1,1} = dir_products_specular_reflectivity;
    
else
    
    pathname_reff_coeff = dir_products_specular_reff_coeff_diel_profiles;
    pathname_reflectivity = dir_products_specular_reflectivity_diel_profiles;
    
end


filename_veg_1 = 'Veg1';
filename_veg_2 = 'Veg2';
filename_veg_01 = 'Veg01';
filename_veg_02 = 'Veg02';
filename_bare_1 = 'Bare1';
filename_bare_2 = 'Bare2';
filename_bare_01 = 'Bare01';
filename_bare_02 = 'Bare02';

for ii = 1 : num_diel_profiles
    
    if gnd_structure_id == Constants.id_gnd_single_layered ...
            || ( gnd_structure_id == Constants.id_gnd_multi_layered && calc_diel_profile_fit_functions(ii, 1) )
    
        % Field: 2 X 1
        if gnd_cover_id == Constants.id_veg_cover

            writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_veg_1, sim_counter, b_coh1v(:,ii) )
            writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_veg_2, sim_counter, b_coh2v(:,ii) )
            writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_veg_01, sim_counter, b0_coh1v(:,ii) )
            writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_veg_02, sim_counter, b0_coh2v(:,ii) )

        end

        writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_bare_1, sim_counter, b_coh1b(:,ii) )
        writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_bare_2, sim_counter, b_coh2b(:,ii) )
        writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_bare_01, sim_counter, b0_coh1b(:,ii) )
        writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_bare_02, sim_counter, b0_coh2b(:,ii) )

        % Power: 4 X 1
        if gnd_cover_id == Constants.id_veg_cover

            writeVarIncremental(pathname_reflectivity{1,ii}, filename_veg_1, sim_counter, P_coh1v(:,ii) );
            writeVarIncremental(pathname_reflectivity{1,ii}, filename_veg_2, sim_counter, P_coh2v(:,ii));
            writeVarIncremental(pathname_reflectivity{1,ii}, filename_veg_01, sim_counter, P0_coh1v(:,ii) );
            writeVarIncremental(pathname_reflectivity{1,ii}, filename_veg_02, sim_counter, P0_coh2v(:,ii) );

        end

        writeVarIncremental(pathname_reflectivity{1,ii}, filename_bare_1, sim_counter, P_coh1b(:,ii) );
        writeVarIncremental(pathname_reflectivity{1,ii}, filename_bare_2, sim_counter, P_coh2b(:,ii) );
        writeVarIncremental(pathname_reflectivity{1,ii}, filename_bare_01, sim_counter, P0_coh1b(:,ii) );
        writeVarIncremental(pathname_reflectivity{1,ii}, filename_bare_02, sim_counter, P0_coh2b(:,ii) );
    
    end
    
end


end



%% Calculate Specular Reflection Matrix(SRM)
function  [R_sv, R_sb, r_sv, r_sb] = CalcSRM()

%% GET GLOBAL DIRECTORIES
dir_afsa = SimulationFolders.getInstance.afsa;


%% GET GLOBAL PARAMETERS
% Simulation Settings
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
% Ground Parameters
gnd_structure_id = GndParams.getInstance.gnd_structure_id;
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


%% CALCULATIONS
% Specular Reflection Matrix
thsd = AngT2S_sf(1, 1) ;

dKz_s = squeeze(dKz(:, ANGDEG == round(thsd), :)) ;

%% Transmissivity Matrix
ArgH = 0 ;
ArgV = 0 ;

% If ground cover is Vegetation-cover
if gnd_cover_id == Constants.id_veg_cover
    
    
    %% GET GLOBAL PARAMETERS
    % Vegetation Parameters
    dim_layers_m = VegParams.getInstance.dim_layers_m;
    num_veg_layers = VegParams.getInstance.num_layers;                  
    
    % Calculate the attenuation
    for ii = 1 : num_veg_layers

        ArgH = ArgH + dKz_s(1, ii) * dim_layers_m(ii, 1) ;
        ArgV = ArgV + dKz_s(2, ii) * dim_layers_m(ii, 1) ;

    end
    
end

% vegetation transmissivity, if any
t_sv = [exp(+1i * ArgV) 0; 0 exp(+1i * ArgH)] ;
% unity transmissivity for bare soil
t_sb = [1 0; 0 1] ;


%% GROUND REFLECTION MATRIX
% If single-layered ground
if gnd_structure_id == Constants.id_gnd_single_layered
    
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
elseif gnd_structure_id == Constants.id_gnd_multi_layered
    
    %% GET GLOBAL PARAMETERS
    calc_diel_profile_fit_functions = GndMLParams.getInstance.calc_diel_profile_fit_functions;

    
    [r_g_cell] = reflectionCoeffsML();
    
    num_diel_profiles = length(r_g_cell);
    
    r_sv = cell(num_diel_profiles, 1);
    r_sb = cell(num_diel_profiles, 1);
    R_sv = cell(num_diel_profiles, 1);
    R_sb = cell(num_diel_profiles, 1);
    
    for ii = 1 : num_diel_profiles
        
        if calc_diel_profile_fit_functions(ii, 1)
            
            r_sv{ii,1} = t_sv * r_g_cell{ii,1} * t_sv ;
            r_sb{ii,1} = t_sb * r_g_cell{ii,1} * t_sb ;
            % 4 X 4
            R_sv{ii,1} = calc_Muller(r_sv{ii,1});
            R_sb{ii,1} = calc_Muller(r_sb{ii,1});
            
        end
        
    end

end


end