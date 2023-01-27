
function specularTerm
% function specularTerm 
%
%   Calculates the specular reflection coefficient (coherent forward 
%   scattering through specular reflection point) and reflectivity. Stores 
%   the calculated values into simulation output folders in an incremental 
%   fashion as the simulation iterations continue.
%
%   See also mainScoBi, directTerm.

%   Copyright � 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



%% GET GLOBAL DIRECTORIES
dir_products = SimulationFolders.getInstance.products;
dir_products_specular = SimulationFolders.getInstance.products_specular;
dir_products_specular_reff_coeff = SimulationFolders.getInstance.products_specular_reff_coeff;
dir_products_specular_reflectivity = SimulationFolders.getInstance.products_specular_reflectivity;
dir_products_specular_reff_coeff_diel_profiles = SimulationFolders.getInstance.products_specular_reff_coeff_diel_profiles;
dir_products_specular_reflectivity_diel_profiles = SimulationFolders.getInstance.products_specular_reflectivity_diel_profiles;
dir_products_specular_pen_dep = SimulationFolders.getInstance.products_specular_pen_dep;
dir_products_specular_pen_dep_diel_profiles = SimulationFolders.getInstance.products_specular_pen_dep_diel_profiles;

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
if gnd_structure_id == Constants.ID_GND_MULTI_LAYERED
    calc_diel_profile_fit_functions = GndMLParams.getInstance.calc_diel_profile_fit_functions;
    calculate_penetration_depth = GndMLParams.getInstance.calculate_penetration_depth;
else
    calculate_penetration_depth = false ;

end
% Bistatic Dynamic Parameters
AllPoints_m = BistaticDynParams.getInstance.AllPoints_m;
AngS2R_rf = BistaticDynParams.getInstance.AngS2R_rf; % SP->Rx Rotation Angle
% Rotation Matrices Dynamic Parameters
u_ts = RotMatDynParams.getInstance.u_ts;
u_sr = RotMatDynParams.getInstance.u_sr;


%% INITIALIZE REQUIRED PARAMETERS
% AllPoints_m = [pos_Tx_m, pos_TxI_m, pos_SP_m, pos_Rx_m, pos_RxI_m, pos_Gnd_m, pos_B_Rx_m, pos_FZ_m] ;
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

f_Hz = f_MHz * Constants.MHZ_TO_HZ ;
lambda_m = Constants.LIGHTSPEED / f_Hz ;     % Wavelength
k0 = 2 * pi * f_Hz / Constants.LIGHTSPEED ;    % Wave number

% Factor Kc
K = 1i * sqrt(EIRP) * sqrt(G0r) * lambda_m / (4 * pi) ;

Kc = K * exp(1i * k0 * (r_st + r_sr)) / (r_st + r_sr) ; 

thrd = AngS2R_rf(1, 1) ;
phrd = AngS2R_rf(2, 1) ;

ant_pat_res_factor = 1 / ant_pat_res_deg;

thd = round( ant_pat_res_factor * rad2deg(th)) / ant_pat_res_factor ; % rounding operation is due to accuracy concerns
phd = round( ant_pat_res_factor * rad2deg(ph)) / ant_pat_res_factor ; % to make the angles multiples of ant_pat_res_deg

% Receiver Antenna values in the specular direction
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

% Transmitter-Receiver Rotation Matrix
% 4 X 4
U_ts = calcMuller(u_ts) ;
% 4 X 4
U_sr = calcMuller(u_sr) ;

% Calculate specular reflection matrices and penetration depth if
% applicable
[R_sv, R_sb, r_sv, r_sb, pd_cell] = calcSRM() ; % r_sv for vegetation, r_sb for bare soil


%% SPECULAR TERM
% For single ground layer
num_diel_profiles = 1;

% For multiple ground layers
if gnd_structure_id == Constants.ID_GND_MULTI_LAYERED
    
    % Number of dielectric profiles depends on number of ground layers.
    [~, num_diel_profiles] = size( Constants.DIEL_PROFILES );
    
end

% P = 4 x 1
% b = 2 x 1

if gnd_structure_id == Constants.ID_GND_SINGLE_LAYERED

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
    
elseif gnd_structure_id == Constants.ID_GND_MULTI_LAYERED 

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
if sim_mode_id == Constants.ID_TIME_SERIES
    
    DoY = DoYs( sim_counter );    
    filename = 'DoYs';
    writeVarIncremental( dir_products, filename, sim_counter, DoY )
    
end

% Kc Factor
filename = 'Kc';
writeComplexVar(dir_products_specular, filename, Kc)


% Determine the output folders w.r.t number of ground layers
if num_diel_profiles == 1
    
    pathname_reff_coeff{1,1} = dir_products_specular_reff_coeff;
    pathname_reflectivity{1,1} = dir_products_specular_reflectivity;
    pathname_pen_dep{1,1} = dir_products_specular_pen_dep ;
    
else
    
    pathname_reff_coeff = dir_products_specular_reff_coeff_diel_profiles;
    pathname_reflectivity = dir_products_specular_reflectivity_diel_profiles;
    pathname_pen_dep = dir_products_specular_pen_dep_diel_profiles ;
    
end


filename_veg_1 = 'Veg1';
filename_veg_2 = 'Veg2';
filename_veg_01 = 'Veg01';
filename_veg_02 = 'Veg02';
filename_bare_1 = 'Bare1';
filename_bare_2 = 'Bare2';
filename_bare_01 = 'Bare01';
filename_bare_02 = 'Bare02';
filename_pendep_1 = 'PenDep1';
filename_pendep_2 = 'PenDep2';

for ii = 1 : num_diel_profiles
    
    if gnd_structure_id == Constants.ID_GND_SINGLE_LAYERED ...
            || ( gnd_structure_id == Constants.ID_GND_MULTI_LAYERED && calc_diel_profile_fit_functions(ii, 1) )
    
        % Reflection coefficient: 2 X 1
        if gnd_cover_id == Constants.ID_VEG_COVER

            writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_veg_1, sim_counter, b_coh1v(:,ii) )
            writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_veg_2, sim_counter, b_coh2v(:,ii) )
            writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_veg_01, sim_counter, b0_coh1v(:,ii) )
            writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_veg_02, sim_counter, b0_coh2v(:,ii) )

        end

        writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_bare_1, sim_counter, b_coh1b(:,ii) )
        writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_bare_2, sim_counter, b_coh2b(:,ii) )
        writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_bare_01, sim_counter, b0_coh1b(:,ii) )
        writeComplexVarIncremental(pathname_reff_coeff{1,ii}, filename_bare_02, sim_counter, b0_coh2b(:,ii) )

        % Reflectivity: 4 X 1
        if gnd_cover_id == Constants.ID_VEG_COVER

            writeVarIncremental(pathname_reflectivity{1,ii}, filename_veg_1, sim_counter, P_coh1v(:,ii) );
            writeVarIncremental(pathname_reflectivity{1,ii}, filename_veg_2, sim_counter, P_coh2v(:,ii));
            writeVarIncremental(pathname_reflectivity{1,ii}, filename_veg_01, sim_counter, P0_coh1v(:,ii) );
            writeVarIncremental(pathname_reflectivity{1,ii}, filename_veg_02, sim_counter, P0_coh2v(:,ii) );

        end

        writeVarIncremental(pathname_reflectivity{1,ii}, filename_bare_1, sim_counter, P_coh1b(:,ii) );
        writeVarIncremental(pathname_reflectivity{1,ii}, filename_bare_2, sim_counter, P_coh2b(:,ii) );
        writeVarIncremental(pathname_reflectivity{1,ii}, filename_bare_01, sim_counter, P0_coh1b(:,ii) );
        writeVarIncremental(pathname_reflectivity{1,ii}, filename_bare_02, sim_counter, P0_coh2b(:,ii) );
        
        % Penetration depth
        if calculate_penetration_depth
             pdval = pd_cell{ii,:};
             writeVarIncremental(pathname_pen_dep{1,ii}, filename_pendep_1, sim_counter, pdval(1) );
             writeVarIncremental(pathname_pen_dep{1,ii}, filename_pendep_2, sim_counter, pdval(2) );
             clear pdval;
        end
    
    end
    
end


end



function  [R_sv, R_sb, r_sv, r_sb, pd_cell] = calcSRM()
% Calculate Specular Reflection Matrix(SRM)


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
if gnd_cover_id == Constants.ID_VEG_COVER
    
    
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
if gnd_structure_id == Constants.ID_GND_SINGLE_LAYERED
    
    ths = deg2rad(thsd) ;
    [RGHIF, RGVIF, ~, ~] = reflectionCoeffSingle(ths, ths, eps_g, h) ;

    r_g = [RGVIF 0; 0 RGHIF] ;


    %% SPECULAR REFLECTION MATRIX
    % 2 X 2
    r_sv{1,1} = t_sv * r_g * t_sv ;
    r_sb{1,1} = t_sb * r_g * t_sb ;
    % 4 X 4
    R_sv{1,1} = calcMuller(r_sv{1,1}) ;
    R_sb{1,1} = calcMuller(r_sb{1,1}) ;
    
    % empty cell; no penetration depth in single layer mode
    pd_cell = {};

% Else if multi-layered ground
elseif gnd_structure_id == Constants.ID_GND_MULTI_LAYERED
    
    %% GET GLOBAL PARAMETERS
    calc_diel_profile_fit_functions = GndMLParams.getInstance.calc_diel_profile_fit_functions;

    
    [r_g_cell, pd_cell] = reflectionCoeffsML();
    
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
            R_sv{ii,1} = calcMuller(r_sv{ii,1});
            R_sb{ii,1} = calcMuller(r_sb{ii,1});
            
        end
        
    end

end


end