
%% Mehmet Kurum
% September 04, 2017


%%

% -----   inf  -------
% -------  air --------
% ----- layer 1 -------
% ----- layer 2 -------
% ..................
% ----- layer N -------
% ---     inf  -------


function [Rp1, Rp2, Rp1_L, Rp2_L, Rp1_2, Rp2_2, Rp1_3, Rp2_3] = multiLayerModel( fig1, fig2, VSM_gcm3, DoY )


%% GET GLOBAL PARAMETERS
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
f_Hz = f_MHz * Constants.MHz2Hz;
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
num_Th = length(th0_Tx_list_deg);
% Ground Parameters
gnd_layer_depth_m = GndParams.getInstance.layer_depth_m; % Layer Measurement Depth [in m]
sand_ratio = GndParams.getInstance.sand_ratio;  % Texture at various depths 
clay_ratio = GndParams.getInstance.clay_ratio;  % Texture at various depths 
rhob_gcm3 = GndParams.getInstance.rhob_gcm3;    % Soil bulk density at various depths 
% Ground MultiLayer Parameters
delZ_m = GndMLParams.getInstance.delZ_m;    % Layer discretization
zA_m = GndMLParams.getInstance.zA_m;    % Air layer
zB_m = GndMLParams.getInstance.zB_m;    % Bottom-most layer


%% CALCULATIONS
% Wavelength
lambda_m = Constants.c / f_Hz ;

% Soil Dielectric Constant
eps_diel_soil = round(10 * dielg(VSM_gcm3, f_Hz, sand_ratio, clay_ratio, rhob_gcm3)) / 10 ;


%% LAYER PROFILE
% Layer bottom
L_bottom = [0; gnd_layer_depth_m(1 : end - 1) + diff(gnd_layer_depth_m) / 2 ] ;

zS = zA_m + L_bottom(end) + zB_m ; % total
z = (0 : delZ_m : zS)' ;


%% GENERATE DIELECTRIC PROFILES
[eps_diel_z2, eps_diel_z3, eps_diel_zL, eps_diel_zS] = generateDielProfiles(eps_diel_soil, z); 


%% PLOT DIELECTRIC PROFILES
plotDielProfiles( fig1, eps_diel_soil, eps_diel_z2, eps_diel_z3, eps_diel_zL, eps_diel_zS, z );


%% CALCULATE REFLECTION COEFFICIENTS
zzb = (zA_m + L_bottom) ;
Lzb = diff(zzb)' ; % / lambda_m ; % complex optical length in units of lambda_m
nA = sqrte(Constants.eps_diel_air) ;
nS = sqrte(eps_diel_soil) ;
nS = nS.';

na = [nA; nA; nA] ;
ns = [nS; nS; nS] ;

% input to multidiel
n = [na, ns] ;

% Initialize
Rp1 = zeros( num_Th, 1) ;
Rp2 = Rp1 ;
Rp1_L = zeros(num_Th, 1) ;
Rp2_L = Rp1_L ;
Rp1_2 = zeros(num_Th, 1) ;
Rp2_2 = Rp1_2 ;
Rp1_3 = zeros(num_Th, 1) ;
Rp2_3 = Rp1_3 ;

for ii = 1 : num_Th
        
    % Set theta index
    ParamsManager.index_Th( ii );
    
    % Initialize the directories depending on dynamic parameters
    SimulationFolders.getInstance.initializeDynamicDirs();
    
    %% GET GLOBAL PARAMETERS
    % Transmitter Parameters
    th0_Tx_deg = th0_Tx_list_deg( ParamsManager.index_Th );
    
    mainSCoBi;
    
    % Reflection Coefficient for 
    rh = multidiel(n, Lzb, 1, th0_Tx_deg, 'te') ;
    rv = multidiel(n, Lzb, 1, th0_Tx_deg, 'th') ;
    r_s = [rv 0; 0 rh] ;
    
    [r0_coh1b, r0_coh2b] = SpecularReflection(r_s) ;

    
    % Reflection Coefficient for Logistic Profile
    [rh_L, rv_L] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z, eps_diel_zL) ;
    r_s_L = [rv_L 0; 0 rh_L] ;
    
    [r0_coh1b_L, r0_coh2b_L] = SpecularReflection(r_s_L) ;
    
    
    % Reflection Coefficient for 2nd Order Profile
    [rh_2, rv_2] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z, eps_diel_z2 ) ;
    r_s_2 = [rv_2 0; 0 rh_2] ;
    
    [r0_coh1b_2, r0_coh2b_2] = SpecularReflection(r_s_2) ;
    
    
    % Reflection Coefficient for 3rd Order Profile
    [rh_3, rv_3] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z, eps_diel_z3 ) ;
    r_s_3 = [rv_3 0; 0 rh_3] ;
    
    [r0_coh1b_3, r0_coh2b_3] = SpecularReflection(r_s_3) ;
 
    
    if pol_Tx == 'X' && pol_Rx == 'X'
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh2b(2)) .^ 2 ; % at ?/?0 = 1
        Rp1_L(ii) = abs(r0_coh1b_L(1)) .^ 2 ; % at ?/?0 = 1
        Rp2_L(ii) = abs(r0_coh2b_L(2)) .^ 2 ; % at ?/?0 = 1
        Rp1_2(ii) = abs(r0_coh1b_2(1)) .^ 2 ; % at ?/?0 = 1
        Rp2_2(ii) = abs(r0_coh2b_2(2)) .^ 2 ; % at ?/?0 = 1
        Rp1_3(ii) = abs(r0_coh1b_3(1)) .^ 2 ; % at ?/?0 = 1
        Rp2_3(ii) = abs(r0_coh2b_3(2)) .^ 2 ; % at ?/?0 = 1
    else
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh1b(2)) .^ 2 ; % at ?/?0 = 1
        Rp1_L(ii) = abs(r0_coh1b_L(1)) .^ 2 ; % at ?/?0 = 1
        Rp2_L(ii) = abs(r0_coh1b_L(2)) .^ 2 ; % at ?/?0 = 1
        Rp1_2(ii) = abs(r0_coh1b_2(1)) .^ 2 ; % at ?/?0 = 1
        Rp2_2(ii) = abs(r0_coh1b_2(2)) .^ 2 ; % at ?/?0 = 1
        Rp1_3(ii) = abs(r0_coh1b_3(1)) .^ 2 ; % at ?/?0 = 1
        Rp2_3(ii) = abs(r0_coh1b_3(2)) .^ 2 ; % at ?/?0 = 1
    end
    
end

plotReflectivityForProfiles( fig1, fig2, DoY, Rp1, Rp2, Rp1_L, Rp2_L, Rp1_2, Rp2_2, Rp1_3, Rp2_3 );

plotReflectivityVsTh( Rp1, Rp2, Rp1_L, Rp2_L, Rp1_2, Rp2_2, Rp1_3, Rp2_3 );

end

function [rh, rv] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z, eps_diel_z)


Lz = diff(z)' / lambda_m ; % complex optical length in units of lambda_m

nAz = sqrte(Constants.eps_diel_air) ;
nmz = sqrte(eps_diel_z(2 : end, :)) ;
nSz = sqrte(eps_diel_z(end, :)) ;

% Air - % isotropic
na = [nAz; nAz; nAz] ;

% Dielectric Profile : isotropic
nm = [nmz(:, 1).'; nmz(:, 1).'; nmz(:, 1).'] ;

% Soil - isotropic
nb = [nSz(:, 1); nSz(:, 1); nSz(:, 1)] ;

%% input to multidiel
n = [na, nm, nb] ;

%% Reflection Coeffficeint
rh = multidiel(n, Lz, 1, th0_Tx_deg, 'te') ;
rv = multidiel(n, Lz, 1, th0_Tx_deg, 'th') ;

end

