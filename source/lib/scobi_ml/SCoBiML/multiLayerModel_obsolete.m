
function multiLayerModel_obsolete


%% GET GLOBAL DIRECTORIES
dir_out_ml_ref = SimulationFolders.getInstance.out_ml_ref;


%% GET GLOBAL PARAMETERS
% Dynamic Parameters
DoYs = DynParams.getInstance.DoYs;
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
num_Th = length(th0_Tx_list_deg);
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;   
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
f_Hz = f_MHz * Constants.MHz2Hz;
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Ground Parameters
sand_ratio = GndParams.getInstance.sand_ratio;  % Texture at various depths 
clay_ratio = GndParams.getInstance.clay_ratio;  % Texture at various depths 
rhob_gcm3 = GndParams.getInstance.rhob_gcm3;    % Soil bulk density at various depths 
% Ground MultiLayer Parameters
layer_bottom_m = GndMLParams.getInstance.layer_bottom_m;
zA_m = GndMLParams.getInstance.zA_m;    % Air layer
z_m = GndMLParams.getInstance.z_m;    % Layer profile


% Wavelength
lambda_m = Constants.c / f_Hz ;


%% PLOT
[fig1, fig2] = plotSMdata();


%% MULTILAYER CALCULATIONS

% Initialize reflectivities for four different profiles
DoYlist = [];
Rp1_list = []; 
Rp2_list = [];
Rp1_L_list = [];
Rp2_L_list = [];
Rp1_2nd_list = [];
Rp2_2nd_list = [];
Rp1_3rd_list = [];
Rp2_3rd_list = [];

Ninterval = length(DoYs);

for ii = 1 : 10 : Ninterval
    
    % Get the corresponding VSM at the current index
    ParamsManager.index_VSM( ii );
    VSM_cm3cm3 = VSM_list_cm3cm3(ii,:)';
    
    % Get the corresponding DoY at the current index
    DoY = DoYs(ii);

%     [Rp1, Rp2, Rp1_L, Rp2_L, Rp1_2nd, Rp2_2nd, Rp1_3rd, Rp2_3rd] = multiLayerModel( fig1, fig2, VSM_cm3cm3, DoY ) ;

    % Soil Dielectric Constant
    eps_diel_soil = round(10 * dielg(VSM_cm3cm3, f_Hz, sand_ratio, clay_ratio, rhob_gcm3)) / 10 ;


    %% GENERATE DIELECTRIC PROFILES
    [eps_diel_z2, eps_diel_z3, eps_diel_zL, eps_diel_zS] = generateDielProfiles(eps_diel_soil, z_m); 


    %% PLOT DIELECTRIC PROFILES
    plotDielProfiles( fig1, eps_diel_soil, eps_diel_z2, eps_diel_z3, eps_diel_zL, eps_diel_zS, z_m );


    %% CALCULATE REFLECTION COEFFICIENTS
    zzb = (zA_m + layer_bottom_m) ;
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
    Rp1_2nd = zeros(num_Th, 1) ;
    Rp2_2nd = Rp1_2nd ;
    Rp1_3rd = zeros(num_Th, 1) ;
    Rp2_3rd = Rp1_3rd ;

    for tt = 1 : num_Th

        % Set theta index
        ParamsManager.index_Th( tt );

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
        [rh_L, rv_L] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z_m, eps_diel_zL) ;
        r_s_L = [rv_L 0; 0 rh_L] ;

        [r0_coh1b_L, r0_coh2b_L] = SpecularReflection(r_s_L) ;


        % Reflection Coefficient for 2nd Order Profile
        [rh_2nd, rv_2nd] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z_m, eps_diel_z2 ) ;
        r_s_2nd = [rv_2nd 0; 0 rh_2nd] ;

        [r0_coh1b_2nd, r0_coh2b_2nd] = SpecularReflection(r_s_2nd) ;


        % Reflection Coefficient for 3rd Order Profile
        [rh_3rd, rv_3rd] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z_m, eps_diel_z3 ) ;
        r_s_3rd = [rv_3rd 0; 0 rh_3rd] ;

        [r0_coh1b_3rd, r0_coh2b_3rd] = SpecularReflection(r_s_3rd) ;


        if pol_Tx == 'X' && pol_Rx == 'X'
            % Reflectivity
            Rp1(tt) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
            Rp2(tt) = abs(r0_coh2b(2)) .^ 2 ; % at ?/?0 = 1
            Rp1_L(tt) = abs(r0_coh1b_L(1)) .^ 2 ; % at ?/?0 = 1
            Rp2_L(tt) = abs(r0_coh2b_L(2)) .^ 2 ; % at ?/?0 = 1
            Rp1_2nd(tt) = abs(r0_coh1b_2nd(1)) .^ 2 ; % at ?/?0 = 1
            Rp2_2nd(tt) = abs(r0_coh2b_2nd(2)) .^ 2 ; % at ?/?0 = 1
            Rp1_3rd(tt) = abs(r0_coh1b_3rd(1)) .^ 2 ; % at ?/?0 = 1
            Rp2_3rd(tt) = abs(r0_coh2b_3rd(2)) .^ 2 ; % at ?/?0 = 1
        else
            % Reflectivity
            Rp1(tt) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
            Rp2(tt) = abs(r0_coh1b(2)) .^ 2 ; % at ?/?0 = 1
            Rp1_L(tt) = abs(r0_coh1b_L(1)) .^ 2 ; % at ?/?0 = 1
            Rp2_L(tt) = abs(r0_coh1b_L(2)) .^ 2 ; % at ?/?0 = 1
            Rp1_2nd(tt) = abs(r0_coh1b_2nd(1)) .^ 2 ; % at ?/?0 = 1
            Rp2_2nd(tt) = abs(r0_coh1b_2nd(2)) .^ 2 ; % at ?/?0 = 1
            Rp1_3rd(tt) = abs(r0_coh1b_3rd(1)) .^ 2 ; % at ?/?0 = 1
            Rp2_3rd(tt) = abs(r0_coh1b_3rd(2)) .^ 2 ; % at ?/?0 = 1
        end

    end

    plotReflectivityForProfiles( fig1, fig2, DoY, Rp1, Rp2, Rp1_L, Rp2_L, Rp1_2nd, Rp2_2nd, Rp1_3rd, Rp2_3rd );

    plotReflectivityVsTh( Rp1, Rp2, Rp1_L, Rp2_L, Rp1_2nd, Rp2_2nd, Rp1_3rd, Rp2_3rd );

    DoYlist = [DoYlist; DoY];
    Rp1_list = [Rp1_list; Rp1]; 
    Rp2_list = [Rp2_list; Rp2];
    Rp1_L_list = [Rp1_L_list; Rp1_L];
    Rp2_L_list = [Rp2_L_list; Rp2_L];
    Rp1_2nd_list = [Rp1_2nd_list; Rp1_2nd];
    Rp2_2nd_list = [Rp2_2nd_list; Rp2_2nd];
    Rp1_3rd_list = [Rp1_3rd_list; Rp1_3rd];
    Rp2_3rd_list = [Rp2_3rd_list; Rp2_3rd];

    [M1(ii), M2(ii)] = plotAddSMpoint( fig1, fig2, DoY, VSM_cm3cm3 );

end

% plotMovie( M1 );

close(fig1);
close(fig2);

%% SAVE
Rp1s = [Rp1_list, Rp1_L_list, Rp1_2nd_list, Rp1_3rd_list];
Rp2s = [Rp2_list, Rp2_L_list, Rp2_2nd_list, Rp2_3rd_list];

filename01 = strcat('ML_Ref01') ;
filename02 = strcat('ML_Ref02') ;
filenameDoY = strcat('DoY') ;

writeVar(dir_out_ml_ref, filename01, Rp1s);
writeVar(dir_out_ml_ref, filename02, Rp2s);
writeVar(dir_out_ml_ref, filenameDoY, DoYlist);

end


function [rh, rv] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z_m, eps_diel_z)


Lz = diff(z_m)' / lambda_m ; % complex optical length in units of lambda_m

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
