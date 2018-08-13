
function plotNBRCS_Periodic_Random_vs_TH


%% ANALYSIS INPUTS
inputFile_sys_tag = 'sysInput-CORN_PAPER-';
inputFile_veg_tag = 'vegVirRowInput-CORN_PAPER-';

cornStagesStruct = struct('s1', 'V1_V9', ...
                            's2', 'V10_VT', ...
                            's3', 'R1_R4', ...
                            's4', 'R5', ...
                            's5', 'R6') ;


%% Choices
% TO-DO: Check here
stage_choice = cornStagesStruct.s3;
FZ_choice = 1 ;
PH_choice = 1;
VSM_choice = 1 ;
RMSH_choice = 1;


    
%% NBRCS - PERIODIC

%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_sys = strcat( inputFile_sys_tag, 'PERIODIC.xml' );
inputFile_veg = strcat( inputFile_veg_tag, 'row0-', stage_choice, '-PERIODIC.xml' );

getInput( inputFile_sys, inputFile_veg );

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();
% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nfz = SimParams.getInstance.Nfz;
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg ;
num_Th = length( th0_Tx_list_deg );

% Initialization
PP1_inc_dB_PERIODIC = zeros(2, 2, num_Th) ; % % PP1_inc1_dB = zeros(2, 2, num_Th) ; 
PP1_inc_dB2_PERIODIC = zeros(2, 2, Nfz, num_Th) ;

KKi_dB_PERIODIC = zeros(num_Th, 1) ;
P1_areas_dB_PERIODIC = zeros(num_Th, 1) ;
P1_areas_dB2_PERIODIC = zeros(Nfz, num_Th) ;


%% READ        
% Assign the index of interest to each, in this analysis
ParamsManager.index_Ph( PH_choice );
ParamsManager.index_VSM( VSM_choice );
ParamsManager.index_RMSH( RMSH_choice );
    
for tt = 1 : num_Th

    % Set theta index
    ParamsManager.index_Th( tt );

    % Initialize the directories depending on theta phi, VSM, and RMSH
    SimulationFolders.getInstance.initializeDynamicDirs();


    %% GET GLOBAL DIRECTORIES
    dir_config = SimulationFolders.getInstance.config;
    dir_freqdiff = SimulationFolders.getInstance.freqdiff;
    dir_out_diffuse_P1_tuple = SimulationFolders.getInstance.out_diffuse_P1_tuple;


    %% NBRCS
    % read Ki
    filename1 = strcat('Ki') ;
    Ki = readComplexVar( dir_freqdiff, filename1) ;
    KKi_dB_PERIODIC(tt) = 10 * log10(abs(Ki) ^ 2 / 4 / pi) ;

    % Fresnel ellipses
    filenamex = 'ellipses_FZ_m' ;
    ellipses_FZ_m = readVar( dir_config, filenamex) ;
    areas_FZ_m = pi * ellipses_FZ_m(:, 1) .* ellipses_FZ_m(:, 2) ;
    P1_areas_dB_PERIODIC(tt) = 10 * log10(areas_FZ_m(FZ_choice)) ;
    P1_areas_dB2_PERIODIC(:, tt) = 10 * log10(areas_FZ_m) ;

    % transmit port 1 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t1_dB') ;
    PP1_inc_dB_PERIODIC(1, :, tt) = squeeze(xx(:, FZ_choice)) ;
    PP1_inc_dB2_PERIODIC(1, :, :, tt) = xx;

    % transmit port 2 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t2_dB') ;
    PP1_inc_dB_PERIODIC(2, :, tt) = squeeze(xx(:, FZ_choice)) ;
    PP1_inc_dB2_PERIODIC(2, :, :, tt) = xx;

end

%% Normalization
NBRCS_dB_PERIODIC = PP1_inc_dB_PERIODIC ...
    - reshape(repmat(KKi_dB_PERIODIC, 4, 1), 2, 2, num_Th) ...
    - reshape(repmat(P1_areas_dB_PERIODIC, 4, 1), 2, 2, num_Th) ; 


    
%% NBRCS - RANDOM

%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_sys = strcat( inputFile_sys_tag, 'RANDOM.xml' );
inputFile_veg = strcat( inputFile_veg_tag, 'RANDOM-', stage_choice, '.xml' );

getInput( inputFile_sys, inputFile_veg );

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();
% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nfz = SimParams.getInstance.Nfz;
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg ;
num_Th = length( th0_Tx_list_deg );

% Initialization
PP1_inc_dB_RANDOM = zeros(2, 2, num_Th) ; % % PP1_inc1_dB = zeros(2, 2, num_Th) ; 
PP1_inc_dB2_RANDOM = zeros(2, 2, Nfz, num_Th) ;

KKi_dB_RANDOM = zeros(num_Th, 1) ;
P1_areas_dB_RANDOM = zeros(num_Th, 1) ;
P1_areas_dB2_RANDOM = zeros(Nfz, num_Th) ;


%% READ        
% Assign the index of interest to each, in this analysis
ParamsManager.index_Ph( PH_choice );
ParamsManager.index_VSM( VSM_choice );
ParamsManager.index_RMSH( RMSH_choice );
    
for tt = 1 : num_Th

    % Set theta index
    ParamsManager.index_Th( tt );

    % Initialize the directories depending on theta phi, VSM, and RMSH
    SimulationFolders.getInstance.initializeDynamicDirs();


    %% GET GLOBAL DIRECTORIES
    dir_config = SimulationFolders.getInstance.config;
    dir_freqdiff = SimulationFolders.getInstance.freqdiff;
    dir_out_diffuse_P1_tuple = SimulationFolders.getInstance.out_diffuse_P1_tuple;


    %% NBRCS
    % read Ki
    filename1 = strcat('Ki') ;
    Ki = readComplexVar( dir_freqdiff, filename1) ;
    KKi_dB_RANDOM(tt) = 10 * log10(abs(Ki) ^ 2 / 4 / pi) ;

    % Fresnel ellipses
    filenamex = 'ellipses_FZ_m' ;
    ellipses_FZ_m = readVar( dir_config, filenamex) ;
    areas_FZ_m = pi * ellipses_FZ_m(:, 1) .* ellipses_FZ_m(:, 2) ;
    P1_areas_dB_RANDOM(tt) = 10 * log10(areas_FZ_m(FZ_choice)) ;
    P1_areas_dB2_RANDOM(:, tt) = 10 * log10(areas_FZ_m) ;

    % transmit port 1 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t1_dB') ;
    PP1_inc_dB_RANDOM(1, :, tt) = squeeze(xx(:, FZ_choice)) ;
    PP1_inc_dB2_RANDOM(1, :, :, tt) = xx;

    % transmit port 2 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t2_dB') ;
    PP1_inc_dB_RANDOM(2, :, tt) = squeeze(xx(:, FZ_choice)) ;
    PP1_inc_dB2_RANDOM(2, :, :, tt) = xx;

end

%% Normalization
NBRCS_dB_RANDOM = PP1_inc_dB_RANDOM ...
    - reshape(repmat(KKi_dB_RANDOM, 4, 1), 2, 2, num_Th) ...
    - reshape(repmat(P1_areas_dB_RANDOM, 4, 1), 2, 2, num_Th) ; 




% TO-DO: Check this

%% GET GLOBAL DIRECTORIES
dir_analysis = SimulationFolders.getInstance.analysis;
if ~exist(dir_analysis, 'dir')
    mkdir(dir_analysis)
end
        
        
%% GET GLOBAL PARAMETERS
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
ph0_Tx_deg = ph0_Tx_list_deg( ParamsManager.index_Ph );
VSM_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3( ParamsManager.index_VSM );
RMSH_cm = DynParams.getInstance.RMSH_list_cm( ParamsManager.index_RMSH );
% Vegetation Parameters
vegetation_stage = VegParams.getInstance.vegetation_stage;


%% FIGURES
figure 
hold on
axis([0 80 -25 -10])
xlabel('\theta_s [\circ]')
ylabel('\sigma^0_e [dB]')
xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
grid

%% PLOT PERIODIC
NBRCS_dB_RR_PERIODIC = squeeze(NBRCS_dB_PERIODIC(1, 1, :)) ;
NBRCS_dB_RL_PERIODIC = squeeze(NBRCS_dB_PERIODIC(1, 2, :)) ;

plot(th0_Tx_list_deg, NBRCS_dB_RL_PERIODIC, '--ob', 'MarkerFaceColor', 'blue', 'MarkerSize', 5)
plot(th0_Tx_list_deg, NBRCS_dB_RR_PERIODIC, '--or', 'MarkerFaceColor', 'red', 'MarkerSize', 5)


%% PLOT RANDOM
NBRCS_dB_RR_RANDOM = squeeze(NBRCS_dB_RANDOM(1, 1, :)) ;
NBRCS_dB_RL_RANDOM = squeeze(NBRCS_dB_RANDOM(1, 2, :)) ;

plot(th0_Tx_list_deg, NBRCS_dB_RL_RANDOM, ':^b', 'MarkerFaceColor', 'blue', 'MarkerSize', 5)
plot(th0_Tx_list_deg, NBRCS_dB_RR_RANDOM, ':^r', 'MarkerFaceColor', 'red', 'MarkerSize', 5)


text(23, -19, strcat( '\phi_s=', num2str( ph0_Tx_deg ), '\circ' ) )
text(23, -20.5, strcat( 'VSM=', num2str( VSM_cm3cm3 ), ' cm^3/cm^3' ))
text(23, -22, strcat( 'RMSH=', num2str( RMSH_cm ), ' cm' ))
text(23, -23.5, strcat('fzone=', num2str(FZ_choice) ) )
text(3, NBRCS_dB_RL_PERIODIC(1,1), strcat('RL'))
text(3, NBRCS_dB_RR_PERIODIC(1,1), strcat('RR'))
title('NBRCS Variations due to Corn Row Structure')

legend('Periodic', 'Periodic', 'Random', 'Random')

fname = strcat('NBRCS(Periodic_Random)_vs_TH-', vegetation_stage, '-', pol_Tx, pol_Rx, '-FZ_', num2str(FZ_choice), '-PH_', num2str( ph0_Tx_deg ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname = strrep( fname, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname), 'tiff')
close


end