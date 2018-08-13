
% Mehmet Kurum
% 11/15/2017
% Diffuse + Specular

function plotPower_VWC_stage


%% ANALYSIS INPUTS
inputFile_sys_tag = 'sysInput-CORN_PAPER-row_';
inputFile_veg_tag = 'vegVirRowInput-Corn-row';

cornStages = { 'V1_V9', 'V10_VT', 'R1_R4', 'R5', 'R6' };
num_stages = length( cornStages );

% Corn field row azimuth angles
rowPhis = [0 45 90] ;


%% Choices
% TO-DO: Check here
rowPhi_choice = 1;
FZ_choice = 1 ;
TH_choice = 4 ;
PH_choice = 1 ;
SM_choice = 3 ;
RMSH_choice = 2 ;
EIRP_choice_dB = 27;
G0r_choice_dB = 0;


%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(rowPhi_choice) ), '.xml' );
inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(rowPhi_choice) ), '-', cornStages{1}, '.xml' );

getInput( inputFile_sys, inputFile_veg );

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();
% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nfz = SimParams.getInstance.Nfz;
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;

    
%% NBRCS
PP1_inc_dB = zeros(2, 2, num_stages) ; % % PP1_inc1_dB = zeros(2, 2, num_stages) ; 
% % PP1_inc2_dB = zeros(2, 2, num_stages) ; PP1_inc3_dB = zeros(2, 2, num_stages) ; PP1_inc4_dB = zeros(2, 2, num_stages) ;
PP1_inc_dB2 = zeros(2, 2, Nfz, num_stages) ;

KKi_dB = zeros(num_stages) ;
P1_areas_dB = zeros(num_stages) ;
P1_areas_dB2 = zeros(Nfz, num_stages) ;      
    
ParamsManager.index_VSM( SM_choice );

for ss = 1 : num_stages

    %% GET INPUT
    inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(rowPhi_choice) ), '.xml' );
    inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(rowPhi_choice) ), '-', cornStages{ ss }, '.xml' );

    getInput( inputFile_sys, inputFile_veg );

    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();

    % Set theta index
    ParamsManager.index_Th( TH_choice );
    ParamsManager.index_Ph( PH_choice );

    % Assign the index of interest to each, in this analysis
    ParamsManager.index_RMSH( RMSH_choice );

    % Initialize the directories depending on theta phi, VSM, and RMSH
    SimulationFolders.getInstance.initializeDynamicDirs();


    %% GET GLOBAL DIRECTORIES
    dir_config = SimulationFolders.getInstance.config;
    dir_freqdiff = SimulationFolders.getInstance.freqdiff;
    dir_out_diffuse_P1_tuple = SimulationFolders.getInstance.out_diffuse_P1_tuple;
    
    %% GET GLOBAL PARAMETERS
    % Virtual Row-Structured Vegetation Parameters
    plugin = VegVirRowParams.getInstance.plugin;

    vwc(ss) = plugin.getVWC;
    [diel_re(ss), diel_im(ss)] = plugin.getDiel;
    diel_abs(ss) = abs( diel_re(ss) + 1i * diel_im(ss) );

    % read Ki
    filename1 = strcat('Ki') ;
    Ki = readComplexVar( dir_freqdiff, filename1) ;
    KKi_dB(ss) = 10 * log10(abs(Ki) ^ 2 / 4 / pi) ;

    % Fresnel ellipses
    filenamex = 'ellipses_FZ_m' ;
    ellipses_FZ_m = readVar( dir_config, filenamex) ;
    areas_FZ_m = pi * ellipses_FZ_m(:, 1) .* ellipses_FZ_m(:, 2) ;
    P1_areas_dB(ss) = 10 * log10(areas_FZ_m(FZ_choice)) ;
    P1_areas_dB2(:, ss) = 10 * log10(areas_FZ_m) ;

    % transmit port 1 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t1_dB') ;
    PP1_inc_dB(1, :, ss) = squeeze(xx(:, FZ_choice)) ;
    PP1_inc_dB2(1, :, :, ss) = xx;

    % transmit port 2 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t2_dB') ;
    PP1_inc_dB(2, :, ss) = squeeze(xx(:, FZ_choice)) ;
    PP1_inc_dB2(2, :, :, ss) = xx;

end



%% RECEIVED POWER - SPECULAR
% Intilization
% over vegetation
P_cohveg_dB = zeros(2, 2, num_stages) ;
% over bare soil
P_cohbare_dB = zeros(2, 2, num_stages) ;

KKc_dB = zeros(num_stages) ;
        
ParamsManager.index_VSM( SM_choice );

for ss = 1 : num_stages

    %% GET INPUT
    inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(rowPhi_choice) ), '.xml' );
    inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(rowPhi_choice) ), '-', cornStages{ ss }, '.xml' );

    getInput( inputFile_sys, inputFile_veg );

    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();

    % Set theta index
    ParamsManager.index_Th( TH_choice );
    ParamsManager.index_Ph( PH_choice );

    % TO-DO: Assign the correct index to each
    ParamsManager.index_RMSH( RMSH_choice );

    % Initialize the directories depending on theta phi, VSM, and RMSH
    SimulationFolders.getInstance.initializeDynamicDirs();


    %% GET GLOBAL DIRECTORIES
    dir_out_specular = SimulationFolders.getInstance.out_specular;

    % read Kc
    filename1 = strcat('Kc') ;
    Kc = readComplexVar( dir_out_specular, filename1) ;
    KKc_dB(ss) = 20 * log10(abs(Kc)) ;

    [~, ~, ~, ~, ~, ~, ~, ~, ...
        P_coh1vegx, P_coh2vegx, ~, ~, ...
        P_coh1barex, P_coh2barex, ~, ~] = readSpecular ;

    P_cohveg_dB(1, :, ss) = 10 * log10( P_coh1vegx(1 : 2) ) + KKc_dB(ss) ;
    P_cohveg_dB(2, :, ss) = 10 * log10( P_coh2vegx(1 : 2) ) + KKc_dB(ss) ;

    P_cohbare_dB(1, :, ss) = 10 * log10( P_coh1barex(1 : 2) ) + KKc_dB(ss) ;
    P_cohbare_dB(2, :, ss) = 10 * log10( P_coh2barex(1 : 2) ) + KKc_dB(ss) ;
    
    if ss == 1
        
        dir_analysis = SimulationFolders.getInstance.analysis;
        if ~exist(dir_analysis, 'dir')
            mkdir(dir_analysis)
        end
        
    end

end


%% Total Power for different Transmitter PhI angles
        
ParamsManager.index_VSM( SM_choice );      

Pvi_RR = squeeze(PP1_inc_dB(1, 1, :)) + EIRP_choice_dB + G0r_choice_dB ;
Pvi_RL = squeeze(PP1_inc_dB(1, 2, :)) + EIRP_choice_dB + G0r_choice_dB ;    
Pvc_RR = squeeze(P_cohveg_dB(1, 1, :)) + EIRP_choice_dB + G0r_choice_dB ;
Pvc_RL = squeeze(P_cohveg_dB(1, 2, :)) + EIRP_choice_dB + G0r_choice_dB ;    
Pbc_RR = squeeze(P_cohbare_dB(1, 1, :)) + EIRP_choice_dB + G0r_choice_dB ;
Pbc_RL = squeeze(P_cohbare_dB(1, 2, :)) + EIRP_choice_dB + G0r_choice_dB ;    


%% LEFT PLOT - X-POL
subplot(1,2,1)
% Specular - Bare
plot(1:5, Pbc_RL, ':or', 'MarkerFaceColor', 'red', 'markersize', 4)
hold
% Specular - Vegetation
plot(1:5, Pvc_RL, '-dr', 'MarkerFaceColor', 'red', 'markersize', 5)
% Diffuse
plot(1:5, Pvi_RL, '-^g', 'MarkerFaceColor', 'green', 'markersize', 5)
% Setup the visualization
axis([0.5, 5.5, -237+EIRP_choice_dB+G0r_choice_dB, -187+EIRP_choice_dB+G0r_choice_dB])
xlabel('Growth Stages')
ylabel('Received Power [dB]')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
xtickangle(30)
grid
text(1.5, Pvi_RL(end) + 13, strcat( '\theta_s=', num2str( th0_Tx_list_deg(TH_choice) ), '\circ' ))
text(1.5, Pvi_RL(end) + 10, strcat( 'VSM=', num2str( VSM_list_cm3cm3(SM_choice) ), ' cm^3/cm^3' ))
text(1.5, Pvi_RL(end) + 7, strcat( 'RMSH=', num2str( RMSH_list_cm(RMSH_choice) ), ' cm' ))
title('X-Pol [RL]')

% Display VWC
yyaxis right
% plot( vwc, ':bo', 'markersize', 5 )
% axis([0.5 5.5 0 10])
% set(gca, 'ytick', 0:3);
plot( vwc, ':bo', 'markersize', 5 )
axis([0.5 5.5 0 10])
set(gca, 'ytick', 0 : 0.5 : 3);
%ylabel('Dielectric Constant','Color','b')
yyaxis left


%% RIGHT PLOT - CO-POL
subplot(1,2,2)
% Specular - Bare
plot(1:5, Pbc_RR, ':or', 'MarkerFaceColor', 'red', 'markersize', 4)
hold
% Specular - Vegetation
plot(1:5, Pvc_RR, '-dr', 'MarkerFaceColor', 'red', 'markersize', 5)
% Diffuse
plot(1:5, Pvi_RR, '-^g', 'MarkerFaceColor', 'green', 'markersize', 5)
% Setup the visualization
axis([0.5, 5.5, -237+EIRP_choice_dB+G0r_choice_dB, -187+EIRP_choice_dB+G0r_choice_dB])
%ylabel('Received Power [dB]')
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
xtickangle(30)
grid
title('Co-Pol [RR]')

% Display VWC
yyaxis right
plot( vwc, ':bo', 'markersize', 5 )
axis([0.5 5.5 0 10])
set(gca, 'ytick', 0 : 0.5 : 3);
ylabel('VWC [kg/m^2]','Color','b')
yyaxis left


legend('Specular - Bare Soil', 'Specular - Vegetation', 'Diffuse', 'VWC')


% TO-DO: Check this
%% GET GLOBAL PARAMETERS
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
VSM_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3( ParamsManager.index_VSM );
RMSH_cm = DynParams.getInstance.RMSH_list_cm( ParamsManager.index_RMSH );


%%
fname1 = strcat('ReceivedPower_VWC_vs_stages-FZ_', num2str(FZ_choice), '-', pol_Tx, pol_Rx, '-TH_', num2str( th0_Tx_list_deg(TH_choice) ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end