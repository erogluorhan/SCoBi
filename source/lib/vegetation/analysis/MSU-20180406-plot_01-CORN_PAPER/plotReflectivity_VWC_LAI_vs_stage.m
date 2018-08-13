
% Mehmet Kurum
% 11/15/2017
% Diffuse + Specular

function plotReflectivity_VWC_LAI_vs_stage


%% ANALYSIS INPUTS
inputFile_sys = 'sysInput-CORN_PAPER-row0.xml';
inputFile_sys_CONST_DIEL = 'sysInput-CORN_PAPER-CONST_DIEL.xml';
inputFile_veg_tag = 'vegVirRowInput-CORN_PAPER-row0-';

cornStages = { 'V1_V9', 'V10_VT', 'R1_R4', 'R5', 'R6' };
num_stages = length( cornStages );


%% Choices
% TO-DO: Check here
FZ_choice = 1 ;
TH_choice = 4 ;
PH_choice = 1 ;
VSM_choice = 2 ;
RMSH_choice = 3 ;
TH_choice_CONST_DIEL = 2 ;
PH_choice_CONST_DIEL = 1 ;
SM_choice_CONST_DIEL = 1 ;
RMSH_choice_CONST_DIEL = 2;


%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_veg = strcat( inputFile_veg_tag, cornStages{1}, '.xml' );

getInput( inputFile_sys, inputFile_veg );

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();
% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();


%% GET GLOBAL PARAMETERS
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm; 


%% RECEIVED POWER - SPECULAR
% Intilization
% over vegetation
P_cohveg_dB_CONST_DIEL = zeros(2, 2, num_stages) ;

KKc_dB_CONST_DIEL = zeros(num_stages,1) ;
        
ParamsManager.index_VSM( SM_choice_CONST_DIEL );

for ss = 1 : num_stages

    %% GET INPUT
    inputFile_veg_CONST_DIEL = strcat( inputFile_veg_tag, cornStages{ss}, '-CONST_DIEL.xml' );

    getInput( inputFile_sys_CONST_DIEL, inputFile_veg_CONST_DIEL );

    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();

    % Set theta index
    ParamsManager.index_Th( TH_choice_CONST_DIEL );
    ParamsManager.index_Ph( PH_choice_CONST_DIEL );

    % TO-DO: Assign the correct index to each
    ParamsManager.index_RMSH( RMSH_choice_CONST_DIEL );

    % Initialize the directories depending on theta phi, VSM, and RMSH
    SimulationFolders.getInstance.initializeDynamicDirs();


    %% GET GLOBAL DIRECTORIES
    dir_out_specular = SimulationFolders.getInstance.out_specular;
    
    % read Kc
    filename1 = strcat('Kc') ;
    Kc_CONST_DIEL = readComplexVar( dir_out_specular, filename1) ;
    KKc_dB_CONST_DIEL(ss) = 20 * log10(abs(Kc_CONST_DIEL)) ;

    [~, ~, ~, ~, ~, ~, ~, ~, ...
        P_coh1vegx_CONST_DIEL, P_coh2vegx_CONST_DIEL, ~, ~, ~, ~, ~, ~] = readSpecular ;

    P_cohveg_dB_CONST_DIEL(1, :, ss) = 10 * log10( P_coh1vegx_CONST_DIEL(1 : 2) ) + KKc_dB_CONST_DIEL(ss) ;
    P_cohveg_dB_CONST_DIEL(2, :, ss) = 10 * log10( P_coh2vegx_CONST_DIEL(1 : 2) ) + KKc_dB_CONST_DIEL(ss) ;
 
end


%% Total (Specular == Total due to dominance) Reflectivity for different Transmitter PhI angles  
REF_vc_RR_CONST_DIEL = squeeze(P_cohveg_dB_CONST_DIEL(1, 1, :)) - KKc_dB_CONST_DIEL;
REF_vc_RL_CONST_DIEL = squeeze(P_cohveg_dB_CONST_DIEL(1, 2, :)) - KKc_dB_CONST_DIEL;    


%% RECEIVED POWER - SPECULAR
% Intilization
% over vegetation
P_cohveg_dB = zeros(2, 2, num_stages) ;
% over bare soil
P_cohbare_dB = zeros(2, 2, num_stages) ;

KKc_dB = zeros(num_stages,1) ;
        
ParamsManager.index_VSM( VSM_choice );

for ss = 1 : num_stages

    %% GET INPUT
    inputFile_veg = strcat( inputFile_veg_tag, cornStages{ ss }, '.xml' );

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
    
    
    %% GET GLOBAL PARAMETERS
    % Virtual Row-Structured Vegetation Parameters
    plugin = VegVirRowParams.getInstance.plugin;

    VWC(ss) = plugin.getVWC;
    LAI(ss) = plugin.getLAI;

    
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


%% Total (Specular == Total due to dominance) Reflectivity for different Transmitter PhI angles
        
ParamsManager.index_VSM( VSM_choice );      
  
REF_vc_RR = squeeze(P_cohveg_dB(1, 1, :)) - KKc_dB;
REF_vc_RL = squeeze(P_cohveg_dB(1, 2, :)) - KKc_dB;    
REF_bc_RR = squeeze(P_cohbare_dB(1, 1, :)) - KKc_dB;
REF_bc_RL = squeeze(P_cohbare_dB(1, 2, :)) - KKc_dB;    


%% LEFT PLOT - X-POL
subplot(1,2,1)
% Reflectivity - Bare
plot(1:5, REF_bc_RL, ':or', 'MarkerFaceColor', 'red', 'markersize', 4)
hold
% Reflectivity - Vegetation
plot(1:5, REF_vc_RL, '-dr', 'MarkerFaceColor', 'red', 'markersize', 5)
plot(1:5, REF_vc_RL_CONST_DIEL, '--dr', 'markersize', 6)

% Setup the visualization
% axis([0.5, 5.5, -222-KKc_dB(1,1), -187-KKc_dB(1,1)])
axis([0.5, 5.5, -40, 0])
xlabel('Growth Stages')
ylabel('Reflectivity [dB]')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
xtickangle(30)
title('X-Pol [RL]')

% Display VWC
yyaxis right
% plot( VWC, ':bo', 'markersize', 5 )
% axis([0.5 5.5 0 10])
% set(gca, 'ytick', 0:3);
plot( VWC, '--bo', 'markersize', 5 )
axis([0.5 5.5 0 8])
set(gca, 'ytick', 0 : 2, 'YColor', 'blue');
ylabel('VWC [kg/m^2]','Color','b')
yyaxis left

set(gca, 'YColor', 'red');
grid
legend('Bare Soil', 'Vegetation', 'Veg. (Const. Diel.)', 'VWC')


%% RIGHT PLOT - CO-POL
subplot(1,2,2)
% Specular - Bare
plot(1:5, REF_bc_RR, ':or', 'MarkerFaceColor', 'red', 'markersize', 4)
hold
% Specular - Vegetation
plot(1:5, REF_vc_RR, '-dr', 'MarkerFaceColor', 'red', 'markersize', 5)
plot(1:5, REF_vc_RR_CONST_DIEL, '--dr', 'markersize', 6)

% Setup the visualization
% axis([0.5, 5.5, -222-KKc_dB(1,1), -187-KKc_dB(1,1)])
axis([0.5, 5.5, -40, 0])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
xtickangle(30)
text(1.5, REF_bc_RL(1) - 1, strcat( '\theta_s=', num2str( th0_Tx_list_deg(TH_choice) ), '\circ' ))
text(1.5, REF_bc_RL(1) - 4, strcat( 'VSM=', num2str( VSM_list_cm3cm3(VSM_choice) ), ' cm^3/cm^3' ))
text(1.5, REF_bc_RL(1) - 7, strcat( 'RMSH=', num2str( RMSH_list_cm(RMSH_choice) ), ' cm' ))
title('Co-Pol [RR]')

% Display LAI
yyaxis right
plot( LAI, '--g^', 'markersize', 5 )
axis([0.5 5.5 0 2])
set(gca, 'ytick', 0 : 0.25 : 0.5, 'YColor', 'green');
ylabel('LAI [m^2/m^2]','Color','g')
yyaxis left


set(gca, 'YColor', 'red');
grid
legend('Bare Soil', 'Vegetation', 'Veg. (Const. Diel.)', 'LAI')


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
fname1 = strcat('Reflectivity_VWC_LAI_vs_stages(with_Const_Diel)-FZ_', num2str(FZ_choice), '-', pol_Tx, pol_Rx, '-TH_', num2str( th0_Tx_list_deg(TH_choice) ), '-PH_', num2str( ph0_Tx_list_deg(PH_choice) ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end