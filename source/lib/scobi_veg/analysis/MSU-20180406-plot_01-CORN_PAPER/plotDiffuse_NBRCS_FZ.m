
% Mehmet Kurum
% 11/12/2017
% NBRCS vs. Fresnel Zones


function plotDiffuse_NBRCS_FZ


%% ANALYSIS INPUTS
inputFile_sys_tag = 'sysInput-CORN_PAPER-row_';
inputFile_veg_tag = 'vegVirRowInput-Corn-row';

cornStagesStruct = struct('s1', 'V1_V9', ...
                            's2', 'V10_VT', ...
                            's3', 'R1_R4', ...
                            's4', 'R5', ...
                            's5', 'R6') ;

% Corn field row azimuth angles
rowPhis = [0 45 90] ;
num_rowPhis =  length(rowPhis) ;


% TO-DO: Check here
%% Choices
stage_choice = cornStagesStruct.s3;
PH_choice = 1;
TH_choice = 4 ;
SM_choice = 3 ;
RMSH_choice = 2 ;


%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(1) ), '.xml' );
inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(1) ), '-', stage_choice, '.xml' );

getInput( inputFile_sys, inputFile_veg );


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nfz = SimParams.getInstance.Nfz;


%% NBRCS
NBRCS_dB = zeros(2, 2, Nfz, num_rowPhis) ;
NBRCS1_dB = zeros(2, 2, Nfz, num_rowPhis) ;
NBRCS2_dB = zeros(2, 2, Nfz, num_rowPhis) ;
NBRCS3_dB = zeros(2, 2, Nfz, num_rowPhis) ;
NBRCS4_dB = zeros(2, 2, Nfz, num_rowPhis) ;

for jj = 1 : num_rowPhis
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(jj) ), '.xml' );
    inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(jj) ), '-', stage_choice, '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    
    % TO-DO: Assign the correct index to each
    ParamsManager.index_Th( TH_choice );
    ParamsManager.index_Ph( PH_choice );
    ParamsManager.index_VSM( SM_choice );
    ParamsManager.index_RMSH( RMSH_choice );
        
    % Initialize the directories depending on theta phi, VSM, and RMSH
    SimulationFolders.getInstance.initializeDynamicDirs();


    %% GET GLOBAL DIRECTORIES
    dir_out_diffuse_NBRCS_tuple = SimulationFolders.getInstance.out_diffuse_NBRCS_tuple;

    % NBRCS   
    % transmit port 1 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS_t1_dB') ;
    NBRCS_dB(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS1_t1_dB') ;
    NBRCS1_dB(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS2_t1_dB') ;
    NBRCS2_dB(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS3_t1_dB') ;
    NBRCS3_dB(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS4_t1_dB') ;
    NBRCS4_dB(1, :, :, jj) = xx ;
    % transmit port 2 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS_t2_dB') ;
    NBRCS_dB(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS1_t2_dB') ;
    NBRCS1_dB(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS2_t2_dB') ;
    NBRCS2_dB(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS3_t2_dB') ;
    NBRCS3_dB(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS4_t2_dB') ;
    NBRCS4_dB(2, :, :, jj) = xx ;
    
   
end


%% plotting as a function of fresnel zones
for jj = 1 : num_rowPhis
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(jj) ), '.xml' );
    inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(jj) ), '-', stage_choice, '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    
    %% GET GLOBAL PARAMETERS
    % Receiver Parameters
    hr_m = RxParams.getInstance.hr_m;
    
    
%     NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, jj)) ;
%     NBRCS_dB_RL = squeeze(NBRCS_dB(1, 2, :, jj)) ;
    
    NBRCS1_dB_RR = squeeze(NBRCS1_dB(1, 1, :, jj)) ;
    NBRCS1_dB_RL = squeeze(NBRCS1_dB(1, 2, :, jj)) ;
    
    NBRCS2_dB_RR = squeeze(NBRCS2_dB(1, 1, :, jj)) ;
    NBRCS2_dB_RL = squeeze(NBRCS2_dB(1, 2, :, jj)) ;
    
%     NBRCS3_dB_RR = squeeze(NBRCS3_dB(1, 1, :, jj)) ;
%     NBRCS3_dB_RL = squeeze(NBRCS3_dB(1, 2, :, jj)) ;
    
    NBRCS4_dB_RR = squeeze(NBRCS4_dB(1, 1, :, jj)) ;
    NBRCS4_dB_RL = squeeze(NBRCS4_dB(1, 2, :, jj)) ;
    
    if jj == 1
        figure
        subplot(1,2,1)
        title('X-Pol [RL]')
        hold
        axis([0 Nfz+2 -50 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS [dB]')
        xticks(1 : Nfz)
        grid
        text( Nfz+0.5, NBRCS2_dB_RL(end) - 3, strcat('\phi_r_o_w [\circ]'))
        
        
        subplot(1,2,2)
        title('Co-Pol [RR]')
        hold
        axis([0 Nfz+2 -50 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS [dB]')
        xticks(1 : Nfz)
        grid
    end
    subplot(1,2,1)
    plot(NBRCS1_dB_RL, ':sb', 'MarkerFaceColor', 'blue')
    plot(NBRCS2_dB_RL, ':dg', 'MarkerFaceColor', 'green')
    plot(NBRCS4_dB_RL, ':or', 'MarkerFaceColor', 'red')
    text(Nfz+0.05, NBRCS2_dB_RL(end), strcat(num2str( rowPhis(jj) )))
    
    
    subplot(1,2,2)
    plot(NBRCS1_dB_RR, ':sb', 'MarkerFaceColor', 'blue')
    plot(NBRCS2_dB_RR, ':dg', 'MarkerFaceColor', 'green')
    plot(NBRCS4_dB_RR, ':or', 'MarkerFaceColor', 'red')
    legend('Single Bounce', 'Double Bounce', 'Triple Bounce')
    
    if jj == 1 
        
        %% GET GLOBAL DIRECTORIES
        dir_analysis = SimulationFolders.getInstance.analysis;
        if ~exist(dir_analysis, 'dir')
            mkdir(dir_analysis)
        end
        
    end

end


%% GET GLOBAL PARAMETERS
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
ph0_Tx_deg = ph0_Tx_list_deg( ParamsManager.index_Ph );
th0_Tx_deg = DynParams.getInstance.th0_Tx_list_deg( ParamsManager.index_Th );
VSM_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3( ParamsManager.index_VSM );
RMSH_cm = DynParams.getInstance.RMSH_list_cm( ParamsManager.index_RMSH );


%%
fname1 = strcat('TH_', num2str( th0_Tx_deg ), '-', pol_Tx, pol_Rx, '-PH_', num2str( ph0_Tx_deg ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close

%%

for jj = 1 : num_rowPhis
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    
    %% GET GLOBAL PARAMETERS
    hr_m = RxParams.getInstance.hr_m;
    
    
    NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, jj)) ;
    NBRCS_dB_RL = squeeze(NBRCS_dB(1, 2, :, jj)) ;
    
   
    if jj == 1
        figure
        hold
        axis([0 Nfz+2 -20 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS [dB]')
        xticks(1 : Nfz)
        grid
    
    end
    plot(NBRCS_dB_RL, ':sb', 'MarkerFaceColor', 'blue')
    plot(NBRCS_dB_RR, ':dg', 'MarkerFaceColor', 'green')
    text( Nfz+0.5, NBRCS_dB_RL(end), strcat( num2str( rowPhis(jj) ) ) )
    
    legend('X-POL [RL]', 'Co-POL [RR]', 'location', 'southwest')
    
    
end

text(Nfz+0.5, - 1, strcat('\phi_r_o_w [\circ]'))


%% GET GLOBAL PARAMETERS
% Dynamic Parameters
th0_Tx_deg = DynParams.getInstance.th0_Tx_list_deg( ParamsManager.index_Th );
VSM_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3( ParamsManager.index_VSM );


%%
fname1 = strcat('TH2_', num2str( th0_Tx_deg ), '-', pol_Tx, pol_Rx, '-PH_', num2str( ph0_Tx_deg ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end