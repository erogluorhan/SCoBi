
% Mehmet Kurum
% 11/12/2017
% NBRCS vs. Observation Angle


function plotDiffuse_NBRCS_ANG


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
FZ_choice = 1 ;
PH_choice = 1;
SM_choice = 3 ;
RMSH_choice = 2 ;
rowPhi_choice = 1;


%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(rowPhi_choice) ), '.xml' );
inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(rowPhi_choice) ), '-', stage_choice, '.xml' );

getInput( inputFile_sys, inputFile_veg );


%% GET GLOBAL PARAMETERS
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg ;
num_Th = length( th0_Tx_list_deg );


%% NBRCS

NBRCS_dB = zeros(2, 2, num_Th, num_rowPhis) ; NBRCS1_dB = zeros(2, 2, num_Th, num_rowPhis) ;
NBRCS2_dB = zeros(2, 2, num_Th, num_rowPhis) ; NBRCS3_dB = zeros(2, 2, num_Th, num_rowPhis) ; NBRCS4_dB = zeros(2, 2, num_Th, num_rowPhis) ;


for jj = 1 : num_rowPhis
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
       
    for ii = 1 : num_Th
        
        % Set theta index
        ParamsManager.index_Th( ii );
        
        % TO-DO: Assign the correct index to each
        ParamsManager.index_Ph( PH_choice );
        ParamsManager.index_VSM( SM_choice );
        ParamsManager.index_RMSH( RMSH_choice );
        
        % Initialize the directories depending on theta phi, VSM, and RMSH
        SimulationFolders.getInstance.initializeDynamicDirs();
        
        
        %% GET GLOBAL DIRECTORIES
        dir_out_diffuse_NBRCS_tuple = SimulationFolders.getInstance.out_diffuse_NBRCS_tuple;
                
        % transmit port 1 / receiver ports 1&2
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS_t1_dB') ;
        NBRCS_dB(1, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS1_t1_dB') ;
        NBRCS1_dB(1, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS2_t1_dB') ;
        NBRCS2_dB(1, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS3_t1_dB') ;
        NBRCS3_dB(1, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS4_t1_dB') ;
        NBRCS4_dB(1, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        
        % transmit port 2 / receiver ports 1&2
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS_t2_dB') ;
        NBRCS_dB(2, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS1_t2_dB') ;
        NBRCS1_dB(2, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS2_t2_dB') ;
        NBRCS2_dB(2, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS3_t2_dB') ;
        NBRCS3_dB(2, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS4_t2_dB') ;
        NBRCS4_dB(2, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        
    end
end


%% reading reflectivity
% Intilization
% over vegetation
R_cohveg_dB = zeros(2, 2, num_Th, num_rowPhis) ;
% over bare soil
R_cohbare_dB = zeros(2, 2, num_Th, num_rowPhis) ;

for jj = 1 : num_rowPhis
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    for ii = 1 : num_Th
        
        % Set theta index
        ParamsManager.index_Th( ii );
        
        % TO-DO: Assign the correct index to each
        ParamsManager.index_Ph( PH_choice );
        ParamsManager.index_VSM( SM_choice );
        ParamsManager.index_RMSH( RMSH_choice );
        
        % Initialize the directories depending on theta phi, VSM, and RMSH
        SimulationFolders.getInstance.initializeDynamicDirs();
        
        
        [~, ~, ~, ~, ~, ~, ~, ~, ...
            R_coh1vegx, R_coh2vegx, ~, ~, ...
            R_coh1barex, R_coh2barex, ~, ~] = readSpecular ;
        
        R_cohveg_dB(1, :, ii, jj) = 10 * log10(R_coh1vegx(1 : 2)) ;
        R_cohveg_dB(2, :, ii, jj) = 10 * log10(R_coh2vegx(1 : 2)) ;
        
        R_cohbare_dB(1, :, ii, jj) = 10 * log10(R_coh1barex(1 : 2)) ;
        R_cohbare_dB(2, :, ii, jj) = 10 * log10(R_coh2barex(1 : 2)) ;
        
    end
    
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
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
VSM_cm3cm3 = VSM_list_cm3cm3( ParamsManager.index_VSM );
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;
RMSH_cm = RMSH_list_cm( ParamsManager.index_RMSH );


%% plotting as a function of angle

NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, rowPhi_choice)) ;
NBRCS_dB_RL = squeeze(NBRCS_dB(1, 2, :, rowPhi_choice)) ;

NBRCS1_dB_RR = squeeze(NBRCS1_dB(1, 1, :, rowPhi_choice)) ;
NBRCS1_dB_RL = squeeze(NBRCS1_dB(1, 2, :, rowPhi_choice)) ;
% NBRCS2_dB_RR = squeeze(NBRCS2_dB(1, 1, :, rowPhi_choice)) ;
% NBRCS2_dB_RL = squeeze(NBRCS2_dB(1, 2, :, rowPhi_choice)) ;
% NBRCS3_dB_RR = squeeze(NBRCS3_dB(1, 1, :, rowPhi_choice)) ;
% NBRCS3_dB_RL = squeeze(NBRCS3_dB(1, 2, :, rowPhi_choice)) ;
NBRCS4_dB_RR = squeeze(NBRCS4_dB(1, 1, :, rowPhi_choice)) ;
NBRCS4_dB_RL = squeeze(NBRCS4_dB(1, 2, :, rowPhi_choice)) ;

Rvc_RR = squeeze(R_cohveg_dB(1, 1, :, rowPhi_choice)) ;
Rvc_RL = squeeze(R_cohveg_dB(1, 2, :, rowPhi_choice)) ;

Rbc_RR = squeeze(R_cohbare_dB(1, 1, :, rowPhi_choice)) ;
Rbc_RL = squeeze(R_cohbare_dB(1, 2, :, rowPhi_choice)) ;


figure
subplot(1,2,1)
plot(th0_Tx_list_deg, NBRCS_dB_RL, '-or', 'MarkerFaceColor', 'red')
hold
plot(th0_Tx_list_deg, NBRCS1_dB_RL, '-ok', 'MarkerFaceColor', 'black')
plot(th0_Tx_list_deg, NBRCS4_dB_RL, '-oc', 'MarkerFaceColor', 'cyan')
plot(th0_Tx_list_deg, NBRCS_dB_RR, '-or', 'MarkerFaceColor', 'white')
plot(th0_Tx_list_deg, NBRCS1_dB_RR, '-ok', 'MarkerFaceColor', 'white')
plot(th0_Tx_list_deg, NBRCS4_dB_RR, '-oc', 'MarkerFaceColor', 'white')

% text(10, -10, strcat())
% text(10, -21, strcat(''))
% text(10, -30, strcat('Triple Bounce'))

legend('Double Bounce', 'Single Bounce', 'Triple Bounce', 'location', 'southeast')


text(-3, NBRCS_dB_RL(1), strcat('RL'))
text(-3, NBRCS_dB_RR(1), strcat('RR'))

text(-3, NBRCS1_dB_RL(1), strcat('RL'))
text(73, NBRCS1_dB_RR(end), strcat('RR'))

text(-3, NBRCS4_dB_RL(1), strcat('RL'))
text(73, NBRCS4_dB_RR(end), strcat('RR'))

axis([-5 85 -50 0])
xlabel('\theta_s [\circ]')
ylabel('NBRCS: \sigma^0_e [dB]')
xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
grid
title('Diffuse [Vegetation]')

%         figure(2)
subplot(1,2,2)
plot(th0_Tx_list_deg, Rbc_RL, '-sb', 'MarkerFaceColor', 'blue')
hold
plot(th0_Tx_list_deg, Rvc_RL, '-dg', 'MarkerFaceColor', 'green')
plot(th0_Tx_list_deg, Rvc_RR, '-dg', 'MarkerFaceColor', 'white')
plot(th0_Tx_list_deg, Rbc_RR, '-sb', 'MarkerFaceColor', 'white')
ylabel('Reflectivity: \Gamma_s [dB]')

text(-3, (Rvc_RL(1) + Rbc_RL(1))/2, strcat('RL'))
text(-3, (Rvc_RR(1) + Rbc_RR(1))/2, strcat('RR'))

axis([-5 85 -50 0])
xlabel('\theta_s [\circ]')
xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
grid
title('Specular')
legend('soil', 'vegetation', 'location', 'southeast')

%%
fname1 = strcat('FZ_', num2str(FZ_choice), '-', pol_Tx, pol_Rx, '-row_', num2str( rowPhis(rowPhi_choice) ), '-PH_', num2str( ph0_Tx_deg ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end