
% Mehmet Kurum
% 11/15/2017
% Diffuse + Specular

function plotNBRCS_TH


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
num_RowPhis =  length(rowPhis) ;


%% Choices
% TO-DO: Check here
stage_choice = cornStagesStruct.s3;
rowPhi_choice = 1;
FZ_choice = 1 ;
PH_choice = 1;
SM_choice = 3 ;
RMSH_choice = 2 ;
EIRP_choice_dB = 27;
G0r_choice_dB = 0;


%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(1) ), '.xml' );
inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(1) ), '-', stage_choice, '.xml' );

getInput( inputFile_sys, inputFile_veg );


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nfz = SimParams.getInstance.Nfz;
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg ;
num_Th = length( th0_Tx_list_deg );

    
%% NBRCS
PP1_inc_dB = zeros(2, 2, num_Th, num_RowPhis) ; % % PP1_inc1_dB = zeros(2, 2, num_Th) ; 
% % PP1_inc2_dB = zeros(2, 2, num_Th) ; PP1_inc3_dB = zeros(2, 2, num_Th) ; PP1_inc4_dB = zeros(2, 2, num_Th) ;
PP1_inc_dB2 = zeros(2, 2, Nfz, num_Th, num_RowPhis) ;

KKi_dB = zeros(num_Th, num_RowPhis) ;
P1_areas_dB = zeros(num_Th, num_RowPhis) ;
P1_areas_dB2 = zeros(Nfz, num_Th, num_RowPhis) ;

for jj = 1 : num_RowPhis
        
    inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(jj) ), '.xml' );
    inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(jj) ), '-', stage_choice, '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    for tt = 1 : num_Th
              
        % Set theta index
        ParamsManager.index_Th( tt );
        
        % Assign the index of interest to each, in this analysis
        ParamsManager.index_Ph( PH_choice );
        ParamsManager.index_VSM( SM_choice );
        ParamsManager.index_RMSH( RMSH_choice );
        
        % Initialize the directories depending on theta phi, VSM, and RMSH
        SimulationFolders.getInstance.initializeDynamicDirs();
        
        
        %% GET GLOBAL DIRECTORIES
        dir_config = SimulationFolders.getInstance.config;
        dir_freqdiff = SimulationFolders.getInstance.freqdiff;
        dir_out_diffuse_P1_tuple = SimulationFolders.getInstance.out_diffuse_P1_tuple;
        
        
        % read Ki
        filename1 = strcat('Ki') ;
        Ki = readComplexVar( dir_freqdiff, filename1) ;
        KKi_dB(tt, jj) = 10 * log10(abs(Ki) ^ 2 / 4 / pi) ;
        
        % Fresnel ellipses
        filenamex = 'ellipses_FZ_m' ;
        ellipses_FZ_m = readVar( dir_config, filenamex) ;
        areas_FZ_m = pi * ellipses_FZ_m(:, 1) .* ellipses_FZ_m(:, 2) ;
        P1_areas_dB(tt, jj) = 10 * log10(areas_FZ_m(FZ_choice)) ;
        P1_areas_dB2(:, tt, jj) = 10 * log10(areas_FZ_m) ;
        
        % transmit port 1 / receiver ports 1&2
        xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t1_dB') ;
        PP1_inc_dB(1, :, tt, jj) = squeeze(xx(:, FZ_choice)) ;
        PP1_inc_dB2(1, :, :, tt, jj) = xx;
        
        % transmit port 2 / receiver ports 1&2
        xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t2_dB') ;
        PP1_inc_dB(2, :, tt, jj) = squeeze(xx(:, FZ_choice)) ;
        PP1_inc_dB2(2, :, :, tt, jj) = xx;
        
    end

end


%% Reflectivity
% Intilization
% over vegetation
P_cohveg_dB = zeros(2, 2, num_Th, num_RowPhis) ;
% over bare soil
P_cohbare_dB = zeros(2, 2, num_Th, num_RowPhis) ;

KKc_dB = zeros(num_Th, num_RowPhis) ;

for jj = 1 : num_RowPhis
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(jj) ), '.xml' );
    inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(jj) ), '-', stage_choice, '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    for tt = 1 : num_Th

        % Set theta index
        ParamsManager.index_Th( tt );
        
        % TO-DO: Assign the correct index to each
        ParamsManager.index_Ph( PH_choice );
        ParamsManager.index_VSM( SM_choice );
        ParamsManager.index_RMSH( RMSH_choice );
        
        % Initialize the directories depending on theta phi, VSM, and RMSH
        SimulationFolders.getInstance.initializeDynamicDirs();
        
        
        %% GET GLOBAL DIRECTORIES
        dir_out_specular = SimulationFolders.getInstance.out_specular;
        
        % read Kc
        filename1 = strcat('Kc') ;
        Kc = readComplexVar( dir_out_specular, filename1) ;
        KKc_dB(tt, jj) = 20 * log10(abs(Kc)) ;
        
        [~, ~, ~, ~, ~, ~, ~, ~, ...
            P_coh1vegx, P_coh2vegx, ~, ~, ...
            P_coh1barex, P_coh2barex, ~, ~] = readSpecular ;
        
        P_cohveg_dB(1, :, tt, jj) = 10 * log10( P_coh1vegx(1 : 2) ) + KKc_dB(tt) ;
        P_cohveg_dB(2, :, tt, jj) = 10 * log10( P_coh2vegx(1 : 2) ) + KKc_dB(tt) ;
        
        P_cohbare_dB(1, :, tt, jj) = 10 * log10( P_coh1barex(1 : 2) ) + KKc_dB(tt) ;
        P_cohbare_dB(2, :, tt, jj) = 10 * log10( P_coh2barex(1 : 2) ) + KKc_dB(tt) ;
        
    end

    if jj == 1
        dir_analysis = SimulationFolders.getInstance.analysis;
        if ~exist(dir_analysis, 'dir')
            mkdir(dir_analysis)
        end
    end
    
end

%% Normalization
NBRCS_dB_TOT = 10 * log10(10.^(P_cohveg_dB/10) + 10.^(PP1_inc_dB/10)) ...
    - reshape(repmat(KKi_dB, 4, 1), 2, 2, num_Th, num_RowPhis) ...
    - reshape(repmat(P1_areas_dB, 4, 1), 2, 2, num_Th, num_RowPhis) ; 
NBRCS_dB = PP1_inc_dB ...
    - reshape(repmat(KKi_dB, 4, 1), 2, 2, num_Th, num_RowPhis) ...
    - reshape(repmat(P1_areas_dB, 4, 1), 2, 2, num_Th, num_RowPhis) ; 

for ii = 1 : Nfz
    NBRCS_dB_TOT2(:, :, ii, :, :) = 10 * log10(10.^(P_cohveg_dB/10) + 10.^(squeeze(PP1_inc_dB2(:, :, ii, :, :))/10)) ...
    - reshape(repmat(KKi_dB, 4, 1), 2, 2, num_Th, num_RowPhis) ...
    - reshape(repmat(P1_areas_dB, 4, 1), 2, 2, num_Th, num_RowPhis) ;

    NBRCS_dB3(:, :, ii, :, :) = squeeze(PP1_inc_dB2(:, :, ii, :, :)) ...
    - reshape(repmat(KKi_dB, 4, 1), 2, 2, num_Th, num_RowPhis) ...
    - reshape(repmat(P1_areas_dB, 4, 1), 2, 2, num_Th, num_RowPhis) ;
end

REFL_dB_TOT = 10 * log10(10.^(P_cohveg_dB/10) + 10.^(PP1_inc_dB/10)) - reshape(repmat(KKc_dB, 4, 1), 2, 2, num_Th, num_RowPhis) ; 
REFL_dB = P_cohveg_dB - reshape(repmat(KKc_dB, 4, 1), 2, 2, num_Th, num_RowPhis) ; 

for ii = 1 : Nfz
    REFL_dB_TOT2(:, :, ii, :, :) = 10 * log10(10.^(P_cohveg_dB/10) + 10.^(squeeze(PP1_inc_dB2(:, :, ii, :, :))/10)) - reshape(repmat(KKc_dB, 4, 1), 2, 2, num_Th, num_RowPhis) ;
end


% TO-DO: Check this
%% GET GLOBAL PARAMETERS
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
ph0_Tx_deg = DynParams.getInstance.ph0_Tx_list_deg( ParamsManager.index_Ph );
VSM_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3( ParamsManager.index_VSM );
RMSH_cm = DynParams.getInstance.RMSH_list_cm( ParamsManager.index_RMSH );


%% FIGURES

%% Total NBRCS vs
for jj = 1 : num_RowPhis
        
    NBRCS_dB_TOT_RR = squeeze(NBRCS_dB_TOT(1, 1, :, jj)) ;
    NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, jj)) ;
    
    if jj == 1
        subplot(1,2,1)
        hold
        axis([0 80 -20 30])
        xlabel('\theta_s [\circ]')
        ylabel('\sigma^0_e [dB]')
        xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
        grid
        text(20, -10, strcat( num2str( rowPhis(1) ), '\circ \leq \phi_r_o_w \leq ', num2str( rowPhis(end) ), '\circ' ), 'Interpreter', 'tex')
        text(45, 20, strcat('fzone = 1'))
        title('Co-POL [RR]')
    end
    subplot(1,2,1)
    plot(th0_Tx_list_deg, NBRCS_dB_RR, '-sb', 'MarkerFaceColor', 'blue', 'MarkerSize', 2)
    
    subplot(1,2,1)
    plot(th0_Tx_list_deg, NBRCS_dB_TOT_RR, ':or', 'MarkerFaceColor', 'red', 'MarkerSize', 2)
    legend('Specular', 'Total')

end


%% Total NBRCS 2
for jj = 1 : Nfz
        
    NBRCS_dB_TOT_RR2 = squeeze(NBRCS_dB_TOT2(1, 1, jj, :, rowPhi_choice)) ;
    NBRCS_dB_RR2 = squeeze(NBRCS_dB3(1, 1, jj, :, rowPhi_choice)) ;    
    
    if jj == 1
        
        subplot(1,2,2)
        hold
        
        axis([0 80 -20 30])
        xlabel('\theta_s [\circ]')
        ylabel('\sigma^0_e [dB]')
        xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
        grid
        text(25, -5, '1\leq fzone\leq 10','Interpreter','tex')
        text(45, 20, strcat('\phi_r_o_w = ', num2str( rowPhis(rowPhi_choice) ), '\circ'))
        title('Co-POL [RR]')
        
    end
    subplot(1,2,2)
    plot(th0_Tx_list_deg, NBRCS_dB_RR2, '-sb', 'MarkerFaceColor', 'blue', 'MarkerSize', 2)
    
    subplot(1,2,2)
    plot(th0_Tx_list_deg, NBRCS_dB_TOT_RR2, ':or', 'MarkerFaceColor', 'red', 'MarkerSize', 2)
    legend('Diffuse', 'Total')

end


fname1 = strcat('NBRCS_vs_TH-RR-FZ_', num2str(FZ_choice), '-', pol_Tx, pol_Rx, '-PH_', num2str( ph0_Tx_deg ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end