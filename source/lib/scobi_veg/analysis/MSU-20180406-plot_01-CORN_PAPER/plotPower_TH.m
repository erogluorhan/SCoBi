
% Mehmet Kurum
% 11/15/2017
% Diffuse + Specular

function plotPower_TH


%% ANALYSIS INPUTS
inputFile_sys_tag = 'sysInput-CORN_PAPER-row0';
inputFile_veg_tag = 'vegVirRowInput-CORN_PAPER-row0';

cornStagesStruct = struct('s1', 'V1_V9', ...
                            's2', 'V10_VT', ...
                            's3', 'R1_R4', ...
                            's4', 'R5', ...
                            's5', 'R6') ;

% Corn field row azimuth angles
rowPhis = [0 30 45 60 90] ;
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
inputFile_sys = strcat( inputFile_sys_tag, '.xml' );
inputFile_veg = strcat( inputFile_veg_tag, '-', stage_choice, '.xml' );

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
        
    inputFile_sys = strcat( inputFile_sys_tag, '.xml' );
    inputFile_veg = strcat( inputFile_veg_tag, '-', stage_choice, '.xml' );
        
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
P0_cohveg_dB = zeros(2, 2, num_Th, num_RowPhis) ;
% over bare soil
P_cohbare_dB = zeros(2, 2, num_Th, num_RowPhis) ;
P0_cohbare_dB = zeros(2, 2, num_Th, num_RowPhis) ;

KKc_dB = zeros(num_Th, num_RowPhis) ;

for jj = 1 : num_RowPhis
    
    inputFile_sys = strcat( inputFile_sys_tag, '.xml' );
    inputFile_veg = strcat( inputFile_veg_tag, '-', stage_choice, '.xml' );
        
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
            P_coh1vegx, P_coh2vegx, P0_coh1vegx, P0_coh2vegx, ...
            P_coh1barex, P_coh2barex, P0_coh1barex, P0_coh2barex] = readSpecular ;
        
        P_cohveg_dB(1, :, tt, jj) = 10 * log10( P_coh1vegx(1 : 2) ) + KKc_dB(tt) ;
        P0_cohveg_dB(1, :, tt, jj) = 10 * log10( P0_coh1vegx(1 : 2) ) + KKc_dB(tt) ;
        P_cohveg_dB(2, :, tt, jj) = 10 * log10( P_coh2vegx(1 : 2) ) + KKc_dB(tt) ;
        P0_cohveg_dB(2, :, tt, jj) = 10 * log10( P0_coh2vegx(1 : 2) ) + KKc_dB(tt) ;
        
        P_cohbare_dB(1, :, tt, jj) = 10 * log10( P_coh1barex(1 : 2) ) + KKc_dB(tt) ;
        P0_cohbare_dB(1, :, tt, jj) = 10 * log10( P0_coh1barex(1 : 2) ) + KKc_dB(tt) ;
        P_cohbare_dB(2, :, tt, jj) = 10 * log10( P_coh2barex(1 : 2) ) + KKc_dB(tt) ;
        P0_cohbare_dB(2, :, tt, jj) = 10 * log10( P0_coh2barex(1 : 2) ) + KKc_dB(tt) ;
        
    end

    if jj == 1
        dir_analysis = SimulationFolders.getInstance.analysis;
        if ~exist(dir_analysis, 'dir')
            mkdir(dir_analysis)
        end
    end
    
end


% TO-DO: Check this
%% GET GLOBAL PARAMETERS
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
PH0_deg = ph0_Tx_list_deg( ParamsManager.index_Ph );
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
VSM_cm3cm3 = VSM_list_cm3cm3( ParamsManager.index_VSM );
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;
RMSH_cm = RMSH_list_cm( ParamsManager.index_RMSH );


%% FIGURE

%% Total Power for different row-phi angles
    
inputFile_sys = strcat( inputFile_sys_tag, '.xml' );
inputFile_veg = strcat( inputFile_veg_tag, '-', stage_choice, '.xml' );

getInput( inputFile_sys, inputFile_veg );
% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();
% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();


Pvi_RR = squeeze(PP1_inc_dB(1, 1, :, rowPhi_choice)) + EIRP_choice_dB + G0r_choice_dB ;
Pvi_RL = squeeze(PP1_inc_dB(1, 2, :, rowPhi_choice)) + EIRP_choice_dB + G0r_choice_dB  ;    
Pvc_RR = squeeze(P_cohveg_dB(1, 1, :, rowPhi_choice)) + EIRP_choice_dB + G0r_choice_dB  ; 
P0vc_RR = squeeze(P0_cohveg_dB(1, 1, :, rowPhi_choice)) + EIRP_choice_dB + G0r_choice_dB  ;
Pvc_RL = squeeze(P_cohveg_dB(1, 2, :, rowPhi_choice)) + EIRP_choice_dB + G0r_choice_dB  ;    
Pbc_RR = squeeze(P_cohbare_dB(1, 1, :, rowPhi_choice)) + EIRP_choice_dB + G0r_choice_dB  ;  
P0bc_RR = squeeze(P0_cohbare_dB(1, 1, :, rowPhi_choice)) + EIRP_choice_dB + G0r_choice_dB  ;
Pbc_RL = squeeze(P_cohbare_dB(1, 2, :, rowPhi_choice)) + EIRP_choice_dB + G0r_choice_dB  ;    

subplot(1,2,1)
plot(th0_Tx_list_deg, Pbc_RL, '-ob', 'MarkerFaceColor', 'blue', 'markersize', 4)
hold
plot(th0_Tx_list_deg, Pvc_RL, '-dr', 'MarkerFaceColor', 'red', 'markersize', 5)

minPow_dB = -237 + EIRP_choice_dB + G0r_choice_dB;
maxPow_dB = -177 + EIRP_choice_dB + G0r_choice_dB;        
axis([0 80 minPow_dB maxPow_dB])
xlabel('\theta_s [\circ]')
ylabel('Received Power [dB]')
xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
yticks(minPow_dB : 5 : maxPow_dB)
grid
%text(72, Pvi_RL(end) + 3, strcat('\phi_r_o_w [\circ]'))
title('X-Pol [RL]')

subplot(1,2,2)
plot(th0_Tx_list_deg, Pbc_RR, '-ob', 'MarkerFaceColor', 'blue', 'markersize', 4)
hold on
plot(th0_Tx_list_deg, P0bc_RR, ':ob', 'markersize', 4)
plot(th0_Tx_list_deg, Pvc_RR, '-dr', 'MarkerFaceColor', 'red', 'markersize', 5)
plot(th0_Tx_list_deg, P0vc_RR, ':dr', 'markersize', 5)

minPow_dB = -237 + EIRP_choice_dB + G0r_choice_dB;
maxPow_dB = -177 + EIRP_choice_dB + G0r_choice_dB;
axis([0, 80, minPow_dB,  maxPow_dB ])
xlabel('\theta_s [\circ]')
xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
yticks(minPow_dB : 5 : maxPow_dB)
grid
%text(72, Pvi_RR(end) + 3, strcat('\phi_r_o_w [\circ]'))
title('Co-Pol [RR]')

subplot(1,2,1)
plot(th0_Tx_list_deg, Pvi_RL, '-^g', 'MarkerFaceColor', 'green', 'markersize', 5)

text(10, -152, strcat( 'VSM=', num2str( VSM_cm3cm3 ), ' cm^3/cm^3' ))
text(10, -156, strcat( 'RMSH=', num2str( RMSH_cm ), ' cm' ))
text(10, -160, strcat('fzone = ', num2str(FZ_choice)))

subplot(1,2,2)
plot(th0_Tx_list_deg, Pvi_RR, '-^g', 'MarkerFaceColor', 'green', 'markersize', 5)
%text( 75, Pvi_RR(end), strcat( num2str( rowPhis(jj) ) ) )
legend('Coherent-Bare Soil', 'Coherent-Bare (Ideal)', 'Coherent-Vegetation', 'Coherent-Veg. (Ideal)', 'Diffuse')


%% GET GLOBAL PARAMETERS
% Vegetation Parameters
vegetation_stage = VegParams.getInstance.vegetation_stage;


%%
fname1 = strcat('ReceivedPower(plusIdeals)_vs_Th-', vegetation_stage, '-FZ_', num2str(FZ_choice), '-', pol_Tx, pol_Rx, '-PH_', num2str( PH0_deg ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end