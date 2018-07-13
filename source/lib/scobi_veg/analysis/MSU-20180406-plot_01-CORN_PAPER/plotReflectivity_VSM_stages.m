
% Mehmet Kurum
% 11/15/2017
% Diffuse + Specular

function plotReflectivity_VSM_stages


%% ANALYSIS INPUTS
inputFile_sys = 'sysInput-CORN_PAPER-row0.xml';
inputFile_veg_tag = 'vegVirRowInput-CORN_PAPER-row0-';


cornStages = { 'V1_V9', 'V10_VT', 'R1_R4', 'R5', 'R6' };
num_stages = length( cornStages );


%% Choices
% TO-DO: Check here
FZ_choice = 1 ;
TH_choice = 3;
PH_choice = 1;
RMSH_choice = 2 ;
stage_choice = 4;


%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_veg = strcat( inputFile_veg_tag, cornStages{1}, '.xml' );

getInput( inputFile_sys, inputFile_veg );


%% GET GLOBAL PARAMETERS
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg ;
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3 ;
num_VSM = length( VSM_list_cm3cm3 );
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm ;
num_RMSH = length( RMSH_list_cm );


%% REFLECTIVITY

% Intilization 
%over vegetation
P_cohveg_dB = zeros(2, 2, num_stages, num_VSM, num_RMSH);
P0_cohveg_dB = zeros(2, 2, num_stages, num_VSM, num_RMSH);
% over bare soil
P_cohbare_dB = zeros(2, 2, num_stages, num_VSM, num_RMSH);
P0_cohbare_dB = zeros(2, 2, num_stages, num_VSM, num_RMSH);

for ss = 1 : num_stages
    
    inputFile_veg = strcat( inputFile_veg_tag, cornStages{ss}, '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
              
    % Assign the index of interest to each, in this analysis
    ParamsManager.index_Th( TH_choice );
    ParamsManager.index_Ph( PH_choice );
    
    for vv = 1 : num_VSM

        % Set VSM index
        ParamsManager.index_VSM( vv );
        
        for rr = 1 : num_RMSH

            % Set VSM index
            ParamsManager.index_RMSH( rr );
        
            % Initialize the directories depending on theta phi, VSM, and RMSH
            SimulationFolders.getInstance.initializeDynamicDirs();


            %% GET GLOBAL DIRECTORIES
            dir_out_specular = SimulationFolders.getInstance.out_specular;

            % read Kc
            filename1 = strcat('Kc') ;
            Kc = readComplexVar( dir_out_specular, filename1) ;
            KKc_dB = 20 * log10(abs(Kc)) ;

            [~, ~, ~, ~, ~, ~, ~, ~, ...
             P_coh1vegx, P_coh2vegx, P0_coh1vegx, P0_coh2vegx, ...
             P_coh1barex, P_coh2barex, P0_coh1barex, P0_coh2barex] = readSpecular ;
            

            P_cohveg_dB(1, :, ss, vv, rr) = 10 * log10( P_coh1vegx(1 : 2) ) + KKc_dB;
            P0_cohveg_dB(1, :, ss, vv, rr) = 10 * log10( P0_coh1vegx(1 : 2) ) + KKc_dB;
            P_cohveg_dB(2, :, ss, vv, rr) = 10 * log10( P_coh2vegx(1 : 2) ) + KKc_dB;
            P0_cohveg_dB(2, :, ss, vv, rr) = 10 * log10( P0_coh2vegx(1 : 2) ) + KKc_dB;

            P_cohbare_dB(1, :, ss, vv, rr) = 10 * log10( P_coh1barex(1 : 2) ) + KKc_dB;
            P0_cohbare_dB(1, :, ss, vv, rr) = 10 * log10( P0_coh1barex(1 : 2) ) + KKc_dB;
            P_cohbare_dB(2, :, ss, vv, rr) = 10 * log10( P_coh2barex(1 : 2) ) + KKc_dB;
            P0_cohbare_dB(2, :, ss, vv, rr) = 10 * log10( P0_coh2barex(1 : 2) ) + KKc_dB;
        
        end
        
    end
    
    if ss == 1
        
        dir_analysis = SimulationFolders.getInstance.analysis;
        if ~exist(dir_analysis, 'dir')
            mkdir(dir_analysis)
        end
    end
    
end


%% NORMALIZATION
REFL_dB = P_cohveg_dB - KKc_dB; 
REFL0_dB = P0_cohveg_dB - KKc_dB; 
REFL_bare_dB = P_cohbare_dB - KKc_dB; 
REFL0_bare_dB = P0_cohbare_dB - KKc_dB; 


% TO-DO: Check this
%% GET GLOBAL PARAMETERS
% Transmiiter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
ph0_Tx_deg = ph0_Tx_list_deg( PH_choice );
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;
RMSH_cm = RMSH_list_cm( RMSH_choice );


%% FIGURES

%% Total Reflectivity
figure
pointSyms = {'-^b', '-db', '-vb', '-ob', '-sb' };
for ss = 1 : num_stages
        
    subplot(1,2,1)
    REFL_dB_RL = squeeze(REFL_dB(2, 1, ss, :, RMSH_choice));
    plot( VSM_list_cm3cm3', REFL_dB_RL, pointSyms{ss} );
    hold on
    
    subplot(1,2,2)
    REFL_dB_RR = squeeze(REFL_dB(1, 1, ss, :, RMSH_choice));
    plot( VSM_list_cm3cm3', REFL_dB_RR, pointSyms{ss} );
    hold on

end

% Add Bare-soil
subplot(1,2,1)
REFL_bare_dB_RL = squeeze(REFL_bare_dB(2, 1, 1, :, RMSH_choice));
plot( VSM_list_cm3cm3', REFL_bare_dB_RL, '-*r' );

subplot(1,2,2)
REFL_bare_dB_RR = squeeze(REFL_bare_dB(1, 1, 1, :, RMSH_choice));
plot( VSM_list_cm3cm3', REFL_bare_dB_RR, '-*r' );
REFL0_bare_dB_RR = squeeze(REFL0_bare_dB(1, 1, 1, :, RMSH_choice));
plot( VSM_list_cm3cm3', REFL0_bare_dB_RR, ':*r' );

% Add labels etc.
subplot(1,2,1)
grid
xlabel('VSM [cm^3/cm^3]')
ylabel('Reflectivity [dB]')
axis( [0, 0.45, -25, 0] )

text(0.05, -3.5, strcat( '\theta_s=', num2str( th0_Tx_list_deg(TH_choice) ), '\circ' ))
text(0.05, -2, strcat( 'RMSH=', num2str(RMSH_list_cm(RMSH_choice)), ' cm'));
title('X-Pol [RL]')

subplot(1,2,2)
grid
xlabel('VSM [cm^3/cm^3]')
axis( [0, 0.45, -50, 0] )
title('Co-Pol [RR]')

legend('V1-V9', 'V10-VT', 'R1-R4', 'R5', 'R6', 'Bare', 'Bare (Ideal)')


%%
fname1 = strcat('Reflectivity_vs_VSM_stages(plusIdeal)-', pol_Tx, pol_Rx, '-FZ_', num2str(FZ_choice), '-TH_', num2str( th0_Tx_list_deg(TH_choice) ), '-PH_', num2str( ph0_Tx_deg ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close



%% COLLORBAR 
REFL_dB_RL = [];
REFL_dB_RR = [];
stages_sct = [];
VSM_sct = [];

%% Total Reflectivity
for ss = 1 : num_stages
    
    stages_sct = [ stages_sct; ss * ones(num_VSM,1) ];
    
    VSM_sct = [ VSM_sct; VSM_list_cm3cm3' ];
    
    REFL_dB_RR = [ REFL_dB_RR; squeeze(REFL_dB(1, 1, ss, :, RMSH_choice)) ] ;
    REFL_dB_RL = [ REFL_dB_RL; squeeze(REFL_dB(2, 1, ss, :, RMSH_choice)) ] ;

end

% Add Bare-soil
stages_sct = [ zeros(num_VSM,1); stages_sct ];
VSM_sct = [ VSM_sct; VSM_list_cm3cm3' ];
REFL_dB_RR = [ squeeze(REFL_bare_dB(1, 1, 1, :, RMSH_choice)); REFL_dB_RR ] ;
REFL_dB_RL = [ squeeze(REFL_bare_dB(2, 1, 1, :, RMSH_choice)); REFL_dB_RL ] ;

figure
subplot(1,2,1)
scatter( VSM_sct, REFL_dB_RL, 100, stages_sct, 'filled' );
xlabel('VSM (cm^3/cm^3)')
ylabel('Reflectivity (dB)')
axis( [0, 0.45, -19, -4] )
% yyaxis right
% cb = colorbar;
% ylabel(cb, 'Growth Stages')
% set(gca, 'ytick', 1:5, 'YtickLabel', strrep( cornStages, '_', '-' ) );
% yyaxis left
text(0.25, -16, strcat( '\theta_s=', num2str( th0_Tx_list_deg(TH_choice) ), '\circ' ))
text(0.25, -17, strcat( 'RMSH=', num2str(RMSH_list_cm(RMSH_choice)), ' cm'));
title('X-Pol [RL]')

subplot(1,2,2)
scatter( VSM_sct, REFL_dB_RR, 100, stages_sct, 'filled' );
xlabel('VSM (cm^3/cm^3)')
axis( [0, 0.45, -48, -18] )
%yyaxis right
cb = colorbar;
ylabel(cb, 'Growth Stages')
cb.Ticks = 0 : 5 ; %Create 8 ticks from zero to 1
cornStages_bare = { 'Bare', 'V1_V9', 'V10_VT', 'R1_R4', 'R5', 'R6' };
cb.TickLabels = strrep( cornStages_bare, '_', '-' ) ;    %Replace the labels of these 8 ticks with the numbers 1 to 8
%set(gca, 'ytick', 1:5, 'YtickLabel', strrep( cornStages, '_', '-' ) );
%yyaxis left
title('Co-Pol [RR]')


%%
fname1 = strcat('Reflectivity_vs_VSM_stages-colorbar-', pol_Tx, pol_Rx, '-FZ_', num2str(FZ_choice), '-TH_', num2str( th0_Tx_list_deg(TH_choice) ), '-PH_', num2str( ph0_Tx_deg ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

% saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


%% ERRORBAR
figure
pointSyms = {'-sg', '-sr', '-sb', '-sc', '-sm' };
markerColors = {'green', 'red', 'blue', 'cyan', 'magenta' };
for ss = 1 : num_stages
    
    REFL_dB_RL_plot = squeeze(REFL_dB(2, 1, ss, :, 2));   
    REFL_dB_RL_max = squeeze(REFL_dB(2, 1, ss, :, 1)); 
    REFL_dB_RL_min = squeeze(REFL_dB(2, 1, ss, :, num_RMSH));
    REFL_dB_RL_pos = REFL_dB_RL_max - REFL_dB_RL_plot;
    REFL_dB_RL_neg = REFL_dB_RL_plot - REFL_dB_RL_min;
    
    REFL_dB_RR_plot = squeeze(REFL_dB(1, 1, ss, :, 2));
    REFL_dB_RR_max = squeeze(REFL_dB(1, 1, ss, :, 1));
    REFL_dB_RR_min = squeeze(REFL_dB(1, 1, ss, :, num_RMSH));
    REFL_dB_RR_pos = REFL_dB_RR_max - REFL_dB_RR_plot;
    REFL_dB_RR_neg = REFL_dB_RR_plot - REFL_dB_RR_min;
    
    subplot(1,2,1)
    hold on
%     errorbar( VSM_list_cm3cm3, REFL_dB_RL_plot, REFL_dB_RL_pos, REFL_dB_RL_neg, pointSyms{ss},'MarkerSize',5,...
%     'MarkerEdgeColor', markerColors{ss},'MarkerFaceColor', markerColors{ss});
    errorbar( VSM_list_cm3cm3, REFL_dB_RL_plot, REFL_dB_RL_pos, REFL_dB_RL_neg);
    
    subplot(1,2,2)
    hold on
    errorbar( VSM_list_cm3cm3, REFL_dB_RR_plot, REFL_dB_RR_pos, REFL_dB_RR_neg);


    if ss == 1
        subplot(1,2,1)
        xlabel('VSM (cm^3/cm^3)')
        ylabel('Reflectivity (dB)')
        axis([0, 0.45, -18, -3])
        title('X-Pol [RL]')
        
        subplot(1,2,2)
        xlabel('VSM (cm^3/cm^3)')
        axis([0, 0.45, -28, -15])
        title('Co-Pol [RR]')
    end

end

subplot(1,2,1)
RMSH_list_str = sprintf('%.1f,' , RMSH_list_cm);
RMSH_list_str = RMSH_list_str( 1 : end-1 );% strip final comma
text( 0.05, -3.5, strcat( 'RMSH = \{', RMSH_list_str, '\} cm' ), 'Interpreter', 'tex')
legend( strrep( cornStages, '_', '-' )  )


%%
fname1 = strcat('Reflectivity_vs_VSM_RMSH_stages-errorbar-', pol_Tx, pol_Rx, '-FZ_', num2str(FZ_choice), '-TH_', num2str( th0_Tx_list_deg(TH_choice) ), '-PH_', num2str( ph0_Tx_deg ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

% saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


%% FOR RMSH
REFL_dB_RL_plot = squeeze(REFL_dB(2, 1, stage_choice, :, 2));   
REFL_dB_RL_max = squeeze(REFL_dB(2, 1, stage_choice, :, 1)); 
REFL_dB_RL_min = squeeze(REFL_dB(2, 1, stage_choice, :, num_RMSH));
REFL_dB_RL_pos = REFL_dB_RL_max - REFL_dB_RL_plot;
REFL_dB_RL_neg = REFL_dB_RL_plot - REFL_dB_RL_min;

REFL_dB_RR_plot = squeeze(REFL_dB(1, 1, stage_choice, :, 2));
REFL_dB_RR_max = squeeze(REFL_dB(1, 1, stage_choice, :, 1));
REFL_dB_RR_min = squeeze(REFL_dB(1, 1, stage_choice, :, num_RMSH));
REFL_dB_RR_pos = REFL_dB_RR_max - REFL_dB_RR_plot;
REFL_dB_RR_neg = REFL_dB_RR_plot - REFL_dB_RR_min;
    
figure
errorbar( VSM_list_cm3cm3, REFL_dB_RL_plot, REFL_dB_RL_pos, REFL_dB_RL_neg, '-sb','MarkerSize',5,...
'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue' );

hold on
errorbar( VSM_list_cm3cm3, REFL_dB_RR_plot, REFL_dB_RR_pos, REFL_dB_RR_neg, '-sr','MarkerSize',5,...
'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red' );

xlabel('VSM [cm^3/cm^3]')
ylabel('Reflectivity [dB]')
axis([0, 0.45, -25, 0])
RMSH_list_str = sprintf('%.1f,' , RMSH_list_cm);
RMSH_list_str = RMSH_list_str( 1 : end-1 );% strip final comma
text( 0.05, -3, strcat( 'RMSH = \{', RMSH_list_str, '\} cm' ), 'Interpreter', 'tex')
text( 0.05, -7, strcat( '\theta_s=', num2str( th0_Tx_list_deg(TH_choice) ), '\circ' ))
text( 0.02, REFL_dB_RL_plot(1), 'RL' )
text( 0.02, REFL_dB_RR_plot(1), 'RR' )
grid
title('Reflectivity Variation due to VSM and RMSH')


inputFile_veg = strcat( inputFile_veg_tag, cornStages{stage_choice}, '.xml' );

getInput( inputFile_sys, inputFile_veg );
% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();

dir_analysis = SimulationFolders.getInstance.analysis;
if ~exist(dir_analysis, 'dir')
    mkdir(dir_analysis)
end


%%
fname1 = strcat('Reflectivity_vs_VSM_RMSH-errorbar-', cornStages{stage_choice}, '-', pol_Tx, pol_Rx, '-FZ_', num2str(FZ_choice), '-TH_', num2str( th0_Tx_list_deg(TH_choice) ), '-PH_', num2str( ph0_Tx_deg ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

% saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close
