
function plotReflectivity_VWC_VSM_TH_RMSH


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
TH_choice = 4;
PH_choice = 1;
RMSH_choice = 2 ;
stage_choice = 2;


%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(1) ), '.xml' );
inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(1) ), '-', cornStages{1}, '.xml' );

getInput( inputFile_sys, inputFile_veg );


%% GET GLOBAL PARAMETERS
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg ;
num_TH = length( th0_Tx_list_deg );
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3 ;
num_VSM = length( VSM_list_cm3cm3 );
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm ;
num_RMSH = length( RMSH_list_cm );


%% REFLECTIVITY

% Intilization 
%over vegetation
P_cohveg_dB = zeros(2, 2, num_stages, num_VSM, num_RMSH);
% over bare soil
P_cohbare_dB = zeros(2, 2, num_stages, num_VSM, num_RMSH);

for ss = 1 : num_stages
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(rowPhi_choice) ), '.xml' );
    inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(rowPhi_choice) ), '-', cornStages{ss}, '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    %% GET GLOBAL PARAMETERS
    % Virtual Row-Structured Vegetation Parameters
    plugin = VegVirRowParams.getInstance.plugin;

    vwc(ss) = plugin.getVWC;
        
    for tt = 1 : num_TH
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
                    P_coh1vegx, P_coh2vegx, ~, ~, ...
                    P_coh1barex, P_coh2barex, ~, ~] = readSpecular ;

    %             % Ideal
    %             [~, ~, ~, ~, ~, ~, ~, ~, ...
    %                 ~, ~, P_coh1vegx, P_coh2vegx, ...
    %                 ~, ~, P_coh1barex, P_coh2barex] = readSpecular ;

                P_cohveg_dB(1, :, ss, tt, vv, rr) = 10 * log10( P_coh1vegx(1 : 2) ) + KKc_dB;
                P_cohveg_dB(2, :, ss, tt, vv, rr) = 10 * log10( P_coh2vegx(1 : 2) ) + KKc_dB;

                P_cohbare_dB(1, :, ss, tt, vv, rr) = 10 * log10( P_coh1barex(1 : 2) ) + KKc_dB;
                P_cohbare_dB(2, :, ss, tt, vv, rr) = 10 * log10( P_coh2barex(1 : 2) ) + KKc_dB;

            end

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
REFL_bare_dB = P_cohbare_dB - KKc_dB; 


% TO-DO: Check this
%% GET GLOBAL PARAMETERS
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Recevier Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
ph0_Tx_deg = ph0_Tx_list_deg( PH_choice );
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;
RMSH_cm = RMSH_list_cm( RMSH_choice );


%% FIGURES

%% COLLORBAR 
% REFL_dB_RL = [];
% REFL_dB_RR = [];
% stages_sct = [];
VSM_sct = [];
VWC_sct = [];

%% Total Reflectivity
% for tt = 1 : num_TH
%     
%     for rr = 1 : num_RMSH
        
        VSM_temp = [];
        for ss = 1 : num_stages
            VSM_temp = [ VSM_temp; VSM_list_cm3cm3' ];
        end

        VSM_sct = [ VSM_sct; VSM_temp ];
        
        for vv = 1 : num_VSM            

%         	VWC_sct = [ VWC_sct; vwc' ];

%             REFL_dB_RR(tt, vv, rr, :) = squeeze(REFL_dB(1, 1, :, tt, vv, rr)) ;
%             REFL_dB_RL(tt, vv, rr, :) = squeeze(REFL_dB(2, 1, :, tt, vv, rr)) ;

            REFL_dB_RR(vv, :) = squeeze(REFL_dB(1, 1, :, TH_choice, vv, RMSH_choice)) ;
            REFL_dB_RL(vv, :) = squeeze(REFL_dB(2, 1, :, TH_choice, vv, RMSH_choice)) ;
        
        end

%     end
% 
% end

vwc = repmat(vwc, num_VSM, 1)
VWC_sct = reshape(vwc, [num_VSM * num_stages, 1]);

% REFL_dB_RR_sct = reshape(REFL_dB_RR, [num_TH * num_VSM * num_RMSH * num_stages, 1]);
% REFL_dB_RL_sct = reshape(REFL_dB_RL, [num_TH * num_VSM * num_RMSH * num_stages, 1]);
REFL_dB_RR_sct = reshape(REFL_dB_RR, [num_VSM * num_stages, 1]);
REFL_dB_RL_sct = reshape(REFL_dB_RL, [num_VSM * num_stages, 1]);

% % Add Bare-soil
% stages_sct = [ zeros(num_VSM,1); stages_sct ];
% VSM_sct = [ VSM_sct; VSM_list_cm3cm3' ];
% REFL_dB_RR = [ squeeze(REFL_bare_dB(1, 1, 1, :, RMSH_choice)); REFL_dB_RR ] ;
% REFL_dB_RL = [ squeeze(REFL_bare_dB(2, 1, 1, :, RMSH_choice)); REFL_dB_RL ] ;

figure
%subplot(1,2,1)
% scatter( VWC_sct, REFL_dB_RL_sct, 100, stages_sct, 'filled' );
scatter( VWC_sct, REFL_dB_RL_sct, 100, VSM_sct, 'filled' );
xlabel('VWC (kg/m^2)')
ylabel('Reflectivity (dB)')
axis( [0, 2.5, -20, 0] )
% yyaxis right
% cb = colorbar;
% ylabel(cb, 'Growth Stages')
% set(gca, 'ytick', 1:5, 'YtickLabel', strrep( cornStages, '_', '-' ) );
% yyaxis left
% text(0.25, -16, strcat( '\theta_s=', num2str( th0_Tx_list_deg(TH_choice) ), '\circ' ))
% text(0.25, -17, strcat( 'RMSH=', num2str(RMSH_list_cm(RMSH_choice)), ' cm'));
title('X-Pol [RL]')
cb = colorbar;
ylabel(cb, 'VSM (cm^3/cm^3)')


%%
fname1 = strcat('Reflectivity_vs_VWC-Scatter2-', pol_Tx, pol_Rx, '-FZ_', num2str(FZ_choice), '-TH_', num2str( th0_Tx_list_deg(TH_choice) ), '-PH_', num2str( ph0_Tx_deg ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close

end
