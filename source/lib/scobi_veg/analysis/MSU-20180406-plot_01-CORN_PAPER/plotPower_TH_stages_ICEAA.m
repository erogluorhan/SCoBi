
% Mehmet Kurum
% 11/15/2017
% Diffuse + Specular

function plotPower_TH_stages_ICEAA


%% ANALYSIS INPUTS
inputFile_sys = 'sysInput-CORN_PAPER-row0.xml';
inputFile_veg_tag = 'vegVirRowInput-CORN_PAPER-row0-';


cornStages = { 'V1_V9', 'V10_VT', 'R1_R4', 'R5', 'R6' };
num_stages = length( cornStages );


%% Choices
% TO-DO: Check here
FZ_choice = 1 ;
ph_choice = 1;
VSM_choice = 4;
RMSH_choice = 2 ;
EIRP_choice_dB = 27;
G0r_choice_dB = 0;


%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_veg = strcat( inputFile_veg_tag, cornStages{1}, '.xml' );

getInput( inputFile_sys, inputFile_veg );


%% GET GLOBAL PARAMETERS
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg ;
num_TH = length( th0_Tx_list_deg );
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm ;
num_RMSH = length( RMSH_list_cm );


%% REFLECTIVITY

% Intilization 
%over vegetation
P_cohveg_dB = zeros(2, 2, num_stages, num_TH, num_RMSH);
P0_cohveg_dB = zeros(2, 2, num_stages, num_TH, num_RMSH);
% over bare soil
P_cohbare_dB = zeros(2, 2, num_stages, num_TH, num_RMSH);
P0_cohbare_dB = zeros(2, 2, num_stages, num_TH, num_RMSH);

for ss = 1 : num_stages
    
    inputFile_veg = strcat( inputFile_veg_tag, cornStages{ss}, '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
              
    % Assign the index of interest to each, in this analysis
    ParamsManager.index_VSM( VSM_choice );
    ParamsManager.index_Ph( ph_choice );
    
    for tt = 1 : num_TH

        % Set VSM index
        ParamsManager.index_Th( tt );
        
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
            

            P_cohveg_dB(1, :, ss, tt, rr) = 10 * log10( P_coh1vegx(1 : 2) ) + KKc_dB + EIRP_choice_dB + G0r_choice_dB;
            P0_cohveg_dB(1, :, ss, tt, rr) = 10 * log10( P0_coh1vegx(1 : 2) ) + KKc_dB + EIRP_choice_dB + G0r_choice_dB;
            P_cohveg_dB(2, :, ss, tt, rr) = 10 * log10( P_coh2vegx(1 : 2) ) + KKc_dB + EIRP_choice_dB + G0r_choice_dB;
            P0_cohveg_dB(2, :, ss, tt, rr) = 10 * log10( P0_coh2vegx(1 : 2) ) + KKc_dB + EIRP_choice_dB + G0r_choice_dB;

            P_cohbare_dB(1, :, ss, tt, rr) = 10 * log10( P_coh1barex(1 : 2) ) + KKc_dB + EIRP_choice_dB + G0r_choice_dB;
            P0_cohbare_dB(1, :, ss, tt, rr) = 10 * log10( P0_coh1barex(1 : 2) ) + KKc_dB + EIRP_choice_dB + G0r_choice_dB;
            P_cohbare_dB(2, :, ss, tt, rr) = 10 * log10( P_coh2barex(1 : 2) ) + KKc_dB + EIRP_choice_dB + G0r_choice_dB;
            P0_cohbare_dB(2, :, ss, tt, rr) = 10 * log10( P0_coh2barex(1 : 2) ) + KKc_dB + EIRP_choice_dB + G0r_choice_dB;
        
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
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
ph0_Tx_deg = ph0_Tx_list_deg( ph_choice );
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3 ;
VSM_cm3cm3 = VSM_list_cm3cm3( VSM_choice );
RMSH_cm = DynParams.getInstance.RMSH_list_cm( RMSH_choice );


%% FIGURES

%% Total Reflectivity
figure
pointSyms = {'-^b', '-db', '-vb', '-ob', '-sb' };
for ss = 1 : num_stages
        
    P_cohveg_dB_RL = squeeze(P_cohveg_dB(2, 1, ss, :, RMSH_choice));
    plot( th0_Tx_list_deg', P_cohveg_dB_RL, pointSyms{ss} );
    hold on

end

% Add Bare-soil
P_cohbare_dB_RL = squeeze(P_cohbare_dB(2, 1, 1, :, RMSH_choice));
plot( th0_Tx_list_deg', P_cohbare_dB_RL, '-*r' );


% Add labels etc.
grid
xlabel('\theta_s [\circ]')
xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
ylabel('Received Power [dB]')
axis( [0, 80, -180, -160] )

title('Received Power (Coherent) due to Growth Stages')


legend('V1-V9', 'V10-VT', 'R1-R4', 'R5', 'R6', 'Bare', 'Bare (Ideal)')


%%
fname1 = strcat('Power_vs_TH_stages(ICEAA)-', pol_Tx, pol_Rx, '-FZ_', num2str(FZ_choice), '-PH_', num2str( ph0_Tx_deg ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close
