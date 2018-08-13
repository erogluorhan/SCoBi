
% Mehmet Kurum
% 11/15/2017
% Diffuse + Specular

function plotVWC_Avg_stage


%% ANALYSIS INPUTS
inputFile_sys = 'sysInput-CORN_PAPER-row0.xml';
inputFile_veg_tag = 'vegVirRowInput-CORN_PAPER-row0-';

cornStages = { 'V1_V9', 'V10_VT', 'R1_R4', 'R5', 'R6' };
num_stages = length( cornStages );


%% Choices
% TO-DO: Check here
FZ_choice = 1 ;
TH_choice = 4 ;
PH_choice = 1 ;
SM_choice = 3 ;
RMSH_choice = 2 ;


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

        
    
ParamsManager.index_VSM( SM_choice );

for ss = 1 : num_stages

    %% GET INPUT
    inputFile_veg = strcat( inputFile_veg_tag, cornStages{ss}, '.xml' );

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
    
    %% GET GLOBAL PARAMETERS
    % Virtual Row-Structured Vegetation Parameters
    plugin = VegVirRowParams.getInstance.plugin;

    vwc_avg(ss) = plugin.getVWC_VolumeNormalized;
    vwc(ss) = plugin.getVWC;
    [diel_re(ss), diel_im(ss)] = plugin.getDiel;
    diel_abs(ss) = abs( diel_re(ss) + 1i * diel_im(ss) );
    
    if ss == 1
        
        dir_analysis = SimulationFolders.getInstance.analysis;
        if ~exist(dir_analysis, 'dir')
            mkdir(dir_analysis)
        end
        
    end

end


%% Total Power for different Transmitter PhI angles
        
ParamsManager.index_VSM( SM_choice );        


% Display VWC
figure
plot( vwc_avg, ':bo', 'markersize', 5 )

%axis([0.5 5.5 0 10])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
ylabel('VWC [kg/m^2]','Color','b')
%set(gca, 'ytick', 0 : 0.5 : 3);
title('Average VWC per volume')


% TO-DO: Check this
%% GET GLOBAL PARAMETERS
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
VSM_cm3cm3 = VSM_list_cm3cm3( ParamsManager.index_VSM );
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;
RMSH_cm = RMSH_list_cm( ParamsManager.index_RMSH );


%%
fname1 = strcat('ReceivedPower_VWC_vs_stages-FZ_', num2str(FZ_choice), '-', pol_Tx, pol_Rx, '-TH_', num2str( th0_Tx_list_deg(TH_choice) ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end