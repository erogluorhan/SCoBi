
% Orhan Eroglu
% 5/4/2018

function plotStage_VWC


%% ANALYSIS INPUTS
inputFile_sys = 'sysInput-CORN_PAPER-row0.xml';
inputFile_veg_tag = 'vegVirRowInput-CORN_PAPER-row0-';

cornStages = { 'V1_V9', 'V10_VT', 'R1_R4', 'R5', 'R6' };
num_stages = length( cornStages );


%% Choices
% TO-DO: Check here
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
    ParamsManager.index_RMSH( RMSH_choice );

    % Initialize the directories depending on theta phi, VSM, and RMSH
    SimulationFolders.getInstance.initializeDynamicDirs();
    
    %% GET GLOBAL PARAMETERS
    % Virtual Row-Structured Vegetation Parameters
    plugin = VegVirRowParams.getInstance.plugin;

    [ vwcStalk_VolNormalized(ss), vwcLeaf_VolNormalized(ss), vwcCob_VolNormalized(ss) ] = plugin.getVWC_VolumeNormalized();
    
%     vwc_avg(ss) = plugin.getVWC_VolumeNormalized;
%     vwc(ss) = plugin.getVWC;
%     [diel_re(ss), diel_im(ss)] = plugin.getDiel;
%     diel_abs(ss) = abs( diel_re(ss) + 1i * diel_im(ss) );
    
    if ss == 1
        
        dir_analysis = SimulationFolders.getInstance.analysis;
        if ~exist(dir_analysis, 'dir')
            mkdir(dir_analysis)
        end
        
    end

end  

vwcCob_VolNormalized( isnan( vwcCob_VolNormalized ) ) = 0;


%% FIGURES

% Display Normalized VWC by Stalk Volume
figure
plot( vwcStalk_VolNormalized, '-b^', 'markersize', 8 )

grid
% axis([0.5 5.5 0 130])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
% set(gca, 'ytick', 0 : 20 : 130 );
ylabel('Normalized VWC (kg/m^3)')
title('Normalized Stalk VWC by Volume')


%%
fname = 'StalkVWCNormalized_vs_stages';

% saveas(gcf, strcat(dir_analysis, '\', fname), 'tiff')
close

% Display Normalized VWC by Leaf Volume
figure
plot( vwcLeaf_VolNormalized, '-b^', 'markersize', 8 )

grid
% axis([0.5 5.5 0 130])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
% set(gca, 'ytick', 0 : 20 : 130 );
ylabel('Normalized VWC (kg/m^3)')
title('Normalized Leaf VWC by Volume')


%%
fname = 'LeafVWCNormalized_vs_stages';

% saveas(gcf, strcat(dir_analysis, '\', fname), 'tiff')
close

% Display Normalized VWC by Cob Volume
figure
plot( vwcCob_VolNormalized, '-b^', 'markersize', 8 )

grid
% axis([0.5 5.5 0 130])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
% set(gca, 'ytick', 0 : 20 : 130 );
ylabel('Normalized VWC (kg/m^3)')
title('Normalized Cob VWC by Volume')


%%
fname = 'CobVWCNormalized_vs_stages';

% saveas(gcf, strcat(dir_analysis, '\', fname), 'tiff')
close


%% ALL TOGETHER
% Display Normalized VWC by Stalk Volume
figure
plot( vwcStalk_VolNormalized, '-bo', 'markersize', 8 )
hold on
plot( vwcLeaf_VolNormalized, '-g^', 'markersize', 8 )
plot( vwcCob_VolNormalized, '-rd', 'markersize', 8 )

grid
axis([0.5 5.5 0 1200])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
% set(gca, 'ytick', 0 : 20 : 130 );
ylabel('Normalized Moisture (kg/m^3)')
title('Constituent Moisture Content Normalized by Volume')
legend('Stalk', 'Leaves', 'Cobs')


%%
fname = 'AllConstituents-VWCNormalized_vs_stages';

saveas(gcf, strcat(dir_analysis, '\', fname), 'tiff')
close

end