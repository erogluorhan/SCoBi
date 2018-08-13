
% Orhan Eroglu
% 5/4/2018

function plotStageAvgDimensions


%% ANALYSIS INPUTS
inputFile_sys = 'sysInput-CORN_PAPER-row0.xml';
inputFile_veg_tag = 'vegVirRowInput-CORN_PAPER-row0-';

cornStages = { 'V1_V9', 'V10_VT', 'R1_R4', 'R5', 'R6' };
num_stages = length( cornStages );

% Corn field row azimuth angles
rowPhis = [0 45 90] ;


%% Choices
% TO-DO: Check here
rowPhi_choice = 1;
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


stalkHeight_mu = [];
stalkBottomDiameter_mu = [];
leafCount_mu = [];
leafLength_mu = [];
leafWidth_mu = [];
cobLength_mu = []; 
cobDiameter_mu = [];

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

    [stalkHeight_mu(ss), stalkHeight_sigma(ss), stalkBottomDiameter_mu(ss), ...
     stalkBottomDiameter_sigma(ss), stalkTopDiameter_mu(ss), ...
     stalkTopDiameter_sigma(ss), leafCount_mu(ss), leafCount_sigma(ss), ...
     leafLength_mu(ss), leafLength_sigma(ss), leafWidth_mu(ss), ...
     leafWidth_sigma(ss), leafThickness_mu(ss), leafThickness_sigma(ss), ...
     cobLength_mu(ss), cobLength_sigma(ss), cobDiameter_mu(ss), ...
     cobDiameter_sigma(ss)] = plugin.getAvgParameters();
    
    if ss == 1
        
        dir_analysis = SimulationFolders.getInstance.analysis;
        if ~exist(dir_analysis, 'dir')
            mkdir(dir_analysis)
        end
        
    end

end      


%% FIGURES

%% STALK DIMENSIONS
% Display Average Stalk Height
stalkHeight_mu_cm = stalkHeight_mu * Constants.m2cm;
stalkHeight_sigma_cm = stalkHeight_sigma * Constants.m2cm;

figure
% plot( stalkHeight_mu_cm, '-b^', 'markersize', 6 )
errorbar( stalkHeight_mu_cm, stalkHeight_sigma_cm, '-b^', 'markersize', 6 )

% Display Average Stalk Diameter
stalkAvgDiameter_mu_mm = ( stalkBottomDiameter_mu + stalkTopDiameter_mu ) / 2 * Constants.m2mm;
stalkAvgDiameter_sigma_mm = ( stalkBottomDiameter_sigma + stalkTopDiameter_sigma ) / 2 * Constants.m2mm;
yyaxis right
% plot( stalkAvgDiameter_mu_mm, '-ro', 'markersize', 6 )
errorbar( stalkAvgDiameter_mu_mm, stalkAvgDiameter_sigma_mm, '-ro', 'markersize', 6 )
axis([0.5 5.5 0 20])
set(gca, 'ytick', 0 : 4 : 20, 'YColor', 'r' );
ylabel('Average Diameter (mm)')
% text(0.8, 13, strcat( 's=', num2str( sprintf( '%.1f', stalkAvgDiameter_sigma_mm(1) ) ), 'mm' ), 'Color', 'r')
% text(1.8, 17.8, strcat( 's=', num2str( sprintf( '%.1f', stalkAvgDiameter_sigma_mm(2) ) ), 'mm' ), 'Color', 'r')
% text(2.7, 11.8, strcat( 's=', num2str( sprintf( '%.1f', stalkAvgDiameter_sigma_mm(3) ) ), 'mm' ), 'Color', 'r')
% text(3.8, 8, strcat( 's=', num2str( sprintf( '%.1f', stalkAvgDiameter_sigma_mm(4) ) ), 'mm' ), 'Color', 'r')
% text(4.8, 9, strcat( 's=', num2str( sprintf( '%.1f', stalkAvgDiameter_sigma_mm(5) ) ), 'mm' ), 'Color', 'r')

yyaxis left
grid
axis([0.5 5.5 0 150])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
set(gca, 'ytick', 0 : 30 : 150, 'YColor', 'blue' );
ylabel('Height (cm)', 'Color', 'b')
title('Average Stalk Dimensions per Growth Stage')
% text(0.8, 17, strcat( 's=', num2str( sprintf( '%.1f', stalkHeight_sigma_cm(1) ) ), 'cm' ), 'Color', 'b')
% text(1.9, 57, strcat( 's=', num2str( sprintf( '%.1f', stalkHeight_sigma_cm(2) ) ), 'cm' ), 'Color', 'b')
% text(2.8, 125, strcat( 's=', num2str( sprintf( '%.1f', stalkHeight_sigma_cm(3) ) ), 'cm' ), 'Color', 'b')
% text(3.8, 133, strcat( 's=', num2str( sprintf( '%.1f', stalkHeight_sigma_cm(4) ) ), 'cm' ), 'Color', 'b')
% text(4.8, 105, strcat( 's=', num2str( sprintf( '%.1f', stalkHeight_sigma_cm(5) ) ), 'cm' ), 'Color', 'b')


%Save
% fname = 'StalkDimensions_vs_stages';
fname = 'StalkDimensions_vs_stages-errorbar';

saveas(gcf, strcat(dir_analysis, '\', fname), 'tiff')
close


%% LEAVE COUNT
% Display Average Number of Leaves
figure
plot( leafCount_mu, '-b*', 'markersize', 8 )

grid
axis([0.5 5.5 8 14])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
set(gca, 'ytick', 8 : 14 );
ylabel('Leaf Count')
title('Average Number of Leaves per Growth Stage')


%Save
fname = 'LeafCount_vs_stages';

% saveas(gcf, strcat(dir_analysis, '\', fname), 'tiff')
close


%% LEAF DIMENSIONS
% Display Average Leaf Length
leafLength_mu_cm = leafLength_mu * Constants.m2cm;
leafLength_sigma_cm = leafLength_sigma * Constants.m2cm;

figure
% plot( leafLength_mu_cm, '-bv', 'markersize', 6 )
errorbar( leafLength_mu_cm, leafLength_sigma_cm, '-bv', 'markersize', 6 )

% Display Average Leaf Width
leafWidth_mu_cm = leafWidth_mu * Constants.m2cm;
leafWidth_sigma_cm = leafWidth_sigma * Constants.m2cm;
yyaxis right
% plot( leafWidth_mu_cm, '-rs', 'markersize', 6 )
errorbar( leafWidth_mu_cm, leafWidth_sigma_cm, '-rs', 'markersize', 6 )
axis([0.5 5.5 0 20])
set(gca, 'ytick', [0 : 4 : 20], 'YColor', 'r' );
ylabel('Width (cm)','Color','r' )
% text(0.8, 4.6, strcat( 's=', num2str( sprintf( '%.1f', leafWidth_sigma_cm(1) ) ), 'cm' ), 'Color', 'r')
% text(1.8, 8.5, strcat( 's=', num2str( sprintf( '%.1f', leafWidth_sigma_cm(2) ) ), 'cm' ), 'Color', 'r')
% text(2.7, 8.7, strcat( 's=', num2str( sprintf( '%.1f', leafWidth_sigma_cm(3) ) ), 'cm' ), 'Color', 'r')
% text(3.8, 7.4, strcat( 's=', num2str( sprintf( '%.1f', leafWidth_sigma_cm(4) ) ), 'cm' ), 'Color', 'r')
% text(4.8, 6.5, strcat( 's=', num2str( sprintf( '%.1f', leafWidth_sigma_cm(5) ) ), 'cm' ), 'Color', 'r')

yyaxis left
grid
axis([0.5 5.5 0 100])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
set(gca, 'ytick', 0 : 20 : 100, 'YColor', 'b' );
ylabel('Length (cm)', 'Color', 'b' )
title('Average Leaf Dimensions per Growth Stage')
% text(0.7, 48, strcat( 's=', num2str( sprintf( '%.1f', leafLength_sigma_cm(1) ) ), 'cm' ), 'Color', 'b')
% text(1.7, 76, strcat( 's=', num2str( sprintf( '%.1f', leafLength_sigma_cm(2) ) ), 'cm' ), 'Color', 'b')
% text(2.8, 75, strcat( 's=', num2str( sprintf( '%.1f', leafLength_sigma_cm(3) ) ), 'cm' ), 'Color', 'b')
% text(3.8, 63, strcat( 's=', num2str( sprintf( '%.1f', leafLength_sigma_cm(4) ) ), 'cm' ), 'Color', 'b')
% text(4.8, 60, strcat( 's=', num2str( sprintf( '%.1f', leafLength_sigma_cm(5) ) ), 'cm' ), 'Color', 'b')


%%
% fname = 'LeafDimensions_vs_stages';
fname = 'LeafDimensions_vs_stages-errorbar';

% saveas(gcf, strcat(dir_analysis, '\', fname), 'tiff')
close


%% COB DIMENSIONS
% Display Average Cob Length
cobLength_mu_cm = cobLength_mu * Constants.m2cm;
cobLength_sigma_cm = cobLength_sigma * Constants.m2cm;

figure
% plot( cobLength_mu_cm, '-bd', 'markersize', 6 )
errorbar( cobLength_mu_cm, cobLength_sigma_cm, '-bd', 'markersize', 6 )

% Display Average Cob Diameter
cobDiameter_mu_mm = cobDiameter_mu * Constants.m2mm;
cobDiameter_sigma_mm = cobDiameter_sigma * Constants.m2mm;
yyaxis right
% plot( cobDiameter_mu_mm, '-rd', 'markersize', 6 )
errorbar( cobDiameter_mu_mm, cobDiameter_sigma_mm, '-rd', 'markersize', 6 )
axis([0.5 5.5 0 60])
set(gca, 'ytick', 0 : 10 : 60, 'YColor', 'r' );
ylabel('Diameter (mm)', 'Color', 'r' )
% text(0.75, 7, strcat( 's=', num2str( sprintf( '%.1f', cobDiameter_sigma_mm(1) ) ), 'mm' ), 'Color', 'r')
% text(2.22, 7, strcat( 's=', num2str( sprintf( '%.1f', cobDiameter_sigma_mm(2) ) ), 'mm' ), 'Color', 'r')
% text(2.9, 27, strcat( 's=', num2str( sprintf( '%.1f', cobDiameter_sigma_mm(3) ) ), 'mm' ), 'Color', 'r')
% text(3.8, 39, strcat( 's=', num2str( sprintf( '%.1f', cobDiameter_sigma_mm(4) ) ), 'mm' ), 'Color', 'r')
% text(4.77, 44, strcat( 's=', num2str( sprintf( '%.1f', cobDiameter_sigma_mm(5) ) ), 'mm' ), 'Color', 'r')

yyaxis left
grid
axis([0.5 5.5 0 30])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
set(gca, 'ytick', 0 : 5 : 30, 'YColor', 'b' );
ylabel('Length (cm)', 'Color', 'b' )
title('Average Cob Dimensions per Growth Stage')
% text(0.75, 2, strcat( 's=', num2str( sprintf( '%.1f', cobLength_sigma_cm(1) ) ), 'cm' ), 'Color', 'b')
% text(2.12, 2, strcat( 's=', num2str( sprintf( '%.1f', cobLength_sigma_cm(2) ) ), 'cm' ), 'Color', 'b')
% text(2.75, 23.5, strcat( 's=', num2str( sprintf( '%.1f', cobLength_sigma_cm(3) ) ), 'cm' ), 'Color', 'b')
% text(3.8, 18.8, strcat( 's=', num2str( sprintf( '%.1f', cobLength_sigma_cm(4) ) ), 'cm' ), 'Color', 'b')
% text(4.78, 11, strcat( 's=', num2str( sprintf( '%.1f', cobLength_sigma_cm(5) ) ), 'cm' ), 'Color', 'b')


%Save
% fname = 'CobDimensions_vs_stages';
fname = 'CobDimensions_vs_stages-errorbar';

saveas(gcf, strcat(dir_analysis, '\', fname), 'tiff')
close


end