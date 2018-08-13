
% Mehmet Kurum
% 11/15/2017
% Diffuse + Specular

function plotDielConst


%% ANALYSIS INPUTS
inputFile_sys_tag = 'sysInput-CORN_PAPER-row';
inputFile_veg_tag = 'vegVirRowInput-CORN_PAPER-row';

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
inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(rowPhi_choice) ), '.xml' );
inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(rowPhi_choice) ), '-', cornStages{1}, '.xml' );

getInput( inputFile_sys, inputFile_veg );

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();
% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    
ParamsManager.index_VSM( SM_choice );

for ss = 1 : num_stages

    %% GET INPUT
    inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(rowPhi_choice) ), '.xml' );
    inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(rowPhi_choice) ), '-', cornStages{ ss }, '.xml' );

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
              
    [dielStalk_RE_mu(ss), dielStalk_RE_sigma(ss), ...
     dielStalk_IM_mu(ss), dielStalk_IM_sigma(ss), ...
     dielLeaf_RE_mu(ss), dielLeaf_RE_sigma(ss), ... 
     dielLeaf_IM_mu(ss), dielLeaf_IM_sigma(ss), ...
     dielCob_RE_mu(ss), dielCob_RE_sigma(ss), ...
     dielCob_IM_mu(ss), dielCob_IM_sigma(ss) ] = plugin.getDiel();
    
    if ss == 1
        
        dir_analysis = SimulationFolders.getInstance.analysis;
        if ~exist(dir_analysis, 'dir')
            mkdir(dir_analysis)
        end
        
    end

end


%% Total Power for different Transmitter PhI angles

% Meeting
dielStalk_RE_mu = [49.6, 46.8, 45.7, 33, 37.9];
dielLeaf_RE_mu = [34.23, 34, 32.8, 25, 1.6];

% Display Dielectric Constants Real Part
figure
plot( dielStalk_RE_mu, '-bo', 'markersize', 6 )
hold on
plot( dielLeaf_RE_mu, '-g^', 'markersize', 6 )
plot( dielCob_RE_mu, '-rd', 'markersize', 6 )
axis([0.5 5.5 0 60])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
% xtickangle(30)
grid
title('Dielectric Contant (Real) vs. Growth Stages')
set(gca, 'ytick', 0:10:60);
ylabel('Dielectric Constant (Real)')

legend('Stalk', 'Leaves', 'Cobs')


%%
fname = strcat('AllConstituents-Diel(RE)_Meeting-vs_stages' ) ;

saveas(gcf, strcat(dir_analysis, '\', fname), 'tiff')
close

% Display Dielectric Constants Imaginary Part
figure
plot( dielStalk_IM_mu, '-bo', 'markersize', 6 )
hold on
plot( dielLeaf_IM_mu, '-g^', 'markersize', 6 )
plot( dielCob_IM_mu, '-rd', 'markersize', 6 )
axis([0.5 5.5 0 15])
xlabel('Growth Stages')
set(gca, 'xtick', 1:5, 'XtickLabel', strrep( cornStages, '_', '-' ) );
xtickangle(30)
grid
title('Dielectric Contant (Imaginary) vs. Growth Stages')
set(gca, 'ytick', 0:2.5:15);
ylabel('Dielectric Constant (Imaginary)')

legend('Stalk', 'Leaves', 'Cobs')


%%
fname = strcat('AllConstituents-Diel(IM)-vs_stages' ) ;

% saveas(gcf, strcat(dir_analysis, '\', fname), 'tiff')
close


end