
% Mehmet Kurum
% 11/12/2017
% NBRCS vs. Fresnel Zones


function plotDiffuse_NBRCS_FZ


%% INPUT FILES
inputFile_sys_tag = 'sysInput-Paulownia-PAPER_PBAND-hr_';
inputFile_veg = 'vegHomInput-Paulownia.xml';


%%
Nfz = 10;
SM_choice = 2 ;
ANG_choice = 4 ;


%% Reciever height
hhr = [20 50 100 500] ;


%% reading
% NBRCS
Nhr =  length(hhr) ;

NBRCS_dB = zeros(2, 2, Nfz, Nhr) ;
NBRCS1_dB = zeros(2, 2, Nfz, Nhr) ;
NBRCS2_dB = zeros(2, 2, Nfz, Nhr) ;
NBRCS3_dB = zeros(2, 2, Nfz, Nhr) ;
NBRCS4_dB = zeros(2, 2, Nfz, Nhr) ;

for jj = 1 : Nhr
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    
    % TO-DO: Assign the correct index to each
    ParamsManager.index_Th( ANG_choice );
    ParamsManager.index_Ph( 1 );
    ParamsManager.index_VSM( SM_choice );
    ParamsManager.index_RMSH( 1 );
        
    % Initialize the directories depending on theta phi, VSM, and RMSH
    SimulationFolders.getInstance.initializeDynamicDirs();


    %% GET GLOBAL DIRECTORIES
    dir_out_diffuse_NBRCS_tuple = SimulationFolders.getInstance.out_diffuse_NBRCS_tuple;

    % NBRCS   
    % transmit port 1 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS_t1_dB') ;
    NBRCS_dB(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS1_t1_dB') ;
    NBRCS1_dB(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS2_t1_dB') ;
    NBRCS2_dB(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS3_t1_dB') ;
    NBRCS3_dB(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS4_t1_dB') ;
    NBRCS4_dB(1, :, :, jj) = xx ;
    % transmit port 2 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS_t2_dB') ;
    NBRCS_dB(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS1_t2_dB') ;
    NBRCS1_dB(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS2_t2_dB') ;
    NBRCS2_dB(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS3_t2_dB') ;
    NBRCS3_dB(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS4_t2_dB') ;
    NBRCS4_dB(2, :, :, jj) = xx ;
    
   
end


%% plotting as a function of fresnel zones
for jj = 1 : Nhr
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    
    %% GET GLOBAL PARAMETERS
    % Receiver Parameters
    hr_m = RecParams.getInstance.hr_m;
    
    
%     NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, jj)) ;
%     NBRCS_dB_RL = squeeze(NBRCS_dB(1, 2, :, jj)) ;
    
    NBRCS1_dB_RR = squeeze(NBRCS1_dB(1, 1, :, jj)) ;
    NBRCS1_dB_RL = squeeze(NBRCS1_dB(1, 2, :, jj)) ;
    
    NBRCS2_dB_RR = squeeze(NBRCS2_dB(1, 1, :, jj)) ;
    NBRCS2_dB_RL = squeeze(NBRCS2_dB(1, 2, :, jj)) ;
    
%     NBRCS3_dB_RR = squeeze(NBRCS3_dB(1, 1, :, jj)) ;
%     NBRCS3_dB_RL = squeeze(NBRCS3_dB(1, 2, :, jj)) ;
    
    NBRCS4_dB_RR = squeeze(NBRCS4_dB(1, 1, :, jj)) ;
    NBRCS4_dB_RL = squeeze(NBRCS4_dB(1, 2, :, jj)) ;
    
    if jj == 1
        figure
        subplot(1,2,1)
        title('X-Pol [RL]')
        hold
        axis([0 Nfz+2 -50 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS [dB]')
        xticks(1 : Nfz)
        grid
        text( Nfz+0.5, NBRCS2_dB_RL(end) - 3, strcat('h_r [m]'))
        
        
        subplot(1,2,2)
        title('Co-Pol [RR]')
        hold
        axis([0 Nfz+2 -50 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS [dB]')
        xticks(1 : Nfz)
        grid
    end
    subplot(1,2,1)
    plot(NBRCS1_dB_RL, ':sb', 'MarkerFaceColor', 'blue')
    plot(NBRCS2_dB_RL, ':dg', 'MarkerFaceColor', 'green')
    plot(NBRCS4_dB_RL, ':or', 'MarkerFaceColor', 'red')
    text(Nfz+0.05, NBRCS2_dB_RL(end), strcat(num2str(hr_m)))
    
    
    subplot(1,2,2)
    plot(NBRCS1_dB_RR, ':sb', 'MarkerFaceColor', 'blue')
    plot(NBRCS2_dB_RR, ':dg', 'MarkerFaceColor', 'green')
    plot(NBRCS4_dB_RR, ':or', 'MarkerFaceColor', 'red')
    legend('Single Bounce', 'Double Bounce', 'Triple Bounce')
    
    


end


%% GET GLOBAL DIRECTORIES
dir_analysis = SimulationFolders.getInstance.analysis;
if ~exist(dir_analysis, 'dir')
    mkdir(dir_analysis)
end


%% GET GLOBAL PARAMETERS
% Satellite Parameters
th0_deg = SatParams.getInstance.th0_list_deg( ParamsManager.index_Th );
polT = SatParams.getInstance.polT;
% Receiver Parameters
polR = RecParams.getInstance.polR;
% Ground Parameters
VSM_cm3cm3 = GndParams.getInstance.VSM_list_cm3cm3( ParamsManager.index_VSM );



%%
fname1 = strcat('TH_', num2str( th0_deg ), '-', polT, polR, '-VSM_', num2str( VSM_cm3cm3 )) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close

%%

for jj = 1 : Nhr
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    
    %% GET GLOBAL PARAMETERS
    hr_m = RecParams.getInstance.hr_m;
    
    
    NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, jj)) ;
    NBRCS_dB_RL = squeeze(NBRCS_dB(1, 2, :, jj)) ;
    
   
    if jj == 1
        figure
        hold
        axis([0 Nfz+2 -20 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS [dB]')
        xticks(1 : Nfz)
        grid
    
    end
    plot(NBRCS_dB_RL, ':sb', 'MarkerFaceColor', 'blue')
    plot(NBRCS_dB_RR, ':dg', 'MarkerFaceColor', 'green')
    text(Nfz+0.5, NBRCS_dB_RL(end), strcat(num2str(hr_m)))
    
    legend('X-POL [RL]', 'Co-POL [RR]', 'location', 'southwest')
    
    
end

text(Nfz+0.5, - 1, strcat('h_r [m]'))


%% GET GLOBAL DIRECTORIES
dir_analysis = SimulationFolders.getInstance.analysis;
if ~exist(dir_analysis, 'dir')
    mkdir(dir_analysis)
end


%% GET GLOBAL PARAMETERS
% Satellite Parameters
th0_deg = SatParams.getInstance.th0_list_deg( ParamsManager.index_Th );
% Ground Parameters
VSM_cm3cm3 = GndParams.getInstance.VSM_list_cm3cm3( ParamsManager.index_VSM );


%%
fname1 = strcat('TH2_', num2str( th0_deg ), '-', polT, polR, '-VSM_', num2str( VSM_cm3cm3 ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end