
% Mehmet Kurum
% 11/12/2017
% NBRCS vs. Fresnel Zones


function plotDiffuse_NBRCS_FZ2


%% INPUT FILES
inputFile_sys_tag = 'sysInput-Paulownia-PAPER_PBAND-hr_';
inputFile_sys_tag_B1 = 'sysInput-Paulownia-PAPER_PBAND_B1-hr_';
inputFile_sys_tag_G1 = 'sysInput-Paulownia-PAPER_PBAND_G1-hr_';
inputFile_sys_tag_B1_G1 = 'sysInput-Paulownia-PAPER_PBAND_B1_G1-hr_';
inputFile_veg = 'vegHomInput-Paulownia.xml';


% TO-DO: Check here
%% Choices
Nfz = 10;
SM_choice = 2 ;
ANG_choice = 4 ;



%% Reciever height
hhr = [20 50 100 500] ;


%% MIGHT BE NEEDED IN THE FUTURE
% %% Polarization
% 
% % t(trans)-r(rec)
% if (polT == 'X') && (polR == 'X')    
%     pols11 = 'XX' ;     pols12 = 'XY' ;
%     pols21 = 'YX' ;     pols22 = 'YY' ;
% elseif (polT == 'R') && (polR == 'R')
%     pols11 = 'RR' ;     pols12 = 'RL' ;
%     pols21 = 'LR' ;     pols22 = 'LL' ;
% elseif (polT == 'R') && (polR == 'X')
%     pols11 = 'RX' ;     pols12 = 'RY' ;
%     pols21 = 'LX' ;     pols22 = 'LY' ;
% elseif (polT == 'L') && (polR == 'L')
%     pols11 = 'LL' ;     pols12 = 'LR' ;
%     pols21 = 'RL' ;     pols22 = 'RR' ;
% elseif (polT == 'L') && (polR == 'X')
%     pols11 = 'LX' ;     pols12 = 'LY' ;
%     pols21 = 'RX' ;     pols22 = 'RY' ;
% end
% 
% pols = {pols11, pols12; pols21 pols22} ; %#ok<NASGU>


%% reading
% NBRCS
Nhr =  length(hhr) ;

NBRCS_dB = zeros(2, 2, Nfz, Nhr) ;
NBRCS1_dB = zeros(2, 2, Nfz, Nhr) ;
NBRCS2_dB = zeros(2, 2, Nfz, Nhr) ;
NBRCS3_dB = zeros(2, 2, Nfz, Nhr) ;
NBRCS4_dB = zeros(2, 2, Nfz, Nhr) ;

NBRCS_dB2 = zeros(2, 2, Nfz, Nhr) ;
NBRCS1_dB2 = zeros(2, 2, Nfz, Nhr) ;
NBRCS2_dB2 = zeros(2, 2, Nfz, Nhr) ;
NBRCS3_dB2 = zeros(2, 2, Nfz, Nhr) ;
NBRCS4_dB2 = zeros(2, 2, Nfz, Nhr) ;

NBRCS_dB3 = zeros(2, 2, Nfz, Nhr) ;
NBRCS1_dB3 = zeros(2, 2, Nfz, Nhr) ;
NBRCS2_dB3 = zeros(2, 2, Nfz, Nhr) ;
NBRCS3_dB3 = zeros(2, 2, Nfz, Nhr) ;
NBRCS4_dB3 = zeros(2, 2, Nfz, Nhr) ;

NBRCS_dB4 = zeros(2, 2, Nfz, Nhr) ;
NBRCS1_dB4 = zeros(2, 2, Nfz, Nhr) ;
NBRCS2_dB4 = zeros(2, 2, Nfz, Nhr) ;
NBRCS3_dB4 = zeros(2, 2, Nfz, Nhr) ;
NBRCS4_dB4 = zeros(2, 2, Nfz, Nhr) ;

for jj = 1 : Nhr
    
    
    %% ORIGINAL 
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


    % NBRCS - Original 
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
    
    
    
    %% B = 1 
    inputFile_sys_B1 = strcat( inputFile_sys_tag_B1, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys_B1, inputFile_veg );
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
    dir_out_diffuse_NBRCS_tuple_B1 = SimulationFolders.getInstance.out_diffuse_NBRCS_tuple;

    
    % NBRCS - B1    
    % transmit port 1 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1, 'NBRCS_t1_dB') ;
    NBRCS_dB2(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1, 'NBRCS1_t1_dB') ;
    NBRCS1_dB2(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1, 'NBRCS2_t1_dB') ;
    NBRCS2_dB2(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1, 'NBRCS3_t1_dB') ;
    NBRCS3_dB2(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1, 'NBRCS4_t1_dB') ;
    NBRCS4_dB2(1, :, :, jj) = xx ;
    % transmit port 2 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1, 'NBRCS_t2_dB') ;
    NBRCS_dB2(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1, 'NBRCS1_t2_dB') ;
    NBRCS1_dB2(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1, 'NBRCS2_t2_dB') ;
    NBRCS2_dB2(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1, 'NBRCS3_t2_dB') ;
    NBRCS3_dB2(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1, 'NBRCS4_t2_dB') ;
    NBRCS4_dB2(2, :, :, jj) = xx ;
    
    
    
    %% G = 1 
    inputFile_sys_G1 = strcat( inputFile_sys_tag_G1, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys_G1, inputFile_veg );
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
    dir_out_diffuse_NBRCS_tuple_G1 = SimulationFolders.getInstance.out_diffuse_NBRCS_tuple;

    
    % NBRCS - G1    
    % transmit port 1 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_NBRCS_tuple_G1, 'NBRCS_t1_dB') ;
    NBRCS_dB3(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_G1, 'NBRCS1_t1_dB') ;
    NBRCS1_dB3(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_G1, 'NBRCS2_t1_dB') ;
    NBRCS2_dB3(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_G1, 'NBRCS3_t1_dB') ;
    NBRCS3_dB3(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_G1, 'NBRCS4_t1_dB') ;
    NBRCS4_dB3(1, :, :, jj) = xx ;
    % transmit port 2 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_NBRCS_tuple_G1, 'NBRCS_t2_dB') ;
    NBRCS_dB3(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_G1, 'NBRCS1_t2_dB') ;
    NBRCS1_dB3(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_G1, 'NBRCS2_t2_dB') ;
    NBRCS2_dB3(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_G1, 'NBRCS3_t2_dB') ;
    NBRCS3_dB3(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_G1, 'NBRCS4_t2_dB') ;
    NBRCS4_dB3(2, :, :, jj) = xx ;
    
    
    
    %% B = 1, G=1 
    inputFile_sys_B1_G1 = strcat( inputFile_sys_tag_B1_G1, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys_B1_G1, inputFile_veg );
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
    dir_out_diffuse_NBRCS_tuple_B1_G1 = SimulationFolders.getInstance.out_diffuse_NBRCS_tuple;

    
    % NBRCS - B1, G1    
    % transmit port 1 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1_G1, 'NBRCS_t1_dB') ;
    NBRCS_dB4(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1_G1, 'NBRCS1_t1_dB') ;
    NBRCS1_dB4(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1_G1, 'NBRCS2_t1_dB') ;
    NBRCS2_dB4(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1_G1, 'NBRCS3_t1_dB') ;
    NBRCS3_dB4(1, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1_G1, 'NBRCS4_t1_dB') ;
    NBRCS4_dB4(1, :, :, jj) = xx ;
    % transmit port 2 / receiver ports 1&2
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1_G1, 'NBRCS_t2_dB') ;
    NBRCS_dB4(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1_G1, 'NBRCS1_t2_dB') ;
    NBRCS1_dB4(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1_G1, 'NBRCS2_t2_dB') ;
    NBRCS2_dB4(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1_G1, 'NBRCS3_t2_dB') ;
    NBRCS3_dB4(2, :, :, jj) = xx ;
    xx = readVar(dir_out_diffuse_NBRCS_tuple_B1_G1, 'NBRCS4_t2_dB') ;
    NBRCS4_dB4(2, :, :, jj) = xx ;
    
   
end

%% plotting as a function of fresnel zones

% % for jj = 1 : Nhr
% %     
% %     hr = hhr(jj) ;
% %     
% % %     NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, jj)) ;
% % %     NBRCS_dB_RL = squeeze(NBRCS_dB(1, 2, :, jj)) ;
% %     
% %     NBRCS1_dB_RR = squeeze(NBRCS1_dB(1, 1, :, jj)) ;
% %     NBRCS1_dB_RL = squeeze(NBRCS1_dB(1, 2, :, jj)) ;
% %     
% %     NBRCS2_dB_RR = squeeze(NBRCS2_dB(1, 1, :, jj)) ;
% %     NBRCS2_dB_RL = squeeze(NBRCS2_dB(1, 2, :, jj)) ;
% %     
% % %     NBRCS3_dB_RR = squeeze(NBRCS3_dB(1, 1, :, jj)) ;
% % %     NBRCS3_dB_RL = squeeze(NBRCS3_dB(1, 2, :, jj)) ;
% %     
% %     NBRCS4_dB_RR = squeeze(NBRCS4_dB(1, 1, :, jj)) ;
% %     NBRCS4_dB_RL = squeeze(NBRCS4_dB(1, 2, :, jj)) ;
% %     
% %     if jj == 1
% %         figure
% %         subplot(1,2,1)
% %         title('X-Pol [RL]')
% %         hold
% %         axis([0 Nfz+2 -50 0])
% %         xlabel('Fresnel Zones')
% %         ylabel('NBRCS [dB]')
% %         xticks(1 : Nfz)
% %         grid
% %         text(Nfz+0.5, NBRCS2_dB_RL(end) - 3, strcat('h_r [m]'))
% %         
% %         
% %         subplot(1,2,2)
% %         title('Co-Pol [RR]')
% %         hold
% %         axis([0 Nfz+2 -50 0])
% %         xlabel('Fresnel Zones')
% %         ylabel('NBRCS [dB]')
% %         xticks(1 : Nfz)
% %         grid
% %     end
% %     subplot(1,2,1)
% %     plot(NBRCS1_dB_RL, ':sb', 'MarkerFaceColor', 'blue')
% %     plot(NBRCS2_dB_RL, ':dg', 'MarkerFaceColor', 'green')
% %     plot(NBRCS4_dB_RL, ':or', 'MarkerFaceColor', 'red')
% %     text(Nfz+0.05, NBRCS2_dB_RL(end), strcat(num2str(hr)))
% %     
% %     
% %     subplot(1,2,2)
% %     plot(NBRCS1_dB_RR, ':sb', 'MarkerFaceColor', 'blue')
% %     plot(NBRCS2_dB_RR, ':dg', 'MarkerFaceColor', 'green')
% %     plot(NBRCS4_dB_RR, ':or', 'MarkerFaceColor', 'red')
% %     legend('Single Bounce', 'Double Bounce', 'Triple Bounce')
% %     
% %     
% % 
% % 
% % end

% % %%
% % fname1 = strcat('TH_', num2str(th0d), '_', polT, polR, '_SM_', num2str(tp(SM_choice))) ;
% % 
% % FolderPath_write = strcat(pwd, '\Figures') ;
% % if ~exist(FolderPath_write, 'dir')
% %     mkdir(FolderPath_write)
% % end
% % 
% % saveas(gcf, strcat(FolderPath_write, '\', fname1), 'jpg')
% % close

%%


%% ORIGINAL TO B1 COMPARISON PLOT
for jj = 1 : Nhr
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    
    %% GET GLOBAL PARAMETERS
    % Receiver Parameters
    hr_m = RecParams.getInstance.hr_m;
    
    %%
    NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, jj)) ;
    NBRCS_dB_RL = squeeze(NBRCS_dB(1, 2, :, jj)) ;
    %%
    NBRCS_dB_RR2 = squeeze(NBRCS_dB2(1, 1, :, jj)) ;
    NBRCS_dB_RL2 = squeeze(NBRCS_dB2(1, 2, :, jj)) ;
    
    %%
   
    if jj == 1
        figure
        subplot(1,2,1)
        hold
        axis([0 Nfz+2 -20 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS: \sigma^0_e [dB]')
        xticks(1 : Nfz)
        text(Nfz+0.5, - 1, strcat('h_r [m]'))
        grid
        
        subplot(1,2,2)
        hold
        axis([0 Nfz+2 -20 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS: \sigma^0_e [dB]')
        xticks(1 : Nfz)
        text(Nfz+0.5, - 8.5, strcat('h_r [m]'))
        grid
                

    end
    subplot(1,2,1)
    
    plot(NBRCS_dB_RL, '-sb', 'MarkerSize', 4, 'MarkerFaceColor', 'blue')
    plot(NBRCS_dB_RL2, ':sb', 'MarkerSize', 4, 'MarkerFaceColor', 'white')
    text(Nfz+0.5, NBRCS_dB_RL(end), strcat(num2str(hr_m)))
    title('X-Pol [RL]')
    
    subplot(1,2,2)
    
    plot(NBRCS_dB_RR, '-or', 'MarkerSize', 4, 'MarkerFaceColor', 'red')
    plot(NBRCS_dB_RR2, ':or', 'MarkerSize', 4, 'MarkerFaceColor', 'white')
    text(Nfz+0.5, NBRCS_dB_RR(end), strcat(num2str(hr_m)))
    title('Co-Pol [RR]')    
    
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


subplot(1,2,1)
        
legend('Actual Antenna', 'Actual Antenna - B1', 'location', 'southwest')

subplot(1,2,2)

legend('Actual Antenna', 'Actual Antenna - B1', 'location', 'northwest')


%%
fname1 = strcat('XTH2_', num2str(th0_deg), '-', polT, polR, '-VSM_', num2str( VSM_cm3cm3 ) );
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close



%% ORIGINAL TO G1 COMPARISON PLOT
for jj = 1 : Nhr
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    
    %% GET GLOBAL PARAMETERS
    % Receiver Parameters
    hr_m = RecParams.getInstance.hr_m;
    
    %%
    NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, jj)) ;
    NBRCS_dB_RL = squeeze(NBRCS_dB(1, 2, :, jj)) ;
    %%
    NBRCS_dB_RR3 = squeeze(NBRCS_dB3(1, 1, :, jj)) ;
    NBRCS_dB_RL3 = squeeze(NBRCS_dB3(1, 2, :, jj)) ;
    
    %%
   
    if jj == 1
        figure
        subplot(1,2,1)
        hold
        axis([0 Nfz+2 -20 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS: \sigma^0_e [dB]')
        xticks(1 : Nfz)
        text(Nfz+0.5, - 1, strcat('h_r [m]'))
        grid
        
        subplot(1,2,2)
        hold
        axis([0 Nfz+2 -20 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS: \sigma^0_e [dB]')
        xticks(1 : Nfz)
        text(Nfz+0.5, - 8.5, strcat('h_r [m]'))
        grid
                

    end
    subplot(1,2,1)
    
    plot(NBRCS_dB_RL, '-sb', 'MarkerSize', 4, 'MarkerFaceColor', 'blue')
    plot(NBRCS_dB_RL3, ':sb', 'MarkerSize', 4, 'MarkerFaceColor', 'white')
    text(Nfz+0.5, NBRCS_dB_RL(end), strcat(num2str(hr_m)))
    title('X-Pol [RL]')
    
    subplot(1,2,2)
    
    plot(NBRCS_dB_RR, '-or', 'MarkerSize', 4, 'MarkerFaceColor', 'red')
    plot(NBRCS_dB_RR3, ':or', 'MarkerSize', 4, 'MarkerFaceColor', 'white')
    text(Nfz+0.5, NBRCS_dB_RR(end), strcat(num2str(hr_m)))
    title('Co-Pol [RR]')
    
end

subplot(1,2,1)
        
legend('Actual Antenna', 'Ideal Antenna', 'location', 'southwest')

subplot(1,2,2)

legend('Actual Antenna', 'Ideal Antenna', 'location', 'northwest')


%%
fname1 = strcat('XTH3_', num2str(th0_deg), '-', polT, polR, '-VSM_', num2str( VSM_cm3cm3 ) );
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


for jj = 1 : Nhr
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    
    %% GET GLOBAL PARAMETERS
    % Receiver Parameters
    hr_m = RecParams.getInstance.hr_m;
    
    %%
    NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, jj)) ;
    NBRCS_dB_RL = squeeze(NBRCS_dB(1, 2, :, jj)) ;
    %%
    NBRCS_dB_RR4 = squeeze(NBRCS_dB4(1, 1, :, jj)) ;
    NBRCS_dB_RL4 = squeeze(NBRCS_dB4(1, 2, :, jj)) ;
    
    %%
   
    if jj == 1
        figure
        subplot(1,2,1)
        hold
        axis([0 Nfz+2 -20 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS: \sigma^0_e [dB]')
        xticks(1 : Nfz)
        text(Nfz+0.5, - 1, strcat('h_r [m]'))
        grid
        
        subplot(1,2,2)
        hold
        axis([0 Nfz+2 -20 0])
        xlabel('Fresnel Zones')
        ylabel('NBRCS: \sigma^0_e [dB]')
        xticks(1 : Nfz)
        text(Nfz+0.5, - 8.5, strcat('h_r [m]'))
        grid
                

    end
    subplot(1,2,1)
    
    plot(NBRCS_dB_RL, '-sb', 'MarkerSize', 4, 'MarkerFaceColor', 'blue')
    plot(NBRCS_dB_RL4, ':sb', 'MarkerSize', 4, 'MarkerFaceColor', 'white')
    text(Nfz+0.5, NBRCS_dB_RL(end), strcat(num2str(hr_m)))
    title('X-Pol [RL]')
    
    subplot(1,2,2)
    
    plot(NBRCS_dB_RR, '-or', 'MarkerSize', 4, 'MarkerFaceColor', 'red')
    plot(NBRCS_dB_RR4, ':or', 'MarkerSize', 4, 'MarkerFaceColor', 'white')
    text(Nfz+0.5, NBRCS_dB_RR(end), strcat(num2str(hr_m)))
    title('Co-Pol [RR]')    
    
end

subplot(1,2,1)
        
legend('Actual Antenna', 'Ideal Antenna - B1', 'location', 'southwest')

subplot(1,2,2)

legend('Actual Antenna', 'Ideal Antenna - B1', 'location', 'northwest')


%%
fname1 = strcat('XTH4_', num2str(th0_deg), '-', polT, polR, '-VSM_', num2str( VSM_cm3cm3 ) );
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end