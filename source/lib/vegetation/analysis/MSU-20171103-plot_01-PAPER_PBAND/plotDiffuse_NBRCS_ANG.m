
% Mehmet Kurum
% 11/12/2017
% NBRCS vs. Observation Angle


function plotDiffuse_NBRCS_ANG


%% INPUT FILES
inputFile_sys_tag = 'sysInput-Paulownia-PAPER_PBAND-hr_';
inputFile_veg = 'vegHomInput-Paulownia.xml';


% TO-DO: Check here
%% Choices
FZ_choice = 1 ;
SM_choice = 2 ;
hr_choice = 4 ;

%% angle 
th0d = 10 : 10 : 70 ;

%% Reciever height
hhr = [20 50 100 500] ;


%% reading
% NBRCS
Nth0d =  length(th0d) ;
Nhr =  length(hhr) ;

NBRCS_dB = zeros(2, 2, Nth0d, Nhr) ; NBRCS1_dB = zeros(2, 2, Nth0d, Nhr) ;
NBRCS2_dB = zeros(2, 2, Nth0d, Nhr) ; NBRCS3_dB = zeros(2, 2, Nth0d, Nhr) ; NBRCS4_dB = zeros(2, 2, Nth0d, Nhr) ;


for jj = 1 : Nhr
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
       
    for ii = 1 :Nth0d
        
        % Set theta index
        ParamsManager.index_Th( ii );
        
        % TO-DO: Assign the correct index to each
        ParamsManager.index_Ph( 1 );
        ParamsManager.index_VSM( SM_choice );
        ParamsManager.index_RMSH( 1 );
        
        % Initialize the directories depending on theta phi, VSM, and RMSH
        SimulationFolders.getInstance.initializeDynamicDirs();
        
        
        %% GET GLOBAL DIRECTORIES
        dir_out_diffuse_NBRCS_tuple = SimulationFolders.getInstance.out_diffuse_NBRCS_tuple;
                
        % transmit port 1 / receiver ports 1&2
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS_t1_dB') ;
        NBRCS_dB(1, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS1_t1_dB') ;
        NBRCS1_dB(1, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS2_t1_dB') ;
        NBRCS2_dB(1, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS3_t1_dB') ;
        NBRCS3_dB(1, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS4_t1_dB') ;
        NBRCS4_dB(1, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        
        % transmit port 2 / receiver ports 1&2
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS_t2_dB') ;
        NBRCS_dB(2, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS1_t2_dB') ;
        NBRCS1_dB(2, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS2_t2_dB') ;
        NBRCS2_dB(2, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS3_t2_dB') ;
        NBRCS3_dB(2, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        xx = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS4_t2_dB') ;
        NBRCS4_dB(2, :, ii, jj) = squeeze(xx(:, FZ_choice)) ;
        
    end
end


%% reading reflectivity
% Intilization
% over vegetation
R_cohveg_dB = zeros(2, 2, Nth0d, Nhr) ;
% over bare soil
R_cohbare_dB = zeros(2, 2, Nth0d, Nhr) ;

for jj = 1 : Nhr
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    for ii = 1 : Nth0d
        
        % Set theta index
        ParamsManager.index_Th( ii );
        
        % TO-DO: Assign the correct index to each
        ParamsManager.index_Ph( 1 );
        ParamsManager.index_VSM( SM_choice );
        ParamsManager.index_RMSH( 1 );
        
        % Initialize the directories depending on theta phi, VSM, and RMSH
        SimulationFolders.getInstance.initializeDynamicDirs();
        
        
        [~, ~, ~, ~, ~, ~, ~, ~, ...
            R_coh1vegx, R_coh2vegx, ~, ~, ...
            R_coh1barex, R_coh2barex, ~, ~] = readSpecular ;
        
        R_cohveg_dB(1, :, ii, jj) = 10 * log10(R_coh1vegx(1 : 2)) ;
        R_cohveg_dB(2, :, ii, jj) = 10 * log10(R_coh2vegx(1 : 2)) ;
        
        R_cohbare_dB(1, :, ii, jj) = 10 * log10(R_coh1barex(1 : 2)) ;
        R_cohbare_dB(2, :, ii, jj) = 10 * log10(R_coh2barex(1 : 2)) ;
        
    end
    
end


%% GET GLOBAL DIRECTORIES
dir_analysis = SimulationFolders.getInstance.analysis;
if ~exist(dir_analysis, 'dir')
    mkdir(dir_analysis)
end


%% GET GLOBAL PARAMETERS
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Ground Parameters
VSM_cm3cm3 = GndParams.getInstance.VSM_list_cm3cm3( ParamsManager.index_VSM );


%% plotting as a function of angle

NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, hr_choice)) ;
NBRCS_dB_RL = squeeze(NBRCS_dB(1, 2, :, hr_choice)) ;

NBRCS1_dB_RR = squeeze(NBRCS1_dB(1, 1, :, hr_choice)) ;
NBRCS1_dB_RL = squeeze(NBRCS1_dB(1, 2, :, hr_choice)) ;
% NBRCS2_dB_RR = squeeze(NBRCS2_dB(1, 1, :, hr_choice)) ;
% NBRCS2_dB_RL = squeeze(NBRCS2_dB(1, 2, :, hr_choice)) ;
% NBRCS3_dB_RR = squeeze(NBRCS3_dB(1, 1, :, hr_choice)) ;
% NBRCS3_dB_RL = squeeze(NBRCS3_dB(1, 2, :, hr_choice)) ;
NBRCS4_dB_RR = squeeze(NBRCS4_dB(1, 1, :, hr_choice)) ;
NBRCS4_dB_RL = squeeze(NBRCS4_dB(1, 2, :, hr_choice)) ;

Rvc_RR = squeeze(R_cohveg_dB(1, 1, :, hr_choice)) ;
Rvc_RL = squeeze(R_cohveg_dB(1, 2, :, hr_choice)) ;

Rbc_RR = squeeze(R_cohbare_dB(1, 1, :, hr_choice)) ;
Rbc_RL = squeeze(R_cohbare_dB(1, 2, :, hr_choice)) ;


figure
subplot(1,2,1)
plot(th0d, NBRCS_dB_RL, '-or', 'MarkerFaceColor', 'red')
hold
plot(th0d, NBRCS1_dB_RL, '-ok', 'MarkerFaceColor', 'black')
plot(th0d, NBRCS4_dB_RL, '-oc', 'MarkerFaceColor', 'cyan')
plot(th0d, NBRCS_dB_RR, '-or', 'MarkerFaceColor', 'white')
plot(th0d, NBRCS1_dB_RR, '-ok', 'MarkerFaceColor', 'white')
plot(th0d, NBRCS4_dB_RR, '-oc', 'MarkerFaceColor', 'white')

% text(10, -10, strcat())
% text(10, -21, strcat(''))
% text(10, -30, strcat('Triple Bounce'))

legend('Double Bounce', 'Single Bounce', 'Triple Bounce', 'location', 'southeast')


text(-3, NBRCS_dB_RL(1), strcat('RL'))
text(-3, NBRCS_dB_RR(1), strcat('RR'))

text(-3, NBRCS1_dB_RL(1), strcat('RL'))
text(73, NBRCS1_dB_RR(end), strcat('RR'))

text(-3, NBRCS4_dB_RL(1), strcat('RL'))
text(73, NBRCS4_dB_RR(end), strcat('RR'))

axis([-5 85 -50 0])
xlabel('\theta_s [\circ]')
ylabel('NBRCS: \sigma^0_e [dB]')
xticks(th0d(1) : 10 : th0d(end))
grid
title('Diffuse [Vegetation]')

%         figure(2)
subplot(1,2,2)
plot(th0d, Rbc_RL, '-sb', 'MarkerFaceColor', 'blue')
hold
plot(th0d, Rvc_RL, '-dg', 'MarkerFaceColor', 'green')
plot(th0d, Rvc_RR, '-dg', 'MarkerFaceColor', 'white')
plot(th0d, Rbc_RR, '-sb', 'MarkerFaceColor', 'white')
ylabel('Reflectivity: \Gamma_s [dB]')

text(-3, (Rvc_RL(1) + Rbc_RL(1))/2, strcat('RL'))
text(-3, (Rvc_RR(1) + Rbc_RR(1))/2, strcat('RR'))

axis([-5 85 -50 0])
xlabel('\theta_s [\circ]')
xticks(th0d(1) : 10 : th0d(end))
grid
title('Specular')
legend('soil', 'vegetation', 'location', 'southeast')

%%
fname1 = strcat('FZ_', num2str(FZ_choice), '-', pol_Tx, pol_Rx, '-VSM_', num2str( VSM_cm3cm3 ), '-hr_', num2str( hhr(hr_choice) )) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end


%% read results
function [b_coh1veg, b_coh2veg, b0_coh1veg, b0_coh2veg, ...
    b_coh1bare, b_coh2bare, b0_coh1bare, b0_coh2bare, ...
    P_coh1veg, P_coh2veg, P0_coh1veg, P0_coh2veg, ...
    P_coh1bare, P_coh2bare, P0_coh1bare, P0_coh2bare] = readSpecular


%% GET GLOBAL DIRECTORIES
dir_out_specular_tuple = SimulationFolders.getInstance.out_specular_tuple;


%% GET GLOBAL PARAMETERS
pol_Tx = TxParams.getInstance.pol_Tx;
pol_Rx = RxParams.getInstance.pol_Rx;


%% READ SAVED OUTPUT

% 2 X 2
filename1 = strcat('Veg1', pol_Tx, pol_Rx) ;
filename2 = strcat('Veg2', pol_Tx, pol_Rx) ;
filename01 = strcat('Veg01', pol_Tx, pol_Rx) ;
filename02 = strcat('Veg02', pol_Tx, pol_Rx) ;
b_coh1veg = readComplexVar(dir_out_specular_tuple, filename1) ;
b_coh2veg = readComplexVar(dir_out_specular_tuple, filename2) ;
b0_coh1veg = readComplexVar(dir_out_specular_tuple, filename01) ;
b0_coh2veg = readComplexVar(dir_out_specular_tuple, filename02) ;

filename1 = strcat('Bare1', pol_Tx, pol_Rx) ;
filename2 = strcat('Bare2', pol_Tx, pol_Rx) ;
filename01 = strcat('Bare01', pol_Tx, pol_Rx) ;
filename02 = strcat('Bare02', pol_Tx, pol_Rx) ;
b_coh1bare = readComplexVar(dir_out_specular_tuple, filename1) ;
b_coh2bare = readComplexVar(dir_out_specular_tuple, filename2) ;
b0_coh1bare = readComplexVar(dir_out_specular_tuple, filename01) ;
b0_coh2bare = readComplexVar(dir_out_specular_tuple, filename02) ;

% 4 X 4
filename1 = strcat('P_Veg1', pol_Tx, pol_Rx) ;
filename2 = strcat('P_Veg2', pol_Tx, pol_Rx) ;
filename01 = strcat('P_Veg01', pol_Tx, pol_Rx) ;
filename02 = strcat('P_Veg02', pol_Tx, pol_Rx) ;
P_coh1veg = readVar(dir_out_specular_tuple, filename1) ;
P_coh2veg = readVar(dir_out_specular_tuple, filename2) ;
P0_coh1veg = readVar(dir_out_specular_tuple, filename01) ;
P0_coh2veg = readVar(dir_out_specular_tuple, filename02) ;

filename1 = strcat('P_Bare1', pol_Tx, pol_Rx) ;
filename2 = strcat('P_Bare2', pol_Tx, pol_Rx) ;
filename01 = strcat('P_Bare01', pol_Tx, pol_Rx) ;
filename02 = strcat('P_Bare02', pol_Tx, pol_Rx) ;
P_coh1bare = readVar(dir_out_specular_tuple, filename1) ;
P_coh2bare = readVar(dir_out_specular_tuple, filename2) ;
P0_coh1bare = readVar(dir_out_specular_tuple, filename01) ;
P0_coh2bare = readVar(dir_out_specular_tuple, filename02) ;

end