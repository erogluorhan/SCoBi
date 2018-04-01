
% Mehmet Kurum
% 11/12/2017
% Received Power vs. Fresnel Zones

function plotDiffuse_ReceivedPower_vsFZ


%% INPUT FILES
inputFile_sys = 'sysInput-Paulownia-PAPER_PBAND-hr_20.xml';
inputFile_veg = 'vegHomInput-Paulownia.xml';


%% GET INPUT FOR INITIAL START
getInput( inputFile_sys, inputFile_veg );

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();

% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();


%% GET GLOBAL PARAMETERS
% Satellite Parameters
th0_list_deg = SatParams.getInstance.th0_list_deg;
PH0_list_deg = SatParams.getInstance.PH0_list_deg;
polT = SatParams.getInstance.polT;
% Receiver Parameters
polR = RecParams.getInstance.polR;
% Ground Parameters
VSM_list_cm3cm3 = GndParams.getInstance.VSM_list_cm3cm3;
RMSH_list_cm = GndParams.getInstance.RMSH_list_cm;

num_Th = length( th0_list_deg );
num_Ph = length( PH0_list_deg );
num_VSM = length( VSM_list_cm3cm3 );
num_RMSH = length( RMSH_list_cm );


%% Choices
% TO-DO: Check here (This can be only valid for SNAPSHOT)
TH_choice = 4;
PH_choice = 1;
VSM_choice = 1;
RMSH_choice = 1;


%% Polarization
% t(trans)-r(rec)
if (polT == 'X') && (polR == 'X')    
    pols11 = 'XX' ;     pols12 = 'XY' ;
    pols21 = 'YX' ;     pols22 = 'YY' ;
elseif (polT == 'R') && (polR == 'R')
    pols11 = 'RR' ;     pols12 = 'RL' ;
    pols21 = 'LR' ;     pols22 = 'LL' ;
elseif (polT == 'R') && (polR == 'X')
    pols11 = 'RX' ;     pols12 = 'RY' ;
    pols21 = 'LX' ;     pols22 = 'LY' ;
elseif (polT == 'L') && (polR == 'L')
    pols11 = 'LL' ;     pols12 = 'LR' ;
    pols21 = 'RL' ;     pols22 = 'RR' ;
elseif (polT == 'L') && (polR == 'X')
    pols11 = 'LX' ;     pols12 = 'LY' ;
    pols21 = 'RX' ;     pols22 = 'RY' ;
end

pols = {pols11, pols12; pols21 pols22} ;


% Set parameters' indices
ParamsManager.index_Th( TH_choice );
ParamsManager.index_Ph( PH_choice );
ParamsManager.index_RMSH( RMSH_choice );

for vv = 1 : num_VSM
    
    ParamsManager.index_VSM( vv );

    % Initialize the directories depending on theta, phi, VSM, and RMSH
    SimulationFolders.getInstance.initializeDynamicDirs();


    %% GET GLOBAL DIRECTORIES
    dir_out_diffuse_P1_tuple = SimulationFolders.getInstance.out_diffuse_P1_tuple;
    % dir_out_diffuse_P2_tuple = SimulationFodlers.getInstance.out_diffuse_P2_tuple;
    % dir_out_diffuse_P3_tuple = SimulationFodlers.getInstance.out_diffuse_P3_tuple;
    % dir_out_diffuse_P4_tuple = SimulationFodlers.getInstance.out_diffuse_P4_tuple;


    %% reading
    % P1
    PP1_inc1_t1_dB(:, vv, :) = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc1_t1_dB');
    PP1_inc2_t1_dB(:, vv, :) = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc2_t1_dB');
    PP1_inc3_t1_dB(:, vv, :) = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc3_t1_dB');
    PP1_inc4_t1_dB(:, vv, :) = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc4_t1_dB');
    PP1_inc_t1_dB(:, vv, :) = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t1_dB');

    PP1_inc1_t2_dB(:, vv, :) = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc1_t2_dB');
    PP1_inc2_t2_dB(:, vv, :) = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc2_t2_dB');
    PP1_inc3_t2_dB(:, vv, :) = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc3_t2_dB');
    PP1_inc4_t2_dB(:, vv, :) = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc4_t2_dB');
    PP1_inc_t2_dB(:, vv, :) = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t2_dB');

end


%% plotting as a function of fresnel zones
P11_inc = squeeze(PP1_inc_t1_dB(1, :, :)) ;
P11_inc1 = squeeze(PP1_inc1_t1_dB(1, :, :)) ;
P11_inc2 = squeeze(PP1_inc2_t1_dB(1, :, :)) ;
P11_inc3 = squeeze(PP1_inc3_t1_dB(1, :, :)) ;
P11_inc4 = squeeze(PP1_inc4_t1_dB(1, :, :)) ;

plotFZ(P11_inc, P11_inc1, P11_inc2, P11_inc3, P11_inc4, pols{1, 1})

P12_inc = squeeze(PP1_inc_t1_dB(2, :, :)) ;
P12_inc1 = squeeze(PP1_inc1_t1_dB(2, :, :)) ;
P12_inc2 = squeeze(PP1_inc2_t1_dB(2, :, :)) ;
P12_inc3 = squeeze(PP1_inc3_t1_dB(2, :, :)) ;
P12_inc4 = squeeze(PP1_inc4_t1_dB(2, :, :)) ;

plotFZ(P12_inc, P12_inc1, P12_inc2, P12_inc3, P12_inc4, pols{1, 2})

P21_inc = squeeze(PP1_inc_t2_dB(1, :, :)) ;
P21_inc1 = squeeze(PP1_inc1_t2_dB(1, :, :)) ;
P21_inc2 = squeeze(PP1_inc2_t2_dB(1, :, :)) ;
P21_inc3 = squeeze(PP1_inc3_t2_dB(1, :, :)) ;
P21_inc4 = squeeze(PP1_inc4_t2_dB(1, :, :)) ;

plotFZ(P21_inc, P21_inc1, P21_inc2, P21_inc3, P21_inc4, pols{2, 1})

P22_inc = squeeze(PP1_inc_t2_dB(2, :, :)) ;
P22_inc1 = squeeze(PP1_inc1_t2_dB(2, :, :)) ;
P22_inc2 = squeeze(PP1_inc2_t2_dB(2, :, :)) ;
P22_inc3 = squeeze(PP1_inc3_t2_dB(2, :, :)) ;
P22_inc4 = squeeze(PP1_inc4_t2_dB(2, :, :)) ;

plotFZ(P22_inc, P22_inc1, P22_inc2, P22_inc3, P22_inc4, pols{2, 2})


% % % P2
% % PP2_inc1_t1_dB = readVar(dir_out_diffuse_P2_tuple, 'PP2_inc1_t1_dB') ;
% % PP2_inc2_t1_dB = readVar(dir_out_diffuse_P2_tuple, 'PP2_inc2_t1_dB') ;
% % PP2_inc3_t1_dB = readVar(dir_out_diffuse_P2_tuple, 'PP2_inc3_t1_dB') ;
% % PP2_inc4_t1_dB = readVar(dir_out_diffuse_P2_tuple, 'PP2_inc4_t1_dB') ;
% % PP2_inc_t1_dB = readVar(dir_out_diffuse_P2_tuple, 'PP2_inc_t1_dB') ;
% % 
% % PP2_inc1_t2_dB = readVar(dir_out_diffuse_P2_tuple, 'PP2_inc1_t2_dB') ;
% % PP2_inc2_t2_dB = readVar(dir_out_diffuse_P2_tuple, 'PP2_inc2_t2_dB') ;
% % PP2_inc3_t2_dB = readVar(dir_out_diffuse_P2_tuple, 'PP2_inc3_t2_dB') ;
% % PP2_inc4_t2_dB = readVar(dir_out_diffuse_P2_tuple, 'PP2_inc4_t2_dB') ;
% % PP2_inc_t2_dB = readVar(dir_out_diffuse_P2_tuple, 'PP2_inc_t2_dB') ;
% % 
% % % P3
% % PP3_inc1_t1_dB = readVar(dir_out_diffuse_P3_tuple, 'PP3_inc1_t1_dB') ;
% % PP3_inc2_t1_dB = readVar(dir_out_diffuse_P3_tuple, 'PP3_inc2_t1_dB') ;
% % PP3_inc3_t1_dB = readVar(dir_out_diffuse_P3_tuple, 'PP3_inc3_t1_dB') ;
% % PP3_inc4_t1_dB = readVar(dir_out_diffuse_P3_tuple, 'PP3_inc4_t1_dB') ;
% % PP3_inc_t1_dB = readVar(dir_out_diffuse_P3_tuple, 'PP3_inc_t1_dB') ;
% % 
% % PP3_inc1_t2_dB = readVar(dir_out_diffuse_P3_tuple, 'PP3_inc1_t2_dB') ;
% % PP3_inc2_t2_dB = readVar(dir_out_diffuse_P3_tuple, 'PP3_inc2_t2_dB') ;
% % PP3_inc3_t2_dB = readVar(dir_out_diffuse_P3_tuple, 'PP3_inc3_t2_dB') ;
% % PP3_inc4_t2_dB = readVar(dir_out_diffuse_P3_tuple, 'PP3_inc4_t2_dB') ;
% % PP3_inc_t2_dB = readVar(dir_out_diffuse_P3_tuple, 'PP3_inc_t2_dB') ;
% % 
% % % P4
% % PP4_inc1_t1_dB = readVar(dir_out_diffuse_P4_tuple, 'PP4_inc1_t1_dB') ;
% % PP4_inc2_t1_dB = readVar(dir_out_diffuse_P4_tuple, 'PP4_inc2_t1_dB') ;
% % PP4_inc3_t1_dB = readVar(dir_out_diffuse_P4_tuple, 'PP4_inc3_t1_dB') ;
% % PP4_inc4_t1_dB = readVar(dir_out_diffuse_P4_tuple, 'PP4_inc4_t1_dB') ;
% % PP4_inc_t1_dB = readVar(dir_out_diffuse_P4_tuple, 'PP4_inc_t1_dB') ;
% % 
% % PP4_inc1_t2_dB = readVar(dir_out_diffuse_P4_tuple, 'PP4_inc1_t2_dB') ;
% % PP4_inc2_t2_dB = readVar(dir_out_diffuse_P4_tuple, 'PP4_inc2_t2_dB') ;
% % PP4_inc3_t2_dB = readVar(dir_out_diffuse_P4_tuple, 'PP4_inc3_t2_dB') ;
% % PP4_inc4_t2_dB = readVar(dir_out_diffuse_P4_tuple, 'PP4_inc4_t2_dB') ;
% % PP4_inc_t2_dB = readVar(dir_out_diffuse_P4_tuple, 'PP4_inc_t2_dB') ;





end

function plotFZ(P_inc, P_inc1, P_inc2, P_inc3, P_inc4, pols)


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nr = SimParams.getInstance.Nr;
Nfz = SimParams.getInstance.Nfz;
% Ground Parameters
RMSH_cm = GndParams.getInstance.RMSH_list_cm( ParamsManager.index_RMSH );
% Satellite Parameters
th0_deg = SatParams.getInstance.th0_list_deg( ParamsManager.index_Th );
PH0_deg = SatParams.getInstance.PH0_list_deg( ParamsManager.index_Ph );
% Receiver Parameters
hr_m = RecParams.getInstance.hr_m;
hpbw_deg = RecParams.getInstance.hpbw_deg;
SLL_dB = RecParams.getInstance.SLL_dB;
XPL_dB = RecParams.getInstance.XPL_dB;


%% GET GLOBAL DIRECTORIES
dir_fig_diffuse_P1_vsFZ = SimulationFolders.getInstance.fig_diffuse_P1_vsFZ;

fig_out_pathname = strcat( dir_fig_diffuse_P1_vsFZ, '\TH_', num2str( th0_deg ) );
if ~exist(fig_out_pathname, 'dir')
    mkdir(fig_out_pathname)
end


FolderNameAnt = strcat('thsd_', num2str(hpbw_deg),...
    '-SLL_', num2str(SLL_dB), '-XPL_', num2str(XPL_dB)) ;

ylabely = 'Received Power [dB]' ;
P_range = 40 ;

ymin1 = floor((min(min(P_inc1))) / 10) * 10 ;
% ymax1 = ceil(max(max(max(P_inc1))) / 10) * 10 ;
ymin2 = floor((min(min(P_inc2))) / 10) * 10 ;
% ymax2 = ceil(max(max(max(B2_inc1))) / 10) * 10 ;
ymin3 = floor((min(min(P_inc3))) / 10) * 10 ;
% ymax3 = ceil(max(max(max(B3_inc1))) / 10) * 10 ;
ymin4 = floor((min(min(P_inc4))) / 10) * 10 ;
% ymax4 = ceil(max(max(max(B4_inc1))) / 10) * 10 ;
ymin = floor((min(min(P_inc))) / 10) * 10 ;
% ymax = ceil(max(max(max(B4_inc))) / 10) * 10 ;

%%
figure
ymini = ymin1 ; ymaxi = ymin1 + P_range ;

subplot(2,2,1)
plot(squeeze(P_inc1(1, :)), ':or')
hold
plot(squeeze(P_inc1(2, :)), ':sb')
plot(squeeze(P_inc1(3, :)), ':db')

title(strcat('(o^+, i^-)'))
text(1, ymaxi - 4, strcat('\theta_0=', num2str(th0_deg), '\circ'))
text(1, ymaxi - 10, strcat('N_r=', num2str(Nr)))
% text(1, ymaxi - 15, strcat(IQUV))
% text(5, ymaxi - 5, strcat('Medium'))
text(5, ymaxi - 4, strcat('h_r=', num2str(hr_m), 'm'))
text(5, ymaxi - 10, strcat(pols, '-pol'))
axis([0 Nfz + 1 ymini ymaxi])
grid
ylabel(ylabely)

%%
ymini = ymin2 ; ymaxi = ymin2 + P_range ;
subplot(2,2,2)
plot(squeeze(P_inc2(1, :)), ':or')
hold
plot(squeeze(P_inc2(2, :)), ':sb')
plot(squeeze(P_inc2(3, :)), ':db')
title(strcat('(o_I^-, i^-)'))
text(1, ymaxi - 5, strrep( strrep( FolderNameAnt, '_', '=' ), '-', ', ' ) )
axis([0 Nfz+1 ymini ymaxi])
grid

%%
ymini = ymin3 ; ymaxi = ymin3 + P_range ;
subplot(2,2,3)
plot(squeeze(P_inc3(1, :)), ':or')
hold
plot(squeeze(P_inc3(2, :)), ':sb')
plot(squeeze(P_inc3(3, :)), ':db')
title(strcat('(o^+, i_I^+)'))
axis([0 Nfz + 1 ymini ymaxi])
grid
xlabel('Fresnel Zones')
ylabel(ylabely)

%%
ymini = ymin4 ; ymaxi = ymin4 + P_range ;
subplot(2,2,4)
plot(squeeze(P_inc4(1, :)), ':or')
hold
plot(squeeze(P_inc4(2, :)), ':sb')
plot(squeeze(P_inc4(3, :)), ':db')
title(strcat('(o_I^-, i_I^+)'))
axis([0 Nfz + 1 ymini ymaxi])
grid
xlabel('Fresnel Zones')

fname = strcat('FZ-ReceivedPower_IND_dB', '-TH_', num2str( th0_deg ), '-PH_', num2str( PH0_deg ), '-RMSH_', num2str( RMSH_cm ) );
fname = strrep( fname, '.', 'dot' );
saveas(gcf, strcat(fig_out_pathname, '\', fname), 'tiff')
close


%%
figure
ymax = ymin + P_range ;
plot(squeeze(P_inc(1, :)), ':or')
hold
plot(squeeze(P_inc(2, :)), ':sb')
plot(squeeze(P_inc(3, :)), ':db')
title(strcat('Medium - All Mechanism *** ', ...
    '\theta_0 = ', num2str(th0_deg), '\circ', ...
    ' *** N_r = ', num2str(Nr), ' ***', pols, '-pol'))
axis([0 Nfz + 1 ymin ymax])
grid
xlabel('Fresnel Zones')
ylabel('Received Power [dB]')
text(1, ymax - 3, strrep( strrep( FolderNameAnt, '_', '=' ), '-', ',   ' ) )
text(1, ymax - 6, strcat('h_r=', num2str(hr_m), 'm'))

fname = strcat('FZ-ReceivedPower_ALL_dB', '-TH_', num2str( th0_deg ), '-PH_', num2str( PH0_deg ), '-RMSH_', num2str( RMSH_cm ) );
fname = strrep( fname, '.', 'dot' );
saveas(gcf, strcat(fig_out_pathname, '\',fname), 'tiff')
close

end