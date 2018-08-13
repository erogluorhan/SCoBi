
% Mehmet Kurum
% 11/12/2017
% NBRCS vs. Fresnel Zones


function plotDiffuse_NBRCS_vsFZ


%% INPUT FILES
inputFile_sys = 'sysInput-Paulownia-PAPER_PBAND-hr_20.xml';
inputFile_veg = 'vegHomInput-Paulownia.xml';


%% GET INPUT FOR INITIAL START
getInput( inputFile_sys, inputFile_veg );

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();

% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();


%% GET GLOBAL PARAMETERS
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;

num_VSM = length( VSM_list_cm3cm3 );


%% Choices
% TO-DO: Check here (This can be only valid for SNAPSHOT)
TH_choice = 4;
PH_choice = 1;
RMSH_choice = 1;


%% Polarization
% t(trans)-r(rec)
if ( pol_Tx == Constants.polarizations{Constants.id_pol_X} ) ...
 && (pol_Rx == Constants.polarizations{Constants.id_pol_X}) 

    pols11 = 'XX' ;     pols12 = 'XY' ;
    pols21 = 'YX' ;     pols22 = 'YY' ;
    
elseif (pol_Tx == Constants.polarizations{Constants.id_pol_R}) ...
    && (pol_Rx == Constants.polarizations{Constants.id_pol_R})

    pols11 = 'RR' ;     pols12 = 'RL' ;
    pols21 = 'LR' ;     pols22 = 'LL' ;
    
elseif (pol_Tx == Constants.polarizations{Constants.id_pol_R}) ...
    && (pol_Rx == Constants.polarizations{Constants.id_pol_X})

    pols11 = 'RX' ;     pols12 = 'RY' ;
    pols21 = 'LX' ;     pols22 = 'LY' ;
    
elseif (pol_Tx == Constants.polarizations{Constants.id_pol_L}) ...
    && (pol_Rx == Constants.polarizations{Constants.id_pol_L})
    
    pols11 = 'LL' ;     pols12 = 'LR' ;
    pols21 = 'RL' ;     pols22 = 'RR' ;
    
elseif (pol_Tx == Constants.polarizations{Constants.id_pol_L}) ...
    && (pol_Rx == Constants.polarizations{Constants.id_pol_X})

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
    dir_out_diffuse_NBRCS_tuple = SimulationFolders.getInstance.out_diffuse_NBRCS_tuple;

    
    %% reading
    % NBRCS
    NBRCS1_t1_dB(:, vv, :) = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS1_t1_dB') ;
    NBRCS2_t1_dB(:, vv, :) = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS2_t1_dB') ;
    NBRCS3_t1_dB(:, vv, :) = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS3_t1_dB') ;
    NBRCS4_t1_dB(:, vv, :) = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS4_t1_dB') ;
    NBRCS_t1_dB(:, vv, :) = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS_t1_dB') ;

    NBRCS1_t2_dB(:, vv, :) = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS1_t2_dB') ;
    NBRCS2_t2_dB(:, vv, :) = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS2_t2_dB') ;
    NBRCS3_t2_dB(:, vv, :) = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS3_t2_dB') ;
    NBRCS4_t2_dB(:, vv, :) = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS4_t2_dB') ;
    NBRCS_t2_dB(:, vv, :) = readVar(dir_out_diffuse_NBRCS_tuple, 'NBRCS_t2_dB') ;

end


%% plotting as a function of fresnel zones
NBRCS11_inc = squeeze(NBRCS_t1_dB(1, :, :)) ;
NBRCS11_inc1 = squeeze(NBRCS1_t1_dB(1, :, :)) ;
NBRCS11_inc2 = squeeze(NBRCS2_t1_dB(1, :, :)) ;
NBRCS11_inc3 = squeeze(NBRCS3_t1_dB(1, :, :)) ;
NBRCS11_inc4 = squeeze(NBRCS4_t1_dB(1, :, :)) ;

plotFZ(NBRCS11_inc, NBRCS11_inc1, NBRCS11_inc2, NBRCS11_inc3, NBRCS11_inc4, pols{1, 1})

NBRCS12_inc = squeeze(NBRCS_t1_dB(2, :, :)) ;
NBRCS12_inc1 = squeeze(NBRCS1_t1_dB(2, :, :)) ;
NBRCS12_inc2 = squeeze(NBRCS2_t1_dB(2, :, :)) ;
NBRCS12_inc3 = squeeze(NBRCS3_t1_dB(2, :, :)) ;
NBRCS12_inc4 = squeeze(NBRCS4_t1_dB(2, :, :)) ;

plotFZ(NBRCS12_inc, NBRCS12_inc1, NBRCS12_inc2, NBRCS12_inc3, NBRCS12_inc4, pols{1, 2})

NBRCS21_inc = squeeze(NBRCS_t2_dB(1, :, :)) ;
NBRCS21_inc1 = squeeze(NBRCS1_t2_dB(1, :, :)) ;
NBRCS21_inc2 = squeeze(NBRCS2_t2_dB(1, :, :)) ;
NBRCS21_inc3 = squeeze(NBRCS3_t2_dB(1, :, :)) ;
NBRCS21_inc4 = squeeze(NBRCS4_t2_dB(1, :, :)) ;

plotFZ(NBRCS21_inc, NBRCS21_inc1, NBRCS21_inc2, NBRCS21_inc3, NBRCS21_inc4, pols{2, 1})

NBRCS22_inc = squeeze(NBRCS_t2_dB(2, :, :)) ;
NBRCS22_inc1 = squeeze(NBRCS1_t2_dB(2, :, :)) ;
NBRCS22_inc2 = squeeze(NBRCS2_t2_dB(2, :, :)) ;
NBRCS22_inc3 = squeeze(NBRCS3_t2_dB(2, :, :)) ;
NBRCS22_inc4 = squeeze(NBRCS4_t2_dB(2, :, :)) ;

plotFZ(NBRCS22_inc, NBRCS22_inc1, NBRCS22_inc2, NBRCS22_inc3, NBRCS22_inc4, pols{2, 2})


end

function plotFZ(NBRCS, NBRCS1, NBRCS2, NBRCS3, NBRCS4, pols)


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nr = SimParams.getInstance.Nr;
Nfz = SimParams.getInstance.Nfz;
% Ground Parameters
RMSH_cm = GndParams.getInstance.RMSH_list_cm( ParamsManager.index_RMSH );
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
th0_Tx_deg = th0_Tx_list_deg( ParamsManager.index_Th );
ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
ph0_Tx_deg = ph0_Tx_list_deg( ParamsManager.index_Ph );
% Receiver Parameters
hr_m = RxParams.getInstance.hr_m;
hpbw_deg = RxParams.getInstance.hpbw_deg;
SLL_dB = RxParams.getInstance.SLL_dB;
XPL_dB = RxParams.getInstance.XPL_dB;


%% GET GLOBAL DIRECTORIES
dir_fig_diffuse_NBRCS_vsFZ = SimulationFolders.getInstance.fig_diffuse_NBRCS_vsFZ;

fig_out_pathname = strcat( dir_fig_diffuse_NBRCS_vsFZ, '\TH_', num2str( th0_Tx_deg ) );
if ~exist(fig_out_pathname, 'dir')
    mkdir(fig_out_pathname)
end


FolderNameAnt = strcat('thsd_', num2str(hpbw_deg),...
    '-SLL_', num2str(SLL_dB), '-XPL_', num2str(XPL_dB)) ;

ylabely = 'NBRCS [dB]' ;
P_range = 40 ;

ymin1 = floor((min(min(NBRCS1))) / 10) * 10 ;
% ymax1 = ceil(max(max(max(NBRCS1))) / 10) * 10 ;
ymin2 = floor((min(min(NBRCS2))) / 10) * 10 ;
% ymax2 = ceil(max(max(max(B2_inc1))) / 10) * 10 ;
ymin3 = floor((min(min(NBRCS3))) / 10) * 10 ;
% ymax3 = ceil(max(max(max(B3_inc1))) / 10) * 10 ;
ymin4 = floor((min(min(NBRCS4))) / 10) * 10 ;
% ymax4 = ceil(max(max(max(B4_inc1))) / 10) * 10 ;
ymin = floor((min(min(NBRCS))) / 10) * 10 ;
% ymax = ceil(max(max(max(B4_inc))) / 10) * 10 ;

%%
figure
ymini = ymin1 ; ymaxi = ymin1 + P_range ;

subplot(2,2,1)
plot(squeeze(NBRCS1(1, :)), ':or')
hold
plot(squeeze(NBRCS1(2, :)), ':sb')
plot(squeeze(NBRCS1(3, :)), ':db')

title(strcat('(o^+, i^-)'))
text(1, ymaxi - 4, strcat('\theta_0=', num2str(th0_Tx_deg), '\circ'))
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
plot(squeeze(NBRCS2(1, :)), ':or')
hold
plot(squeeze(NBRCS2(2, :)), ':sb')
plot(squeeze(NBRCS2(3, :)), ':db')
title(strcat('(o_I^-, i^-)'))
text(1, ymaxi - 5, strrep( strrep( FolderNameAnt, '_', '=' ), '-', ',   ' ) )
axis([0 Nfz+1 ymini ymaxi])
grid

%%
ymini = ymin3 ; ymaxi = ymin3 + P_range ;
subplot(2,2,3)
plot(squeeze(NBRCS3(1, :)), ':or')
hold
plot(squeeze(NBRCS3(2, :)), ':sb')
plot(squeeze(NBRCS3(3, :)), ':db')
title(strcat('(o^+, i_I^+)'))
axis([0 Nfz + 1 ymini ymaxi])
grid
xlabel('Fresnel Zones')
ylabel(ylabely)

%%
ymini = ymin4 ; ymaxi = ymin4 + P_range ;
subplot(2,2,4)
plot(squeeze(NBRCS4(1, :)), ':or')
hold
plot(squeeze(NBRCS4(2, :)), ':sb')
plot(squeeze(NBRCS4(3, :)), ':db')
title(strcat('(o_I^-, i_I^+)'))
axis([0 Nfz + 1 ymini ymaxi])
grid
xlabel('Fresnel Zones')

fname = strcat('FZ-NBRCS_IND_dB', '-TH_', num2str( th0_Tx_deg ), '-PH_', num2str( ph0_Tx_deg ), '-RMSH_', num2str( RMSH_cm ) );
fname = strrep( fname, '.', 'dot' );
saveas(gcf, strcat(fig_out_pathname, '\', fname), 'tiff')
close


%%
figure
ymax = ymin + P_range ;
plot(squeeze(NBRCS(1, :)), ':or')
hold
plot(squeeze(NBRCS(2, :)), ':sb')
plot(squeeze(NBRCS(3, :)), ':db')
title(strcat('Medium - All Mechanism *** ', ...
    '\theta_0 = ', num2str(th0_Tx_deg), '\circ', ...
    ' *** N_r = ', num2str(Nr), ' ***', pols, '-pol'))
axis([0 Nfz + 1 ymin ymax])
grid
xlabel('Fresnel Zones')
ylabel('NBRCS [dB]')
text(1, ymax - 3, strrep( strrep( FolderNameAnt, '_', '=' ), '-', ',   ' ) )
text(1, ymax - 6, strcat('h_r = ', num2str(hr_m), 'm'))

fname = strcat('FZ-NBRCS_ALL_dB', '-TH_', num2str( th0_Tx_deg ), '-PH_', num2str( ph0_Tx_deg ), '-RMSH_', num2str( RMSH_cm ) );
fname = strrep( fname, '.', 'dot' );
saveas(gcf, strcat(FolderPath_write, '\',fname), 'tiff')
close

end