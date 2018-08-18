
% Mehmet Kurum
% 11/12/2017
% Received Power vs. Observation Angle


function plotDiffuse_ReceivedPower_vsTH


%% INPUT FILES
inputFile_sys = 'sysInput-Paulownia-PAPER_PBAND-hr_20.xml';
inputFile_veg = 'vegHomInput-Paulownia.xml';


%% GET INPUT FOR INITIAL START
getInput( inputFile_sys, inputFile_veg );

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();

% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();


%% GET GLOBAL PARAMETERS
% Simulation Parameters
sim_mode_id = SimSettings.getInstance.sim_mode_id;
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Configuration Parameters
th0_Tx_list_deg = ConfigParams.getInstance.th0_Tx_list_deg;
ph0_Tx_list_deg = ConfigParams.getInstance.ph0_Tx_list_deg;
VSM_list_cm3cm3 = ConfigParams.getInstance.VSM_list_cm3cm3;
RMSH_list_cm = ConfigParams.getInstance.RMSH_list_cm;

num_Th = length( th0_Tx_list_deg );
num_Ph = length( ph0_Tx_list_deg );
num_VSM = length( VSM_list_cm3cm3 );
num_RMSH = length( RMSH_list_cm );


%% Choices
% TO-DO: Check here
FZ_choice = 10 ;
SM_choice = 2 ;


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


%% reading
% Received Power
PP1_inc_dB = zeros(2, 2, num_Th) ; PP1_inc1_dB = zeros(2, 2, num_Th) ; 
PP1_inc2_dB = zeros(2, 2, num_Th) ; PP1_inc3_dB = zeros(2, 2, num_Th) ; PP1_inc4_dB = zeros(2, 2, num_Th) ;


if sim_mode_id == Constants.id_snapshot
    
    % For each phi (azimuth angle)
    for pp = 1 : num_Ph

        % Set phi index
        ParamsManager.index_Ph( pp );

        % For each VSM (volumetric soil moisture)
        for vv = 1 : num_VSM

            % Set VSM index
            ParamsManager.index_VSM( vv );

            % For each RMSH (root mean square height - roughness)
            for rr = 1 : num_RMSH

                % Set RMSH index
                ParamsManager.index_RMSH( rr );
    
                for tt = 1 : num_Th

                    % Set theta index
                    ParamsManager.index_Th( tt );

                    % Initialize the directories depending on theta,
                    % phi, VSM, and RMSH
                    SimulationFolders.getInstance.initializeDynamicDirs();


                    %% GET GLOBAL DIRECTORIES
                    dir_out_diffuse_P1_tuple = SimulationFolders.getInstance.out_diffuse_P1_tuple;

                    % transmit port 1 / receiver ports 1&2
                    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t1_dB') ;
                    PP1_inc_dB(1, :, tt) = squeeze(xx(:, FZ_choice)) ;    
                    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc1_t1_dB') ;
                    PP1_inc1_dB(1, :, tt) = squeeze(xx(:, FZ_choice)) ;    
                    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc2_t1_dB') ;
                    PP1_inc2_dB(1, :, tt) = squeeze(xx(:, FZ_choice)) ;    
                    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc3_t1_dB') ;
                    PP1_inc3_dB(1, :, tt) = squeeze(xx(:, FZ_choice)) ;    
                    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc4_t1_dB') ;
                    PP1_inc4_dB(1, :, tt) = squeeze(xx(:, FZ_choice)) ;

                    % transmit port 2 / receiver ports 1&2
                    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t2_dB') ;
                    PP1_inc_dB(2, :, tt) = squeeze(xx(:, FZ_choice)) ;    
                    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc1_t2_dB') ;
                    PP1_inc1_dB(2, :, tt) = squeeze(xx(:, FZ_choice)) ;    
                    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc2_t2_dB') ;
                    PP1_inc2_dB(2, :, tt) = squeeze(xx(:, FZ_choice)) ;    
                    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc3_t2_dB') ;
                    PP1_inc3_dB(2, :, tt) = squeeze(xx(:, FZ_choice)) ;    
                    xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc4_t2_dB') ;
                    PP1_inc4_dB(2, :, tt) = squeeze(xx(:, FZ_choice)) ;

                end
                
                %% Plot as a function of TH
                plotANG(FZ_choice, pols, PP1_inc_dB, PP1_inc1_dB, PP1_inc2_dB, PP1_inc3_dB, PP1_inc4_dB)

            end

        end
        
    end               
    
% Time-series simulation
else
    
    % For each corresponding tuple of theta (looking angle), 
    % phi (azimuth angle), VSM (volumetric soil moisture), and 
    % RMSH (root mean square height - roughness)
    for ii = 1 : num_Th  % The length of each is equal

        % Set theta, phi, VSM, ad RMSH index the same
        ParamsManager.index_Th( ii );
        ParamsManager.index_Ph( ii );
        ParamsManager.index_VSM( ii );
        ParamsManager.index_RMSH( ii );

        % Initialize the directories depending on theta, phi, VSM, and 
        % RMSH
        SimulationFolders.getInstance.initializeDynamicDirs();

        %% TO-DO: PLOT

    end

end


end

function  plotANG( FZ_choice, pols,...
    PP1_inc_dB, PP1_inc1_dB, PP1_inc2_dB, PP1_inc3_dB, PP1_inc4_dB)


%% GET GLOBAL DIRECTORIES
dir_fig_diffuse_P1_vsTH = SimulationFolders.getInstance.fig_diffuse_P1_vsTH;

fig_out_pathname = strcat( dir_fig_diffuse_P1_vsTH, '\FZ_', num2str(FZ_choice) );
if ~exist(fig_out_pathname, 'dir')
    mkdir(fig_out_pathname)
end

            
%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nr = SimParams.getInstance.Nr;
% Configuration Parameters
VSM_cm3cm3 = ConfigParams.getInstance.VSM_list_cm3cm3( ParamsManager.index_VSM );
RMSH_cm = ConfigParams.getInstance.RMSH_list_cm( ParamsManager.index_RMSH );
th0_Tx_list_deg = ConfigParams.getInstance.th0_Tx_list_deg;
ph0_Tx_list_deg = ConfigParams.getInstance.ph0_Tx_list_deg;
ph0_Tx_deg = ph0_Tx_list_deg( ParamsManager.index_Ph );
% Receiver Parameters
hr_m = RxParams.getInstance.hr_m;
hpbw_deg = RxParams.getInstance.hpbw_deg;
SLL_dB = RxParams.getInstance.SLL_dB;
XPL_dB = RxParams.getInstance.XPL_dB;


FolderNameAnt = strcat('thsd_', num2str( hpbw_deg ),...
                '-SLL_', num2str( SLL_dB ), '-XPL_', num2str( XPL_dB )) ;

%%
P_range = 40 ;
ylabely = 'Received Power [dB]' ;

%%

figure
ymin1 = floor(min(min(min(PP1_inc_dB(1, :, :)))) / 10) * 10 ;
ymin = ymin1 - 10 ;
ymax = ymin1 + P_range ;
plot(th0_Tx_list_deg, squeeze(PP1_inc_dB(1, 1, :)), '-or')
hold
plot(th0_Tx_list_deg, squeeze(PP1_inc_dB(1, 2, :)), ':or')
plot(th0_Tx_list_deg, squeeze(PP1_inc_dB(2, 1, :)), '-sb')
plot(th0_Tx_list_deg, squeeze(PP1_inc_dB(2, 2, :)), ':sb')

axis([th0_Tx_list_deg(1) - 10 th0_Tx_list_deg(end) + 10 ymin ymax])
% text(1, yavg1, pols{1})
% text(1, yavg2, pols{2})
text(10, ymax - 3, strcat('FZ=', num2str(FZ_choice), ',   N_r=', num2str(Nr), ',   h_r=', num2str(hr_m), 'm'))
% text(30, ymax - 5, strcat('N_r = ', num2str(Nr)))
text(10, ymax - 6, strrep( strrep( FolderNameAnt, '_', '=' ), '-', ',   ' ) )
text(10, ymax - 9, strcat('VSM=', num2str( VSM_cm3cm3 ), 'cm^3/cm^3'))
% text(20, ymax - 5, strcat('h_r = ', num2str(hr_m), 'm'))
% text(30, ymax - 10, strcat('Medium'))
grid
xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
yticks(ymin:10:ymax)

xlabel('Incidence Angle [degrees]')
ylabel(ylabely)
ttlstr = 'Diffuse - All Mechanisms - Medium' ;
title(ttlstr)
legend(pols{1}, pols{2}, pols{3}, pols{4})

fname = strcat('TH-ReceivedPower_ALL_dB', '-PH_', num2str( ph0_Tx_deg ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) );
fname = strrep( fname, '.', 'dot' );

saveas(gcf, strcat(fig_out_pathname, '\', fname), 'tiff')
close


%%
figure
ymin1 = floor(min(min(min(PP1_inc1_dB(1, :, :)))) / 10) * 10 ;
ymin2 = floor(min(min(min(PP1_inc2_dB(1, :, :)))) / 10) * 10 ;
ymin3 = floor(min(min(min(PP1_inc3_dB(1, :, :)))) / 10) * 10 ;
ymin4 = floor(min(min(min(PP1_inc4_dB(1, :, :)))) / 10) * 10 ;

ymin = ymin1 - 10 ;
ymax = ymin1 + P_range ;

%
subplot(2, 2, 1)
plot(th0_Tx_list_deg, squeeze(PP1_inc1_dB(1, 1, :)), '-or')
hold
plot(th0_Tx_list_deg, squeeze(PP1_inc1_dB(1, 2, :)), ':or')
plot(th0_Tx_list_deg, squeeze(PP1_inc1_dB(2, 1, :)), '-sb')
plot(th0_Tx_list_deg, squeeze(PP1_inc1_dB(2, 2, :)), ':sb')

axis([th0_Tx_list_deg(1) - 10 th0_Tx_list_deg(end) + 10 ymin ymax])
% text(1, yavg1, pols{1})
% text(1, yavg2, pols{2})
text(5, ymax - 5, strcat('FZ=', num2str(FZ_choice), ',  h_r=', num2str(hr_m), 'm', ',  N_r=', num2str(Nr)))
% text(10, ymax - 10, strcat('N_r = ', num2str(Nr)))
text(5, ymax - 12, strrep( strrep( FolderNameAnt, '_', '=' ), '-', ', ') )
% text(30, ymax - 5, strcat('h_r = ', num2str(hr_m), 'm'))
% text(30, ymax - 10, strcat('Medium'))
grid
xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
yticks(ymin:10:ymax)

% xlabel('Observation Angle [degrees]')
ylabel(ylabely)
ttlstr = 'Diffuse - (o^+, i^-)' ;
title(ttlstr)
% legend(pols{1}, pols{2}, pols{3}, pols{4})

%

ymin = ymin2 - 10 ;
ymax = ymin2 + P_range ;

subplot(2, 2, 2)
plot(th0_Tx_list_deg, squeeze(PP1_inc2_dB(1, 1, :)), '-or')
hold
plot(th0_Tx_list_deg, squeeze(PP1_inc2_dB(1, 2, :)), ':or')
plot(th0_Tx_list_deg, squeeze(PP1_inc2_dB(2, 1, :)), '-sb')
plot(th0_Tx_list_deg, squeeze(PP1_inc2_dB(2, 2, :)), ':sb')

axis([th0_Tx_list_deg(1) - 10 th0_Tx_list_deg(end) + 10 ymin ymax])
text(5, ymax - 5, strcat('VSM=', num2str( VSM_cm3cm3 ), 'cm^3/cm^3'))
grid
xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
yticks(ymin:10:ymax)
% xlabel('Observation Angle [degrees]')
% ylabel(ylabely)
ttlstr = 'Diffuse - (o_I^-, i^-)' ;
title(ttlstr)
% legend(pols{1}, pols{2}, pols{3}, pols{4})

%
ymin = ymin3 - 10 ;
ymax = ymin3 + P_range ;
subplot(2, 2, 3)
plot(th0_Tx_list_deg, squeeze(PP1_inc3_dB(1, 1, :)), '-or')
hold
plot(th0_Tx_list_deg, squeeze(PP1_inc3_dB(1, 2, :)), ':or')
plot(th0_Tx_list_deg, squeeze(PP1_inc3_dB(2, 1, :)), '-sb')
plot(th0_Tx_list_deg, squeeze(PP1_inc3_dB(2, 2, :)), ':sb')

axis([th0_Tx_list_deg(1) - 10 th0_Tx_list_deg(end) + 10 ymin ymax])
grid
xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
yticks(ymin:10:ymax)
xlabel('Incidence Angle [degrees]')
ylabel(ylabely)
ttlstr = 'Diffuse - (o^+, i_I^+)' ;
title(ttlstr)
% legend(pols{1}, pols{2}, pols{3}, pols{4})

%
ymin = ymin4 - 10 ;
ymax = ymin4 + P_range ;
subplot(2, 2, 4)
plot(th0_Tx_list_deg, squeeze(PP1_inc4_dB(1, 1, :)), '-or')
hold
plot(th0_Tx_list_deg, squeeze(PP1_inc4_dB(1, 2, :)), ':or')
plot(th0_Tx_list_deg, squeeze(PP1_inc4_dB(2, 1, :)), '-sb')
plot(th0_Tx_list_deg, squeeze(PP1_inc4_dB(2, 2, :)), ':sb')

axis([th0_Tx_list_deg(1) - 10 th0_Tx_list_deg(end) + 10 ymin ymax])
grid
xticks(th0_Tx_list_deg(1) : 10 : th0_Tx_list_deg(end))
yticks(ymin:10:ymax)
xlabel('Incidence Angle [degrees]')
% ylabel(ylabely)
ttlstr = 'Diffuse - (o_I^-, i_I^+)' ;
title(ttlstr)
hleg1 = legend(pols{1}, pols{2}, pols{3}, pols{4}) ;
set(hleg1,'position',[0.4 0.4 0.2 0.2])

%%
fname = strcat('TH-ReceivedPower_IND_dB', '-PH_', num2str( ph0_Tx_deg ), '-VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str( RMSH_cm ) );
fname = strrep( fname, '.', 'dot' );

saveas(gcf, strcat( fig_out_pathname, '\', fname), 'tiff')
close

end