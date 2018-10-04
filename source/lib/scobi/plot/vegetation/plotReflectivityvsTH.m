% Mehmet Kurum
% modified 11/13/2017
% Reflectivity vs. Observation Angle

function plotReflectivityvsTH

resetWS();


%% GET GLOBAL DIRECTORIES
dir_sims_main = Directories.getInstance.sims_main;


%% GET USER INPUTS
% Get the simulation folder
selpath = uigetdir( dir_sims_main );

% Fixed parameter values for plot
VSM_cm3cm3_fixed = 0.15;
RMSH_cm_fixed = 1.0;
ph0_Tx_deg_fixed = 15;

% Find the input folder of the selected simulation 
dir_sims_input = strcat( selpath, '\input');
inputParamStructFile = strcat(dir_sims_input, '\', ConstantNames.inputParamsStruct_filename );
load(inputParamStructFile, 'inputParamsStruct')

% Get all parameters used within this simulation
ParamsManager.initAllInputParamsFromInputParamsStruct( inputParamsStruct );


%% GET GLOBAL DIRECTORIES
% Reflectivity path
dir_products_specular_reflectivity = strcat(selpath, '\products\specular\reflectivity' );

% Figure path
dir_figure_specular_reflectivity_vsTH = strcat(selpath, '\figure\specular\reflectivity\vs_TH' );


%% GET GLOBAL PARAMETERS
% Simulation Parameters
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Configuration Parameters
th0_Tx_list_deg = ConfigParams.getInstance.th0_Tx_list_deg;
ph0_Tx_list_deg = ConfigParams.getInstance.ph0_Tx_list_deg;
VSM_list_cm3cm3 = ConfigParams.getInstance.VSM_list_cm3cm3;
RMSH_list_cm = ConfigParams.getInstance.RMSH_list_cm;


%% GET CORRESPONDING DATA
% Find unique incidence angles
th0_Tx_list_deg = unique(th0_Tx_list_deg);
th0_Tx_list_deg = th0_Tx_list_deg';

% Find the coresponding parameter indices
VSM_indices = VSM_list_cm3cm3 == VSM_cm3cm3_fixed;
RMHSH_indices = RMSH_list_cm == RMSH_cm_fixed;
ph0_indices = ph0_Tx_list_deg == ph0_Tx_deg_fixed
common_indices = VSM_indices & RMHSH_indices & ph0_indices;
common_indices = common_indices';

% % read Kc (KKc_dB should be added to relectivity when received power is considered)
% Kc = readComplexVar(dir_products_specular, 'Kc') ;
% Kc_dB = 20 * log10(abs(Kc))] ;


%% READ REFLECTIVITY AND FIND MATCHES
% Read the reflectivity outputs of the simulation
P_bare1 = readVar( dir_products_specular_reflectivity, 'Bare1' );
P_bare2 = readVar( dir_products_specular_reflectivity, 'Bare2' );
P_bare01 = readVar( dir_products_specular_reflectivity, 'Bare01' );
P_bare02 = readVar( dir_products_specular_reflectivity, 'Bare02' );

% Get the corresponding reflectivity values
P_bare1 = P_bare1(1:2, common_indices);
P_bare2 = P_bare2(1:2, common_indices);
P_bare01 = P_bare01(1:2, common_indices);
P_bare02 = P_bare02(1:2, common_indices);

% If ground coverf is Vegetation
if gnd_cover_id == Constants.id_veg_cover
    
    % Read the reflectivity outputs of the simulation
    P_veg1 = readVar( dir_products_specular_reflectivity, 'Veg1' );
    P_veg2 = readVar( dir_products_specular_reflectivity, 'Veg2' );
    P_veg01 = readVar( dir_products_specular_reflectivity, 'Veg01' );
    P_veg02 = readVar( dir_products_specular_reflectivity, 'Veg02' );
    
    % Get the corresponding reflectivity values
    P_veg1 = P_veg1(1:2, common_indices);
    P_veg2 = P_veg2(1:2, common_indices);
    P_veg01 = P_veg01(1:2, common_indices);
    P_veg02 = P_veg02(1:2, common_indices);

end


%% PLOT
% Reflectivity in dB 
figure

if (pol_Tx == 'X') && (pol_Tx == 'X')

    bareVV_dB = 10 * log10(squeeze(P_bare1(1, :))) ;
    bareHH_dB = 10 * log10(squeeze(P_bare2(2, :))) ;

    plot(th0_Tx_list_deg, bareVV_dB, ':or') % co-pol VV
    hold
    plot(th0_Tx_list_deg, bareHH_dB, '-or') % co-pol HH
    
    % If ground coverf is Vegetation
    if gnd_cover_id == Constants.id_veg_cover
        vegVV_dB = 10 * log10(squeeze(P_veg1(1, :))) ;
        vegHH_dB = 10 * log10(squeeze(P_veg2(2, :))) ;

        plot(th0_Tx_list_deg, vegVV_dB, ':sg') % co-pol VV
        plot(th0_Tx_list_deg, vegHH_dB, '-sg') % co-pol HH
    end
    
else
    
    BareCO_dB = 10 * log10(squeeze(P_bare1(1, :))) ;
    BareX_dB = 10 * log10(squeeze(P_bare1(2, :))) ;

    plot(th0_Tx_list_deg, BareCO_dB, ':or') % co-pol
    hold
    plot(th0_Tx_list_deg, BareX_dB, '-or') % x-pol
    
    % If ground coverf is Vegetation
    if gnd_cover_id == Constants.id_veg_cover
        VegCO_dB = 10 * log10(squeeze(P_veg1(1, :))) ;
        VegX_dB = 10 * log10(squeeze(P_veg1(2, :))) ;

        plot(th0_Tx_list_deg, VegCO_dB, ':sg') % co-pol
        plot(th0_Tx_list_deg, VegX_dB, '-sg') % x-pol
    end
end

grid
xlabel('Incidence Angle [degrees]')
ylabel('Reflectivity [dB]')

if (pol_Tx == 'R') && (pol_Rx == 'R')

    legend('bare: RR-pol', 'bare: RL-pol', 'veg: RR-pol', 'veg: RL-pol', 'Location', 'southeast')
    title(strcat('Reflectivity - Tx: RHCP, Rx: RHCP', '- PHI:', num2str( ph0_Tx_deg_fixed ), '\circ', '- VSM:', num2str( VSM_cm3cm3_fixed ), '%', '- RMSH:', num2str( RMSH_cm_fixed ), 'cm'))

elseif (pol_Tx == 'R') && (pol_Rx == 'X')

    legend('bare: RX-pol', 'bare: RY-pol', 'veg: RX-pol', 'veg: RY-pol', 'Location', 'southeast')
    title(strcat('Reflectivity - Tx: RHCP, Rx: LINEAR', '- PHI:', num2str( ph0_Tx_deg_fixed ), '\circ', '- VSM:', num2str( VSM_cm3cm3_fixed ), '%', '- RMSH:', num2str( RMSH_cm_fixed ), 'cm'))

elseif (pol_Tx == 'X') && (pol_Tx == 'X')

    legend('bare: XX-pol', 'bare: YY-pol', 'veg: XX-pol', 'veg: YY-pol', 'Location', 'southeast')
    title(strcat('Reflectivity - Tx: LINEAR, Rx: LINEAR', '- PHI:', num2str( ph0_Tx_deg_fixed ), '\circ', '- VSM:', num2str( VSM_cm3cm3_fixed ), '%', '- RMSH:', num2str( RMSH_cm_fixed ), 'cm'))
end

axis([0 80 -50 0])

fname = strcat('Reflectivity_vs_TH', '-PH_', num2str( ph0_Tx_deg_fixed ), '-VSM_', num2str( VSM_cm3cm3_fixed ), '-RMSH_', num2str( RMSH_cm_fixed ) );
fname = strrep( fname, '.', 'dot' );

if ~exist(dir_figure_specular_reflectivity_vsTH, 'dir')        
    mkdir(dir_figure_specular_reflectivity_vsTH)
end

saveas(gcf, strcat(dir_figure_specular_reflectivity_vsTH, '\',fname), 'tiff')
close


end



% Reset the workspace
function resetWS

% Restore search path to defaults
restoredefaultpath

% Add the common "scobi" directory to the path to start running SCoBi
addpath( genpath( strrep(pwd, '\plot\vegetation', '') ) );


% Store current debug breakpoints before doing clear all
myBreakpoints = dbstatus;
save('myBreakpoints.mat', 'myBreakpoints');

% Clear all the workspace
clear all;
clc ;

% Restore debug breakpoints
load('myBreakpoints.mat');
dbstop(myBreakpoints);
clear myBreakpoints;

if ( exist('myBreakpoints.mat','file') ) 
    delete('myBreakpoints.mat'); 
end

end