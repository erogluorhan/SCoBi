
function plotReflectivityVsEL
% function plotReflectivityVsEL 
%
%   Plots the reflectivity as a function of transmitter elevation angle 
%   (for the other variable parameters, i.e. transmitter azimuth, VSM, RMSH,
%   fixed at specific values). 
%   This function is not called within the SCoBi simulation flow. The user
%   can run it after a full simulation for plotting purposes. 
%
%   It is applicable (as is) to the default Agriculture inputs only, which
%   is provided with this version's release (v1.0.0). Modifying this
%   function to apply for other scenarios is upon the developer's
%   initiative.
%
%   See also plotReflectivityVsVSM.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



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
inputParamStructFile = strcat(dir_sims_input, '\', ConstantNames.INPUT_PARAMS_STRUCT_FILENAME );
load(inputParamStructFile, 'inputParamsStruct')

% Get all parameters used within this simulation
ParamsManager.initAllInputParamsFromInputParamsStruct( inputParamsStruct );


%% GET GLOBAL DIRECTORIES
% Reflectivity path
dir_products_specular_reflectivity = strcat(selpath, '\products\specular\reflectivity' );

% Figure path
dir_figure_specular_reflectivity_vsEL = strcat(selpath, '\figure\specular\reflectivity\vs_EL' );


%% GET GLOBAL PARAMETERS
% Simulation Parameters
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
% Configuration Parameters
el0_Tx_list_deg = ConfigParams.getInstance.el0_Tx_list_deg;
ph0_Tx_list_deg = ConfigParams.getInstance.ph0_Tx_list_deg;
VSM_list_cm3cm3 = ConfigParams.getInstance.VSM_list_cm3cm3;
RMSH_list_cm = ConfigParams.getInstance.RMSH_list_cm;


%% GET CORRESPONDING DATA
% Find unique elevation angles
el0_Tx_list_deg = unique(el0_Tx_list_deg);
el0_Tx_list_deg = el0_Tx_list_deg';

% Find the coresponding parameter indices
VSM_indices = VSM_list_cm3cm3 == VSM_cm3cm3_fixed;
RMHSH_indices = RMSH_list_cm == RMSH_cm_fixed;
ph0_indices = ph0_Tx_list_deg == ph0_Tx_deg_fixed;
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
if gnd_cover_id == Constants.ID_VEG_COVER
    
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

    plot(el0_Tx_list_deg, bareVV_dB, ':or') % co-pol VV
    hold
    plot(el0_Tx_list_deg, bareHH_dB, '-or') % co-pol HH
    
    % If ground coverf is Vegetation
    if gnd_cover_id == Constants.ID_VEG_COVER
        vegVV_dB = 10 * log10(squeeze(P_veg1(1, :))) ;
        vegHH_dB = 10 * log10(squeeze(P_veg2(2, :))) ;

        plot(el0_Tx_list_deg, vegVV_dB, ':sg') % co-pol VV
        plot(el0_Tx_list_deg, vegHH_dB, '-sg') % co-pol HH
    end
    
else
    
    BareCO_dB = 10 * log10(squeeze(P_bare1(1, :))) ;
    BareX_dB = 10 * log10(squeeze(P_bare1(2, :))) ;

    plot(el0_Tx_list_deg, BareCO_dB, ':or') % co-pol
    hold
    plot(el0_Tx_list_deg, BareX_dB, '-or') % x-pol
    
    % If ground coverf is Vegetation
    if gnd_cover_id == Constants.ID_VEG_COVER
        VegCO_dB = 10 * log10(squeeze(P_veg1(1, :))) ;
        VegX_dB = 10 * log10(squeeze(P_veg1(2, :))) ;

        plot(el0_Tx_list_deg, VegCO_dB, ':sg') % co-pol
        plot(el0_Tx_list_deg, VegX_dB, '-sg') % x-pol
    end
end

grid
xlabel('Elevation Angle [degrees]')
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

fname = strcat('Reflectivity_vs_EL', '-PH_', num2str( ph0_Tx_deg_fixed ), '-VSM_', num2str( VSM_cm3cm3_fixed ), '-RMSH_', num2str( RMSH_cm_fixed ) );
fname = strrep( fname, '.', 'dot' );

if ~exist(dir_figure_specular_reflectivity_vsEL, 'dir')        
    mkdir(dir_figure_specular_reflectivity_vsEL)
end

saveas(gcf, strcat(dir_figure_specular_reflectivity_vsEL, '\',fname), 'tiff')
close


end



% Reset the workspace
function resetWS

% Restore search path to defaults
restoredefaultpath

% Add the common "scobi" directory to the path to start running SCoBi
addpath( genpath( strrep( pwd, '\plot\scobi', '' ) ) );

% Add "input" directory to the path
addpath( genpath( Directories.getInstance.input ) );


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