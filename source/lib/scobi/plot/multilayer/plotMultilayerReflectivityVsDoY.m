function plotMultilayerReflectivityVsDoY
% function plotMultilayer
%
%   Plots the mutli-layer reflectivity as a function of DoY (for the other  
%   variable parameters, i.e. transmitter incidence and azimuth angles, and 
%   RMSH fixed at specific values). 
%   This function is not called within the SCoBi simulation flow. The user
%   can run it after a full simulation for plotting purposes. 
%
%   See also plotDielProfiles, plotSMdata, plotReflectivityForProfiles.

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

% Find the input folder of the selected simulation 
dir_sims_input = strcat( selpath, '\input');
inputParamStructFile = strcat(dir_sims_input, '\', ConstantNames.INPUT_PARAMS_STRUCT_FILENAME );
load(inputParamStructFile, 'inputParamsStruct')

% Get all parameters used within this simulation
ParamsManager.initAllInputParamsFromInputParamsStruct( inputParamsStruct );


% %% GET GLOBAL DIRECTORIES
% % Reflectivity path
% dir_products_specular_reflectivity = strcat(selpath, '\products\specular\reflectivity' );
% 
% % Figure path
% dir_figure_specular_reflectivity_vsDoY = strcat(selpath, '\figure\specular\reflectivity\vs_DoY' );
% 
% 
% %% GET GLOBAL PARAMETERS
% % Simulation Parameters
% gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
% % Transmitter Parameters
% pol_Tx = TxParams.getInstance.pol_Tx;
% % Receiver Parameters
% pol_Rx = RxParams.getInstance.pol_Rx;
% % Configuration Parameters
% el0_Tx_list_deg = ConfigParams.getInstance.el0_Tx_list_deg;
% ph0_Tx_list_deg = ConfigParams.getInstance.ph0_Tx_list_deg;
% VSM_list_cm3cm3 = ConfigParams.getInstance.VSM_list_cm3cm3;
% RMSH_list_cm = ConfigParams.getInstance.RMSH_list_cm;
% 
% 
% %% GET CORRESPONDING DATA
% % % Find unique VSM values
% % VSM_list_cm3cm3 = unique(VSM_list_cm3cm3);
% % VSM_list_cm3cm3 = VSM_list_cm3cm3';
% 
% % % Find the coresponding parameter indices
% % RMHSH_indices = RMSH_list_cm == RMSH_cm_fixed;
% % th0_indices = el0_Tx_list_deg == el0_Tx_deg_fixed;
% % ph0_indices = ph0_Tx_list_deg == ph0_Tx_deg_fixed
% % common_indices = th0_indices & RMHSH_indices & ph0_indices;
% % common_indices = common_indices';
% 
% 
% % Mey be needed if power is analyzed rather than reflectivity
% % % read Kc (KKc_dB should be added to relectivity when received power is considered)
% % Kc = readComplexVar(dir_products_specular, 'Kc') ;
% % Kc_dB = 20 * log10(abs(Kc))] ;
% 
% 
% %% READ REFLECTIVITY AND FIND MATCHES
% % Read the reflectivity outputs of the simulation
% P_bare1 = readVar( dir_products_specular_reflectivity, 'Bare1' );
% P_bare2 = readVar( dir_products_specular_reflectivity, 'Bare2' );
% P_bare01 = readVar( dir_products_specular_reflectivity, 'Bare01' );
% P_bare02 = readVar( dir_products_specular_reflectivity, 'Bare02' );
% 
% % Get the corresponding reflectivity values
% P_bare1 = P_bare1(1:2, :);
% P_bare2 = P_bare2(1:2, :);
% P_bare01 = P_bare01(1:2, :);
% P_bare02 = P_bare02(1:2, :);
% 
% % If ground coverf is Vegetation
% if gnd_cover_id == Constants.ID_VEG_COVER
%     
%     % Read the reflectivity outputs of the simulation
%     P_veg1 = readVar( dir_products_specular_reflectivity, 'Veg1' );
%     P_veg2 = readVar( dir_products_specular_reflectivity, 'Veg2' );
%     P_veg01 = readVar( dir_products_specular_reflectivity, 'Veg01' );
%     P_veg02 = readVar( dir_products_specular_reflectivity, 'Veg02' );
%     
%     % Get the corresponding reflectivity values
%     P_veg1 = P_veg1(1:2, :);
%     P_veg2 = P_veg2(1:2, :);
%     P_veg01 = P_veg01(1:2, :);
%     P_veg02 = P_veg02(1:2, :);
% 
% end


%% PLOT
fig1 = figure(1) ;
fig2 = figure(2);


%% PLOT VSM data first
plotSMdata(fig1, fig2);
   
   
%% PLOT DIELECTRIC PROFILES    
plotDielProfiles( fig1, selpath );

    
% PLOT REFLECTIVITY
plotReflectivityForProfiles( Rp1, Rp2 )


end



% Reset the workspace
function resetWS

% Restore search path to defaults
restoredefaultpath

% Add the common "scobi" directory to the path to start running SCoBi
addpath( genpath( strrep( pwd, '\plot\multilayer', '' ) ) );

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

