%-e-------------------------*--. --- --. .--. ...*--------------------------
%
%                    %%%%%  %%%%%   %%%%%  %%%%% %%%
%                    %      %       %   %  %   %  %
%                    %%%%%  %       %   %  %%%%   %
%                        %  %       %   %  %   %  %
%                    %%%%%  %%%%%   %%%%%  %%%%  %%%
%
%
%--------------------------------------------------------------------------
%                         SCoBi-Veg v1.0
%
%    Copyright (C) 2018-2023 Mehmet Kurum, Orhan Eroglu
%--------------------------------------------------------------------------
% dx
%    This program is free software: You can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implid warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

%To-DO: Ensure about Copyright and GNU

function runSCoBi


%% WORKSPACE MANAGEMENT
% Reset workspace
resetWS();

% Add the common path to start running SCoBi
addpath( genpath( strcat(pwd, '/common') ) );


%% MAIN SCOBI WINDOW (SIMULATOR SELECTION)
% Open the main SCoBi GUI for simulator selection
simulator_id = gui_SCoBiMain();

% If no output from GUI, terminate program
if (isempty(simulator_id))
    
    fprintf('SCoBi terminated by the user\n');
    
    return
end

% Initialize workspace for the selected simulator
initSCoBiWS( simulator_id );


%% START GUI FOR THE SELECTED SIMULATOR
inputStruct = startSelectedGUI( simulator_id );

% If no output from GUI, terminate program
if (isempty(inputStruct))
    
    fprintf('SCoBi terminated by the user\n');
    
    return
end

if simulator_id == Constants.id_multi_layer
    
    runSCoBiML( simulator_id, inputStruct );
    
elseif simulator_id == Constants.id_veg_agr || ...
       simulator_id == Constants.id_veg_for
   
   runSCoBiVeg( simulator_id, inputStruct );
   
end


end 


function resetWS

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


function runSCoBiML( simulator_id, inputStruct )
    
    % Get input and check validity
    getInput(simulator_id, inputStruct);

    % TO-DO: Input Validity check
    isInputValid = initWithInputs();
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % Create simulations' directories
    SimulationFolders.getInstance.makeStaticDirs();
    isInputValid = 1;
        
    
    % If input is valid
    if isInputValid

        % To-DO: Find a more structural way
        % Write all inputs to a text file and show it
        writeToFile( simulator_id, inputStruct );


        %% GET GLOBAL PARAMETERS
        % Simulation Settings
        sim_mode_id = SimSettings.getInstance.sim_mode_id;
        % Dynamic Parameters
        DoYs = DynParams.getInstance.DoYs;
        th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
        ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
        VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
        RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;

        
        num_DoYs = length( DoYs );
        num_Th = length( th0_Tx_list_deg );
        num_Ph = length( ph0_Tx_list_deg );
        num_VSM = length( VSM_list_cm3cm3 );
        num_RMSH = length( RMSH_list_cm );

        % Snapshot simulation
        if sim_mode_id == Constants.id_snapshot

            % For each theta (looking angle)
            for tt = 1 : num_Th

                % Set theta index
                ParamsManager.index_Th( tt );

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

                            % Initialize the directories depending on theta,
                            % phi, VSM, and RMSH
                            SimulationFolders.getInstance.initializeDynamicDirs();

                            % Call SCoBi main flow
                            mainSCoBi();

                        end

                    end

                end

            end


        % Time-series simulation
        else

            % For each corresponding tuple of theta (looking angle), 
            % phi (azimuth angle), VSM (volumetric soil moisture), and 
            % RMSH (root mean square height - roughness)
            for ii = 1 : 10 : num_DoYs  % The length of each is equal

                % Set theta, phi, VSM, ad RMSH index the same
                ParamsManager.index_DoY( ii );
                ParamsManager.index_Th( ii );
                ParamsManager.index_Ph( ii );
                ParamsManager.index_VSM( ii );
                ParamsManager.index_RMSH( ii );

                % Initialize the directories depending on theta, phi, VSM, and 
                % RMSH
                SimulationFolders.getInstance.initializeDynamicDirs();

                % Call SCoBi main flow
                mainSCoBi();

            end

        end
        

%         % TO-DO: Check for multiple theta or phi values
%         ParamsManager.index_Th(1);
%         ParamsManager.index_Ph(1);
%         ParamsManager.index_VSM(1);
%         ParamsManager.index_RMSH(1);
%         % Initialize the directories depending on theta,
%         % phi, VSM, and RMSH
%         SimulationFolders.getInstance.initializeDynamicDirs();
% 
%         multiLayerModel();

    % Else if input is NOT valid
    else

        return

    end
    
end


function runSCoBiVeg( simulator_id, inputStruct )

    % Get input and check validity
    getInput(simulator_id, inputStruct);

    isInputValid = initWithInputs();

    % If input is valid
    if isInputValid

        % To-DO: Find a more structural way
        % Write all inputs to a text file and show it
        writeToFile( simulator_id, inputStruct );


        %% GET GLOBAL PARAMETERS
        % Simulation Settings
        sim_mode_id = SimSettings.getInstance.sim_mode_id;
        % Dynamic Parameters
        DoYs = DynParams.getInstance.DoYs;
        th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
        ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
        VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
        RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;

        
        num_DoYs = length( DoYs );
        num_Th = length( th0_Tx_list_deg );
        num_Ph = length( ph0_Tx_list_deg );
        num_VSM = length( VSM_list_cm3cm3 );
        num_RMSH = length( RMSH_list_cm );

        % Snapshot simulation
        if sim_mode_id == Constants.id_snapshot

            % For each theta (looking angle)
            for tt = 1 : num_Th

                % Set theta index
                ParamsManager.index_Th( tt );

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

                            % Initialize the directories depending on theta,
                            % phi, VSM, and RMSH
                            SimulationFolders.getInstance.initializeDynamicDirs();

                            % Call SCoBi main flow
                            mainSCoBi;

                        end

                    end

                end

            end


        % Time-series simulation
        else

            % For each corresponding tuple of theta (looking angle), 
            % phi (azimuth angle), VSM (volumetric soil moisture), and 
            % RMSH (root mean square height - roughness)
            for ii = 1 : num_DoYs  % The length of each is equal

                % Set theta, phi, VSM, ad RMSH index the same
                ParamsManager.index_DoY( ii );
                ParamsManager.index_Th( ii );
                ParamsManager.index_Ph( ii );
                ParamsManager.index_VSM( ii );
                ParamsManager.index_RMSH( ii );

                % Initialize the directories depending on theta, phi, VSM, and 
                % RMSH
                SimulationFolders.getInstance.initializeDynamicDirs();

                % Call SCoBi main flow
                mainSCoBi;

            end

        end

    % Else if input is NOT valid
    else

        return

    end
    
end



function inputStruct = startSelectedGUI( simulator_id )

inputStruct = [];

% TO-DO: Check this!
% If the OS is not UNIX OR it is MAC and Matlab version below 7.14
if (~isunix || (ismac && verLessThan('matlab', '7.14')))
    
    % If simulator is SCoBi-Veg (Agriculture)
    if simulator_id == Constants.id_veg_agr
        
        inputStruct = gui_SCoBi_Veg;
    
    % Else if simulator is SCoBi-Veg (Forest)
    elseif simulator_id == Constants.id_veg_for
        
        inputStruct = gui_SCoBi_Veg;
    
    % Else if simulator is SCoBi-ML (MultiLayer)
    elseif simulator_id == Constants.id_multi_layer
        
        inputStruct = gui_SCoBi_ML;
        
    else
        
        % GUIs for new simulators should be added here.
     
    end
    
end

end

function writeToFile( simulator_id, inputStruct )


%% GET GLOBAL DIRECTORIES
dir_simulations = SimulationFolders.getInstance.simulations;
inputFile = strcat(dir_simulations, '\', 'input_report.txt');


%% GET GLOBAL PARAMETERS
% Simulation Settings
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
sim_mode_id = SimSettings.getInstance.sim_mode_id;
% Simulation Parameters
veg_method_id = SimParams.getInstance.veg_method_id;
% Receiver Parameters
orientation_Rx_id = RxParams.getInstance.orientation_Rx_id;
ant_pat_Rx_id = RxParams.getInstance.ant_pat_Rx_id;


% Open file to write inputs
fileID = fopen(inputFile,'w');


% Write the selected simulator name
simulators = Constants.simulators;
simulatorString = simulators{ 1, simulator_id };
fprintf(fileID, sprintf( strcat('++++++++++++\t\t', simulatorString, '\t\t++++++++++++\n' ) ) );

% Write the current date time
t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')
fprintf(fileID, strcat( string(t), '\n\n' ) );


%%SIMULATION PARAMETERS
fprintf(fileID, strcat( 'Simulation name:\t', inputStruct.sim_name, '\n' ) );
fprintf(fileID, strcat( 'Campaign:\t\t\t', inputStruct.campaign, '\n' ) );
fprintf(fileID, strcat( 'Campaign Date:\t\t', inputStruct.campaign_date, '\n' ) );
fprintf(fileID, strcat( 'Campaign Plot:\t\t', inputStruct.plot, '\n' ) );

% If gnd_cover is Vegetation, then look at sim_mode
if gnd_cover_id == Constants.id_veg_cover
    
    % If sim_mode is Snapshot, then write veg_method
    if sim_mode_id == Constants.id_snapshot

        fprintf(fileID, strcat( 'Vegetation method:\t', inputStruct.veg_method, '\n' ) );
    
    end
    
end  

% If veg_method is 'Virtual' only, then write virtual vegetation
% orientation
if veg_method_id == Constants.id_veg_vir
    fprintf(fileID, strcat( 'Virtual Vegetation orientation:\t', inputStruct.veg_vir_orientation, '\n' ) );
end

% If gnd_cover is Vegetation, then write veg_plant
if gnd_cover_id == Constants.id_veg_cover
    
    fprintf(fileID, strcat( 'Vegetation plant:\t', inputStruct.veg_plant, '\n' ) );   
    
end

if simulator_id == Constants.id_veg_agr ...
        || simulator_id == Constants.id_veg_for
    % Write number of realizations
    fprintf(fileID, strcat( 'Number of realizations:\t', num2str(inputStruct.Nr), '\n' ) );
    % Write number of Fresnel zones
    fprintf(fileID, strcat( 'Number of Fresnel zones:\t', num2str(inputStruct.Nfz), '\n' ) );
end



%% SIMULATION SETTINGS
% Write sim_mode
fprintf(fileID, strcat( 'Simulation mode:\t', inputStruct.sim_mode, '\n' ) );
% Write gnd_cover
fprintf(fileID, strcat( 'Ground cover:\t\t', inputStruct.gnd_cover, '\n' ) );

if simulator_id == Constants.id_veg_agr ...
        || simulator_id == Constants.id_veg_for
    
    % Write user preferences
    fprintf(fileID, strcat( 'Write attenuation?:\t\t', num2str(inputStruct.write_attenuation), '\n' ) );
    fprintf(fileID, strcat( 'Calculate Direct Term?:\t\t', num2str(inputStruct.calc_direct_term), '\n' ) );
    fprintf(fileID, strcat( 'Calculate Specular Term?:\t\t', num2str(inputStruct.calc_specular_term), '\n' ) );
    if sim_mode_id == Constants.id_snapshot ...
            && gnd_cover_id == Constants.id_veg_cover
        fprintf(fileID, strcat( 'Calculate Diffuse Term?:\t\t', num2str(inputStruct.calc_diffuse_term), '\n' ) );
    end
    
end

%fprintf(fileID, strcat( 'Draw live plots?:\t\t', num2str(inputStruct.draw_live_plots), '\n' ) );


%% TRANSMITTER PARAMETERS
fprintf(fileID, strcat( '\n\nTRANSMITTER PARAMETERS\n' ) );
fprintf(fileID, strcat( 'Operating frequency (MHz):\t', num2str(inputStruct.f_MHz), '\n' ) );
fprintf(fileID, strcat( 'Range to Earth center (km):\t', num2str(inputStruct.r_Tx_km), '\n' ) );
fprintf(fileID, strcat( 'EIRP (dB):\t\t\t\t\t', num2str(inputStruct.EIRP_dB), '\n' ) );
fprintf(fileID, strcat( 'Polarization:\t\t\t\t', inputStruct.pol_Tx, '\n' ) );


%% RECEIVER PARAMETERS
fprintf(fileID, strcat( '\n\nRECEIVER PARAMETERS\n' ) );
fprintf(fileID, strcat( 'Altitude (m):\t\t\t', num2str(inputStruct.hr_m), '\n' ) );
fprintf(fileID, strcat( 'Gain (dB):\t\t\t\t', num2str(inputStruct.G0r_dB), '\n' ) );
fprintf(fileID, strcat( 'Polarization:\t\t\t', inputStruct.pol_Rx, '\n' ) );
fprintf(fileID, strcat( 'Orientation:\t\t\t', inputStruct.orientation_Rx, '\n' ) );

if orientation_Rx_id == Constants.id_Rx_fixed

    fprintf(fileID, strcat( 'Incidence angle (deg):\t', num2str(inputStruct.th0_Rx_deg), '\n' ) );
    fprintf(fileID, strcat( 'Azimuth angle (deg):\t', num2str(inputStruct.ph0_Rx_deg), '\n' ) );

end

fprintf(fileID, strcat( 'Antenna Pattern:\t\t', inputStruct.ant_pat_Rx, '\n' ) );

% If receiver antenna pattern is Generalized-Gaussian
if ant_pat_Rx_id == Constants.id_Rx_GG
    
    fprintf(fileID, strcat( 'Antenna Pattern Resolution (deg):\t', num2str(inputStruct.ant_pat_res_deg_Rx), '\n' ) );
    
    fprintf(fileID, strcat( 'Half-power beamwidth (deg):\t\t\t', num2str(inputStruct.hpbw_deg), '\n' ) );
    fprintf(fileID, strcat( 'Side-lobe level (dB):\t\t\t\t', num2str(inputStruct.SLL_dB), '\n' ) );
    fprintf(fileID, strcat( 'Cross-polarization level (dB):\t\t', num2str(inputStruct.XPL_dB), '\n' ) );

% Else if receiver antenna pattern is User-defined
elseif ant_pat_Rx_id == Constants.id_Rx_user_defined
    
    filename = strrep(inputStruct.ant_pat_Rx_file, '\', '/');
    fprintf(fileID, strcat( 'Antenna pattern file:\t', filename, '\n' ) );
    
end


%% GROUND PARAMETERS
fprintf(fileID, strcat( '\n\nGROUND PARAMETERS\n' ) );
    fprintf(fileID, strcat( 'Dielectric model:\t', inputStruct.diel_model, '\n' ) );
% If the simulator is SCoBi-Veg (Agriculture OR Forest)
if simulator_id == Constants.id_veg_agr ...
        || simulator_id == Constants.id_veg_for
        
    fprintf(fileID, strcat( 'Sand ratio:\t\t', num2str(inputStruct.sand_ratio), '\n' ) );
    fprintf(fileID, strcat( 'Clay ratio:\t\t', num2str(inputStruct.clay_ratio), '\n' ) );
    fprintf(fileID, strcat( 'Bulk density (g/cm3):\t', num2str(inputStruct.rhob_gcm3), '\n' ) );
    
% Else if simulator is SCoBi-ML (Multilayer)
elseif simulator_id == Constants.id_multi_layer
    
    filename = strrep(inputStruct.dyn_inputs_file, '\', '/');
    fprintf(fileID, strcat( 'Ground input file:\t', filename, '\n' ) );
    
end


%% VEGETATION PARAMETERS
fprintf(fileID, strcat( '\n\nVEGETETATION PARAMETERS\n' ) );
filename = strrep(inputStruct.veg_inputs_file, '\', '/');
fprintf(fileID, strcat( 'Vegetation input file:\t', filename, '\n' ) );



fclose(fileID);

winopen( inputFile );

end