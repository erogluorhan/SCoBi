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
%     isInputValid = initWithInputs();
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % Create but not creating simulations' directories
    SimulationFolders.getInstance.makeStaticDirs();
    
    % TO-DO: Check for multiple theta or phi values
    ParamsManager.index_Th(1);
    ParamsManager.index_Ph(1);
    ParamsManager.index_VSM(1);
    ParamsManager.index_RMSH(1);
    % Initialize the directories depending on theta,
    % phi, VSM, and RMSH
    SimulationFolders.getInstance.initializeDynamicDirs();
    
    mainSCoBiML();
    
end


function runSCoBiVeg( simulator_id, inputStruct )

    % Get input and check validity
    getInput(simulator_id, inputStruct);

    isInputValid = initWithInputs();

    % If input is valid
    if isInputValid


        %% GET GLOBAL PARAMETERS
        % Simulation Settings
        sim_mode_id = SimSettings.getInstance.sim_mode_id;
        % Dynamic Parameters
        th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
        ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
        VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
        RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;

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
            for ii = 1 : num_Th  % The length of each is equal

                % Set theta, phi, VSM, ad RMSH index the same
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