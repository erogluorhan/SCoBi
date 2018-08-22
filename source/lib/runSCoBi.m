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
%                         SCoBi v1.0
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

%TO-DO: Ensure about Copyright and GNU

function runSCoBi


%% WORKSPACE MANAGEMENT
% Add the common "scobi" directory to the path to start running SCoBi
addpath( genpath( strcat(pwd, '/scobi') ) );

% Add "input" directory to the path
addpath( genpath( Directories.getInstance.input ) );

% Reset workspace
resetWS();


%% GUI: MAIN SCOBI WINDOW (SIMULATOR SELECTION)
% Open the main SCoBi GUI for simulator selection
simulator_id = gui_SCoBiMain();

% If no output from GUI, terminate program
if (isempty(simulator_id))
    
    fprintf('SCoBi terminated by the user\n');
    
    return
    
end


%% GUI: START THE SELECTED SIMULATOR'S GUI
inputStruct = startSelectedGUI( simulator_id );

% If no output from GUI, terminate program
if (isempty(inputStruct))
    
    fprintf('SCoBi terminated by the user\n');
    
    return

end

% Initialize all input parameters by using the inputs from GUI
varout = initAllInputParams(simulator_id, inputStruct);

% Check input validity
[isInputValid, sim_counter_start] = initWithInputs( varout );


%% SIMULATIONS
% If input is valid
if isInputValid
    
    
    %% GET GLOBAL PARAMETERS
    num_sims = ParamsManager.num_sims;
    

    % Write all inputs to a text file and show it
    ParamsManager.writeInputReports( inputStruct );
    
    % Run the simulations
    for ii = sim_counter_start : num_sims

        ParamsManager.sim_counter( ii );

        % Call SCoBi main flow
        mainSCoBi;

    end

% Else if input is NOT valid
else

    return

end

end 