%--------------------------------------------------------------------------
% function runSCoBi
%
%   Simulation engine function. 
%
%   runSCoBi starts the simulation, calls GUI classes to get the user 
%   inputs, performs input validity check by using ParamsManager class and 
%   runs the simulation iterations for the required number of simulations.
%
%   See also mainSCoBi.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program (SCoBi) is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.3



%                    %%%%%  %%%%%   %%%%%  %%%%% %%%
%                    %      %       %   %  %   %  %
%                    %%%%%  %       %   %  %%%%   %
%                        %  %       %   %  %   %  %
%                    %%%%%  %%%%%   %%%%%  %%%%  %%%
%
%
%--------------------------------------------------------------------------
%                         SCoBi v1.0.3
%
%    Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd
%--------------------------------------------------------------------------
% 
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
%    along with this program (COPYING.txt).  If not, 
%    see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
%  %%%%%%%%%%%%%%%%%%%%%%%%%%  UPDATE HISTORY  %%%%%%%%%%%%%%%%%%%%%%%%%  %
%
%   Version 1.0.3
%
%   June 6, 2020
%
% Added functionality for backwards compatibility for inputStruct files made
% with older versions of SCoBi.

function runSCoBi


%% WORKSPACE MANAGEMENT
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
inputStruct = [];

% If the OS is not UNIX OR it is MAC and Matlab version below 7.14
if (~isunix || (ismac && verLessThan('matlab', '7.14')))
        
    [simulator_id, inputStruct ] = gui_SCoBi(simulator_id);
    
end

% If no output from GUI, terminate program
if (isempty(inputStruct))
    
    fprintf('SCoBi terminated by the user\n');
    
    return

end

% Initialize all input parameters by using the inputs from GUI
[varout, inputStruct] = initAllInputParams(simulator_id, inputStruct);

% Check input validity
[isInputValid, sim_counter_start] = initWithInputs( varout );


%% SIMULATIONS
% If input is valid
if isInputValid
    
    
    %% GET GLOBAL PARAMETERS
    num_sims = ParamsManager.num_sims;
    

    % Write all inputs to a text file and show it
    ParamsManager.writeInputReports( inputStruct );

    %% START SIMULATIONS
    disp('++++++++++++++++++++++++++++   START SIMULATIONS   ++++++++++++++++++++++++++++++++')
    
    sim_start = datetime('now');
    disp( strcat('Simulation Start: ', string(sim_start) ) )
    
    % Run the simulations
    for ii = sim_counter_start : num_sims

        ParamsManager.sim_counter( ii );

        % Call SCoBi main flow
        mainSCoBi();

    end


    %% END SIMULATIONS
    disp('-------------------------------   END SIMULATIONS   -------------------------------')
    
    
    disp( strcat('Simulation Start: ', string(sim_start) ) )
    sim_stop = datetime('now');
    disp( strcat('Simulation End: ', string(sim_stop) ) )
    duration = sim_stop - sim_start          
    disp( strcat('Duration: ', string(duration) ) )

% Else if input is NOT valid
else

    return

end

end 


% Reset the workspace
function resetWS

% Restore search path to defaults
restoredefaultpath

% Add the common "scobi" directory to the path to start running SCoBi
addpath( genpath( strcat(pwd, '\scobi') ) );

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