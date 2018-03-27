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
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

%To-DO: Ensure about Copyright and GNU

function runSCoBi

% Reset the current state of the workspace
resetCurrentState;

%startGUI;

% Get input and check validity
getInput;

isInputValid = initWithInputs();

% If input is valid
if isInputValid
    
    num_Th = length( SatParams.getInstance.th0_deg );
    num_Ph = length( SatParams.getInstance.PH0_deg );
    num_VSM = length( GndParams.getInstance.VSM_cm3cm3 );
    num_RMSH = length( GndParams.getInstance.RMSH_cm );
    
    % Snapshot simulation
    if SimSettings.getInstance.sim_mode == Constants.sim_mode.SNAPSHOT
        
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



function resetCurrentState

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
if (exist('myBreakpoints.mat','file')) delete('myBreakpoints.mat'); end

% Add all subdirectories to the path
addpath( genpath( strcat(pwd, '/constants') ) ); % This is required to first addpath
addpath( genpath( Directories.getInstance.input ) );
addpath( genpath( Directories.getInstance.bistatic ) );
addpath( genpath( Directories.getInstance.gui ) );
addpath( genpath( Directories.getInstance.init ) );
addpath( genpath( Directories.getInstance.monte_carlo ) );
addpath( genpath( Directories.getInstance.param ) );
addpath( genpath( Directories.getInstance.products ) );
addpath( genpath( Directories.getInstance.util ) );

end



function startGUI

if (~isunix || (ismac && verLessThan('matlab', '7.14')))
    [mode, mode_vinc, mode_data, mode_ref, flag_ms_pos, flag_ms, flag_ge, flag_cov, flag_NTRIP, flag_amb, ...
        flag_skyplot, flag_plotproc, flag_var_dyn_model, flag_stopGOstop, flag_SP3, flag_SBAS, flag_IAR, ...
        filerootIN, filerootOUT, filename_R_obs, filename_M_obs, ...
        filename_nav, filename_ref, filename_pco, pos_M_man, protocol_idx, multi_antenna_rf, iono_model, tropo_model,fsep_char] = gui_goGPS;
end

end