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

getInput;

isInputValid = initWithInputs();

if isInputValid
    
    num_Th = length( SatParams.getInstance.th0_deg );
    num_Ph = length( SatParams.getInstance.PH0_deg );
    num_VSM = length( GndParams.getInstance.VSM );
    num_RMSH = length( GndParams.getInstance.RMSH );
    
    % Cross-simulation
    if SimSettings.getInstance.sim_mode == Constants.sim_mode.CROSS
        
        for tt = 1 : num_Th
            
            ParamsManager.index_Th( tt );
            
            for pp = 1 : num_Ph
            
                ParamsManager.index_Ph( pp );
                
                for vv = 1 : num_VSM
            
                    ParamsManager.index_VSM( vv );
                    
                    for rr = 1 : num_RMSH
                        
                        ParamsManager.index_RMSH( rr );
                        
                        SimulationFolders.getInstance.initializeDynamicDirs();
                        
                        mainSCoBi;
                        
                    end
                    
                end
                
            end
            
        end
        
        
    % Time-series simulation
    else
        
        for ii = 1 : num_Th  % The length of each is equal
            
            ParamsManager.index_Th( ii );
            ParamsManager.index_Ph( ii );
            ParamsManager.index_VSM( ii );
            ParamsManager.index_RMSH( ii );
            
            SimulationFolders.getInstance.initializeDynamicDirs();
            
            mainSCoBi;
        
        end
        
    end
    
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

% Include all subdirectories
addpath( genpath( strcat(pwd, '/dir') ) ); % This is required to first addpath
addpath( genpath( Directories.getInstance.input ) );
addpath( genpath( Directories.getInstance.gui ) );
addpath( genpath( Directories.getInstance.init ) );
addpath( genpath( Directories.getInstance.monte_carlo ) );
addpath( genpath( Directories.getInstance.param ) );
addpath( genpath( Directories.getInstance.products ) );
addpath( genpath( Directories.getInstance.SCoBi ) );
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