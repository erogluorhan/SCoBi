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

simulator_id = 2;

initSCoBiWS( simulator_id );

% startGUI;

% if (~isunix || (ismac && verLessThan('matlab', '7.14')))
%     [mode, mode_vinc, mode_data, mode_ref, flag_ms_pos, flag_ms, flag_ge, flag_cov, flag_NTRIP, flag_amb, ...
%         flag_skyplot, flag_plotproc, flag_var_dyn_model, flag_stopGOstop, flag_SP3, flag_SBAS, flag_IAR, ...
%         filerootIN, filerootOUT, filename_R_obs, filename_M_obs, ...
%         filename_nav, filename_ref, filename_pco, pos_M_man, protocol_idx, multi_antenna_rf, iono_model, tropo_model,fsep_char] = gui_SCoBi;
% end
% 
% if (isempty(mode))
%     
%     fprintf('SCoBi terminated by the user\n');
%     
%     return
% end

if simulator_id == Constants.id_multi_layer
    
    inputFile_sys = 'default_input-scobi_ml.xml';

    inputFile_veg = 'vegHomInput-Paulownia.xml';
% 
    inputFile_dyn = 'default_input-ml.xlsx';

    % Get input and check validity
    getInput(inputFile_sys, inputFile_veg, inputFile_dyn);

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

    Main_MSU_ML_process
    
elseif simulator_id == Constants.id_veg

    inputFile_sys = 'default_input-scobi_veg.xml';
    % inputFile_sys = 'sysInput-Paulownia-PAPER_PBAND-hr_20.xml';
    % inputFile_sys = 'sysInput-Corn.xml';

    inputFile_veg = 'vegHomInput-Paulownia.xml';
    % inputFile_veg = 'vegVirRowInput-CORN_PAPER-row0-R1_R4.xml';

    inputFile_dyn = 'default_input-snapshot.xlsx';
    % inputFile_dyn = 'default_input-time_series.xlsx';

    % Get input and check validity
    getInput(inputFile_sys, inputFile_veg, inputFile_dyn);

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


end 



function startGUI

if (~isunix || (ismac && verLessThan('matlab', '7.14')))
    [mode, mode_vinc, mode_data, mode_ref, flag_ms_pos, flag_ms, flag_ge, flag_cov, flag_NTRIP, flag_amb, ...
        flag_skyplot, flag_plotproc, flag_var_dyn_model, flag_stopGOstop, flag_SP3, flag_SBAS, flag_IAR, ...
        filerootIN, filerootOUT, filename_R_obs, filename_M_obs, ...
        filename_nav, filename_ref, filename_pco, pos_M_man, protocol_idx, multi_antenna_rf, iono_model, tropo_model,fsep_char] = gui_SCoBi;
end

if (isempty(mode))
    
    fprintf('SCoBi terminated by the user\n');
    
    return
    
end

end