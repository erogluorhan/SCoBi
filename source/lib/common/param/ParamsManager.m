classdef ParamsManager < handle
    %% PARAMSMANAGER CLASS - Check the validity of input parameters
    %
    
    
    methods (Static)
                
        function result = index_Th( value )
           
            persistent index_Th
            
            % Set function
            if nargin ~= 0
                
                %% GET GLOBAL PARAMETERS
                % Dynamic Parameters
                th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
                % Receiver Parameters
                orientation_Rx_id = RxParams.getInstance.orientation_Rx_id;
                
                index_Th = value;
                
                th0_Tx_deg = th0_Tx_list_deg( index_Th );
                
                % If receiver orientation is Specular-facing, then set
                % receiver th0_Rx_deg equal to transmitter th0_Tx_deg
                if orientation_Rx_id == Constants.id_Rx_specular_facing
                    RxParams.getInstance.set_th0_Rx_deg( th0_Tx_deg );
                end
            
            end
            % Get function
            
            result = index_Th;
            
        end
                
        function result = index_Ph( value )
           
            persistent index_Ph
            
            % Set function
            if nargin ~= 0
                
                %% GET GLOBAL PARAMETERS
                % Dynamic Parameters
                ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
                % Receiver Parameters
                orientation_Rx_id = RxParams.getInstance.orientation_Rx_id;
                
                index_Ph = value;
                
                ph0_Tx_deg = ph0_Tx_list_deg( index_Ph );
                
                % If receiver orientation is Specular-facing, then set
                % receiver ph0_Rx_deg equal to transmitter ph0_Tx_deg
                if orientation_Rx_id == Constants.id_Rx_specular_facing
                    RxParams.getInstance.set_ph0_Rx_deg( ph0_Tx_deg );
                end
                
            end
            % Get function
            
            result = index_Ph;
            
        end
                
        function result = index_VSM( value )
           
            persistent index_VSM
            
            if nargin ~= 0
                index_VSM = value;
            end
            
            result = index_VSM;
            
        end
                
        function result = index_RMSH( value )
           
            persistent index_RMSH
            
            if nargin ~= 0
                index_RMSH = value;
            end
            
            result = index_RMSH;
            
        end
        
        function [validityResult, terminateResult, terminateMsg] = isInputValid()
                        
            %% GET GLOBAL DIRECTORIES
            dir_hr = SimulationFolders.getInstance.hr;
            
            
            %% GET GLOBAL PARAMETERS
            sim_name = SimParams.getInstance.sim_name;
            
            
            %% DEFINE PERSISTENT VARIABLES
            persistent isValid 
            persistent isTerminate
            persistent msg
            
            isValid = 0;
            isTerminate = 0;
            
            
            if exist( strcat( dir_hr, '\', ConstantNames.input_params_filename), 'file') == 2 
                        
                isSysInputEqual = ParamsManager.isSysInputEqual();
                isVegInputEqual = ParamsManager.isVegInputEqual();
                isInputValidForSimMode = ParamsManager.isInputValidForSimMode();
                
                isValid = isSysInputEqual & isVegInputEqual & isInputValidForSimMode;

                if ~isValid
                    
                    % If the reason for invalidity is the simulation mode
                    % (Cross simulation OR Time-series)
                    if ~isInputValidForSimMode
                        
                        isTerminate = 1;
                        msg = 'Input for dynamic system parameters are not valid for simulation mode. SCoBi simulator will be terminated for user check!';
                        
                    % Else if the reason is input inequality with the
                    % existing simulation folders
                    else

                        % TO-DO: Design SCoBi-s own GUI
                        %newSimName;

                        prompt = {'Either one or more inputs conflict with the existing simulations. Enter new simulation name:'};
                        dlg_title = 'Input Conflict!';
                        num_lines = 1;
                        defaultans = {sim_name};
                        new_sim_name = char( inputdlg(prompt,dlg_title,num_lines,defaultans) );  

                        if ~isempty(new_sim_name)
                            SimParams.getInstance.updateSimName( new_sim_name );

                            % Initializing simulations' directories for the updated simulation name
                            SimulationFolders.getInstance.initialize();

                            ParamsManager.isInputValid();
                        else
                            isTerminate = 1;
                            msg = 'A simulation with the same name exists, but input_params file is missing. SCoBi simulator will be terminated for user check!';
                        end
                    end
                end
            
            else
                
                % The first time this simulation is being created
                if ~exist( dir_hr)
                
                    isValid = ParamsManager.isInputValidForSimMode();
                
                % Simulation exists, but input_params_filename looks 
                % to be deleted. That is an error. User should check!
                else
                    isValid = 0;
                    isTerminate = 1;
                    msg = 'Input_params file is missing from the simulation directory; it might have been deleted. SCoBi simulator will be terminated for user check!'
                end
                
            end
            
            validityResult = isValid;
            terminateResult = isTerminate;
            terminateMsg = msg;
            
        end
        
        
        function isValid = isInputValidForSimMode
            
            %% GET GLOBAL PARAMETERS
            % Simulation Parameters
            sim_mode_id = SimSettings.getInstance.sim_mode_id;
            % Dynamic Parameters
            th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
            num_Th = length( th0_Tx_list_deg );
            ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
            num_Ph = length( ph0_Tx_list_deg );
            VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
            num_VSM = length( VSM_list_cm3cm3 );
            RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;
            num_RMSH = length( RMSH_list_cm );
            
            
            isValid = 0;
                            
            maxnum = max( [num_Th, num_Ph, num_VSM, num_RMSH] );
            
            % Simulation mode is Time-series
            if sim_mode_id == Constants.id_time_series
                
                % Add a control for timestamps as well
                if (num_Th == maxnum || num_Th == 1) && ...
                   (num_Ph == maxnum || num_Ph == 1) && ...
                   (num_VSM == maxnum || num_VSM == 1 ) && ...
                   (num_RMSH == maxnum || num_RMSH == 1)
                    
                    isValid = 1;
                    
                    if num_Th == 1
                        DynParams.getInstance.rep_Th(maxnum);
                    end
                    
                    if num_Ph == 1
                        DynParams.getInstance.rep_Ph(maxnum);
                    end
                    
                    if num_VSM == 1
                        DynParams.getInstance.rep_VSM(maxnum);
                    end                    
                    
                    if num_RMSH == 1
                        DynParams.getInstance.rep_RMSH(maxnum);
                    end

                end
                
            % Simulation mode is Snapshot
            elseif sim_mode_id == Constants.id_snapshot
                
                if ~( num_Th == 0 || num_Ph == 0 || num_VSM == 0 || ...
                        num_RMSH == 0 )
                
                    isValid = 1;
                
                end
                
            end
            
        end
        
        
        function isEqual = isSysInputEqual
                        
            %% GET GLOBAL DIRECTORIES
            dir_hr = SimulationFolders.getInstance.hr;            
            
            
            %% GET GLOBAL PARAMETERS
            % Simulation Parameters
            Nfz = SimParams.getInstance.Nfz;
            % Ground Parameters
            sand_ratio = GndParams.getInstance.sand_ratio;
            clay_ratio = GndParams.getInstance.clay_ratio;
            rhob_gcm3 = GndParams.getInstance.rhob_gcm3;
            % Transmitter Parameters
            f_MHz = TxParams.getInstance.f_MHz;
            r_Tx_m = TxParams.getInstance.r_Tx_m;
            EIRP_dB = TxParams.getInstance.EIRP_dB;
            
            
            %% LOAD META-DATA           
            load( strcat( dir_hr, '\', ConstantNames.input_params_filename), ConstantNames.input_params_var );
            
            isEqual = 0;

            if input_params(ConstantNames.sim_Nfz) == Nfz
                if input_params(ConstantNames.gnd_sand_ratio) == sand_ratio
                    if input_params(ConstantNames.gnd_clay_ratio) == clay_ratio
                        if input_params(ConstantNames.gnd_rhob_gcm3) == rhob_gcm3
                            if input_params(ConstantNames.Tx_f_MHz) == f_MHz
                                if input_params(ConstantNames.r_Tx_m) == r_Tx_m
                                    if input_params(ConstantNames.Tx_EIRP_dB) == EIRP_dB
                                        isEqual = 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
        end
        
        function isEqual = isVegInputEqual
                        
            %% GET GLOBAL DIRECTORIES
            dir_hr = SimulationFolders.getInstance.hr;  
            
            
            %% GET GLOBAL PARAMETERS
            % Simulation Parameters
            veg_method_id = SimParams.getInstance.veg_method_id;
            veg_vir_orientation_id = SimParams.getInstance.veg_vir_orientation_id;
            
            
            %% LOAD INPUT PARAMS
            load( strcat( dir_hr, '\', ConstantNames.input_params_filename), ConstantNames.input_params_var );
            
            
            isEqual = 0;
            
            % If vegetation method is homogenous
            if veg_method_id == Constants.id_veg_hom
            
                %% GET GLOBAL PARAMETERS
                % Vegetation Parameters
                vegetation_stage = VegParams.getInstance.vegetation_stage;
                TYPES = VegParams.getInstance.TYPES;
                dim_layers_m = VegParams.getInstance.dim_layers_m;
                TYPKND = VegParams.getInstance.TYPKND;
                scat_cal_veg = VegParams.getInstance.scat_cal_veg;
                LTK = VegParams.getInstance.LTK;
                dsty = VegParams.getInstance.dsty;
                dim1_m = VegParams.getInstance.dim1_m;
                dim2_m = VegParams.getInstance.dim2_m;
                dim3_m = VegParams.getInstance.dim3_m;
                epsr = VegParams.getInstance.epsr;
                parm1_deg = VegParams.getInstance.parm1_deg;
                parm2_deg = VegParams.getInstance.parm2_deg;

                if strcmp( input_params(ConstantNames.veg_vegetationStage), vegetation_stage )
                    if isequal( input_params(ConstantNames.veg_hom_TYPES), TYPES )
                        if isequal( input_params(ConstantNames.veg_hom_dimLayers_m), dim_layers_m )
                            if isequal( input_params(ConstantNames.veg_hom_TYPKND), TYPKND )                                          
                                if isequal( input_params(ConstantNames.veg_hom_scatCalVeg), scat_cal_veg )
                                    if isequal( input_params(ConstantNames.veg_hom_LTK), LTK )
                                        if isequal( input_params(ConstantNames.veg_hom_dsty), dsty )
                                            if isequal( input_params(ConstantNames.veg_hom_dim1_m), dim1_m )
                                                if isequal( input_params(ConstantNames.veg_hom_dim2_m), dim2_m )
                                                    if isequal( input_params(ConstantNames.veg_hom_dim3_m), dim3_m )
                                                        if isequal( input_params(ConstantNames.veg_hom_epsr), epsr )
                                                            if isequal( input_params(ConstantNames.veg_hom_parm1_deg), parm1_deg )
                                                                if isequal( input_params(ConstantNames.veg_hom_parm2_deg), parm2_deg )
                                                                    isEqual = 1;
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
            % Else if vegetation method is virtual    
            else
                
                % If veg_method is Virtual and veg_vir_orientation is Row-crop
                if veg_vir_orientation_id == Constants.id_veg_vir_row_crop
                    %% GET GLOBAL PARAMETERS
                    % Virtual Row-Structured Vegetation Parameters
                    vegetation_stage = VegParams.getInstance.vegetation_stage;
                    plugin = VegVirRowParams.getInstance.plugin;
                    row_space_m = VegVirRowParams.getInstance.row_space_m;
                    col_space_m = VegVirRowParams.getInstance.col_space_m;
                    phi_row_deg = VegVirRowParams.getInstance.phi_row_deg;
                    seed_fluctuation_m = VegVirRowParams.getInstance.seed_fluctuation_m;

                    if strcmp( input_params(ConstantNames.veg_vegetationStage), vegetation_stage )
                        if input_params(ConstantNames.veg_vir_row_plugin).isTheSame( plugin )
                            if input_params(ConstantNames.veg_vir_row_rowSpace_m) == row_space_m
                                if input_params(ConstantNames.veg_vir_row_colSpace_m) == col_space_m
                                    if input_params(ConstantNames.veg_vir_row_phiRow_deg) == phi_row_deg
                                        if input_params(ConstantNames.veg_vir_row_seedFluctuation_m) == seed_fluctuation_m
                                            isEqual = 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                  
                % Else if veg_method is virtual random
                elseif veg_vir_orientation_id == Constants.id_veg_vir_random_spread
                   
                    % TO-DO: Implement for virtual random vegetation 
                end
            end
            
        end
        
        function [result, Nr_current, dispMsg] = isToGenerateScatPos()
            
            %% GET GLOBAL DIRECTORIES
            dir_position = SimulationFolders.getInstance.position;
            
            %% GET GLOBAL PARAMETERS
            % Simulation Settings
            calc_specular_term = SimSettings.getInstance.calc_specular_term;
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
            % Simulation Parameters
            Nr = SimParams.getInstance.Nr;
            veg_method_id = SimParams.getInstance.veg_method_id;
            % Vegetation Parameters            
            readExistingVegParams;
            scat_cal_veg = VegParams.getInstance.scat_cal_veg;
            num_scat_cal = sum(sum(sum(scat_cal_veg))) ;
            
            % First check the user preferences
            if veg_method_id == Constants.id_veg_hom
                if ~calc_diffuse_term
                    dispMsg = 'Generate Scatterer Positions - SKIPPED (User Preferences - No Diffuse Term)';
                    Nr_current = NaN;
                    result = Constants.need_for_run.NO;
                    return
                end
            
            elseif veg_method_id == Constants.id_veg_vir
                if ~calc_diffuse_term && ~calc_specular_term
                    dispMsg = 'Generate Scatterer Positions - SKIPPED (User Preferences - No Specular and Diffuse Term)';
                    Nr_current = NaN;
                    result = Constants.need_for_run.NO;
                    return
                end
            end
                                    
            % If passes user preferences, check the existence of files
            all_files = dir(dir_position);
            num_files = numel(all_files) - 2;
            Nr_current = num_files / num_scat_cal ;

            if isnan(Nr_current) || isinf(Nr_current), Nr_current = 0; end
            
            if Nr_current >= Nr
                dispMsg = 'Generate Scatterer Positions - SKIPPED - Already exists!';
                result = Constants.need_for_run.NO;
            else
                if Nr_current > 0
                    dispMsg = 'Generate Scatterer Positions - Partially exists!';
                    result = Constants.need_for_run.PARTIAL;
                else
                    dispMsg = 'Generate Scatterer Positions';
                    result = Constants.need_for_run.FULL;
                end
            end
            
        end
        
        function [result, Nr_current, dispMsg] = isToCalculateFScatAmp()
                        
            %% GET GLOBAL DIRECTORIES
            dir_fscat = SimulationFolders.getInstance.fscat;
            
            %% GET GLOBAL PARAMETERS
            % Simulation Settings
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
            % Simulation Parameters
            Nr = SimParams.getInstance.Nr;
            % Vegetation Parameters
            scat_cal_veg = VegParams.getInstance.scat_cal_veg;
            num_scat_cal = sum(sum(sum(scat_cal_veg))) ;
                        
            % First check the user preferences
            if ~calc_diffuse_term
                dispMsg = 'Scattering Amplitudes - SKIPPED (User Preferences - No Diffuse Term)';
                Nr_current = NaN;
                result = Constants.need_for_run.NO;
                return
            end
            
            % If passes user preferences, check the existence of files
            all_files = dir(dir_fscat);
            num_files = numel(all_files) - 2;
            Nr_current = num_files / num_scat_cal / Constants.factor_fscat ;

            if isnan(Nr_current), Nr_current = 0; end
            
            if Nr_current >= Nr
                dispMsg = 'Scattering Amplitudes - SKIPPED - Already exists!';
                result = Constants.need_for_run.NO;
            else
                if Nr_current > 0
                    dispMsg = 'Scattering Amplitudes - Partially exists!';
                    result = Constants.need_for_run.PARTIAL;
                else
                    dispMsg = 'Scattering Amplitudes';
                    result = Constants.need_for_run.FULL;
                end
            end
            
        end
        
        function [result, write_attenuation, dispMsg] = isToCalcPropagation()
            
            %% GET GLOBAL DIRECTORIES
            dir_afsa = SimulationFolders.getInstance.afsa;
            
            %% GET GLOBAL PARAMETERS
            % Simulation Settings
            calc_specular_term = SimSettings.getInstance.calc_specular_term;
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
            write_attenuation = SimSettings.getInstance.write_attenuation;
                        
            % First check the user preferences
            if ~( calc_specular_term || calc_diffuse_term )
                dispMsg = 'Propagation - SKIPPED (User Preferences - No Specular and Diffuse Term)';
                result = Constants.need_for_run.NO;
                return
            end
            
            % If passes user preferences, check the existence of files
            all_files = dir(dir_afsa);
            num_files = numel(all_files) - 2;
            
            if num_files >= Constants.num_afsa
                dispMsg = 'Propagation - SKIPPED - Already exists!';
                result = Constants.need_for_run.NO;
            else
                dispMsg = 'Propagation';
                result = Constants.need_for_run.FULL;
            end
            
        end
        
        function [result, dispMsg] = isToCalcRxAntPatMatrix()
                       
            %% GET GLOBAL DIRECTORIES
            dir_ant_lookup = SimulationFolders.getInstance.ant_lookup;
            
            % If passes user preferences, check the existence of files
            all_files = dir(dir_ant_lookup);
            num_files = numel(all_files) - 2;
            
            if num_files >= Constants.num_ant_lookup
                dispMsg = 'Antenna Pattern Matrix - SKIPPED - Already exists!';
                result = Constants.need_for_run.NO;
            else
                dispMsg = 'Antenna Pattern Matrix';
                result = Constants.need_for_run.FULL;
            end
            
        end
        
        function [result, Nr_current, dispMsg] = isToRealizeAntennaPattern()
            
            %% GET GLOBAL DIRECTORIES
            dir_ant_real = SimulationFolders.getInstance.ant_real;
            
            %% GET GLOBAL PARAMETERS
            % Simulation Settings
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
            % Simulation Parameters
            Nr = SimParams.getInstance.Nr;
            % Vegetation Parameters
            scat_cal_veg = VegParams.getInstance.scat_cal_veg;
            num_scat_cal = sum(sum(sum(scat_cal_veg))) ;
                        
            % First check the user preferences
            if ~calc_diffuse_term
                dispMsg = 'Antenna Pattern Realizations - SKIPPED (User Preferences - No Diffuse Term)';
                Nr_current = NaN;
                result = Constants.need_for_run.NO;
                return
            end
            
            % If passes user preferences, check the existence of files
            all_files = dir(dir_ant_real);
            num_files = numel(all_files) - 2;
            Nr_current = num_files / num_scat_cal / Constants.factor_ant_real ;

            if isnan(Nr_current), Nr_current = 0; end
            
            if Nr_current >= Nr
                dispMsg = 'Antenna Pattern Realizations - SKIPPED - Already exists!';
                result = Constants.need_for_run.NO;
            else
                if Nr_current > 0
                    dispMsg = 'Antenna Pattern Realizations - Partially exists!';
                    result = Constants.need_for_run.PARTIAL;
                else
                    dispMsg = 'Antenna Pattern Realizations';
                    result = Constants.need_for_run.FULL;
                end
            end
            
        end
        
        function [result, dispMsg] = isToCalcRotationMatrices()
            
            %% GET GLOBAL DIRECTORIES
            dir_rot_lookup = SimulationFolders.getInstance.rot_lookup;
            
            % If passes user preferences, check the existence of files
            all_files = dir(dir_rot_lookup);
            num_files = numel(all_files) - 2;
            
            if num_files >= Constants.num_rot_lookup
                dispMsg = 'Rotation Matrices - SKIPPED - Already exists!';
                result = Constants.need_for_run.NO;
            else
                dispMsg = 'Rotation Matrices';
                result = Constants.need_for_run.FULL;
            end
            
        end
        
        function [result, Nr_current, dispMsg] = isToRealizeRotations()
            
            %% GET GLOBAL DIRECTORIES
            dir_rot_real = SimulationFolders.getInstance.rot_real;
            
            %% GET GLOBAL PARAMETERS
            % Simulation Settings
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
            % Simulation Parameters
            Nr = SimParams.getInstance.Nr;
            % Vegetation Parameters
            scat_cal_veg = VegParams.getInstance.scat_cal_veg;
            num_scat_cal = sum(sum(sum(scat_cal_veg))) ;
                        
            % First check the user preferences
            if ~calc_diffuse_term
                dispMsg = 'Rotation Realizations - SKIPPED (User Preferences - No Diffuse Term)';
                Nr_current = NaN;
                result = Constants.need_for_run.NO;
                return
            end
            
            % If passes user preferences, check the existence of files
            all_files = dir(dir_rot_real);
            num_files = numel(all_files) - 2;
            Nr_current = num_files / num_scat_cal / Constants.factor_rot_real ;

            if isnan(Nr_current), Nr_current = 0; end
            
            if Nr_current >= Nr
                dispMsg = 'Rotation Realizations - SKIPPED - Already exists!';
                result = Constants.need_for_run.NO;
            else
                if Nr_current > 0
                    dispMsg = 'Rotation Realizations - Partially exists!';
                    result = Constants.need_for_run.PARTIAL;
                else
                    dispMsg = 'Rotation Realizations';
                    result = Constants.need_for_run.FULL;
                end
            end
            
        end
        
        function [result, dispMsg] = isToCalculateSpecularTerm()
            
            %% GET GLOBAL DIRECTORIES
            dir_out_specular_tuple = SimulationFolders.getInstance.out_specular_tuple;
            
            %% GET GLOBAL PARAMETERS
            % Simulation Settings
            calc_specular_term = SimSettings.getInstance.calc_specular_term;
                        
            % First check the user preferences
            if ~calc_specular_term
                dispMsg = 'Specular Term - SKIPPED (User Preferences - No Specular Term)';
                result = Constants.need_for_run.NO;
                return
            end            
            
            % If passes user preferences, check the existence of files           
            num_files = 0;
            
            if ( exist(dir_out_specular_tuple,'dir') == 7 )
                all_files = dir(dir_out_specular_tuple);
                num_files = numel(all_files) - 2;
            end
            
            if num_files == Constants.num_out_specular
                dispMsg = 'Specular Term - SKIPPED - Already exists!';
                result = Constants.need_for_run.NO;
            else
                dispMsg = 'Specular Term';
                result = Constants.need_for_run.FULL;
            end
            
        end
        
        function [result, dispMsg] = isToCalculateDiffuseTerm()
            
            %% GET GLOBAL DIRECTORIES
            dir_freqdiff_b1_tuple = SimulationFolders.getInstance.freqdiff_b1_tuple;
            
            %% GET GLOBAL PARAMETERS
            % Simulation Settings
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
            % Simulation Parameters
            Nr = SimParams.getInstance.Nr;
                        
            % First check the user preferences
            if ~calc_diffuse_term
                dispMsg = 'Diffuse Term - SKIPPED (User Preferences - No Diffuse Term)';
                result = Constants.need_for_run.NO;
                return
            end            
            
            % If passes user preferences, check the existence of files            
            Nr_current = 0;
            
            if ( exist(dir_freqdiff_b1_tuple,'dir') == 7 )
                all_files = dir(dir_freqdiff_b1_tuple);
                num_files = numel(all_files) - 2;
                Nr_current = num_files / Constants.factor_frediff_b1;
            end
            
            if Nr_current >= Nr
                dispMsg = 'Diffuse Term - SKIPPED - Already exists!';
                result = Constants.need_for_run.NO;
            else
                dispMsg = 'Diffuse Term';
                result = Constants.need_for_run.FULL;
            end
            
        end
        
        function saveSimParams()
            
            %% GET GLOBAL DIRECTORIES
            dir_hr = SimulationFolders.getInstance.hr;
            
            
            %% GET GLOBAL PARAMETERS
            % Simulation Parameters
            version = SimParams.getInstance.version;
            Nfz = SimParams.getInstance.Nfz;
            veg_vir_orientation_id = SimParams.getInstance.veg_vir_orientation_id;
            veg_method_id = SimParams.getInstance.veg_method_id;
            % Ground Parameters
            sand_ratio = GndParams.getInstance.sand_ratio;
            clay_ratio = GndParams.getInstance.clay_ratio;
            rhob_gcm3 = GndParams.getInstance.rhob_gcm3;
            % Transmitter Parameters
            f_MHz = TxParams.getInstance.f_MHz;
            r_Tx_m = TxParams.getInstance.r_Tx_m;
            EIRP_dB = TxParams.getInstance.EIRP_dB;
                    
            
            keySet = {ConstantNames.version, ...
                      ConstantNames.sim_Nfz, ...
                      ConstantNames.gnd_sand_ratio, ...
                      ConstantNames.gnd_clay_ratio, ...
                      ConstantNames.gnd_rhob_gcm3, ...
                      ConstantNames.Tx_f_MHz, ...
                      ConstantNames.r_Tx_m, ...
                      ConstantNames.Tx_EIRP_dB};
                  
            valueSet = {version, ...
                        Nfz, ...
                        sand_ratio, ...
                        clay_ratio, ...
                        rhob_gcm3, ...
                        f_MHz, ...
                        r_Tx_m, ...
                        EIRP_dB};
                    
            %% Put vegetation parameters to the map
            % If vegetation method is homogenous
            if veg_method_id == Constants.id_veg_hom
            
                keySet{end+1} = ConstantNames.veg_vegetationStage;
                keySet{end+1} = ConstantNames.veg_hom_TYPES;
                keySet{end+1} = ConstantNames.veg_hom_dimLayers_m;
                keySet{end+1} = ConstantNames.veg_hom_TYPKND;
                keySet{end+1} = ConstantNames.veg_hom_scatCalVeg;
                keySet{end+1} = ConstantNames.veg_hom_LTK;
                keySet{end+1} = ConstantNames.veg_hom_dsty;
                keySet{end+1} = ConstantNames.veg_hom_dim1_m;
                keySet{end+1} = ConstantNames.veg_hom_dim2_m;
                keySet{end+1} = ConstantNames.veg_hom_dim3_m;
                keySet{end+1} = ConstantNames.veg_hom_epsr;
                keySet{end+1} = ConstantNames.veg_hom_parm1_deg;
                keySet{end+1} = ConstantNames.veg_hom_parm2_deg;
                
                
                %% GET GLOBAL PARAMETERS
                % Vegetation Parameters
                % Assign to map values
                valueSet{end+1} = VegParams.getInstance.vegetation_stage;
                valueSet{end+1} = VegParams.getInstance.TYPES;
                valueSet{end+1} = VegParams.getInstance.dim_layers_m;
                valueSet{end+1} = VegParams.getInstance.TYPKND;
                valueSet{end+1} = VegParams.getInstance.scat_cal_veg;
                valueSet{end+1} = VegParams.getInstance.LTK;
                valueSet{end+1} = VegParams.getInstance.dsty;
                valueSet{end+1} = VegParams.getInstance.dim1_m;
                valueSet{end+1} = VegParams.getInstance.dim2_m;
                valueSet{end+1} = VegParams.getInstance.dim3_m;
                valueSet{end+1} = VegParams.getInstance.epsr;
                valueSet{end+1} = VegParams.getInstance.parm1_deg;
                valueSet{end+1} = VegParams.getInstance.parm2_deg;
            
            % If vegetation method is virtual
            else
                % If veg_method is Virtual and veg_vir_orientation is
                % Row-crop
                if veg_vir_orientation_id == Constants.id_veg_vir_row_crop
                    keySet{end+1} = ConstantNames.veg_vegetationStage;
                    keySet{end+1} = ConstantNames.veg_vir_row_plugin;
                    keySet{end+1} = ConstantNames.veg_vir_row_rowSpace_m;
                    keySet{end+1} = ConstantNames.veg_vir_row_colSpace_m;
                    keySet{end+1} = ConstantNames.veg_vir_row_phiRow_deg;
                    keySet{end+1} = ConstantNames.veg_vir_row_seedFluctuation_m;


                    %% GET GLOBAL PARAMETERS
                    % Virtual Row-Structured Vegetation Parameters
                    % Assign to map values
                    valueSet{end+1} = VegParams.getInstance.vegetation_stage;
                    valueSet{end+1} = VegVirRowParams.getInstance.plugin;
                    valueSet{end+1} = VegVirRowParams.getInstance.row_space_m;
                    valueSet{end+1} = VegVirRowParams.getInstance.col_space_m;
                    valueSet{end+1} = VegVirRowParams.getInstance.phi_row_deg;
                    valueSet{end+1} = VegVirRowParams.getInstance.seed_fluctuation_m;
                
                % TO-DO: If veg_method is Virtual, and veg_vir_orientation
                % is Random_spread
                elseif veg_vir_orientation_id == Constants.id_veg_vir_random_spread
                    
                    % Implement for virtual random vegetation
                    
                end
            end
                    
            input_params = containers.Map(keySet,valueSet);

            save( strcat( dir_hr, '\', ConstantNames.input_params_filename), ConstantNames.input_params_var);
            
        end
    end
    
end

