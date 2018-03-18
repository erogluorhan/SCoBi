classdef ParamsManager
    %% PARAMSMANAGER CLASS - Check the validity of input parameters
    %
    
    
    methods (Static)
                
        function result = index_Th( value )
           
            persistent index_Th
            
            if nargin ~= 0
                index_Th = value;
            end
            
            result = index_Th;
            
        end
                
        function result = index_Ph( value )
           
            persistent index_Ph
            
            if nargin ~= 0
                index_Ph = value;
            end
            
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
            
            persistent isValid 
            persistent isTerminate
            persistent msg
            
            isValid = 0;
            isTerminate = 0;
            
            if exist( strcat( SimulationFolders.getInstance.hr, '\', ConstantNames.input_params_filename), 'file') == 2 
                        
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
                        defaultans = {SimParams.getInstance.sim_name};
                        new_sim_name = char( inputdlg(prompt,dlg_title,num_lines,defaultans) );  

                        if ~isempty(new_sim_name)
                            SimParams.getInstance.updateSimName(new_sim_name);

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
                if ~exist( SimulationFolders.getInstance.hr)
                
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
            
            isValid = 0;
                
            num_Th = length( SatParams.getInstance.th0_deg );
            num_Ph = length( SatParams.getInstance.PH0_deg );
            num_VSM = length( GndParams.getInstance.VSM );
            num_RMSH = length( GndParams.getInstance.RMSH );
            
            maxnum = max( [num_Th, num_Ph, num_VSM, num_RMSH] );
            
            if SimSettings.getInstance.sim_mode == Constants.sim_mode.TIME_SERIES
                
                % Add a control for timestamps as well
                if (num_Th == maxnum || num_Th == 1) && ...
                   (num_Ph == maxnum || num_Ph == 1) && ...
                   (num_VSM == maxnum || num_VSM == 1 ) && ...
                   (num_RMSH == maxnum || num_RMSH == 1)
                    
                    isValid = 1;
                    
                    if num_Th == 1
                        SatParams.getInstance.rep_Th(maxnum);
                    end
                    
                    if num_Ph == 1
                        SatParams.getInstance.rep_Ph(maxnum);
                    end
                    
                    if num_VSM == 1
                        GndParams.getInstance.rep_VSM(maxnum);
                    end                    
                    
                    if num_RMSH == 1
                        GndParams.getInstance.rep_RMSH(maxnum);
                    end

                end
                
            else
                
                if ~( num_Th == 0 || num_Ph == 0 || num_VSM == 0 || ...
                        num_RMSH == 0 )
                
                    isValid = 1;
                
                end
                
            end
            
        end
        
        
        function isEqual = isSysInputEqual
            
            isEqual = 0;
                       
            load( strcat( SimulationFolders.getInstance.hr, '\', ConstantNames.input_params_filename), ConstantNames.input_params_var );
            
            sim_Nfz = SimParams.getInstance.Nfz;
            gnd_sand = GndParams.getInstance.sand;
            gnd_clay = GndParams.getInstance.clay;
            gnd_rhoB = GndParams.getInstance.rho_b;
            sat_fMHz = SatParams.getInstance.fMHz;
            sat_rsatKm = SatParams.getInstance.rsat;
            sat_EIRPdB = SatParams.getInstance.EIRP_dB;

            if input_params(ConstantNames.sim_Nfz) == sim_Nfz
                if input_params(ConstantNames.gnd_sand) == gnd_sand
                    if input_params(ConstantNames.gnd_clay) == gnd_clay
                        if input_params(ConstantNames.gnd_rhoB) == gnd_rhoB
                            if input_params(ConstantNames.sat_fMHz) == sat_fMHz
                                if input_params(ConstantNames.sat_rsatKm) == sat_rsatKm
                                    if input_params(ConstantNames.sat_EIRPdB) == sat_EIRPdB
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
            
            isEqual = 0;
                       
            load( strcat( SimulationFolders.getInstance.hr, '\', ConstantNames.input_params_filename), ConstantNames.input_params_var );
            
            if SimParams.getInstance.vegetation_method == Constants.veg_methods.HOMOGENOUS
            
                vegetation_stage = VegParams.getInstance.vegetation_stage;
                TYPES = VegParams.getInstance.TYPES;
                dim_layers = VegParams.getInstance.dim_layers;
                TYPKND = VegParams.getInstance.TYPKND;
                scat_cal_veg = VegParams.getInstance.scat_cal_veg;
                LTK = VegParams.getInstance.LTK;
                dsty = VegParams.getInstance.dsty;
                dim1 = VegParams.getInstance.dim1;
                dim2 = VegParams.getInstance.dim2;
                dim3 = VegParams.getInstance.dim3;
                epsr = VegParams.getInstance.epsr;
                parm1 = VegParams.getInstance.parm1;
                parm2 = VegParams.getInstance.parm2;

                if strcmp( input_params(ConstantNames.veg_hom_vegetationStage), vegetation_stage )
                    if isequal( input_params(ConstantNames.veg_hom_TYPES), TYPES )
                        if isequal( input_params(ConstantNames.veg_hom_dimLayers), dim_layers )
                            if isequal( input_params(ConstantNames.veg_hom_TYPKND), TYPKND )                                          
                                if isequal( input_params(ConstantNames.veg_hom_scatCalVeg), scat_cal_veg )
                                    if isequal( input_params(ConstantNames.veg_hom_LTK), LTK )
                                        if isequal( input_params(ConstantNames.veg_hom_dsty), dsty )
                                            if isequal( input_params(ConstantNames.veg_hom_dim1), dim1 )
                                                if isequal( input_params(ConstantNames.veg_hom_dim2), dim2 )
                                                    if isequal( input_params(ConstantNames.veg_hom_dim3), dim3 )
                                                        if isequal( input_params(ConstantNames.veg_hom_epsr), epsr )
                                                            if isequal( input_params(ConstantNames.veg_hom_parm1), parm1 )
                                                                if isequal( input_params(ConstantNames.veg_hom_parm2), parm2 )
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
                
            else
                
                plugin = VegVirRowParams.getInstance.plugin;
                row_space = VegVirRowParams.getInstance.row_space;
                col_space = VegVirRowParams.getInstance.col_space;
                phi_row = VegVirRowParams.getInstance.phi_row;
                plant_row_spread = VegVirRowParams.getInstance.plant_row_spread;
                plant_col_spread = VegVirRowParams.getInstance.plant_col_spread;

                if input_params(ConstantNames.veg_vir_row_plugin).isTheSame( plugin )
                    if input_params(ConstantNames.veg_vir_row_rowSpace) == row_space
                        if input_params(ConstantNames.veg_vir_row_colSpace) == col_space
                            if input_params(ConstantNames.veg_vir_row_phiRow) == phi_row
                                if input_params(ConstantNames.veg_vir_row_plantRowSpread) == plant_row_spread
                                    if input_params(ConstantNames.veg_vir_row_plantColSpread) == plant_col_spread
                                        isEqual = 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
        end
        
        function [result, Nr_current, dispMsg] = isToGenerateScatPos()
            
            calc_specular_term = SimSettings.getInstance.calc_specular_term;
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
            
            Nr = SimParams.getInstance.Nr;
            vegetation_method = SimParams.getInstance.vegetation_method;
            
            % First check the user preferences
            if vegetation_method == Constants.veg_methods.HOMOGENOUS
                if ~calc_diffuse_term
                    dispMsg = 'Generate Scatterer Positions - SKIPPED (User Preferences - No Diffuse Term)';
                    Nr_current = NaN;
                    result = Constants.need_for_run.NO;
                    return
                end
            
            else vegetation_method == Constants.veg_methods.VIRTUAL
                if ~calc_diffuse_term && ~calc_specular_term
                    dispMsg = 'Generate Scatterer Positions - SKIPPED (User Preferences - No Specular and Diffuse Term)';
                    Nr_current = NaN;
                    result = Constants.need_for_run.NO;
                    return
                end
            end
            
            
            readExistingVegParams;
            
            scat_cal_veg = VegParams.getInstance.scat_cal_veg;
            num_scat_cal = sum(sum(sum(scat_cal_veg))) ;
            
            % If passes user preferences, check the existence of files
            all_files = dir(SimulationFolders.getInstance.position);
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
        
        function [result, Nr_current, dispMsg] = isToCalculateFwdScatAmp()
            
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
            
            Nr = SimParams.getInstance.Nr;
            
            scat_cal_veg = VegParams.getInstance.scat_cal_veg;
            num_scat_cal = sum(sum(sum(scat_cal_veg))) ;
                        
            % First check the user preferences
            if ~calc_diffuse_term
                dispMsg = 'Forward Scatterer Amplitudes - SKIPPED (User Preferences - No Diffuse Term)';
                Nr_current = NaN;
                result = Constants.need_for_run.NO;
                return
            end
            
            % If passes user preferences, check the existence of files
            all_files = dir(SimulationFolders.getInstance.fscat);
            num_files = numel(all_files) - 2;
            Nr_current = num_files / num_scat_cal / Constants.factor_fscat ;

            if isnan(Nr_current), Nr_current = 0; end
            
            if Nr_current >= Nr
                dispMsg = 'Forward Scatterer Amplitudes - SKIPPED - Already exists!';
                result = Constants.need_for_run.NO;
            else
                if Nr_current > 0
                    dispMsg = 'Forward Scatterer Amplitudes - Partially exists!';
                    result = Constants.need_for_run.PARTIAL;
                else
                    dispMsg = 'Forward Scatterer Amplitudes';
                    result = Constants.need_for_run.FULL;
                end
            end
            
        end
        
        function [result, dispMsg] = isToCalcPropagation()
            
            calc_specular_term = SimSettings.getInstance.calc_specular_term;
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
                        
            % First check the user preferences
            if ~( calc_specular_term || calc_diffuse_term )
                dispMsg = 'Propagation - SKIPPED (User Preferences - No Specular and Diffuse Term)';
                result = Constants.need_for_run.NO;
                return
            end
            
            % If passes user preferences, check the existence of files
            all_files = dir(SimulationFolders.getInstance.afsa);
            num_files = numel(all_files) - 2;
            
            if num_files >= Constants.num_afsa
                dispMsg = 'Propagation - SKIPPED - Already exists!';
                result = Constants.need_for_run.NO;
            else
                dispMsg = 'Propagation';
                result = Constants.need_for_run.FULL;
            end
            
        end
        
        function [result, dispMsg] = isToAntennaPatternMatrix()
                       
            % If passes user preferences, check the existence of files
            all_files = dir(SimulationFolders.getInstance.ant_lookup);
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
            
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
            
            Nr = SimParams.getInstance.Nr;
            
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
            all_files = dir(SimulationFolders.getInstance.ant_real);
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
                       
            % If passes user preferences, check the existence of files
            all_files = dir(SimulationFolders.getInstance.rot_lookup);
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
            
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
            
            Nr = SimParams.getInstance.Nr;
            
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
            all_files = dir(SimulationFolders.getInstance.rot_real);
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
            
            calc_specular_term = SimSettings.getInstance.calc_specular_term;
                        
            % First check the user preferences
            if ~calc_specular_term
                dispMsg = 'Specular Term - SKIPPED (User Preferences - No Specular Term)';
                result = Constants.need_for_run.NO;
                return
            end            
            
            % If passes user preferences, check the existence of files
            outName = strcat(SimulationFolders.getInstance.out_specular, '\VSM_', num2str( GndParams.getInstance.VSM( ParamsManager.index_VSM ) ), '-RMSH_', num2str(GndParams.getInstance.RMSH( ParamsManager.index_RMSH )) );
            
            num_files = 0;
            
            if ( exist(outName,'dir') == 7 )
                all_files = dir(outName);
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
            
            calc_diffuse_term = SimSettings.getInstance.calc_diffuse_term;
            
            Nr = SimParams.getInstance.Nr;
                        
            % First check the user preferences
            if ~calc_diffuse_term
                dispMsg = 'Diffuse Term - SKIPPED (User Preferences - No Diffuse Term)';
                result = Constants.need_for_run.NO;
                return
            end            
            
            % If passes user preferences, check the existence of files
            outName = strcat(SimulationFolders.getInstance.freqdiff_b1, '\', 'VSM_', num2str( GndParams.getInstance.VSM( ParamsManager.index_VSM ) ), '-RMSH_', num2str(GndParams.getInstance.RMSH( ParamsManager.index_RMSH )));
            
            Nr_current = 0;
            
            if ( exist(outName,'dir') == 7 )
                all_files = dir(outName);
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
            
            keySet = {ConstantNames.version, ...
                      ConstantNames.sim_Nfz, ...
                      ConstantNames.gnd_sand, ...
                      ConstantNames.gnd_clay, ...
                      ConstantNames.gnd_rhoB, ...
                      ConstantNames.sat_fMHz, ...
                      ConstantNames.sat_rsatKm, ...
                      ConstantNames.sat_EIRPdB};
                  
            valueSet = {SimParams.getInstance.version, ...
                        SimParams.getInstance.Nfz, ...
                        GndParams.getInstance.sand, ...
                        GndParams.getInstance.clay, ...
                        GndParams.getInstance.rho_b, ...
                        SatParams.getInstance.fMHz, ...
                        SatParams.getInstance.rsat, ...
                        SatParams.getInstance.EIRP_dB};
                    
            % Add vegetation parameters        
            if SimParams.getInstance.vegetation_method == Constants.veg_methods.HOMOGENOUS
            
                keySet{end+1} = ConstantNames.veg_hom_vegetationStage;
                keySet{end+1} = ConstantNames.veg_hom_TYPES;
                keySet{end+1} = ConstantNames.veg_hom_dimLayers;
                keySet{end+1} = ConstantNames.veg_hom_TYPKND;
                keySet{end+1} = ConstantNames.veg_hom_scatCalVeg;
                keySet{end+1} = ConstantNames.veg_hom_LTK;
                keySet{end+1} = ConstantNames.veg_hom_dsty;
                keySet{end+1} = ConstantNames.veg_hom_dim1;
                keySet{end+1} = ConstantNames.veg_hom_dim2;
                keySet{end+1} = ConstantNames.veg_hom_dim3;
                keySet{end+1} = ConstantNames.veg_hom_epsr;
                keySet{end+1} = ConstantNames.veg_hom_parm1;
                keySet{end+1} = ConstantNames.veg_hom_parm2;
                           
                valueSet{end+1} = VegParams.getInstance.vegetation_stage;
                valueSet{end+1} = VegParams.getInstance.TYPES;
                valueSet{end+1} = VegParams.getInstance.dim_layers;
                valueSet{end+1} = VegParams.getInstance.TYPKND;
                valueSet{end+1} = VegParams.getInstance.scat_cal_veg;
                valueSet{end+1} = VegParams.getInstance.LTK;
                valueSet{end+1} = VegParams.getInstance.dsty;
                valueSet{end+1} = VegParams.getInstance.dim1;
                valueSet{end+1} = VegParams.getInstance.dim2;
                valueSet{end+1} = VegParams.getInstance.dim3;
                valueSet{end+1} = VegParams.getInstance.epsr;
                valueSet{end+1} = VegParams.getInstance.parm1;
                valueSet{end+1} = VegParams.getInstance.parm2;
            
            % For virtual vegetetation 
            else
                keySet{end+1} = ConstantNames.veg_vir_row_plugin;
                keySet{end+1} = ConstantNames.veg_vir_row_rowSpace;
                keySet{end+1} = ConstantNames.veg_vir_row_colSpace;
                keySet{end+1} = ConstantNames.veg_vir_row_phiRow;
                keySet{end+1} = ConstantNames.veg_vir_row_plantRowSpread;
                keySet{end+1} = ConstantNames.veg_vir_row_plantColSpread;
                               
                valueSet{end+1} = VegVirRowParams.getInstance.plugin;
                valueSet{end+1} = VegVirRowParams.getInstance.row_space;
                valueSet{end+1} = VegVirRowParams.getInstance.col_space;
                valueSet{end+1} = VegVirRowParams.getInstance.phi_row;
                valueSet{end+1} = VegVirRowParams.getInstance.plant_row_spread;
                valueSet{end+1} = VegVirRowParams.getInstance.plant_col_spread;
            end
                    
            input_params = containers.Map(keySet,valueSet);

            save( strcat( SimulationFolders.getInstance.hr, '\', ConstantNames.input_params_filename), ConstantNames.input_params_var);
            
        end
    end
    
end

