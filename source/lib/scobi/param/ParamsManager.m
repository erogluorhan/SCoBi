classdef ParamsManager < handle
    %% PARAMSMANAGER CLASS - Manage the parameters through the whole sim.
    % Static class.
    %
    % Following parameters are managed as persistent variables within
    % methods, but they act as global simulation parameters
    %
    % num_sims:    Total number of simulations (based on sim_mode)
    % sim_counter: Current simulation index within num_sims 
    
    
    methods (Static)
          
        % num_sims: Total number of simulations (based on sim_mode)
        % This function is both a getter and a setter
        function result = num_sims( value )
           
            persistent num_sims
            
            % Set value
            if nargin ~= 0
                num_sims = value;
            end
            
            % Get value
            result = num_sims;
            
        end
        
                
        function result = sim_counter( value )
           
            persistent sim_counter
            
            % Set function
            if nargin ~= 0
                
                %% GET GLOBAL PARAMETERS
                % Configuration Parameters
                th0_Tx_list_deg = ConfigParams.getInstance.th0_Tx_list_deg;
                ph0_Tx_list_deg = ConfigParams.getInstance.ph0_Tx_list_deg;
                % Receiver Parameters
                orientation_Rx_id = RxParams.getInstance.orientation_Rx_id;
                
                sim_counter = value;
                
                th0_Tx_deg = th0_Tx_list_deg( sim_counter );
                ph0_Tx_deg = ph0_Tx_list_deg( sim_counter );
                
                % If receiver orientation is Specular-facing, then set
                % receiver th0_Rx_deg equal to transmitter th0_Tx_deg
                if orientation_Rx_id == Constants.id_Rx_specular_facing
                    
                    RxParams.getInstance.set_th0_Rx_deg( th0_Tx_deg );
                    RxParams.getInstance.set_ph0_Rx_deg( ph0_Tx_deg );
                
                end
            
            end
            % Get function
            
            result = sim_counter;
            
        end
        
        
        function copyInputFiles( varin )
            
            
            %% GET GLOBAL DIRECTORIES
            dir_sim_input_used_files = SimulationFolders.getInstance.sim_input_used_files;
            
            if ~isempty( varin{1,1} )
                
                copyfile( varin{1,1}, dir_sim_input_used_files );
                
            end
            
            if ~isempty( varin{2,1} )
                
                copyfile( varin{2,1}, dir_sim_input_used_files );
                
            end
            
            if ~isempty( varin{3,1} )
                
                copyfile( varin{3,1}, dir_sim_input_used_files );
                
            end
            
        end
        
        
        function sim_counter_start = determineSimCounterStart


            %% GET GLOBAL DIRECTORIES
            dir_products_specular_reflectivity = SimulationFolders.getInstance.products_specular_reflectivity;
            
                                   
            % First read the existing specular variable, if any    
            % Choose the last bare-soil power item that is stored in
            % specularTerm (to make sure the other direct and specular 
            % items should have been stored if this is stored)
            filename02 = 'Bare02';
            currentVar = readVar(dir_products_specular_reflectivity, filename02 );

            % If no current variable, write var as initial
            if isnan( currentVar )

                sim_counter_start = 1;

            % Else if current variable exists, append var to the end
            else
                % First read the current variable's size
                [~, M] = size(currentVar);

                % Jump to the next index for starting the sim_counter
                sim_counter_start = M + 1;

            end
            
        end
        
        
        function inputParamsStruct = generateInputParamsStruct

            
        %% GET GLOBAL PARAMETERS
        % Simulation Settings
        gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
        
            
        %% SIMULATION SETTINGS        
        simSettings = SimSettings.getInstance;
        
        inputParamsStruct.version = simSettings.version;
        inputParamsStruct.campaign = simSettings.campaign;
        inputParamsStruct.sim_name = simSettings.sim_name;
        inputParamsStruct.simulator_id = simSettings.simulator_id;
        inputParamsStruct.sim_mode_id = simSettings.sim_mode_id;
        inputParamsStruct.gnd_cover_id = simSettings.gnd_cover_id;
        
        if gnd_cover_id == Constants.id_veg_cover
            
            inputParamsStruct.write_attenuation = simSettings.write_attenuation;
            
        end
        
        inputParamsStruct.include_in_master_sim_file = simSettings.include_in_master_sim_file;
        inputParamsStruct.draw_live_plots = simSettings.draw_live_plots;
        
            
        %% TRASNMITTER PARAMETERS
        txParams = TxParams.getInstance;
        
        inputParamsStruct.f_MHz = txParams.f_MHz;
        inputParamsStruct.r_Tx_m = txParams.r_Tx_m;
        inputParamsStruct.EIRP_dB = txParams.EIRP_dB;
        inputParamsStruct.pol_Tx = txParams.pol_Tx;
        
            
        %% RECEIVER PARAMETERS
        rxParams = RxParams.getInstance;
        
        inputParamsStruct.hr_m = rxParams.hr_m;
        inputParamsStruct.G0r_dB = rxParams.G0r_dB;
        inputParamsStruct.pol_Rx = rxParams.pol_Rx;
        inputParamsStruct.ant_pat_Rx_id = rxParams.ant_pat_Rx_id;   
        inputParamsStruct.ant_pat_struct_Rx = rxParams.ant_pat_struct_Rx;
        inputParamsStruct.ant_pat_res_deg = rxParams.ant_pat_res_deg;
        inputParamsStruct.orientation_Rx_id = rxParams.orientation_Rx_id;
        
        if rxParams.orientation_Rx_id == Constants.id_Rx_fixed
            
            inputParamsStruct.th0_Rx_deg = rxParams.th0_Rx_deg;
            inputParamsStruct.ph0_Rx_deg = rxParams.ph0_Rx_deg;  
            
        end

        % If receiver antenna pattern is Generalized-Gaussian
        if rxParams.ant_pat_Rx_id == Constants.id_Rx_GG

            rxGGParams = RxGGParams.getInstance;
            
            inputParamsStruct.hpbw_deg = rxGGParams.hpbw_deg;
            inputParamsStruct.SLL_dB = rxGGParams.SLL_dB;
            inputParamsStruct.XPL_dB = rxGGParams.XPL_dB;

        elseif rxParams.ant_pat_Rx_id == Constants.id_Rx_user_defined

            rxUserDefinedParams = RxUserDefinedParams.getInstance;
            
            inputParamsStruct.ant_pat_fullfilename = rxUserDefinedParams.ant_pat_fullfilename;

        % Else if Antenna pattern is Cosine to the power n
        elseif rxParams.ant_pat_Rx_id == Constants.id_Rx_cos_pow_n 

            % Should be implemented when this pattern added

        end

        
        %% GROUND PARAMETERS
        gndParams = GndParams.getInstance;
        
        inputParamsStruct.sand_ratio = gndParams.sand_ratio;
        inputParamsStruct.clay_ratio = gndParams.clay_ratio;
        inputParamsStruct.rhob_gcm3 = gndParams.rhob_gcm3;
        inputParamsStruct.diel_model_id = gndParams.diel_model_id;
        inputParamsStruct.gnd_structure_id = gndParams.gnd_structure_id;
                    
        % If ground is multi-layered
        if gndParams.gnd_structure_id == Constants.id_gnd_multi_layered

            gndMLParams = GndMLParams.getInstance;
        
            inputParamsStruct.layer_depth_m = gndMLParams.layer_depth_m;
            inputParamsStruct.delZ_m = gndMLParams.delZ_m;
            inputParamsStruct.zA_m = gndMLParams.zA_m;
            inputParamsStruct.zB_m = gndMLParams.zB_m;
            inputParamsStruct.calc_diel_profile_fit_functions = gndMLParams.calc_diel_profile_fit_functions;

        end

        
        %% CONFIGURATION PARAMETERS
        configParams = ConfigParams.getInstance;
        
        if simSettings.sim_mode_id == Constants.id_time_series
        
            inputParamsStruct.DoYs = configParams.DoYs;
            
        end
        
        inputParamsStruct.th0_Tx_list_deg = configParams.th0_Tx_list_deg;
        inputParamsStruct.ph0_Tx_list_deg = configParams.ph0_Tx_list_deg;
        inputParamsStruct.VSM_list_cm3cm3 = configParams.VSM_list_cm3cm3;
        inputParamsStruct.RMSH_list_cm = configParams.RMSH_list_cm;

        
        %% VEGETATION PARAMETERS
        % If ground cover is Vegetation, then Vegetation Parameters are set
        if gnd_cover_id == Constants.id_veg_cover
            
            vegParams = VegParams.getInstance;

            inputParamsStruct.dim_layers_m = vegParams.dim_layers_m;
            inputParamsStruct.num_layers = vegParams.num_layers;
            inputParamsStruct.TYPKND = vegParams.TYPKND;
            inputParamsStruct.num_types = vegParams.num_types;
            inputParamsStruct.num_kinds = vegParams.num_kinds;
            inputParamsStruct.LTK = vegParams.LTK;
            inputParamsStruct.dsty = vegParams.dsty;
            inputParamsStruct.dim1_m = vegParams.dim1_m;
            inputParamsStruct.dim2_m = vegParams.dim2_m;
            inputParamsStruct.dim3_m = vegParams.dim3_m;
            inputParamsStruct.epsr = vegParams.epsr;
            inputParamsStruct.parm1_deg = vegParams.parm1_deg;
            inputParamsStruct.parm2_deg = vegParams.parm2_deg;

        end
            
        end  
        
        
        function isValid = isConfigInputsValid
            
            %% GET GLOBAL PARAMETERS
            % Simulation 
            sim_mode_id = SimSettings.getInstance.sim_mode_id;
            % Configuration Parameters
            DoYs = ConfigParams.getInstance.DoYs;
            num_DoY = length( DoYs );
            th0_Tx_list_deg = ConfigParams.getInstance.th0_Tx_list_deg;
            num_Th = length( th0_Tx_list_deg );
            ph0_Tx_list_deg = ConfigParams.getInstance.ph0_Tx_list_deg;
            num_Ph = length( ph0_Tx_list_deg );
            VSM_list_cm3cm3 = ConfigParams.getInstance.VSM_list_cm3cm3;
            num_VSM = length( VSM_list_cm3cm3 );
            RMSH_list_cm = ConfigParams.getInstance.RMSH_list_cm;
            num_RMSH = length( RMSH_list_cm );
            
            
            isValid = 0;
                
            if sim_mode_id == Constants.id_snapshot
                
                if ~( num_Th == 0 || num_Ph == 0 || num_VSM == 0 || ...
                        num_RMSH == 0 )

                    isValid = 1;
                    
                end
                
            elseif sim_mode_id == Constants.id_time_series

                if ~( num_DoY == 0 || num_Th == 0 || num_Ph == 0 || num_VSM == 0 || ...
                        num_RMSH == 0 )

                    isValid = 1;
                    
                end
                
            end
            
        end
        
        
        function [validityResult, terminateResult, terminateMsg] = isInputValid
                   
            
            %% GET GLOBAL DIRECTORIES           
            dir_sim_input = SimulationFolders.getInstance.sim_input; 
            
            
            %% DEFINE PERSISTENT VARIABLES
            persistent isValid 
            persistent isTerminate
            persistent msg
            
            isValid = 0;
            isTerminate = 0;
            

            full_file_name = strcat( dir_sim_input, '\', ConstantNames.inputParamsStruct_filename);
            
            % Get whether configuration inputs are valid for this simulation
            isConfigInputsValid = ParamsManager.isConfigInputsValid();
            
            % First check if a simulation exists with the same name and it contains inputParamsStruct
            if exist( strcat( dir_sim_input, '\', ConstantNames.inputParamsStruct_filename), 'file') == 2 
            
                % Load the previously saved inputParamsStruct
                load( full_file_name, ConstantNames.inputParamsStruct);
                
                % Generate inputParamsStruct for the current simulation
                newInputParamsStruct = ParamsManager.generateInputParamsStruct();
                
                % Check if current and the latest inputParamsStruct are
                % equal
                isSysInputEqual = isequaln( newInputParamsStruct, inputParamsStruct );
                
                
                % Simulation validity is based on input equality (between 
                % current and latest) and and configuration input validity 
                isValid = isSysInputEqual & isConfigInputsValid;
                
                % If configuration inputs are not valid
                if ~isConfigInputsValid

                    isTerminate = 1;
                    msg = 'Configuration inputs are not valid for simulation mode. SCoBi simulator will be terminated for user check!';

                % Else if the reason is input inequality with the
                % existing simulation folders
                elseif ~isSysInputEqual

                    isTerminate = 1;
                    msg = 'Either one or more inputs conflict with an existing simulation. Create a simulation with a new name or remove the conflicting simulation. SCoBi simulator will be terminated for user check!';

                end
            
            else
                
                % The first time this simulation is being created
                if ~exist( dir_sim_input )
                
                    isValid = isConfigInputsValid;
                
                % Simulation exists, but inputStruct file is missing. 
                % That is an error. User should check that!
                else
                    isTerminate = 1;
                    msg = 'inputParamsStruct file is missing from the simulation directory; it might have been deleted. SCoBi simulator will be terminated for user check!'
                end
                
            end
            
            validityResult = isValid;
            terminateResult = isTerminate;
            terminateMsg = msg;
            
        end
        
        
        function [result, dispMsg] = isToCalculateMLReflectivities()
            
            %% GET GLOBAL DIRECTORIES
            % TO-DO: Add folder and file existence controls
            %dir_rot_real = SimulationFolders.getInstance.rot_real;
            
            %% GET GLOBAL PARAMETERS
            % Ground Parameters
            gnd_structure_id = GndParams.getInstance.gnd_structure_id;
            

            % If ground is single-layered
            if gnd_structure_id == Constants.id_gnd_single_layered
                
                dispMsg = '';
                result = Constants.need_for_run.NO;
                return
            
            % Else, it is multi-layered
            elseif gnd_structure_id == Constants.id_gnd_multi_layered
            
            dispMsg = 'Multilayer Reflectivities';
            result = Constants.need_for_run.FULL;

            
            end            
            
        end
        
        
        function [result, write_attenuation, dispMsg] = isToCalcPropagation()
            
            
            %% GET GLOBAL DIRECTORIES
            dir_afsa = SimulationFolders.getInstance.afsa;
            
            
            %% GET GLOBAL PARAMETERS
            % Simulation Settings
            gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
            write_attenuation = SimSettings.getInstance.write_attenuation;            
            
            
            % If ground cover is Bare-soil, then no need for scatterer
            % positions
            if gnd_cover_id == Constants.id_bare_soil
                
                dispMsg = 'Propagation - SKIPPED (Ground cover - Bare-soil)';
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
        
        
        % to-DO: Determine the final version
        function [result, dispMsg] = isToGenerateDielMLProfiles()
            
            
            %% GET GLOBAL PARAMETERS
            % Ground Parameters
            gnd_structure_id = GndParams.getInstance.gnd_structure_id;
            

            % If ground is single-layered, skip
            if gnd_structure_id == Constants.id_gnd_single_layered
                
                dispMsg = 'MultiLayer Dielectric Profiles - SKIPPED - Single-layered ground!';
                result = Constants.need_for_run.NO;
                return
            
            % Else, it is SCoBi-ML
            elseif gnd_structure_id == Constants.id_gnd_multi_layered

            dispMsg = 'MultiLayer Dielectric Profiles';
            result = Constants.need_for_run.FULL;
            
            end            
            
        end
        
        
        function initAllInputParamsFromInputParamsStruct( inputParamsStruct )

            
        %% SIMULATION SETTINGS
        % Simulation Settings        
        version = inputParamsStruct.version;
        sim_name = inputParamsStruct.sim_name;
        campaign = inputParamsStruct.campaign;
        simulator_id = inputParamsStruct.simulator_id;
        sim_mode_id = inputParamsStruct.sim_mode_id;
        gnd_cover_id = inputParamsStruct.gnd_cover_id;
        
        if gnd_cover_id == Constants.id_veg_cover
            
            write_attenuation = inputParamsStruct.write_attenuation;
            
        end
        
        include_in_master_sim_file = inputParamsStruct.include_in_master_sim_file;
        draw_live_plots = inputParamsStruct.draw_live_plots;
        
        % Initialize
        SimSettings.getInstance.initialize( campaign, simulator_id, sim_mode_id, gnd_cover_id, write_attenuation, ...
                include_in_master_sim_file, draw_live_plots );
        
        
        %% TRASNMITTER PARAMETERS
        f_MHz = inputParamsStruct.f_MHz;
        r_Tx_m = inputParamsStruct.r_Tx_m;
        EIRP_dB = inputParamsStruct.EIRP_dB;
        pol_Tx = inputParamsStruct.pol_Tx;
        
        % Initialize
        TxParams.getInstance.initialize( f_MHz, r_Tx_m, EIRP_dB, pol_Tx );
            
        
        %% RECEIVER PARAMETERS        
        hr_m = inputParamsStruct.hr_m;
        G0r_dB = inputParamsStruct.G0r_dB;
        pol_Rx = inputParamsStruct.pol_Rx;
        ant_pat_Rx_id = inputParamsStruct.ant_pat_Rx_id;   
        ant_pat_struct_Rx = inputParamsStruct.ant_pat_struct_Rx;
        ant_pat_res_deg = inputParamsStruct.ant_pat_res_deg;
        orientation_Rx_id = inputParamsStruct.orientation_Rx_id;
        
        th0_Rx_deg = [];
        ph0_Rx_deg = [];
        if orientation_Rx_id == Constants.id_Rx_fixed
            
            th0_Rx_deg = inputParamsStruct.th0_Rx_deg;
            ph0_Rx_deg = inputParamsStruct.ph0_Rx_deg;  
            
        end

        % If receiver antenna pattern is Generalized-Gaussian
        if ant_pat_Rx_id == Constants.id_Rx_GG
            
            hpbw_deg = inputParamsStruct.hpbw_deg;
            SLL_dB = inputParamsStruct.SLL_dB;
            XPL_dB = inputParamsStruct.XPL_dB;
            
            % Initialize
            RxGGParams.getInstance.initialize( hpbw_deg, SLL_dB, XPL_dB );

        elseif ant_pat_Rx_id == Constants.id_Rx_user_defined
            
            ant_pat_fullfilename = inputParamsStruct.ant_pat_fullfilename;
            
            % Initialize
            RxUserDefinedParams.getInstance.initialize( ant_pat_fullfilename );

        % Else if Antenna pattern is Cosine to the power n
        elseif ant_pat_Rx_id == Constants.id_Rx_cos_pow_n 

            % Should be implemented when this pattern added

        end
        
        % Initialize
        RxParams.getInstance.initialize( hr_m, G0r_dB, pol_Rx, ant_pat_Rx_id, ...
                ant_pat_struct_Rx, ant_pat_res_deg, orientation_Rx_id, ...
                th0_Rx_deg, ph0_Rx_deg );

        
        %% GROUND PARAMETERS        
        sand_ratio = inputParamsStruct.sand_ratio;
        clay_ratio = inputParamsStruct.clay_ratio;
        rhob_gcm3 = inputParamsStruct.rhob_gcm3;
        diel_model_id = inputParamsStruct.diel_model_id;
        gnd_structure_id = inputParamsStruct.gnd_structure_id;
                    
        % If ground is multi-layered
        if gnd_structure_id == Constants.id_gnd_multi_layered
            
            addpath( genpath( Directories.getInstance.multi_layer ) );
        
            layer_depth_m = inputParamsStruct.layer_depth_m;
            delZ_m = inputParamsStruct.delZ_m;
            zA_m = inputParamsStruct.zA_m;
            zB_m = inputParamsStruct.zB_m;
            calc_diel_profile_fit_functions = inputParamsStruct.calc_diel_profile_fit_functions;
            
            % Initialize
            GndMLParams.getInstance.initialize( layer_depth_m, delZ_m, zA_m, zB_m, calc_diel_profile_fit_functions );

        end
        
        % Initialize
        GndParams.getInstance.initialize( gnd_structure_id, sand_ratio, clay_ratio, rhob_gcm3, diel_model_id );

        
        %% CONFIGURATION PARAMETERS 
        DoYs = [];
        
        if sim_mode_id == Constants.id_time_series
        
            DoYs = inputParamsStruct.DoYs;
            
        end
        
        th0_Tx_list_deg = inputParamsStruct.th0_Tx_list_deg;
        ph0_Tx_list_deg = inputParamsStruct.ph0_Tx_list_deg;
        VSM_list_cm3cm3 = inputParamsStruct.VSM_list_cm3cm3;
        RMSH_list_cm = inputParamsStruct.RMSH_list_cm;

        % Initialize
        ConfigParams.getInstance.initialize( DoYs, th0_Tx_list_deg, ph0_Tx_list_deg, VSM_list_cm3cm3, RMSH_list_cm );
        
        %% VEGETATION PARAMETERS
        % If ground cover is Vegetation, then Vegetation Parameters are set
        if gnd_cover_id == Constants.id_veg_cover
            
            addpath( genpath( Directories.getInstance.vegetation ) );
            
            dim_layers_m = inputParamsStruct.dim_layers_m;
            num_layers = inputParamsStruct.num_layers;
            TYPKND = inputParamsStruct.TYPKND;
            num_types = inputParamsStruct.num_types;
            num_kinds = inputParamsStruct.num_kinds;
            LTK = inputParamsStruct.LTK;
            dsty = inputParamsStruct.dsty;
            dim1_m = inputParamsStruct.dim1_m;
            dim2_m = inputParamsStruct.dim2_m;
            dim3_m = inputParamsStruct.dim3_m;
            epsr = inputParamsStruct.epsr;
            parm1_deg = inputParamsStruct.parm1_deg;
            parm2_deg = inputParamsStruct.parm2_deg;
                        
            % Initialize
            VegParams.getInstance.initialize( dim_layers_m, num_layers, ...
                       TYPKND, num_types, num_kinds, LTK, dsty, dim1_m, ...
                       dim2_m, dim3_m, epsr, parm1_deg, parm2_deg );

        end
            
        end  
        
        
        function saveInputParams
            
            %% GET GLOBAL DIRECTORIES
            dir_sim_input = SimulationFolders.getInstance.sim_input;
            
            inputParamsStruct = ParamsManager.generateInputParamsStruct();

            full_file_name = strcat( dir_sim_input, '\', ConstantNames.inputParamsStruct_filename);
            
            save( full_file_name, ConstantNames.inputParamsStruct);
            
        end
        
        
        function writeInputReports( inputStruct )
            
            
            %% GET GLOBAL PARAMETERS
            % Simulation Settings
            include_in_master_sim_file = SimSettings.getInstance.include_in_master_sim_file;
            
            
            ParamsManager.writeInputStructToReportFile( inputStruct );
            
            % If selected to be included in the master simulation file
            if include_in_master_sim_file                
                
                ParamsManager.writeSimIntoMasterSimFile();
                
            end
            
        end
        
        
        function writeInputStructToReportFile( inputStruct )


        %% GET GLOBAL DIRECTORIES
        dir_sim_input = SimulationFolders.getInstance.sim_input;
        inputFile = strcat(dir_sim_input, '\', 'input_report.txt');


        %% GET GLOBAL PARAMETERS
        % Simulation Settings
        simulator_id = SimSettings.getInstance.simulator_id;
        gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
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
        fprintf(fileID, strcat( char(t), '\n\n' ) );


        %% SIMULATION SETTINGS
        % Write the simulation software version
        fprintf(fileID, strcat( 'SCoBi Version:\t\t\t', Constants.version, '\n' ) );
        % Write the simulation campaign
        fprintf(fileID, strcat( 'Campaign:\t\t\t', inputStruct.campaign, '\n' ) );
        % Write sim_mode
        fprintf(fileID, strcat( 'Simulation mode:\t', inputStruct.sim_mode, '\n' ) );
        % Write gnd_cover
        fprintf(fileID, strcat( 'Ground cover:\t\t', inputStruct.gnd_cover, '\n' ) );

        if gnd_cover_id == Constants.id_veg_cover

            % Write user preferences
            fprintf(fileID, strcat( 'Write attenuation?:\t\t', num2str(inputStruct.write_attenuation), '\n' ) );
            
        end
        
        fprintf(fileID, strcat( 'Include in Master Simulation File?:\t\t', num2str(inputStruct.include_in_master_sim_file), '\n' ) );

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
        
        filename = strrep(inputStruct.config_inputs_file, '\', '/');
        fprintf(fileID, strcat( 'Ground input file:\t', filename, '\n' ) );

        %% VEGETATION PARAMETERS
        if gnd_cover_id == Constants.id_veg_cover
            
            fprintf(fileID, strcat( '\n\nVEGETATION PARAMETERS\n' ) );
            filename = strrep(inputStruct.veg_inputs_file, '\', '/');
            fprintf(fileID, strcat( 'Vegetation input file:\t', filename, '\n' ) );
            
        end


        %% CONFIGURATION PARAMETERS
        fprintf(fileID, strcat( '\n\nCONFIGURATION PARAMETERS\n' ) );
        filename = strrep(inputStruct.config_inputs_file, '\', '/');
        fprintf(fileID, strcat( 'Configuration input file:\t', filename, '\n' ) );


        fclose(fileID);

        winopen( inputFile );

        end
        
        
        function writeSimIntoMasterSimFile


        %% GET GLOBAL DIRECTORIES
        dir_sims_main_dir = SimulationFolders.getInstance.sims_main_dir;
        masterSimInputFullFile = strcat(dir_sims_main_dir, '\', ConstantNames.master_sim_input_filename );

            
        %% GET GLOBAL PARAMETERS
        % Simulation Settings
        sim_name = SimSettings.getInstance.sim_name;
        gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
        gnd_cover = Constants.gnd_covers{ 1, gnd_cover_id };
        sim_mode_id = SimSettings.getInstance.sim_mode_id;
        sim_mode = Constants.sim_modes{ 1, sim_mode_id };
        
        % Transmitter Parameters
        f_MHz = TxParams.getInstance.f_MHz;
        pol_Tx = TxParams.getInstance.pol_Tx;
        % Receiver Parameters
        hr_m = RxParams.getInstance.hr_m;
        pol_Rx = RxParams.getInstance.pol_Rx;
        % Ground Parameters
        gnd_structure_id = GndParams.getInstance.gnd_structure_id;
        gnd_structure = Constants.gnd_structures{ 1, gnd_structure_id };
        diel_model_id = GndParams.getInstance.diel_model_id;
        diel_model = Constants.diel_models{ 1, diel_model_id };
        
        
        varNames = {'SimName', 'GndStructure', 'SimMode', 'GroundCover', 'TxFreq_MHz', 'TxPol', 'RxAltitude_m', 'RxPol', 'DielModel'};
                
        % Open master inputs file to write this sim into it
        % If this is the first time master inputs file is created
        if ~exist( masterSimInputFullFile, 'file' )
            
            T = table( {char(sim_name)}, {char(gnd_structure)}, {char(sim_mode)}, {char(gnd_cover)}, f_MHz, {char(pol_Tx)}, hr_m, {char(pol_Rx)}, {char(diel_model)}, 'VariableNames', varNames );

        % Else, the master input file exists
        else
            
            T = readtable( masterSimInputFullFile );
            
            if isempty( find(strcmp(sim_name, table2cell(T))) )
                
                [numRows, ~] = size(T);
                
                newRow = table( {char(sim_name)}, {char(gnd_structure)}, {char(sim_mode)}, {char(gnd_cover)}, f_MHz, {char(pol_Tx)}, hr_m, {char(pol_Rx)}, {char(diel_model)}, 'VariableNames', varNames  );
            
                T = [ T; newRow ];
                
            end
        
        end

        writetable( T, masterSimInputFullFile);
            
        end  

end
    
end

