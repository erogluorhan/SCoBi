classdef ConstantNames
    %CONSTANTNAMES Global constant naming strings
    %   Keeps the global constants for the names of anything (e.g. files, 
    %   variables, xml tags, etc.) in the simulation
    
    properties (Constant)
        
        version = 'version'; % String constant for sim version
        
        input_params_var = 'input_params';
        input_params_filename = 'input_params.mat';
        
        attenuation_out_filename = '\Attenuation.xls';
             
        sys_input = 'sys_input'; % String constant for system inputs
        
        veg_input = 'veg_input'; % String constant for vegetation input
        
        % TO-DO: There are changes! for every kind of inputs
        
        % String constants for Ground Inputs and Parameters
        dyn_th0_Tx_list_deg = 'th0_Tx_list_deg';
        dyn_ph0_Tx_list_deg = 'ph0_Tx_list_deg';
        dyn_VSM_list_cm3cm3 = 'VSM_list_cm3cm3';
        dyn_RMSH_list_cm = 'RMSH_list_cm';
        
        % String constants for Ground Inputs and Params
        gnd_sand_ratio = 'sand_ratio';
        gnd_clay_ratio = 'clay_ratio';
        gnd_rhob_gcm3 = 'rhob_gcm3';
        gnd_diel_model = 'diel_model';
        
        % String constants for Receiver Inputs and Params
        Rx_hr_m = 'hr_m';
        Rx_G0r_dB = 'G0r_dB';
        Rx_pol_Rx = 'pol_Rx';
        Rx_orientation_Rx = 'orientation_Rx';
        Rx_th0_Rx_deg = 'th0_Rx_deg';
        Rx_ph0_Rx_deg = 'ph0_Rx_deg';
        Rx_ant_pat_Rx = 'ant_pat_Rx';
        Rx_ant_pat_res_deg = 'ant_pat_res_deg';
        Rx_GG_hpbw_deg = 'hpbw_deg';
        Rx_GG_SLL_dB = 'SLL_dB'
        Rx_GG_XPL_dB = 'XPL_dB';
        Rx_user_def_ant_pat_file = 'ant_pat_file';
        
        % String constants for Transmitter Inputs and Params
        Tx_f_MHz = 'f_MHz';
        r_Tx_m = 'r_Tx_m';
        Tx_rsat_km = 'rsat_km';
        Tx_EIRP_dB = 'EIRP_dB';
        Tx_pol_Tx = 'pol_Tx';
        
        % String constants for Simulation Inputs and Params
        sim_simName = 'sim_name';
        sim_campaign = 'campaign';
        sim_campaignDate = 'campaign_date';
        sim_plot = 'plot';
        sim_vegMethod = 'veg_method';
        sim_vegVirOrientation = 'veg_vir_orientation';
        sim_vegetationPlant = 'vegetation_plant';
        sim_Nr = 'Nr';
        sim_Nfz = 'Nfz';
        
        % String constants for Settings
        set_simulator = 'simulator'
        set_simMode = 'sim_mode';
        set_simMode_snapshot = 'snapshot';
        set_simMode_timeSeries = 'time_series';
        set_gndCover = 'gnd_cover';
        set_writeAttenuation = 'write_attenuation';
        set_calcDirectTerm = 'calc_direct_term';
        set_calcSpecularTerm = 'calc_specular_term';
        set_calcDiffuseTerm = 'calc_diffuse_term';
        
        veg_file = 'veg_input_file';
        
        veg_vegetationStage = 'vegetation_stage';
        
        veg_hom_TYPES = 'TYPES';
        veg_hom_TYPKND = 'TYPKND';
        veg_hom_dimLayers_m = 'dim_layers_m';
        veg_hom_scatCalVeg = 'scat_cal_veg';
        veg_hom_LTK = 'LTK';
        veg_hom_dsty = 'dsty';
        veg_hom_dim1_m = 'dim1_m';
        veg_hom_dim2_m = 'dim2_m';
        veg_hom_dim3_m = 'dim3_m';
        veg_hom_epsr = 'epsr';
        veg_hom_parm1_deg = 'parm1_deg';
        veg_hom_parm2_deg = 'parm2_deg';
        veg_vir_row_plugin = 'plugin';
        veg_vir_row_rowSpace_m = 'row_space_m';
        veg_vir_row_colSpace_m = 'col_space_m';
        veg_vir_row_phiRow_deg = 'phi_row_deg';
        veg_vir_row_seedFluctuation_m = 'seed_fluctuation_m';
        
        
        dyn_input_file = 'dyn_input_file';
        
    end
    
end

