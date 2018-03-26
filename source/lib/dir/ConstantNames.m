classdef ConstantNames
    %CONSTANTNAMES Global constant naming strings
    %   Keeps the global constants for the names of anything (e.g. files, 
    %   variables, xml tags, etc.) in the simulation
    
    properties (Constant)
        
        version = 'version'; % String constant for sim version
        
        input_params_var = 'input_params';
        input_params_filename = 'input_params.mat';
             
        sys_input = 'sys_input'; % String constant for system inputs
        
        veg_input = 'veg_input'; % String constant for vegetation input
        
        % String constants for Ground Inputs and Params
        gnd_VSM_cm3cm3 = 'VSM_cm3cm3';
        gnd_sand_ratio = 'sand_ratio';
        gnd_clay_ratio = 'clay_ratio';
        gnd_rhob_gcm3 = 'rhob_gcm3';
        gnd_RMSH_cm = 'RMSH_cm';
        
        % String constants for Receiver Inputs and Params
        rec_hr_m = 'hr_m';
        rec_G0r_dB = 'G0r_dB';
        rec_hpbw_deg = 'hpbw_deg';
        rec_SLL_dB = 'SLL_dB'
        rec_XPL_dB = 'XPL_dB';
        rec_polR = 'polR';
        
        % String constants for Satellite Inputs and Params
        sat_f_MHz = 'f_MHz';
        sat_rsat_m = 'rsat_m';
        sat_rsat_km = 'rsat_km';
        sat_th0_deg = 'th0_deg';
        sat_PH0_deg = 'PH0_deg';
        sat_EIRP_dB = 'EIRP_dB';
        sat_polT = 'polT';
        
        % String constants for Simulation Inputs and Params
        sim_simName = 'sim_name';
        sim_campaign = 'campaign';
        sim_campaignDate = 'campaign_date';
        sim_plot = 'plot';
        sim_vegetationMethod = 'vegetation_method';
        sim_vegetationIsRow = 'vegetation_isRow';
        sim_vegetationPlant = 'vegetation_plant';
        sim_Nr = 'Nr';
        sim_Nfz = 'Nfz';
        
        % String constants for Settings
        set_simMode = 'sim_mode';
        set_simMode_snapshot = 'snapshot';
        set_simMode_timeSeries = 'time_series';
        set_groundCover = 'ground_cover';
        set_writeAttenuation = 'write_attenuation';
        set_calcDirectTerm = 'calc_direct_term';
        set_calcSpecularTerm = 'calc_specular_term';
        set_calcDiffuseTerm = 'calc_diffuse_term';
        
        veg_hom_vegetationStage = 'vegetation_stage';
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
        veg_vir_row_plantRowSpread_m = 'plant_row_spread_m';
        veg_vir_row_plantColSpread_m = 'plant_col_spread_m';
        
    end
    
end

