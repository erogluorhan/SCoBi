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
        gnd_VSM = 'VSM';
        gnd_sand = 'sand';
        gnd_clay = 'clay';
        gnd_rhoB = 'rho_b';
        gnd_RMSH = 'RMSH';
        
        % String constants for Receiver Inputs and Params
        rec_hr = 'hr';
        rec_G0rdB = 'G0r_dB';
        rec_hpbwDeg = 'hpbw_deg';
        rec_SLLdB = 'SLL_dB'
        rec_XPLdB = 'XPL_dB';
        rec_polR = 'polR';
        
        % String constants for Satellite Inputs and Params
        sat_fMHz = 'fMHz';
        sat_rsatKm = 'rsat_km';
        sat_th0Deg = 'th0_deg';
        sat_PH0Deg = 'PH0_deg';
        sat_EIRPdB = 'EIRP_dB';
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
        set_simMode_cross = 'cross_sim';
        set_simMode_timeSeries = 'time_series';
        set_groundCover = 'ground_cover';
        set_calcMetaData = 'calc_meta_data';
        set_calcDirectTerm = 'calc_direct_term';
        set_calcSpecularTerm = 'calc_specular_term';
        set_calcDiffuseTerm = 'calc_diffuse_term';
        
        veg_hom_vegetationStage = 'vegetation_stage';
        veg_hom_TYPES = 'TYPES';
        veg_hom_TYPKND = 'TYPKND';
        veg_hom_dimLayers = 'dim_layers';
        veg_hom_scatCalVeg = 'scat_cal_veg';
        veg_hom_LTK = 'LTK';
        veg_hom_dsty = 'dsty';
        veg_hom_dim1 = 'dim1';
        veg_hom_dim2 = 'dim2';
        veg_hom_dim3 = 'dim3';
        veg_hom_epsr = 'epsr';
        veg_hom_parm1 = 'parm1';
        veg_hom_parm2 = 'parm2';
        veg_vir_row_plugin = 'plugin';
        veg_vir_row_rowSpace = 'row_space';
        veg_vir_row_colSpace = 'col_space';
        veg_vir_row_phiRow = 'phi_row';
        veg_vir_row_plantRowSpread = 'plant_row_spread';
        veg_vir_row_plantColSpread = 'plant_col_spread';
        
    end
    
end

