classdef ConstantNames
    %CONSTANTNAMES Global constant naming strings
    %   Keeps the global constants for the names of anything (e.g. files, 
    %   variables, xml tags, etc.) in the simulation
    
    properties (Constant)
        
        version = 'version'; % String constant for sim version
        
        
        
        %% INPUT FILES  
        % Master simulations inputs file name
        master_sim_input_filename = 'master_sims.xlsx';
        
        % Default input file that should be under \SCoBi\source\input\sys
        defaultScoBiVegForInputFileName = 'default-veg_for-ss-single.mat';
        defaultScoBiRootZoneInputFileName = 'default-root_zone-ts-ml.mat';
        
        % Full path and name of the input file that is studied in the last
        % execution of SCoBi is kept in the following file, if any 
        lastInputFile = 'last_input.mat';
        
        % Input parameters file that is generated from all inputted static 
        % parameters 
        inputParamsStruct = 'inputParamsStruct';
        inputParamsStruct_filename = 'inputParamsStruct.mat';
        
        attenuation_out_filename = '\Attenuation.xls';
             
        sys_input = 'sys_input'; % String constant for system inputs
        
        veg_input = 'veg_input'; % String constant for vegetation input
               
        
        % String constants for simulators
        scobi_veg_agr = 'SCoBi-Veg(Agriculture)';
        scobi_veg_for = 'SCoBi-Veg(Forest)';
        scobi_ml = 'SCoBi-ML';
        
    end
    
end

