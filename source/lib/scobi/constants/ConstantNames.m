classdef ConstantNames
    %CONSTANTNAMES Global constant naming strings
    %   Keeps the global constants for the names of anything (e.g. files, 
    %   variables, xml tags, etc.) in the simulation
    
    properties (Constant)        
        
        
        %% INPUT FILES  
        % Master simulations inputs file name
        master_sim_input_filename = 'master_sims.xlsx';
        
        % Full path and name of the input file that is studied in the last
        % execution of SCoBi is kept in the following file, if any 
        lastInputFile = 'last_input.mat';
        
        % Input parameters file that is generated from all inputted static 
        % parameters 
        inputParamsStruct = 'inputParamsStruct';
        inputParamsStruct_filename = 'inputParamsStruct.mat';
        
        attenuation_out_filename = '\Attenuation.xls';
                     
    end
    
end

