classdef ConstantNames
% class ConstantNames 
%
%   Global constant naming strings (char arrays). This class only consists 
%   of constant (final) properties that keep the global names for any 
%   purpose (e.g. files, variables, xml tags, etc.) in the simulation.   
%
%   See also Constants, Directories.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0

    
    properties (Constant)        
        
        
        %% INPUT FILES  
        % Master simulations inputs file name
        MASTER_SIM_INPUT_FILENAME = 'master_sims.xlsx';
        
        % Input parameters file that is generated from all inputted static 
        % parameters 
        INPUT_PARAMS_STRUCT = 'inputParamsStruct';
        INPUT_PARAMS_STRUCT_FILENAME = 'inputParamsStruct.mat';
        
        % Excel file name to write the vegetation attenauation 
        ATTENUATION_OUT_FILENAME = '\Attenuation.xlsx';
                     
    end
    
end

