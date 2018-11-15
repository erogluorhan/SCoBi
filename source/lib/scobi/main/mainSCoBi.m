
function mainSCoBi
% function mainSCoBi
% 
%   Simulation iteration function. 
%
%   mainSCoBi is called by runSCoBi for every single simulation iteration. 
%   This function is common for any analysis type such as vegetation, 
%   bare-soil, or root-zone. The need for analysis type-specific function 
%   calls are determined by ParamManager class.
%
%   In each iteration, mainSCoBi calls the following functions:
%   - updateBistaticDynParams
%   - updateGndDynParams
%   - generateDielMLProfiles (if needed)
%   - calcPropagation (if needed)
%   - writeAttenuation (if needed)
%   - updateRotMatDynParams
%   - directTerm
%   - specularTerm
%
%   See also runSCoBi, ParamsManager, updateBistaticDynParams, 
%   updateGndDynParams, updateRotMatDynParams, directTerm, specularTerm.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



%% START OF THE CURRENT ITERATION
disp('+++++++++++++++++++++   START: mainSCoBi   ++++++++++++++++++++++++')



%% CALCULATE AND UPDATE BISTATIC CONFIGURATION AND GROUND DYNAMIC PARAMS
disp('Update Bistatic Configuration')

% Calculate the dynamic bistatic configuration and transmitter parameters
% and update the class BistaticDynParams 
updateBistaticDynParams();


disp('Update Dynamic Ground Parameters')

% Calculate the dynamic ground parameters and update the class GndDynParams
updateGndDynParams();



%% GENERATE DIELECTRIC PROFILES
% Use ParamsManager to determine if multilayer dielectric profiles needed
[needForDielProfiles, dispMsg] = ParamsManager.isToGenerateDielMLProfiles();

disp( dispMsg );
    
% If multilayer dielectric profiles needed
if ~needForDielProfiles == Constants.NEED_TO_RUN_STRUCT.NO

    generateDielMLProfiles(); 
    
end



%% CALCULATE INCREMENTAL PROPAGATION CONSTANT FOR EACH LAYER
% Use ParamsManager to determine if vegetation propagation needed
[needForPropagation, needForWriteAttenuation, dispMsg] = ParamsManager.isToCalcPropagation();

disp( dispMsg );

% If multilayer dielectric profles needed
if needForPropagation == Constants.NEED_TO_RUN_STRUCT.FULL
    
    % Calculate vegetation propagation
    calcPropagation();
    
    % If writing attenuation values to Excel file needed
    if needForWriteAttenuation
        
        disp('VEGETATION - Write Propagation Values to Excel File')
        
        writeAttenuation();
        
    end
    
end



%% CALCULATE AND UPDATE ROTATION MATRICES
disp('Update Polarization Rotation Matrices')     

% Calculate the dynamic polarizaztion rotation matrices and update the 
% class RotMatDynParams
updateRotMatDynParams();



%% CALCULATE DIRECT CONTRIBUTION
disp('Calculate Direct Contribution')
    
% Calculate the direct (line-of-sight) contribution in the received signals
directTerm();



%% CALCULATE SPECULAR CONTRIBUTION
disp('Calculate Specular Contribution')
    
% Calculate the coherent (through specular reflection point) contribution 
% in the received signals
specularTerm();


%% END OF THE CURRENT ITERATION
disp('------------------------   END: mainSCoBi   -----------------------')


end