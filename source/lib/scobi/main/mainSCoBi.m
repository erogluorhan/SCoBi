
function mainSCoBi
%mainSCoBi Simulation iteration function. 
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

%    Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%    This program is free software: You can redistribute it and/or 
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation, either version 3 of the 
%    License, or (at your option) any later version.

%   Version: 1.0.0

%% start of the program
disp('++++++++++++++++   START: mainSCoBi   ++++++++++++++++++++')


%% CALCULATE AND UPDATE BISTATIC CONFIGURATION AND GROUND DYNAMIC PARAMS
disp('++++++++   UPDATE BISTATIC CONFIGURATION AND GROUND   +++++++++')

% Update Bistaic Configuration parameters
updateBistaticDynParams();

% Update Ground Dynamic Parameters
updateGndDynParams();



%% GENERATE DIELECTRIC PROFILES
% TO-DO: Test ParamsManager controls
[needForDielProfiles, dispMsg] = ParamsManager.isToGenerateDielMLProfiles();

disp( dispMsg );
    
if ~needForDielProfiles == Constants.need_for_run.NO

    generateDielMLProfiles(); 
    
end


%% CALCULATE INCREMENTAL PROPAGATION CONSTANT FOR EACH LAYER
disp('++   CALCULATE INCREMENTAL PROPAGATION CONSTANT FOR EACH LAYER   ++')

[needForPropagation, needForWriteAttenuation, dispMsg] = ParamsManager.isToCalcPropagation();

disp( dispMsg );

if needForPropagation == Constants.need_for_run.FULL
    
    calcPropagation();
    
    disp('++++++   WRITE ATTENUATION VALUES TO OUTPUT EXCEL FILE   ++++++')
    
    if needForWriteAttenuation
        
        writeAttenuation();
        
    end
end


%% CALCULATE AND UPDATE ROTATION MATRICES
disp('+++++++++++   CALCULATE AND UPDATE ROTATION MATRICES   ++++++++++++')     

updateRotMatDynParams();


%% CALCULATE DIRECT CONTRIBUTION
disp('+++++++++++++   CALCULATE DIRECT CONTRIBUTION   +++++++++++++++++++')
    
directTerm();


%% CALCULATE SPECULAR CONTRIBUTION
disp('+++++++++++++   CALCULATE SPECULAR CONTRIBUTION   +++++++++++++++++')
    
specularTerm();


%% end of the program
disp('++++++++++++++++   END: mainSCoBi   ++++++++++++++++++++')


end