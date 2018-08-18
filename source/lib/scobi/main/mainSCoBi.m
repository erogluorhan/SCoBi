
% Mehmet Kurum
% April 6, 2017

function mainSCoBi

%% START: SCoBi MAIN PROGRAM
disp('++++++++++++++++   START: SCoBi MAIN PROGRAM   ++++++++++++++++++++')
tstart = datetime('now') %#ok<NOPRT,*NASGU>


%% INITIALIZE BISTATIC CONFIGURATION AND GROUND
disp('++++++++   UPDATE BISTATIC CONFIGURATION AND GROUND   +++++++++')
t = datetime('now') %#ok<NOPRT>

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


%% GENERATE SCATTERING POSITIONS
disp('+++++++++++++++   GENERATE SCATTERING POSITIONS   +++++++++++++++++')
t = datetime('now') %#ok<NOPRT>

[needForScatPos, Nr_current, dispMsg] = ParamsManager.isToGenerateScatPos();

disp( dispMsg );
    
if ~needForScatPos == Constants.need_for_run.NO
    
        generateScatPos( Nr_current );
        
end


%% CALCULATE INCREMENTAL PROPAGATION CONSTANT FOR EACH LAYER
disp('++   CALCULATE INCREMENTAL PROPAGATION CONSTANT FOR EACH LAYER   ++')
t = datetime('now') %#ok<NOPRT>

[needForPropagation, needForWriteAttenuation, dispMsg] = ParamsManager.isToCalcPropagation();

disp( dispMsg );

if needForPropagation == Constants.need_for_run.FULL
    
    calcPropagation ;
    
    disp('++++++++   WRITE ATTENUATION VALUES TO OUTPUT EXCEL FILE   ++++++++')
    
    if needForWriteAttenuation
        
        writeAttenuation ;
        
    end
end


%% CALCULATE ROTATION MATRICES
disp('+++++++++++++++++   CALCULATE ROTATION MATRICES   +++++++++++++++++') 
t = datetime('now') %#ok<NOPRT>

[needForCalcRotationMatrices, dispMsg] = ParamsManager.isToCalcRotationMatrices();

disp( dispMsg );

% TO-DO: may be shrinked more if diffuse & specular will not be calculated
if needForCalcRotationMatrices == Constants.need_for_run.FULL
    
    calcRotationMatrices ;
    
end


%% CALCULATE DIRECT CONTRIBUTION
disp('+++++++++++++   CALCULATE DIRECT CONTRIBUTION   +++++++++++++++++++')
t = datetime('now')
    
directTerm;


%% CALCULATE SPECULAR CONTRIBUTION
disp('+++++++++++++   CALCULATE SPECULAR CONTRIBUTION   +++++++++++++++++')
t = datetime('now')
    
specularTerm ;


%% end of the program
disp('++++++++++++++++   START: SCoBi MAIN PROGRAM   ++++++++++++++++++++')
tstart %#ok<NOPRT>
tstop = datetime('now') %#ok<NOPRT,*NASGU>
duration = tstop - tstart  %#ok<NOPRT>


end