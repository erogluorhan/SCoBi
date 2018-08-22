
% Mehmet Kurum
% April 6, 2017

function mainSCoBi

%% START: SCoBi MAIN PROGRAM
disp('++++++++++++++++   START: SCoBi MAIN PROGRAM   ++++++++++++++++++++')
tstart = datetime('now') %#ok<NOPRT,*NASGU>


%% CALCULATE AND UPDATE BISTATIC CONFIGURATION AND GROUND DYNAMIC PARAMS
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


%% CALCULATE INCREMENTAL PROPAGATION CONSTANT FOR EACH LAYER
disp('++   CALCULATE INCREMENTAL PROPAGATION CONSTANT FOR EACH LAYER   ++')
t = datetime('now') %#ok<NOPRT>

[needForPropagation, needForWriteAttenuation, dispMsg] = ParamsManager.isToCalcPropagation();

disp( dispMsg );

if needForPropagation == Constants.need_for_run.FULL
    
    calcPropagation ;
    
    disp('++++++   WRITE ATTENUATION VALUES TO OUTPUT EXCEL FILE   ++++++')
    
    if needForWriteAttenuation
        
        writeAttenuation ;
        
    end
end


%% CALCULATE AND UPDATE ROTATION MATRICES
disp('+++++++++++   CALCULATE AND UPDATE ROTATION MATRICES   ++++++++++++')     
t = datetime('now')

updateRotMatDynParams();


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