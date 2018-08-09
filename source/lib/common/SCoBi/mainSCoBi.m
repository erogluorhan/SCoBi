
% Mehmet Kurum
% April 6, 2017

function mainSCoBi

%% START: SCoBi MAIN PROGRAM
disp('++++++++++++++++   START: SCoBi MAIN PROGRAM   ++++++++++++++++++++')
tstart = datetime('now') %#ok<NOPRT,*NASGU>


%% INITIALIZE BISTATIC CONFIGURATION AND GROUND
disp('++++++++   INITIALIZE BISTATIC CONFIGURATION AND GROUND   +++++++++')
t = datetime('now') %#ok<NOPRT>

% Setting Configuration parameters
initConfiguration;

% Initialize Surface Dynamics (Diel)
initSurfaceDynamics;



%% GENERATE DIELECTRIC PROFILES
% TO-DO: Test ParamsManager controls
[needForDielProfiles, dispMsg] = ParamsManager.isToGenerateDielProfiles();

disp( dispMsg );
    
if ~needForDielProfiles == Constants.need_for_run.NO

    generateDielProfiles; 
    
end


%% GENERATE SCATTERING POSITIONS
disp('+++++++++++++++   GENERATE SCATTERING POSITIONS   +++++++++++++++++')
t = datetime('now') %#ok<NOPRT>

[needForScatPos, Nr_current, dispMsg] = ParamsManager.isToGenerateScatPos();

disp( dispMsg );
    
if ~needForScatPos == Constants.need_for_run.NO
    
        generateScatPos( Nr_current );
        
end

%% CALCULATE SCATTERING AMPLITUDES
disp('+++++++++++++   CALCULATE SCATTERING AMPLITUDES   +++++++++++++++++')
t = datetime('now') %#ok<NOPRT,*NASGU>

[needForFScatAmp, Nr_current, dispMsg] = ParamsManager.isToCalculateFScatAmp();

disp( dispMsg );

if needForFScatAmp ~= Constants.need_for_run.NO

    calcFScatAmp( Nr_current );

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


%% CALCULATE RECEIVER ANTENNA PATTERN REALIZATIONS
disp('+++++++   CALCULATE RECEIVER ANTENNA PATTERN REALIZATIONS   +++++++') 
t = datetime('now') %#ok<NOPRT>

[needForRealizeAntennaPattern, Nr_current, dispMsg] = ParamsManager.isToRealizeAntennaPattern();

disp( dispMsg );

if needForRealizeAntennaPattern ~= Constants.need_for_run.NO

    realizeAntennaPattern( Nr_current );

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


%% CALCULATE ROTATION REALIZATIONS
disp('+++++++++++++++   CALCULATE ROTATION REALIZATIONS   +++++++++++++++') 
t = datetime('now') %#ok<NOPRT>

[needForRealizeRotations, Nr_current, dispMsg] = ParamsManager.isToRealizeRotations();

disp( dispMsg );

if needForRealizeRotations ~= Constants.need_for_run.NO

    realizeRotation( Nr_current );

end


%% CALCULATE DIRECT CONTRIBUTION
disp('+++++++++++++   CALCULATE DIRECT CONTRIBUTION   +++++++++++++++++++')
t = datetime('now') %#ok<NOPRT,*NASGU>

[needForDirectTerm, dispMsg] = ParamsManager.isToCalculateDirectTerm();

disp( dispMsg );

if needForDirectTerm == Constants.need_for_run.FULL
    
    directTerm ;
    
end



%% CALCULATE SPECULAR CONTRIBUTION
disp('+++++++++++++   CALCULATE SPECULAR CONTRIBUTION   +++++++++++++++++')
t = datetime('now') %#ok<NOPRT,*NASGU>
 
[needForSpecularTerm, dispMsg] = ParamsManager.isToCalculateSpecularTerm();

disp( dispMsg );

if needForSpecularTerm == Constants.need_for_run.FULL
    
    specularTerm ;
    
end


%% CALCULATE DIFFUSE CONTRIBUTION
disp('+++++++++++++++   CALCULATE DIFFUSE CONTRIBUTION   ++++++++++++++++')
t = datetime('now') %#ok<NOPRT,*NASGU>
 
[needForDiffuseTerm, dispMsg] = ParamsManager.isToCalculateDiffuseTerm();

disp( dispMsg );

if needForDiffuseTerm == Constants.need_for_run.FULL

    diffuseTerm();

    % Averaging over realizations
    avgDiffuseTerm ;
    
end

%% end of the program
disp('++++++++++++++++   START: SCoBi MAIN PROGRAM   ++++++++++++++++++++')
tstart %#ok<NOPRT>
tstop = datetime('now') %#ok<NOPRT,*NASGU>
duration = tstop - tstart  %#ok<NOPRT>


end