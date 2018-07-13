
% Mehmet Kurum
% April 6, 2017


function BiE2ES_main_specular_only

%% Global
global FolderPath_afsa     
global FolderPath_Ant FolderPath_rot 

%% Start
% disp('+++++++++++++++++++Start:BiE2ES Main Program++++++++++++++++++++++++')
% tstart = datetime('now') %#ok<NOPRT,*NASGU>

%% Set Input variables and Configuration
% disp('+++++++++++++++++++++++Setting input parameters+++++++++++++++++++++')
% Creating directories
CreateDir3 ;
% Setting Configuration parameters
SetConfig3 ;
% Vegetation Parameters
SetNlayerVeg3
% Ground Parameters
SetGround3 ;

%% Propagation Constant
% disp('++++++++Determine Incremental Propagation Contant for each layer++++')
MainFolder = pwd ;
cd(FolderPath_afsa) ;

if ~isempty(dir('*.dat'))  % Calculate only once
%     disp('Already done!')
    cd(MainFolder)
else
    cd(MainFolder)
%     t = datetime('now') %#ok<NOPRT>
    % Incremental vegetation constants
    Propagation3 ;
    disp('+++++++++++++++++++Writing attenuation values++++++++++++++++++++++++')
%     t = datetime('now') %#ok<NOPRT>
    WriteAttenuation3 ;
end

%% Receive Antenna Pattern Matrix
% disp('+++++++++++++++Calculating Receive Antenna Pattern Matrix+++++++++++') 
% t = datetime('now') %#ok<NOPRT>

% look-up
cd(strcat(FolderPath_Ant, '\', 'LOOK-UP')) ;
if ~isempty(dir('*.mat'))  % Calculate only once
%     disp('Already done!')
    cd(MainFolder)
else
    cd(MainFolder)
    AntennaPatternMatrix3 ;
end


%% Calculate Rotation Matrices
% disp('+++++++++++++++Calculating Rotation Matrices+++++++++++') 
% t = datetime('now') %#ok<NOPRT>

% look-up
cd(strcat(FolderPath_rot, '\', 'LOOK-UP')) ;
if ~isempty(dir('*.mat'))  % Calculate only once
%     disp('Already done!')
    cd(MainFolder)
else
    cd(MainFolder)
    CalcRotMats3 ;
end

%% Calculate Direct Contribution
% disp('++++++++++++++++Calculate Direct Contribution+++++++++++++++++++++++')
% t = datetime('now') %#ok<NOPRT,*NASGU>

% DirectTerm3 ;

%% Calculate Specular Contribution
% disp('++++++++++++++++Calculate Specular Contribution+++++++++++++++++++++')
% t = datetime('now') %#ok<NOPRT,*NASGU>

SpecularTerm3 ;

%% end of the program
% disp('++++++++++++++++++++End: BiE2ES Main Program++++++++++++++++++++++++')
% tstart %#ok<NOPRT>
% tstop = datetime('now') %#ok<NOPRT,*NASGU>
% duration = tstop - tstart  %#ok<NOPRT>


end