
% Mehmet Kurum
% April 6 22, 2017

function CreateDir3


global hr th0d hpbw SLL XPL polT polR typeOfCanopy 
global FolderPath_Config FolderPath_Veg FolderPath_gnd FolderPath_afsa
global FolderPath_hr FolderPath_th0d FolderPath_pol 
global  FolderPath_Ant FolderPath_rot FolderPath_out FolderPath_fig



%% Receiver Height
Folder_hr = strcat('hr-', num2str(hr)) ;
FolderPath_hr = strcat('SIM\', typeOfCanopy, '\', Folder_hr) ;
if ~exist(FolderPath_hr, 'dir')
    mkdir(FolderPath_hr)
end

%% Angle of incidence
Folder_th0d = strcat('th0d-', num2str(th0d)) ;
FolderPath_th0d = strcat(FolderPath_hr, '\', Folder_th0d) ;
if ~exist(FolderPath_th0d, 'dir')
    mkdir(FolderPath_th0d)
end

%% Polarization

FolderPath_pol = strcat(FolderPath_th0d, '\', polT, polR) ;
if ~exist(FolderPath_pol, 'dir')
    mkdir(FolderPath_pol)
end

%% Average Forward Scattering Amplitude
FolderPath_afsa = strcat(FolderPath_hr, '\', 'AFSA') ;
if ~exist(FolderPath_afsa, 'dir')
    mkdir(FolderPath_afsa)
end

% % %% Bistatic Scattering Amplitudes
% % FolderPath_fscat = strcat(FolderPath_th0d, '\', 'FSCAT') ;
% % if ~exist(FolderPath_fscat, 'dir')
% %     mkdir(FolderPath_fscat)
% % end

%% Vegetation

FolderPath_Veg = strcat(FolderPath_hr, '\', 'VEG') ;
if ~exist(FolderPath_Veg, 'dir')
    mkdir(FolderPath_Veg)
end

%% Ground

FolderPath_gnd = strcat(FolderPath_hr, '\', 'GND') ;
if ~exist(FolderPath_gnd, 'dir')
    mkdir(FolderPath_gnd)
end

%% Configuration

FolderPath_Config = strcat(FolderPath_th0d, '\', 'CONFIG') ;
if ~exist(FolderPath_Config, 'dir')
    mkdir(FolderPath_Config)
end

% % %% Geometry
% % FolderPath_geo = strcat(FolderPath_th0d, '\', 'GEO') ;
% % if ~exist(FolderPath_geo, 'dir')
% %     mkdir(FolderPath_geo)
% % end
% % 
% % dirName = strcat(FolderPath_geo, '\', 'POSITION') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_geo, '\', 'FZONES') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_geo, '\', 'DISTANCE') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_geo, '\', 'INCIDENCE') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_geo, '\', 'SCATTERING') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_geo, '\', 'OBSERVATION') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end

%% Antenna

FolderNameAnt = strcat('thsd-', num2str(hpbw),...
    '_SLL-', num2str(SLL), '_XPL-', num2str(XPL)) ;

FolderPath_Ant = strcat(FolderPath_th0d, '\', 'ANT\', FolderNameAnt) ;
if ~exist(FolderPath_Ant, 'dir')
    mkdir(FolderPath_Ant)
end

dirName = strcat(FolderPath_Ant, '\', 'LOOK-UP') ;
if ~exist(dirName, 'dir')
    mkdir(dirName)
end

% % dirName = strcat(FolderPath_Ant, '\', 'REALIZATION') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end

%% Rotation

FolderPath_rot = strcat(FolderPath_pol, '\', 'ROT') ;
if ~exist(FolderPath_rot, 'dir')
    mkdir(FolderPath_rot)
end

dirName = strcat(FolderPath_rot, '\', 'LOOK-UP') ;
if ~exist(dirName, 'dir')
    mkdir(dirName)
end

% % dirName = strcat(FolderPath_rot, '\', 'REALIZATION') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end

% %% Frequncy Response
% dirName = strcat(FolderPath_pol, '\', 'FREQ') ;
% if ~exist(dirName, 'dir')
%     mkdir(dirName)
% end

% % FolderPath_freqdiff = strcat(FolderPath_pol, '\FREQ\DIFFUSE\', FolderNameAnt) ;
% % if ~exist(FolderPath_freqdiff, 'dir')
% %     mkdir(FolderPath_freqdiff)
% % end
% % 
% % dirName = strcat(FolderPath_freqdiff, '\', 'b1') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_freqdiff, '\', 'b2') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_freqdiff, '\', 'b3') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_freqdiff, '\', 'b4') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_freqdiff, '\', 'P1') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_freqdiff, '\', 'P2') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_freqdiff, '\', 'P3') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_freqdiff, '\', 'P4') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end

%% Output
FolderPath_out = strcat(FolderPath_pol, '\OUTPUT\', FolderNameAnt) ;
if ~exist(FolderPath_out, 'dir')
    mkdir(FolderPath_out)
end

dirName = strcat(FolderPath_out, '\', 'DIRECT') ;
if ~exist(dirName, 'dir')
    mkdir(dirName)
end

dirName = strcat(FolderPath_out, '\', 'SPECULAR') ;
if ~exist(dirName, 'dir')
    mkdir(dirName)
end

% % dirName = strcat(FolderPath_out, '\', 'DIFFUSE') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_out, '\', 'DIFFUSE', '\', 'b1') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_out, '\', 'DIFFUSE', '\', 'b2') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_out, '\', 'DIFFUSE', '\', 'b3') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_out, '\', 'DIFFUSE', '\', 'b4') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end

%% Figure

FolderPath_fig = strcat(FolderPath_hr, '\', 'FIGURE\', FolderNameAnt) ;
if ~exist(FolderPath_fig, 'dir')
    mkdir(FolderPath_fig)
end

dirName = strcat(FolderPath_fig, '\', 'DIRECT') ;
if ~exist(dirName, 'dir')
    mkdir(dirName)
end

dirName = strcat(FolderPath_fig, '\', 'SPECULAR') ;
if ~exist(dirName, 'dir')
    mkdir(dirName)
end

% % dirName = strcat(FolderPath_fig, '\', 'DIFFUSE') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_fig, '\', 'DIFFUSE', '\', 'b1') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_fig, '\', 'DIFFUSE', '\', 'b2') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_fig, '\', 'DIFFUSE', '\', 'b3') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % dirName = strcat(FolderPath_fig, '\', 'DIFFUSE', '\', 'b4') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end

%%
% % 
% % dirName = strcat('SIM\', typeOfCanopy, '\', 'FOLDER_NAMES') ;
% % if ~exist(dirName, 'dir')
% %     mkdir(dirName)
% % end
% % 
% % % location of permanent folder names
% % filename = strcat(pwd, '\SIM\FOLDER_NAMES\', Folder_hr, '.txt') ;
% % fid = fopen(filename, 'w') ;
% % fprintf(fid, '%s', Folder_hr) ;
% % fclose(fid) ;
% % 
% % % location of current folder names
% % filename = 'FolderName' ;
% % filename = strcat(pwd, '\', filename, '.txt') ;
% % fid = fopen(filename, 'w') ;
% % fprintf(fid, '%s', Folder_hr) ;
% % fclose(fid) ;

end