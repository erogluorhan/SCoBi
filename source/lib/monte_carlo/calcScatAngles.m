
function calcScatAngles(filename, pP)
% CALCSCATANGLES: Calculates  incident and scattering angles for given
% particle data
% pP: Particle center of gravity positions in ground frame with absolute z-values


%% GET GLOBAL DIRECTORIES
dir_config = SimulationFolders.getInstance.config;
dir_incidence = SimulationFolders.getInstance.incidence;
dir_scattering = SimulationFolders.getInstance.scattering ;
dir_observation = SimulationFolders.getInstance.observation ;
dir_distance = SimulationFolders.getInstance.distance ;


%% READ META-DATA (TRANSFORMATIONS AND CONFIGURATION)
% All Positions
filenamex = 'AllPoints_m' ;
AllPoints_m = readVar(dir_config, filenamex) ;

pT_m = AllPoints_m(:, 1);        % Transmitter
pTI_m = AllPoints_m(:, 2);        % Transmitter Image

pR_m = AllPoints_m(:, 4);        % Receiver
pRI_m = AllPoints_m(:, 5);       % Receiver Image

% Incident unit vectors (only in one direction - transmitter is far away)
filenamex = 'isn' ;
isn = readVar(dir_config, filenamex) ;   % propagation vector (i_s^-)
filenamex = 'osp' ;
isp = readVar(dir_config, filenamex) ;   % propagation vector (i_s^+=o_s^+)

% Transformations
filenamex = 'Tgs' ;
Tgs = readVar(dir_config, filenamex) ;   % G -> S
filenamex = 'Tgr' ;
Tgr = readVar(dir_config, filenamex) ;   % G -> R
filenamex = 'TgrI' ;
TgrI = readVar(dir_config, filenamex) ;  % G -> RI


%% CALCULATIONS
% Number of Particles
[Npart, ~] = size(pP) ; 

% repmat
pR_m = repmat(pR_m', Npart, 1) ;         % Receiver
pRI_m = repmat(pRI_m', Npart, 1) ;       % Image Receiver

pT_m = repmat(pT_m', Npart, 1) ;         % Transmitter
pTI_m = repmat(pTI_m', Npart, 1) ;       % Image Transmitter


% Scatteting directions and distances to Transmitter / Image Transmitter
ki = pT_m - pP ;                  % vector from particle to antenna
ri = sqrt(sum(ki' .^ 2))' ;     % distance from particle to antenna

kiI = pTI_m - pP ;                    % vector from particle to image antenna
riI = sqrt(sum(kiI' .^ 2))' ;       % distance from particle to image antenna


% Scatteting directions and distances to receiver/image reciever
ko = pR_m - pP ;                  % vector from particle to antenna
ro = sqrt(sum(ko' .^ 2))' ;     % distance from particle to antenna
ko = ko ./ repmat(ro, 1, 3) ;   % unit vector

koI = pRI_m - pP ;                    % vector from particle to image antenna
roI = sqrt(sum(koI' .^ 2))' ;       % distance from particle to image antenna
koI = koI ./ repmat(roI, 1, 3) ;    % unit vector

% TO-DO: Decide what to do with this
% figure
% riroI = ri + roI ;
%  riIro = riI + ro ;
% plot(riIro, 'or')
% hold
% plot(riroI, 'ob')

% figure
% subplot(2,2,1)
% plot(ro, ':o')
% title('r_o')
% subplot(2,2,2)
% plot(roI, ':o')
% title('r_o_I')
% subplot(2,2,3)
% plot(ri, ':o')
% title('r_i')
% subplot(2,2,4)
% plot(riI, ':o')
% title('r_i_I')
% 
% figure
% subplot(2,2,1)
% plot(ro+ri, ':o')
% title('dd : r_o + r_i')
% subplot(2,2,2)
% plot(roI+ri, ':o')
% title('rd : r_o_I + r_i')
% subplot(2,2,3)
% plot(ro+riI, ':o')
% title('dr : r_o + r_i_I')
% subplot(2,2,4)
% plot(riI+roI, ':o')
% title('rr : r_o_I + r_i_I')

% Scattering angles from particles 

% propagation vector to receiver in local (specular) frame
ko_sf = (Tgs * ko')' ;

% FROM PARTICLE TO REC
% incidence angle
thsd = acos(ko_sf(:, 3)) * 180 / pi ;
% azimuth angle
phsd = atan2(ko_sf(:, 2), ko_sf(:, 1)) * 180 / pi ;

% propagation vector to image of the receiver in local (specular) frame
koI_sf = (Tgs * koI')' ;

% FROM PARTICLE TO IM of REC
% incidence angle
thsdI = acos(koI_sf(:, 3)) * 180 / pi ;
% azimuth angle
phsdI = atan2(koI_sf(:, 2), koI_sf(:, 1)) * 180 / pi ;


%% Incident angles on particles

% propagation vector TRANS TO PARTICLE in local (specular) frame
ki_sf = (Tgs * isn)' ;

% incidence angle
thid = acos(-ki_sf(:, 3)) * 180 / pi ;
% azimuth angle
phid = atan2(-ki_sf(:, 2), -ki_sf(:, 1)) * 180 / pi ;

% propagation vector IM OF TRANS TO PARTICLE in local (specular) frame
kiI_sf = (Tgs * isp)' ;

% incidence angle
thidI = acos(-kiI_sf(:, 3)) * 180 / pi ;
% azimuth angle
phidI = atan2(-kiI_sf(:, 2), -kiI_sf(:, 1)) * 180 / pi ;

%% incident angles to receiver and image of receiver
% propagation vector in receiver antenna system
ko_rf = (Tgr * ko')' ;

% FROM PARTICLE TO REC
% off-axis angle of zr towards particle
thrd = mod(acos(-ko_rf(:, 3)) * 180 / pi + 360, 360) ;
% Reciever orientation - azimuth
phrd = mod(atan2(-ko_rf(:, 2), -ko_rf(:, 1)) * 180 / pi + 360, 360) ;
% we take mod to make sure the angle is between 0 and 360

% propagation vector in image receiver antenna system
koI_rIf = (TgrI * koI')' ;

% FROM PARTICLE TO IM of REC
% off-axis angle of zr towards particle
thrdI = mod(acos(-koI_rIf(:, 3)) * 180 / pi + 360, 360) ;
% Reciever orientation - azimuth
phrdI = mod(atan2(-koI_rIf(:, 2), -koI_rIf(:, 1)) * 180 / pi + 360, 360) ; 
% we take mod to make sure the angle is between 0 and 360


%% SAVE ALL 
% Incidence Angles
disp('Saving incidence angles...')
tic;
filenamex = strcat('thid_', filename) ;
writeVar(dir_incidence, filenamex, thid)
filenamex = strcat('thidI_', filename) ;
writeVar(dir_incidence, filenamex, thidI)
filenamex = strcat('phid_', filename) ;
writeVar(dir_incidence, filenamex, phid)
filenamex = strcat('phidI_', filename) ;
writeVar(dir_incidence, filenamex, phidI)
toc

% Scattering Angles
disp('Saving scattering angles...')
tic;
filenamex = strcat('thsd_', filename) ;
writeVar(dir_scattering, filenamex, thsd)
filenamex = strcat('thsdI_', filename) ;
writeVar(dir_scattering, filenamex, thsdI)
filenamex = strcat('phsd_', filename) ;
writeVar(dir_scattering, filenamex, phsd)
filenamex = strcat('phsdI_', filename) ;
writeVar(dir_scattering, filenamex, phsdI)
toc

% Receiver Observation Angles
disp('Saving observation angles...')
tic;
filenamex = strcat('thrd_', filename) ;
writeVar(dir_observation, filenamex, thrd)
filenamex = strcat('thrdI_', filename) ;
writeVar(dir_observation, filenamex, thrdI)
filenamex = strcat('phrd_', filename) ;
writeVar(dir_observation, filenamex, phrd)
filenamex = strcat('phrdI_', filename) ;
writeVar(dir_observation, filenamex, phrdI)
toc

%  Particle distance to antennas and image antennas
disp('Saving distances...')
tic;
filenamex = strcat('ri_', filename) ;
writeVar(dir_distance, filenamex, ri)
filenamex = strcat('riI_', filename) ;
writeVar(dir_distance, filenamex, riI)
filenamex = strcat('ro_', filename) ;
writeVar(dir_distance, filenamex, ro)
filenamex = strcat('roI_', filename) ;
writeVar(dir_distance, filenamex, roI)
toc


end