% Mehmet Kurum
% April 6, 2017

function realizeRotation(indRealization)

%% GET GLOBAL PARAMETERS
% Vegetation Parameters
scat_cal_veg = VegParams.getInstance.scat_cal_veg;
TYPKND = VegParams.getInstance.TYPKND;


%% CALCULATIONS
% Layer parameters
[Nlayer, Ntype] = size(TYPKND) ;

for ii = 1 : Nlayer
    
    for jj = 1 : Ntype
        
        Nkind = TYPKND(ii, jj) ;
        
        for kk = 1 : Nkind
            
            filename = strcat('R', num2str(indRealization), '_L', num2str(ii), '_T',...
                num2str(jj), '_K', num2str(kk)) ;
            disp(filename)           
                      
            % If the particle is a scatterer
            if scat_cal_veg(kk, jj, ii) == 1
                
                tic ;
                
                disp('calculating...')

                calc(filename)

                disp('done...')
                
                toc
                
            else
                
                % do not calculate
                disp('skiped...')
                
            end              
            
        end % Nkind
        
    end % Ntype
    
end % Nlayer


end


function calc(filename)

%% GET GLOBAL DIRECTORY
dir_observation = SimulationFolders.getInstance.observation;
dir_rot_lookup = SimulationFolders.getInstance.rot_lookup ;
dir_rot_real = SimulationFolders.getInstance.rot_real ;

%% GET GLOBAL PARAMETERS
% Receiver Parameters
ant_pat_Rx_id = RxParams.getInstance.ant_pat_Rx_id;
ant_pat_res_deg = RxParams.getInstance.ant_pat_res_deg;


%% READ OR LOAD META-DATA
% Receiver Observation Angles
disp('Reading observation angles...')
tic;
filenamex = strcat('thrd_', filename) ;
thrd = readVar(dir_observation, filenamex) ;
filenamex = strcat('thrdI_', filename) ;
thrdI = readVar(dir_observation, filenamex) ;
filenamex = strcat('phrd_', filename) ;
phrd = readVar(dir_observation, filenamex) ;
filenamex = strcat('phrdI_', filename) ;
phrdI = readVar(dir_observation, filenamex) ;
toc

Npart = length(thrd) ;

% Antenna Pattern and Lookup Angles (th and ph)
load([SimulationFolders.getInstance.ant_lookup '\AntPat.mat'], 'th', 'ph')

% Realization of Rotation matrices
disp('loading Rotation Matrices . . .')
tic ;
% 2 X 2
load([dir_rot_lookup '\u_gar.mat'], 'u_gar', 'th', 'ph')
load([dir_rot_lookup '\u_garI.mat'], 'u_garI', 'th', 'ph')
load([dir_rot_lookup '\u_sr.mat'], 'u_sr')
load([dir_rot_lookup '\u_ts.mat'], 'u_ts')
load([dir_rot_lookup '\u_tr.mat'], 'u_tr')
% % % 4 X 4
% % load([dir_rot_lookup '\U_gar.mat'], 'U_gar', 'th', 'ph')
% % load([dir_rot_lookup '\U_garI.mat'], 'U_garI', 'th', 'ph')
% % load([dir_rot_lookup '\U_sr.mat'], 'U_sr')
% % load([dir_rot_lookup '\U_ts.mat'], 'U_ts')
% % load([dir_rot_lookup '\U_tr.mat'], 'U_tr')
toc


if ant_pat_Rx_id == Constants.id_Rx_user_defined
    
    [~, Nth] = size(th);

    % Calculate the antenna pattern resolution in degrees
    ant_pat_res_deg = Constants.ant_pat_th_range_deg / (Nth - 1);

end

ant_pat_res_factor = 1 / ant_pat_res_deg;

thd = round( ant_pat_res_factor * th * Constants.rad2deg ) / ant_pat_res_factor ; % rounding operation is due to accuracy concerns
phd = round( ant_pat_res_factor * ph * Constants.rad2deg ) / ant_pat_res_factor ;

% 2 X 2
u_gar2 = zeros(Npart, 2, 2) ;
u_garI2 = zeros(Npart, 2, 2) ;
% % % 4 X 4
% % U_gar2 = zeros(Npart, 4, 4) ;
% % U_garI2 = zeros(Npart, 4, 4) ;


%%
disp('Rotation Matrices towards receiver. . .')
tic ;
for ii = 1 : Npart

    ind_th = thd == round( ant_pat_res_factor * thrd(ii, 1)) / ant_pat_res_factor ;
    ind_ph = phd == round( ant_pat_res_factor * phrd(ii, 1)) / ant_pat_res_factor ;
    % 2 X 2
    u_gar2(ii, :, :) = cell2mat(u_gar(ind_th & ind_ph)) ;
% %     % 4 X 4
% %     U_gar2(ii, :, :) = cell2mat(U_gar(ind_th & ind_ph)) ;
    
end
toc

%%
disp('Rotation Matrices towards image of receiver. . .')
tic ;
for ii = 1 : Npart
    
    ind_th = thd == round( ant_pat_res_factor * thrdI(ii, 1)) / ant_pat_res_factor ;
    ind_ph = phd == round( ant_pat_res_factor * phrdI(ii, 1)) / ant_pat_res_factor ;
    % 2 X 2
    u_garI2(ii, :, :) = cell2mat(u_garI(ind_th & ind_ph)) ;
% %     % 4 X 4
% %     U_garI2(ii, :, :) = cell2mat(U_garI(ind_th & ind_ph)) ;
    
end
toc


%% SAVE
% Rotation Matrix values
disp('saving rotation matrices...')
tic;

% 2 X 2
filenamex = strcat('u_gar_', filename) ;
writeComplexVar(dir_rot_real, filenamex, u_gar2)
filenamex = strcat('u_garI_', filename) ;
writeComplexVar(dir_rot_real, filenamex, u_garI2)
% % % 4 X 4
% % filenamex = strcat('U_gar_', filename) ;
% % writeComplexVar(dir_rot_real, filenamex, U_gar2)
% % filenamex = strcat('U_garI_', filename) ;
% % writeComplexVar(dir_rot_real, filenamex, U_garI2)

toc

end

