% Mehmet Kurum
% April 6, 2017

function realizeRotation(indRealization)


% indRealization : Nth Realization
% indRealization = 1 ;


%% Reading vegetation parameters...
scat_cal_veg = VegParams.getInstance.scat_cal_veg;

TYPKND = VegParams.getInstance.TYPKND;


%% Layer parameters
[Nlayer, Ntype] = size(TYPKND) ;

%% Calculations...

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

                % +++++++++++++
                calc(filename)
                % +++++++++++++

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

% Npart
% thrd
% phrd
% thrdI
% phrdI

%% Reading Receiver Observation Angles

disp('Reading observation angles...')
tic;
pathname = SimulationFolders.getInstance.observation ;
filenamex = strcat('thrd_', filename) ;
thrd = readVar(pathname, filenamex) ;
filenamex = strcat('thrdI_', filename) ;
thrdI = readVar(pathname, filenamex) ;
filenamex = strcat('phrd_', filename) ;
phrd = readVar(pathname, filenamex) ;
filenamex = strcat('phrdI_', filename) ;
phrdI = readVar(pathname, filenamex) ;
toc

Npart = length(thrd) ;
%% Gennting th and ph
load([SimulationFolders.getInstance.ant_lookup '\AntPat.mat'], 'th', 'ph')


%% Realization of Rotation matrices
disp('loading Rotation Matrices . . .')
tic ;
pathname = SimulationFolders.getInstance.rot_lookup ;
% 2 X 2
load([pathname '\u_gar.mat'], 'u_gar', 'th', 'ph')
load([pathname '\u_garI.mat'], 'u_garI', 'th', 'ph')
load([pathname '\u_sr.mat'], 'u_sr')
load([pathname '\u_ts.mat'], 'u_ts')
load([pathname '\u_tr.mat'], 'u_tr')
% % % 4 X 4
% % load([pathname '\U_gar.mat'], 'U_gar', 'th', 'ph')
% % load([pathname '\U_garI.mat'], 'U_garI', 'th', 'ph')
% % load([pathname '\U_sr.mat'], 'U_sr')
% % load([pathname '\U_ts.mat'], 'U_ts')
% % load([pathname '\U_tr.mat'], 'U_tr')
toc

thd = round( 2 * th * Constants.rad2deg ) / 2 ; % rounding operation is due to accuracy concerns
phd = round( 2 * ph * Constants.rad2deg ) / 2 ;

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

    ind_th = thd == round(2 * thrd(ii, 1)) / 2 ;
    ind_ph = phd == round(2 * phrd(ii, 1)) / 2 ;
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
    
    ind_th = thd == round(2 * thrdI(ii, 1)) / 2 ;
    ind_ph = phd == round(2 * phrdI(ii, 1)) / 2 ;
    % 2 X 2
    u_garI2(ii, :, :) = cell2mat(u_garI(ind_th & ind_ph)) ;
% %     % 4 X 4
% %     U_garI2(ii, :, :) = cell2mat(U_garI(ind_th & ind_ph)) ;
    
end
toc

%% Saving Rotation Matrix values
disp('saving rotation matrices...')
tic;
pathname = SimulationFolders.getInstance.rot_real ;

% 2 X 2
filenamex = strcat('u_gar_', filename) ;
writeComplexVar(pathname, filenamex, u_gar2)
filenamex = strcat('u_garI_', filename) ;
writeComplexVar(pathname, filenamex, u_garI2)
% % % 4 X 4
% % filenamex = strcat('U_gar_', filename) ;
% % writeComplexVar(pathname, filenamex, U_gar2)
% % filenamex = strcat('U_garI_', filename) ;
% % writeComplexVar(pathname, filenamex, U_garI2)

toc





end

