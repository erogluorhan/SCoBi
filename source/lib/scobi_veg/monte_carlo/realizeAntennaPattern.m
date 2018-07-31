% Mehmet Kurum
% April 6, 2017

function realizeAntennaPattern(indRealization)


%% GET GLOBAL PARAMETERS
% Vegetation Parameters
scat_cal_veg = VegParams.getInstance.scat_cal_veg;
TYPKND = VegParams.getInstance.TYPKND;


% Layer parameters
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
                    
                tic
                
                disp('calculating...')

                % +++++++++++++
                calc(filename)
                % +++++++++++++

                disp('done...')
                
                toc
                
            else
                
                disp('skiped...')
                
            end                      
            
        end % Nkind
        
    end % Ntype
    
end % Nlayer




end

function calc(filename)

%% GET GLOBAL DIRECTORIES
dir_observation = SimulationFolders.getInstance.observation;
dir_ant_real = SimulationFolders.getInstance.ant_real;

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

% Antenna Pattern
disp('loading Antenna Pattern Matrix . . .')
tic ;
load([SimulationFolders.getInstance.ant_lookup '\AntPat.mat'], 'g', 'th', 'ph')
toc


if ant_pat_Rx_id == Constants.id_Rx_user_defined
    
    [~, Nth] = size(th);

    % Calculate the antenna pattern resolution in degrees
    ant_pat_res_deg = Constants.ant_pat_th_range_deg / (Nth - 1);

end

ant_pat_res_factor = 1 / ant_pat_res_deg;

thd = round( ant_pat_res_factor * th * Constants.rad2deg) / ant_pat_res_factor ; % rounding operation is due to accuracy concerns
phd = round( ant_pat_res_factor * ph * Constants.rad2deg) / ant_pat_res_factor ;
gr = zeros(Npart, 2, 2) ;
grI = zeros(Npart, 2, 2) ;

%% ANTENNA PATTERN REALIZATIONS
% Receiver
disp('Receiver Antenna values in the scattering directions. . .')
tic ;
for ii = 1 : Npart
    
    ind_th = thd == round( ant_pat_res_factor * thrd(ii, 1) ) / ant_pat_res_factor ; % round is to make it a multiple of ant_pat_res_deg
    ind_ph = phd == round( ant_pat_res_factor * phrd(ii, 1) ) / ant_pat_res_factor ;
    
    %%
    g11 = cell2mat(g(1,1)) ; g12 = cell2mat(g(1,2)) ;
    g21 = cell2mat(g(2,1)) ; g22 = cell2mat(g(2,2)) ;
    
    gg = [g11(ind_th & ind_ph), g12(ind_th & ind_ph); ...
        g21(ind_th & ind_ph), g22(ind_th & ind_ph)] ;
    
    gr(ii, :, :) = gg ;
% %     %%
% %     G11 = cell2mat(G(1,1)) ; G12 = cell2mat(G(1,2)) ; G13 = cell2mat(G(1,3)) ; G14 = cell2mat(G(1,4)) ;
% %     G21 = cell2mat(G(2,1)) ; G22 = cell2mat(G(2,2)) ; G23 = cell2mat(G(2,3)) ; G24 = cell2mat(G(2,4)) ;
% %     G31 = cell2mat(G(3,1)) ; G32 = cell2mat(G(3,2)) ; G33 = cell2mat(G(3,3)) ; G34 = cell2mat(G(1,4)) ;
% %     G41 = cell2mat(G(4,1)) ; G42 = cell2mat(G(4,2)) ; G43 = cell2mat(G(4,3)) ; G44 = cell2mat(G(4,4)) ;
% %     
% %     GG = [G11(ind_th & ind_ph), G12(ind_th & ind_ph), G13(ind_th & ind_ph), G14(ind_th & ind_ph); ...
% %         G21(ind_th & ind_ph), G22(ind_th & ind_ph), G23(ind_th & ind_ph), G24(ind_th & ind_ph); ...
% %         G31(ind_th & ind_ph), G32(ind_th & ind_ph), G33(ind_th & ind_ph), G34(ind_th & ind_ph); ...
% %         G41(ind_th & ind_ph), G42(ind_th & ind_ph), G43(ind_th & ind_ph), G44(ind_th & ind_ph)] ;
% %     
% %     Gr(ii, :, :) = GG ;
    
end
toc

%Receiver Image
disp('Image Receiver Antenna values in the scattering directions. . .')
tic ;
for ii = 1 : Npart

    ind_th = thd == round( ant_pat_res_factor * thrdI(ii, 1)) / ant_pat_res_factor ;  % round is to make it a multiple of ant_pat_res_deg
    ind_ph = phd == round( ant_pat_res_factor * phrdI(ii, 1)) / ant_pat_res_factor ;
    
        %%
    g11 = cell2mat(g(1,1)) ; g12 = cell2mat(g(1,2)) ;
    g21 = cell2mat(g(2,1)) ; g22 = cell2mat(g(2,2)) ;
    
    gg = [g11(ind_th & ind_ph), g12(ind_th & ind_ph); ...
        g21(ind_th & ind_ph), g22(ind_th & ind_ph)] ;
    
    grI(ii, :, :) = gg ;
% %     %%
% %     G11 = cell2mat(G(1,1)) ; G12 = cell2mat(G(1,2)) ; G13 = cell2mat(G(1,3)) ; G14 = cell2mat(G(1,4)) ;
% %     G21 = cell2mat(G(2,1)) ; G22 = cell2mat(G(2,2)) ; G23 = cell2mat(G(2,3)) ; G24 = cell2mat(G(2,4)) ;
% %     G31 = cell2mat(G(3,1)) ; G32 = cell2mat(G(3,2)) ; G33 = cell2mat(G(3,3)) ; G34 = cell2mat(G(1,4)) ;
% %     G41 = cell2mat(G(4,1)) ; G42 = cell2mat(G(4,2)) ; G43 = cell2mat(G(4,3)) ; G44 = cell2mat(G(4,4)) ;
% %     
% %     GG = [G11(ind_th & ind_ph), G12(ind_th & ind_ph), G13(ind_th & ind_ph), G14(ind_th & ind_ph); ...
% %         G21(ind_th & ind_ph), G22(ind_th & ind_ph), G23(ind_th & ind_ph), G24(ind_th & ind_ph); ...
% %         G31(ind_th & ind_ph), G32(ind_th & ind_ph), G33(ind_th & ind_ph), G34(ind_th & ind_ph); ...
% %         G41(ind_th & ind_ph), G42(ind_th & ind_ph), G43(ind_th & ind_ph), G44(ind_th & ind_ph)] ;
% %     
% %     GrI(ii, :, :) = GG ;
   
end
toc


%% SAVE
disp('Saving antenna values...')
tic;
filenamex = strcat('gr_', filename) ;
writeComplexVar(dir_ant_real, filenamex, gr)
filenamex = strcat('grI_', filename) ;
writeComplexVar(dir_ant_real, filenamex, grI)
% % filenamex = strcat('Gr_', filename) ;
% % write_cplxvar(dir_observation, filenamex, Gr)
% % filenamex = strcat('GrI_', filename) ;
% % write_cplxvar(dir_observation, filenamex, GrI)
toc


end

