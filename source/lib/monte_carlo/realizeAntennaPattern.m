% Mehmet Kurum
% April 6, 2017

function realizeAntennaPattern(indRealization)


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
                      
            % +++++++++++++++++++++++++
            if jj == 1  % Leaf
                % do not calculate (leaves contribute absrobtion mostly)
                disp('skiped...')
            else % cylinder
                tic ;
                if scat_cal_veg(kk, jj, ii) == 1
                    disp('calculating...')
                    
                    % +++++++++++++
                    calc(filename)
                    % +++++++++++++
                    
                    disp('done...')
                    toc
                else
                    disp('skiped...')
                end
                
            end
            % +++++++++++++++++++++++++                        
            
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
%% Realization of Gain functions
% thrd, phrd, thrdI, phrdI
disp('loading Antenna Pattern Matrix . . .')
tic ;
load([SimulationFolders.getInstance.ant_lookup '\AntPat.mat'], 'g', 'th', 'ph')
toc

thd = round(2 * th * Constants.rad2deg) / 2 ; % rounding operation is due to accuracy concerns
phd = round(2 * ph * Constants.rad2deg) / 2 ;
gr = zeros(Npart, 2, 2) ;
grI = zeros(Npart, 2, 2) ;

% Gr = zeros(Npart, 4, 4) ;
% GrI = zeros(Npart, 4, 4) ;

%%
disp('Receiver Antenna values in the scattering directions. . .')
tic ;
for ii = 1 : Npart
    
    ind_th = thd == round(2 * thrd(ii, 1)) / 2 ; % round is to make it a multiple of 0.5
    ind_ph = phd == round(2 * phrd(ii, 1)) / 2 ;
    
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

%%
disp('Image Receiver Antenna values in the scattering directions. . .')
tic ;
for ii = 1 : Npart

    ind_th = thd == round(2 * thrdI(ii, 1)) / 2 ;
    ind_ph = phd == round(2 * phrdI(ii, 1)) / 2 ;
    
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

%% Saving Antenna gain values
disp('Saving antenna values...')
tic;
filenamex = strcat('gr_', filename) ;
writeComplexVar(SimulationFolders.getInstance.ant_real, filenamex, gr)
filenamex = strcat('grI_', filename) ;
writeComplexVar(SimulationFolders.getInstance.ant_real, filenamex, grI)
% % filenamex = strcat('Gr_', filename) ;
% % write_cplxvar(pathname, filenamex, Gr)
% % filenamex = strcat('GrI_', filename) ;
% % write_cplxvar(pathname, filenamex, GrI)
toc




end

