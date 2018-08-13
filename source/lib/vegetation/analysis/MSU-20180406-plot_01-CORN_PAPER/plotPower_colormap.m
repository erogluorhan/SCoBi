
% Orhan Eroglu
% 4/13/2018
% Diffuse + Specular

function plotPower_colormap


%% ANALYSIS INPUTS
inputFile_sys_tag = 'sysInput-CORN_PAPER-row_';
inputFile_veg_tag = 'vegVirRowInput-Corn-row';

cornStagesStruct = struct('s1', 'V1_V9', ...
                            's2', 'V10_VT', ...
                            's3', 'R1_R4', ...
                            's4', 'R5', ...
                            's5', 'R6') ;

% Corn field row azimuth angles
rowPhis = [0 45 90] ;


%% Choices
% TO-DO: Check here
stage_choice = cornStagesStruct.s3;
rowPhi_choice = 1;
FZ_choice = 1 ;
TH_choice = 4;
PH_choice = 1 ;
RMSH_choice = 2 ;


%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_sys = strcat( inputFile_sys_tag, num2str( rowPhis(rowPhi_choice) ), '.xml' );
inputFile_veg = strcat( inputFile_veg_tag, num2str( rowPhis(rowPhi_choice) ), '-', stage_choice, '.xml' );

getInput( inputFile_sys, inputFile_veg );

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();
% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();


%% GET GLOBAL PARAMETERS
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg ;
num_Th = length( th0_Tx_list_deg );
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3 ;
num_VSM = length( VSM_list_cm3cm3 );
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm ;
num_RMSH = length( RMSH_list_cm );


%% SPECULAR TERM
% Intilization
% over vegetation
P_cohveg_dB = zeros(2, 2, num_Th, num_VSM, num_RMSH) ;

KKc_dB = zeros(num_Th, num_VSM, num_RMSH) ;

for tt = 1 : num_Th        
    
    ParamsManager.index_Th( tt );
    ParamsManager.index_Ph( PH_choice );
    
    for vv = 1 : num_VSM

        % Set theta index
        ParamsManager.index_VSM( vv );
    
        for rr = 1 : num_RMSH

            % Set theta index
            ParamsManager.index_RMSH( rr );

            % Initialize the directories depending on theta phi, VSM, and RMSH
            SimulationFolders.getInstance.initializeDynamicDirs();


            %% GET GLOBAL DIRECTORIES
            dir_out_specular = SimulationFolders.getInstance.out_specular;

            % read Kc
            filename1 = strcat('Kc') ;
            Kc = readComplexVar( dir_out_specular, filename1) ;
            KKc_dB(tt, vv, rr) = 20 * log10(abs(Kc)) ;

            [~, ~, ~, ~, ~, ~, ~, ~, ...
                P_coh1vegx, P_coh2vegx, ~, ~, ~, ~, ~, ~] = readSpecular ;

            P_cohveg_dB(1, :, tt, vv, rr) = 10 * log10( P_coh1vegx(1 : 2) ) + KKc_dB(tt, vv, rr) ;
            P_cohveg_dB(2, :, tt, vv, rr) = 10 * log10( P_coh2vegx(1 : 2) ) + KKc_dB(tt, vv, rr) ;

        end
        
    end

    if tt == 1
        dir_analysis = SimulationFolders.getInstance.analysis;
        if ~exist(dir_analysis, 'dir')
            mkdir(dir_analysis)
        end
    end
    
end

len = num_Th * num_VSM * num_RMSH;
indices = [];
indices(:,3) = repmat( (1:num_RMSH)', num_Th*num_VSM, 1);
ind = 1;
for vv = 1 : len / num_RMSH
    for rr = 1 : num_RMSH
        val = mod( vv, num_VSM);
        if val == 0
            val = num_VSM;
        end
        indices(ind,2) = val;
        ind = ind + 1;
    end
end
ind = 1;
for vv = 1 : len / num_RMSH / num_VSM
    for rr = 1 : num_VSM * num_RMSH
        indices(ind,1) = vv;
        ind = ind + 1;
    end
end

for ii = 1 : len
    
    P(ii,1) = P_cohveg_dB(1, 2, indices(ii,1), indices(ii,2), indices(ii,3) );
    
end


th0_scatter3 = th0_Tx_list_deg(indices(:,1))';
VSM_scatter3 = VSM_list_cm3cm3(indices(:,2))';
RMSH_scatter3 = RMSH_list_cm(indices(:,3))';
figure
scatter3( th0_scatter3, VSM_scatter3, RMSH_scatter3, 100, P, 'filled' );
xlabel('TH (\circ)')
ylabel('VSM (cm^3/cm^3)')
zlabel('RMSH (cm)')
cb = colorbar;
ylabel(cb, 'Received Power (dB)')



[VSM_RMSH_TH, RMSH_TH_VSM, TH_VSM_RMSH] =  meshgrid( VSM_list_cm3cm3 , RMSH_list_cm, th0_Tx_list_deg );

for rr = 1 : num_RMSH
    for vv = 1 : num_VSM
        for tt = 1 : num_Th
            P_TH_RMSH_VSM(rr,vv,tt) = P_cohveg_dB(1, 2, tt, vv, rr);        
        end
    end
end

figure
slice( VSM_RMSH_TH, RMSH_TH_VSM, TH_VSM_RMSH, P_TH_RMSH_VSM, VSM_list_cm3cm3, RMSH_list_cm, th0_Tx_list_deg );
xlabel('VSM (cm^3/cm^3)')
ylabel('RMSH (cm)')
zlabel('TH (\circ)')
cb = colorbar;
ylabel(cb, 'Received Power (dB)')



[VSM_RMSH, RMSH_VSM] =  meshgrid( VSM_list_cm3cm3 , RMSH_list_cm );

for rr = 1 : num_RMSH
    for vv = 1 : num_VSM
        P_RMSH_VSM(rr,vv) = P_cohveg_dB(1, 2, TH_choice, vv, rr);        
    end
end

figure
surf( VSM_RMSH, RMSH_VSM, P_RMSH_VSM,'FaceColor','interp')
xlabel('VSM (cm^3/cm^3)')
ylabel('RMSH (cm)')
zlabel('Received Power (dB)')
cb = colorbar;
ylabel(cb, 'Received Power (dB)')

figure
imagesc( P_RMSH_VSM)
xlabel('VSM (cm^3/cm^3)')
ylabel('RMSH (cm)')
cb = colorbar;
ylabel(cb, 'Received Power (dB)')


[VSM_TH, TH_VSM] =  meshgrid( VSM_list_cm3cm3, th0_Tx_list_deg );

for tt = 1 : num_Th
    for vv = 1 : num_VSM         
        P_TH_VSM(tt,vv) = P_cohveg_dB(1, 2, tt, vv, RMSH_choice );        
    end
end

figure
surf( VSM_TH, TH_VSM, P_TH_VSM)
xlabel('VSM (cm^3/cm^3)')
ylabel('TH (\circ)')
zlabel('Received Power (dB)')
cb = colorbar;
ylabel(cb, 'Received Power (dB)')

figure
imagesc( P_TH_VSM)
xlabel('VSM (cm^3/cm^3)')
ylabel('TH (\circ)')
cb = colorbar;
ylabel(cb, 'Received Power (dB)')


% TO-DO: Check this
%% GET GLOBAL PARAMETERS
pol_Tx = TxParams.getInstance.pol_Tx;
pol_Rx = RxParams.getInstance.pol_Rx;


%%
fname1 = strcat('3D-FZ_', num2str(FZ_choice), '-', pol_Tx, pol_Rx ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end