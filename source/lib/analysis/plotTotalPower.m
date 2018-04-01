
% Mehmet Kurum
% 11/15/2017
% Diffuse + Specular

function plotTotalPower


%% INPUT FILES
inputFile_sys_tag = 'sysInput-Paulownia-PAPER_PBAND-hr_';
inputFile_veg = 'vegHomInput-Paulownia.xml';


% TO-DO: Check here
%% Choices
Nfz = 10;
FZ_choice = 1 ;
SM_choice = 2 ;  % [5 15 25]

%% angle 
th0d = 10 : 10 : 70 ;

%% Reciever height
hhr = [20 50 100 500] ;

    
%% reading
% NBRCS
Nth0d =  length(th0d) ;
Nhr =  length(hhr) ;

PP1_inc_dB = zeros(2, 2, Nth0d, Nhr) ; % % PP1_inc1_dB = zeros(2, 2, Nth0d) ; 
% % PP1_inc2_dB = zeros(2, 2, Nth0d) ; PP1_inc3_dB = zeros(2, 2, Nth0d) ; PP1_inc4_dB = zeros(2, 2, Nth0d) ;
PP1_inc_dB2 = zeros(2, 2, Nfz, Nth0d, Nhr) ;

KKi_dB = zeros(Nth0d, Nhr) ;
P1_areas_dB = zeros(Nth0d, Nhr) ;
P1_areas_dB2 = zeros(Nfz, Nth0d, Nhr) ;

for jj = 1 : Nhr
        
    inputFile_sys = strcat( inputFile_sys_tag, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    for tt = 1 : Nth0d
              
        % Set theta index
        ParamsManager.index_Th( tt );
        
        % TO-DO: Assign the correct index to each
        ParamsManager.index_Ph( 1 );
        ParamsManager.index_VSM( SM_choice );
        ParamsManager.index_RMSH( 1 );
        
        % Initialize the directories depending on theta phi, VSM, and RMSH
        SimulationFolders.getInstance.initializeDynamicDirs();
        
        
        %% GET GLOBAL DIRECTORIES
        dir_config = SimulationFolders.getInstance.config;
        dir_freqdiff = SimulationFolders.getInstance.freqdiff;
        dir_out_diffuse_P1_tuple = SimulationFolders.getInstance.out_diffuse_P1_tuple;
        
        
        % read Ki
        filename1 = strcat('Ki') ;
        Ki = readComplexVar( dir_freqdiff, filename1) ;
        KKi_dB(tt, jj) = 10 * log10(abs(Ki) ^ 2 / 4 / pi) ;
        
        % Fresnel ellipses
        filenamex = 'ellipse_s_m' ;
        ellipse_s = readVar( dir_config, filenamex) ;
        area_s = pi * ellipse_s(:, 1) .* ellipse_s(:, 2) ;
        P1_areas_dB(tt, jj) = 10 * log10(area_s(FZ_choice)) ;
        P1_areas_dB2(:, tt, jj) = 10 * log10(area_s) ;
        
        % transmit port 1 / receiver ports 1&2
        xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t1_dB') ;
        PP1_inc_dB(1, :, tt, jj) = squeeze(xx(:, FZ_choice)) ;
        PP1_inc_dB2(1, :, :, tt, jj) = xx;
        
        % transmit port 2 / receiver ports 1&2
        xx = readVar(dir_out_diffuse_P1_tuple, 'PP1_inc_t2_dB') ;
        PP1_inc_dB(2, :, tt, jj) = squeeze(xx(:, FZ_choice)) ;
        PP1_inc_dB2(2, :, :, tt, jj) = xx;
        
    end

end


%% reading reflectivity

% Intilization
% over vegetation
P_cohveg_dB = zeros(2, 2, Nth0d, Nhr) ;
% over bare soil
P_cohbare_dB = zeros(2, 2, Nth0d, Nhr) ;

KKc_dB = zeros(Nth0d, Nhr) ;

for jj = 1 : Nhr
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    for tt = 1 : Nth0d

        % Set theta index
        ParamsManager.index_Th( tt );
        
        % TO-DO: Assign the correct index to each
        ParamsManager.index_Ph( 1 );
        ParamsManager.index_VSM( SM_choice );
        ParamsManager.index_RMSH( 1 );
        
        % Initialize the directories depending on theta phi, VSM, and RMSH
        SimulationFolders.getInstance.initializeDynamicDirs();
        
        
        %% GET GLOBAL DIRECTORIES
        dir_out_specular = SimulationFolders.getInstance.out_specular;
        
        % read Kc
        filename1 = strcat('Kc') ;
        Kc = readComplexVar( dir_out_specular, filename1) ;
        KKc_dB(tt, jj) = 20 * log10(abs(Kc)) ;
        
        [~, ~, ~, ~, ~, ~, ~, ~, ...
            P_coh1vegx, P_coh2vegx, ~, ~, ...
            P_coh1barex, P_coh2barex, ~, ~] = readSpecular ;
        
        P_cohveg_dB(1, :, tt, jj) = 10 * log10( P_coh1vegx(1 : 2) ) + KKc_dB(tt) ;
        P_cohveg_dB(2, :, tt, jj) = 10 * log10( P_coh2vegx(1 : 2) ) + KKc_dB(tt) ;
        
        P_cohbare_dB(1, :, tt, jj) = 10 * log10( P_coh1barex(1 : 2) ) + KKc_dB(tt) ;
        P_cohbare_dB(2, :, tt, jj) = 10 * log10( P_coh2barex(1 : 2) ) + KKc_dB(tt) ;
        
    end
    
end

%% Normalization
NBRCS_dB_TOT = 10 * log10(10.^(P_cohveg_dB/10) + 10.^(PP1_inc_dB/10)) ...
    - reshape(repmat(KKi_dB, 4, 1), 2, 2, Nth0d, Nhr) ...
    - reshape(repmat(P1_areas_dB, 4, 1), 2, 2, Nth0d, Nhr) ; 
NBRCS_dB = PP1_inc_dB ...
    - reshape(repmat(KKi_dB, 4, 1), 2, 2, Nth0d, Nhr) ...
    - reshape(repmat(P1_areas_dB, 4, 1), 2, 2, Nth0d, Nhr) ; 

for ii = 1 : Nfz
    NBRCS_dB_TOT2(:, :, ii, :, :) = 10 * log10(10.^(P_cohveg_dB/10) + 10.^(squeeze(PP1_inc_dB2(:, :, ii, :, :))/10)) ...
    - reshape(repmat(KKi_dB, 4, 1), 2, 2, Nth0d, Nhr) ...
    - reshape(repmat(P1_areas_dB, 4, 1), 2, 2, Nth0d, Nhr) ;

    NBRCS_dB3(:, :, ii, :, :) = squeeze(PP1_inc_dB2(:, :, ii, :, :)) ...
    - reshape(repmat(KKi_dB, 4, 1), 2, 2, Nth0d, Nhr) ...
    - reshape(repmat(P1_areas_dB, 4, 1), 2, 2, Nth0d, Nhr) ;
end

REFL_dB_TOT = 10 * log10(10.^(P_cohveg_dB/10) + 10.^(PP1_inc_dB/10)) - reshape(repmat(KKc_dB, 4, 1), 2, 2, Nth0d, Nhr) ; 
REFL_dB = P_cohveg_dB - reshape(repmat(KKc_dB, 4, 1), 2, 2, Nth0d, Nhr) ; 

for ii = 1 : Nfz
    REFL_dB_TOT2(:, :, ii, :, :) = 10 * log10(10.^(P_cohveg_dB/10) + 10.^(squeeze(PP1_inc_dB2(:, :, ii, :, :))/10)) - reshape(repmat(KKc_dB, 4, 1), 2, 2, Nth0d, Nhr) ;
end

%% Total NBRCS
for jj = 1 : Nhr
        
    NBRCS_dB_TOT_RR = squeeze(NBRCS_dB_TOT(1, 1, :, jj)) ;
    NBRCS_dB_RR = squeeze(NBRCS_dB(1, 1, :, jj)) ;
    
    if jj == 1
        subplot(1,2,1)
        hold
        axis([0 80 -20 30])
        xlabel('\theta_s [\circ]')
        ylabel('\sigma^0_e [dB]')
        xticks(th0d(1) : 10 : th0d(end))
        grid
        text(20, -15, '20m\leq h_r\leq 500m','Interpreter','tex')
        text(45, 20, strcat('fzone = 1'))
        title('Co-POL [RR]')
    end
    subplot(1,2,1)
    plot(th0d, NBRCS_dB_RR, '-sb', 'MarkerFaceColor', 'blue', 'MarkerSize', 2)
    
    subplot(1,2,1)
    plot(th0d, NBRCS_dB_TOT_RR, ':or', 'MarkerFaceColor', 'red', 'MarkerSize', 2)
    legend('Specular', 'Total')

end

%% Total NBRCS 2
for jj = 1 : Nfz
        
    NBRCS_dB_TOT_RR2 = squeeze(NBRCS_dB_TOT2(1, 1, jj, :, 4)) ;
    NBRCS_dB_RR2 = squeeze(NBRCS_dB3(1, 1, jj, :, 4)) ;
    
    
    if jj == 1
        
        subplot(1,2,2)
        hold
        
        axis([0 80 -20 30])
        xlabel('\theta_s [\circ]')
        ylabel('\sigma^0_e [dB]')
        xticks(th0d(1) : 10 : th0d(end))
        grid
        text(25, -15, '1\leq fzone\leq 10','Interpreter','tex')
        text(45, 20, strcat('h_r = 500m'))
        title('Co-POL [RR]')
        
    end
    subplot(1,2,2)
    plot(th0d, NBRCS_dB_RR2, '-sb', 'MarkerFaceColor', 'blue', 'MarkerSize', 2)
    
    subplot(1,2,2)
    plot(th0d, NBRCS_dB_TOT_RR2, ':or', 'MarkerFaceColor', 'red', 'MarkerSize', 2)
    legend('Specular', 'Total')

end


%% GET GLOBAL DIRECTORIES
dir_analysis = SimulationFolders.getInstance.analysis;
if ~exist(dir_analysis, 'dir')
    mkdir(dir_analysis)
end


% TO-DO: Check this
%% GET GLOBAL PARAMETERS
polT = SatParams.getInstance.polT;
polR = RecParams.getInstance.polR;
VSM_cm3cm3 = GndParams.getInstance.VSM_list_cm3cm3( ParamsManager.index_VSM );


fname1 = strcat('sig_RR_FZ_', num2str(FZ_choice), '-', polT, polR, '-VSM_', num2str( VSM_cm3cm3 )) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


%% Total Reflectivity
for jj = 1 : Nhr
    
    REFL_dB_TOT_RR = squeeze(REFL_dB_TOT(1, 1, :, jj)) ;
    REFL_dB_RR = squeeze(REFL_dB(1, 1, :, jj)) ;
    
    if jj == 1

        subplot(1,2,1)
        plot(th0d, REFL_dB_RR, '-sb', 'MarkerFaceColor', 'blue')
        hold
       
        axis([0 80 -50 0])
        xlabel('\theta_s [\circ]')
        ylabel('Reflectivity [dB]')
        xticks(th0d(1) : 10 : th0d(end))
        grid
        text(20, -45, '20m\leq h_r\leq 500m','Interpreter','tex')
        text(45, -10, strcat('fzone = 1'))
        title('Co-POL [RR]')
        
        x = [0.2 0.2];
        y = [0.5 0.3];
        annotation('textarrow',x,y)

    end
    
    subplot(1,2,1)
    plot(th0d, REFL_dB_TOT_RR, ':or', 'MarkerFaceColor', 'red', 'MarkerSize', 3)
    legend('Specular', 'Total')

end

%% Total Reflectivity 2
for jj = 1 : Nfz
    
    REFL_dB_TOT_RR2 = squeeze(REFL_dB_TOT2(1, 1, jj, :, 1)) ;
    REFL_dB_RR2 = squeeze(REFL_dB(1, 1, :, 1)) ;
    
    if jj == 1

        subplot(1,2,2)
        plot(th0d, REFL_dB_RR2, '-sb', 'MarkerFaceColor', 'blue')
        hold
       
        axis([0 80 -50 0])
        xlabel('\theta_s [\circ]')
        ylabel('Reflectivity [dB]')
        xticks(th0d(1) : 10 : th0d(end))
        grid
        text(25, -45, '1\leq fzone\leq 10','Interpreter','tex')
        text(45, -10, strcat('h_r = 20m'))
        title('Co-POL [RR]')
        
        x = [0.65 0.65];
        y = [0.4 0.6];
        annotation('textarrow',x,y)
        
    end
    
    subplot(1,2,2)
    plot(th0d, REFL_dB_TOT_RR2, ':or', 'MarkerFaceColor', 'red', 'MarkerSize', 3)
    legend('Specular', 'Total')

end


%%
fname1 = strcat('RR_FZ_', num2str(FZ_choice), '-', polT, polR, '-VSM_', num2str( VSM_cm3cm3 )) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


%% Total Power
for jj = 1 : Nhr
    
    inputFile_sys = strcat( inputFile_sys_tag, num2str( hhr(jj) ), '.xml' );
        
    getInput( inputFile_sys, inputFile_veg );
    % Initialize but not create simulations' directories
    SimulationFolders.getInstance.initializeStaticDirs();
    % [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();
    
    
    % TO-DO: Check this
    %% GET GLOBAL PARAMETERS
    hr_m = RecParams.getInstance.hr_m;
    
    
    Pvi_RR = squeeze(PP1_inc_dB(1, 1, :, jj)) ;
    Pvi_RL = squeeze(PP1_inc_dB(1, 2, :, jj)) ;    
    Pvc_RR = squeeze(P_cohveg_dB(1, 1, :, jj)) ;
    Pvc_RL = squeeze(P_cohveg_dB(1, 2, :, jj)) ;    
    Pbc_RR = squeeze(P_cohbare_dB(1, 1, :, jj)) ;
    Pbc_RL = squeeze(P_cohbare_dB(1, 2, :, jj)) ;    
    
    if jj == 1

        subplot(1,2,1)
        plot(th0d, Pbc_RL, '-sb', 'MarkerFaceColor', 'blue')
        hold
        plot(th0d, Pvc_RL, '-dg', 'MarkerFaceColor', 'green')
        
        axis([0 90 -220 -170])
        xlabel('\theta_s [\circ]')
        ylabel('Received Power [dB]')
        xticks(th0d(1) : 10 : th0d(end))
        grid
        text(72, Pvi_RL(end) + 3, strcat('h_r [m]'))
        title('X-Pol [RL]')

        subplot(1,2,2)
        plot(th0d, Pbc_RR, '-sb', 'MarkerFaceColor', 'blue')
        hold
        plot(th0d, Pvc_RR, '-dg', 'MarkerFaceColor', 'green')
        
        axis([0 90 -220 -170])
        xlabel('\theta_s [\circ]')
        xticks(th0d(1) : 10 : th0d(end))
        grid
        text(72, Pvi_RR(end) + 3, strcat('h_r [m]'))
        title('Co-Pol [RR]')

    end
    
    subplot(1,2,1)
    plot(th0d, Pvi_RL, ':or', 'MarkerFaceColor', 'red', 'markersize', 3)
    text(75, Pvi_RL(end), strcat( num2str(hr_m) ))

    subplot(1,2,2)
    plot(th0d, Pvi_RR, ':or', 'MarkerFaceColor', 'red', 'markersize', 3)
    text( 75, Pvi_RR(end), strcat( num2str(hr_m) ) )
    legend('Specular - Bare Soil', 'Specular - Vegetation', 'Diffuse - Vegetation')

end


%% GET GLOBAL DIRECTORIES
dir_analysis = SimulationFolders.getInstance.analysis;
if ~exist(dir_analysis, 'dir')
    mkdir(dir_analysis)
end


%%
fname1 = strcat('FZ_', num2str(FZ_choice), '-', polT, polR, '-VSM_', num2str( VSM_cm3cm3 ) ) ;
fname1 = strrep( fname1, '.', 'dot' );

saveas(gcf, strcat(dir_analysis, '\', fname1), 'tiff')
close


end



%% read results
function [b_coh1veg, b_coh2veg, b0_coh1veg, b0_coh2veg, ...
    b_coh1bare, b_coh2bare, b0_coh1bare, b0_coh2bare, ...
    P_coh1veg, P_coh2veg, P0_coh1veg, P0_coh2veg, ...
    P_coh1bare, P_coh2bare, P0_coh1bare, P0_coh2bare] = readSpecular


%% GET GLOBAL DIRECTORIES
dir_out_specular_tuple = SimulationFolders.getInstance.out_specular_tuple;


%% GET GLOBAL PARAMETERS
polT = SatParams.getInstance.polT;
polR = RecParams.getInstance.polR;


%% READ SAVED OUTPUT

% 2 X 2
filename1 = strcat('Veg1', polT, polR) ;
filename2 = strcat('Veg2', polT, polR) ;
filename01 = strcat('Veg01', polT, polR) ;
filename02 = strcat('Veg02', polT, polR) ;
b_coh1veg = readComplexVar(dir_out_specular_tuple, filename1) ;
b_coh2veg = readComplexVar(dir_out_specular_tuple, filename2) ;
b0_coh1veg = readComplexVar(dir_out_specular_tuple, filename01) ;
b0_coh2veg = readComplexVar(dir_out_specular_tuple, filename02) ;

filename1 = strcat('Bare1', polT, polR) ;
filename2 = strcat('Bare2', polT, polR) ;
filename01 = strcat('Bare01', polT, polR) ;
filename02 = strcat('Bare02', polT, polR) ;
b_coh1bare = readComplexVar(dir_out_specular_tuple, filename1) ;
b_coh2bare = readComplexVar(dir_out_specular_tuple, filename2) ;
b0_coh1bare = readComplexVar(dir_out_specular_tuple, filename01) ;
b0_coh2bare = readComplexVar(dir_out_specular_tuple, filename02) ;

% 4 X 4
filename1 = strcat('P_Veg1', polT, polR) ;
filename2 = strcat('P_Veg2', polT, polR) ;
filename01 = strcat('P_Veg01', polT, polR) ;
filename02 = strcat('P_Veg02', polT, polR) ;
P_coh1veg = readVar(dir_out_specular_tuple, filename1) ;
P_coh2veg = readVar(dir_out_specular_tuple, filename2) ;
P0_coh1veg = readVar(dir_out_specular_tuple, filename01) ;
P0_coh2veg = readVar(dir_out_specular_tuple, filename02) ;

filename1 = strcat('P_Bare1', polT, polR) ;
filename2 = strcat('P_Bare2', polT, polR) ;
filename01 = strcat('P_Bare01', polT, polR) ;
filename02 = strcat('P_Bare02', polT, polR) ;
P_coh1bare = readVar(dir_out_specular_tuple, filename1) ;
P_coh2bare = readVar(dir_out_specular_tuple, filename2) ;
P0_coh1bare = readVar(dir_out_specular_tuple, filename01) ;
P0_coh2bare = readVar(dir_out_specular_tuple, filename02) ;

end