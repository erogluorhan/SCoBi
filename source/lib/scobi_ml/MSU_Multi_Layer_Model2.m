
%% Mehmet Kurum
% September 04, 2017


%%

% -----   inf  -------
% -------  air --------
% ----- layer 1 -------
% ----- layer 2 -------
% ..................
% ----- layer N -------
% ---     inf  -------


function MSU_Multi_Layer_Model2(vsm, ci_doy_SM)

%% GET GLOBAL PARAMETERS
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
f_Hz = f_MHz * Constants.MHz2Hz;
pol_Tx = TxParams.getInstance.pol_Tx;
%% RECEIVER PARAMETERS
pol_Rx = RxParams.getInstance.pol_Rx;
% Dynamic Parameters
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
num_Th = length(th0_Tx_list_deg);


%% System Parameters

% Operating Frequency [Hz]
lambda_m = Constants.c / f_Hz ;


%% Soil Parameters

% Layer Measurement Depth [in m]
L_depth = [5 10 20 40] * 1e-2 ;

% standard variation of layer boundaries: 10%
sL_depth = 0.1 * L_depth ;

% Texture at various depths (L_depth)
sand = [0.10 0.10 0.10 0.10]  ;
clay = [0.31 0.31 0.31 0.34] ;

% Soil bulk density at various depths (L_depth)
rho_b = [1.4 1.4 1.4 1.5] ;

% Mausred volumetric soil moisture at various depths (L_depth)
% vsm = [0.21 0.34 0.38 0.44]  ; % cm^3/cm^3
% vsm = [0.21 0.31 0.35 0.44]  ; % cm^3/cm^3
% vsm = [0.35 0.33 0.36 0.44]  ; % cm^3/cm^3

% number of layers
N_layer = length(rho_b) ;

% Soil Dielectric Constant
eps_diel_soil = round(10 * dielg(vsm, f_Hz, sand, clay, rho_b)) / 10 ;

% Air Dielectric Constant
eps_diel_air = 1.0 - 0.0 * 1i ;


% Roughness
lambda_cm = lambda_m * 1e2 ;       % in cm
ko = 2 * pi / lambda_cm ;
sig = 0.5 ;                   % in cm
rmsH = (2 * sig * ko) ^ 2 ;    % effective roughness parameter
h = rmsH ;

%% Layer Profile

% Layer bottom
L_bottom = [0 L_depth(1 : end - 1) + diff(L_depth) / 2 ] ;

% Layer thickness
L_thickness = diff(L_bottom) ;

% Layer discretization
Delz = 1e-3 ; % interval in m

% total thickness considered
zA = 10e-2 ;     % air layer
zB = 30e-2 ;    % the bottom most layer
zS = zA + L_bottom(end) + zB ; % total
z = (0 : Delz : zS)' ;

%% 2nd order polyfit
zz = zA + L_depth ;
p2 = polyfit(zz, eps_diel_soil, 2) ;
eps_diel_z2 = polyval(p2, z) ;
eps_diel_z2(z <= zA) = eps_diel_air ;

FigON = 1 ;
if FigON == 3
    figure(1)
    subplot(3,4,1)
    plot(real(eps_diel_z2), (z-zA)*1e2, 'k', 'linewidth', 2)
    hold on
    plot(real(eps_diel_soil), L_depth*1e2, 'o')
    grid
    axis([0 30 0 z(end)])
    set(gca,'YDir','reverse')
    set(gca,'YLim',[0-zA*1e2 z(end, :) * 1e2-zA*1e2])
    % set(gca,'YTick',[0 zz z(end)] * 1e2)
    % aa = sort([zz, z(end), zA + L_depth]) - zA ;
    aa = [0, L_depth] ;
    set(gca,'YTick',aa * 1e2)
    bb = real([eps_diel_air, eps_diel_soil]) ;
    bb = sort(unique(bb)) ;
    set(gca,'XTick',bb)
    xlabel('\epsilon\prime - real part')
    ylabel('z [cm]')
    title('2nd order polynomial')
    hold off
%     saveas(gcf, strcat(pwd, '\', 'Re2'), 'jpg')
%     close(gcf)
end

FigON = 0 ;
if FigON == 1
    figure(10)
    subplot(2,2,1)
    plot(real(eps_diel_z2), (z-zA)*1e2, 'k', 'linewidth', 2)
    hold on
    plot(real(eps_diel_soil), L_depth*1e2, 'o')
    grid
    axis([0 30 0 z(end)])
    set(gca,'YDir','reverse')
    set(gca,'YLim',[0-zA*1e2 z(end, :) * 1e2-zA*1e2])
    % set(gca,'YTick',[0 zz z(end)] * 1e2)
    % aa = sort([zz, z(end), zA + L_depth]) - zA ;
    aa = [0, L_depth] ;
    set(gca,'YTick',aa * 1e2, 'FontSize', 6)
    bb = real([eps_diel_air, eps_diel_soil]) ;
    bb = sort(unique(bb)) ;
    set(gca,'XTick',bb)
%     xlabel('\epsilon\prime - real part')
    ylabel('z [cm]')
    title('2nd order polynomial')
    hold off
%     saveas(gcf, strcat(pwd, '\', 'Re2'), 'jpg')
%     close(gcf)
end


%% 3rd order polyfit
zz = zA + L_depth ;
p3 = polyfit(zz, eps_diel_soil, 3) ;
eps_diel_z3 = polyval(p3, z) ;
eps_diel_z3(z <= zA) = eps_diel_air ;
eps_diel_z3(real(eps_diel_z3) < 1) = 1 ;

FigON = 1 ;
if FigON == 1
    figure(1)
    subplot(3,4,2)
    plot(real(eps_diel_z3), (z-zA)*1e2, 'c', 'linewidth', 2)
    hold on
    plot(real(eps_diel_soil), L_depth*1e2, 'o')
    grid
    axis([0 30 0 z(end)])
    set(gca,'YDir','reverse')
    set(gca,'YLim',[0-zA*1e2 z(end, :) * 1e2-zA*1e2])
    % set(gca,'YTick',[0 zz z(end)] * 1e2)
    % aa = sort([zz, z(end), zA + L_depth]) - zA ;
    aa = [0, L_depth] ;
    set(gca,'YTick',aa * 1e2)
    bb = real([eps_diel_air, eps_diel_soil]) ;
    bb = sort(unique(bb)) ;
    set(gca,'XTick',bb)
    xlabel('\epsilon\prime - real part')
    ylabel('z [cm]')
    title('3rd order polynomial')
    hold off

%     saveas(gcf, strcat(pwd, '\', 'Re3'), 'jpg')
%     close(gcf)
end

FigON = 0 ;
if FigON == 1
    figure(10)
    subplot(2,2,2)
    plot(real(eps_diel_z3), (z-zA)*1e2, 'c', 'linewidth', 2)
    hold on
    plot(real(eps_diel_soil), L_depth*1e2, 'o')
    grid
    axis([0 30 0 z(end)])
    set(gca,'YDir','reverse')
    set(gca,'YLim',[0-zA*1e2 z(end, :) * 1e2-zA*1e2])
    % set(gca,'YTick',[0 zz z(end)] * 1e2)
    % aa = sort([zz, z(end), zA + L_depth]) - zA ;
    aa = [0, L_depth] ;
    set(gca,'YTick',aa * 1e2, 'FontSize', 6)
    bb = real([eps_diel_air, eps_diel_soil]) ;
    bb = sort(unique(bb)) ;
    set(gca,'XTick',bb)
%     xlabel('\epsilon\prime - real part')
%     ylabel('z [cm]')
    title('3rd order plynomial')
    hold off

%     saveas(gcf, strcat(pwd, '\', 'Re3'), 'jpg')
%     close(gcf)
end


%% smoother transition
[z, eps_diel_z0] ...
    = DielPrf0(eps_diel_air, eps_diel_soil, L_thickness, L_depth, sL_depth, zA, z) ;

%% discrete slab

eps_diel_zs = eps_diel_z0 ; 
eps_diel_zs(z < zA) = eps_diel_air ; 

for ii = 1 : length(L_bottom)-1
    eps_diel_zs(z > (zA + L_bottom(ii)) & z < (zA + L_bottom(ii+1))) = eps_diel_soil(ii) ;
end

eps_diel_zs(z > (zA + L_bottom(end))) = eps_diel_soil(end) ; 


FigON = 1 ;
if FigON == 1
    figure(1)
    subplot(3,4,4)
    plot(real(eps_diel_zs), (z-zA)*1e2, 'r', 'linewidth', 2)
    hold on
    plot(real(eps_diel_soil), L_depth*1e2, 'o')
    grid
    axis([0 30 0 z(end)])
    set(gca,'YDir','reverse')
    set(gca,'YLim',[0-zA*1e2 z(end, :) * 1e2-zA*1e2])
    % set(gca,'YTick',[0 zz z(end)] * 1e2)
    % aa = sort([zz, z(end), zA + L_depth]) - zA ;
    aa = [0, L_depth] ;
    set(gca,'YTick',aa * 1e2)
    bb = real([eps_diel_air, eps_diel_soil]) ;
    bb = sort(unique(bb)) ;
    set(gca,'XTick',bb)
    xlabel('\epsilon\prime - real part')
    ylabel('z [cm]')
        title('Slab')
 hold off
 %     saveas(gcf, strcat(pwd, '\', 'Re3'), 'jpg')
%     close(gcf)
end

FigON = 0 ;
if FigON == 1
    figure(10)
    subplot(2,2,4)
    plot(real(eps_diel_zs), (z-zA)*1e2, 'r', 'linewidth', 2)
    hold on
    plot(real(eps_diel_soil), L_depth*1e2, 'o')
    grid
    axis([0 30 0 z(end)])
    set(gca,'YDir','reverse')
    set(gca,'YLim',[0-zA*1e2 z(end, :) * 1e2-zA*1e2])
    % set(gca,'YTick',[0 zz z(end)] * 1e2)
    % aa = sort([zz, z(end), zA + L_depth]) - zA ;
    aa = [0, L_depth] ;
    set(gca,'YTick',aa * 1e2, 'FontSize', 6)
    bb = real([eps_diel_air, eps_diel_soil]) ;
    bb = sort(unique(bb)) ;
    set(gca,'XTick',bb)
    xlabel('\epsilon\prime - real part')
%     ylabel('z [cm]')
        title('Slab')
 hold off
     saveas(gcf, strcat(pwd, '\', 'profiles'), 'jpg')
    close(gcf)
end

%% calculate reflection coefficients

%% Discrete Slab

zzb = (zA + L_bottom)' ;
Lzb = diff(zzb)' ; % / lambda_m ; % complex optical length in units of lambda_m
nA = sqrte(eps_diel_air) ;
nS = sqrte(eps_diel_soil) ;

na = [nA; nA; nA] ;
ns = [nS; nS; nS] ;

% input to multidiel
n = [na, ns] ;

Rp1 = zeros( num_Th, 1) ;
Rp2 = Rp1 ;

for ii = 1 : num_Th
        
    % Set theta index
    ParamsManager.index_Th( ii );
    
    % Initialize the directories depending on dynamic parameters
    SimulationFolders.getInstance.initializeDynamicDirs();
    
    %% GET GLOBAL PARAMETERS
    % Transmitter Parameters
    th0_Tx_deg = th0_Tx_list_deg( ParamsManager.index_Th );
                
    % Reflection Coefficient
    rh = multidiel(n, Lzb, 1, th0_Tx_deg, 'te') ;
    rv = multidiel(n, Lzb, 1, th0_Tx_deg, 'th') ;
    r_s = [rv 0; 0 rh] ;
    [r0_coh1b, r0_coh2b] = AntennaPolVeg(r_s) ;
    
    if pol_Tx == 'X' && pol_Rx == 'X'
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh2b(2)) .^ 2 ; % at ?/?0 = 1
    else
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh1b(2)) .^ 2 ; % at ?/?0 = 1
    end
    
end

figure(1)
subplot(3,4,5:8)
plot(ci_doy_SM, Rp1, '-or', 'MarkerFaceColor', 'r', 'markers', 3)
plot(ci_doy_SM, Rp2, '-or', 'MarkerFaceColor', 'w', 'markers', 3)
ylabel('Reflectivity')
title(strcat('Transmit/Receive Polarization:',  pol_Tx, pol_Rx))

figure(2)
subplot(2,1,1)
plot(ci_doy_SM, Rp1, '-or', 'MarkerFaceColor', 'r', 'markers', 1)
plot(ci_doy_SM, Rp2, '-or', 'MarkerFaceColor', 'w', 'markers', 1)
ylabel('Reflectivity')
% title(strcat('Transmit/Receive Polarization:',  pol_Tx, pol_Rx))

% FigON = 1 ;
% if FigON == 1
%     figure
%     plot(th0_Tx_list_deg, Rp1, '-or', th0_Tx_list_deg, Rp2, '-sb')
%     axis([0 90 0 1])
%     xlabel('Angle of Observation')
%     ylabel('Reflectivity')
%     title('Discrete Slab')
%     
%     saveas(gcf, strcat(pwd, '\', 'R_slab_', pol_Tx, pol_Rx), 'jpg')
%     close(gcf)
% end

%% calculate reflection coefficients


Rp1 = zeros(num_Th, 1) ;
Rp2 = Rp1 ;

for ii = 1 : num_Th
    
    
    % Set theta index
    ParamsManager.index_Th( ii );
    
    % Initialize the directories depending on dynamic parameters
    SimulationFolders.getInstance.initializeDynamicDirs();
    
    %% GET GLOBAL PARAMETERS
    % Transmitter Parameters
    th0_Tx_deg = th0_Tx_list_deg( ParamsManager.index_Th );
    
    [rh, rv] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z, eps_diel_z0, eps_diel_air) ;
    r_s = [rv 0; 0 rh] ;
    [r0_coh1b, r0_coh2b] = AntennaPolVeg(r_s) ;
    
    if pol_Tx == 'X' && pol_Rx == 'X'
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh2b(2)) .^ 2 ; % at ?/?0 = 1
    else
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh1b(2)) .^ 2 ; % at ?/?0 = 1
    end
    
    
end

figure(1)
subplot(3,4,5:8)
plot(ci_doy_SM, Rp1, '-ob', 'MarkerFaceColor', 'b', 'markers', 3)
plot(ci_doy_SM, Rp2, '-ob', 'MarkerFaceColor', 'w', 'markers', 3)
ylabel('Reflectivity')
title(strcat('Transmit/Receive Polarization:',  pol_Tx, pol_Rx))


figure(2)
subplot(2,1,1)
plot(ci_doy_SM, Rp1, '-ob', 'MarkerFaceColor', 'b', 'markers', 1)
plot(ci_doy_SM, Rp2, '-ob', 'MarkerFaceColor', 'w', 'markers', 1)
ylabel('Reflectivity')
% title(strcat('Transmit/Receive Polarization:',  pol_Tx, pol_Rx))


% FigON = 1 ;
% if FigON == 1
%     
%     figure
%     plot(th0_Tx_list_deg, Rp1, '-or', th0_Tx_list_deg, Rp2, '-sb')
%     axis([0 90 0 1])
%     xlabel('Angle of Observation')
%     ylabel('Reflectivity')
%     title('Logistic function')
%     saveas(gcf, strcat(pwd, '\', 'R_smooth_', pol_Tx, pol_Rx), 'jpg')
%     close(gcf)
% end


%%
Rp1 = zeros(num_Th, 1) ;
Rp2 = Rp1 ;

for ii = 1 : num_Th
    
    % Set theta index
    ParamsManager.index_Th( ii );
    
    % Initialize the directories depending on dynamic parameters
    SimulationFolders.getInstance.initializeDynamicDirs();
    
    %% GET GLOBAL PARAMETERS
    % Transmitter Parameters
    th0_Tx_deg = th0_Tx_list_deg( ParamsManager.index_Th );
        
    [rh, rv] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z, eps_diel_z2, eps_diel_air) ;
    r_s = [rv 0; 0 rh] ;
    [r0_coh1b, r0_coh2b] = AntennaPolVeg(r_s) ;
    
    if pol_Tx == 'X' && pol_Rx == 'X'
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh2b(2)) .^ 2 ; % at ?/?0 = 1
    else
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh1b(2)) .^ 2 ; % at ?/?0 = 1
    end
    
    
end

figure(1)
subplot(3,4,5:8)
plot(ci_doy_SM, Rp1, '-ok', 'MarkerFaceColor', 'k', 'markers', 3)
plot(ci_doy_SM, Rp2, '-ok', 'MarkerFaceColor', 'w', 'markers', 3)
ylabel('Reflectivity')
title(strcat('Transmit/Receive Polarization:',  pol_Tx, pol_Rx))

figure(2)
subplot(2,1,1)
plot(ci_doy_SM, Rp1, '-ok', 'MarkerFaceColor', 'k', 'markers', 1)
plot(ci_doy_SM, Rp2, '-ok', 'MarkerFaceColor', 'w', 'markers', 1)
ylabel('Reflectivity')
% title(strcat('Transmit/Receive Polarization:',  pol_Tx, pol_Rx))

% FigON = 1 ;
% if FigON == 1
%     
%     figure
%     plot(th0_Tx_list_deg, Rp1, '-or', th0_Tx_list_deg, Rp2, '-sb')
%     axis([0 90 0 1])
%     xlabel('Angle of Observation')
%     ylabel('Reflectivity')
%     title('polynomial 2nd order')
%     saveas(gcf, strcat(pwd, '\', 'R_2ndorder_', pol_Tx, pol_Rx), 'jpg')
%     close(gcf)
% end


%%
Rp1 = zeros(num_Th, 1) ;
Rp2 = Rp1 ;

for ii = 1 : num_Th
    
    % Set theta index
    ParamsManager.index_Th( ii );
    
    % Initialize the directories depending on dynamic parameters
    SimulationFolders.getInstance.initializeDynamicDirs();
    
    %% GET GLOBAL PARAMETERS
    % Transmitter Parameters
    th0_Tx_deg = th0_Tx_list_deg( ParamsManager.index_Th );
    
    [rh, rv] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z, eps_diel_z3, eps_diel_air) ;
    r_s = [rv 0; 0 rh] ;
    [r0_coh1b, r0_coh2b] = AntennaPolVeg(r_s) ;
    
    if pol_Tx == 'X' && pol_Rx == 'X'
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh2b(2)) .^ 2 ; % at ?/?0 = 1
    else
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh1b(2)) .^ 2 ; % at ?/?0 = 1
    end
    
    
end

figure(1)
subplot(3,4,5:8)
plot(ci_doy_SM, Rp1, '-oc', 'MarkerFaceColor', 'c', 'markers', 3)
plot(ci_doy_SM, Rp2, '-oc', 'MarkerFaceColor', 'w', 'markers', 3)
ylabel('Reflectivity')
title(strcat('Transmit/Receive Polarization:',  pol_Tx, pol_Rx))

figure(2)
subplot(2,1,1)
plot(ci_doy_SM, Rp1, '-oc', 'MarkerFaceColor', 'c', 'markers', 1)
plot(ci_doy_SM, Rp2, '-oc', 'MarkerFaceColor', 'w', 'markers', 1)
ylabel('Reflectivity')
title(strcat('Transmit/Receive Polarization: RV')) %,  pol_Tx, pol_Rx))


% FigON = 1 ;
% if FigON == 1
%     
%     figure
%     plot(th0_Tx_list_deg, Rp1, '-or', th0_Tx_list_deg, Rp2, '-sb')
%     axis([0 90 0 1])
%     xlabel('Angle of Observation')
%     ylabel('Reflectivity')
%     title('polynomial 3rd order')
%     saveas(gcf, strcat(pwd, '\', 'R_3rdorder_', pol_Tx, pol_Rx), 'jpg')
%     close(gcf)
% end


end

function [rh, rv] = Calc_RelectionCoef(lambda_m, th0_Tx_deg, z, eps_diel_z, eps_diel_air)


Lz = diff(z)' / lambda_m ; % complex optical length in units of lambda_m

nAz = sqrte(eps_diel_air) ;
nmz = sqrte(eps_diel_z(2 : end, :)) ;
nSz = sqrte(eps_diel_z(end, :)) ;

% Air - % isotropic
na = [nAz; nAz; nAz] ;

% Dielectric Profile : isotropic
nm = [nmz(:, 1).'; nmz(:, 1).'; nmz(:, 1).'] ;

% Soil - isotropic
nb = [nSz(:, 1); nSz(:, 1); nSz(:, 1)] ;

%% input to multidiel
n = [na, nm, nb] ;

%% Reflection Coeffficeint
rh = multidiel(n, Lz, 1, th0_Tx_deg, 'te') ;
rv = multidiel(n, Lz, 1, th0_Tx_deg, 'th') ;

end


%%  Dielectric Profiles

function [z, eps_diel_z] ...
    = DielPrf0(eps_diel_air, eps_diel_soil, L_thickness, L_depth, sL_depth, zA, z)


N_layer = length(sL_depth) ;

eps1 = eps_diel_air ;
eps_diel_z = eps1  ;
zz = zeros(1, N_layer) ;
for ii = 1 : N_layer
    
    sdel = sL_depth(ii) ;
    if ii == 1
        zz(ii) = zA ;
    else
        zz(ii) = zz(ii - 1) + L_thickness(ii - 1) ;
    end
    eps2 = eps_diel_soil(ii) ;
    
    eps_diel_z = eps_diel_z ...
        + (eps2 - eps1) ./ (1 + exp(-2.197 * (z - zz(ii)) / sdel)) ;
    
    eps1 = eps2 ;
    
end

% to exclude roughness
eps_diel_z(z <= zA) = eps_diel_air ;

FigON = 1 ;
if FigON == 1
    figure(1)
    subplot(3,4,3)
    plot(real(eps_diel_z), z*1e2-zA*1e2, 'linewidth', 2)
    hold on
    plot(real(eps_diel_soil), L_depth*1e2, 'o')
    grid
    axis([0 30 0 z(end)])
    set(gca,'YDir','reverse')
    set(gca,'YLim',[0-zA*1e2 z(end, :) * 1e2-zA*1e2])
    % set(gca,'YTick',[0 zz z(end)] * 1e2)
    % aa = sort([zz, z(end), zA + L_depth]) - zA ;
    aa = [0, L_depth] ;
    set(gca,'YTick',aa * 1e2)
    bb = real([eps_diel_air, eps_diel_soil]) ;
    bb = sort(unique(bb)) ;
    set(gca,'XTick',bb)
    xlabel('\epsilon\prime - real part')
    ylabel('z [cm]')
    hold off
    title('Logistic function fit')
%     saveas(gcf, strcat(pwd, '\', 'Re0'), 'jpg')
%     close(gcf)
end

FigON = 0 ;
if FigON == 1
    figure(10)
    subplot(2,2,3)
    plot(real(eps_diel_z), z*1e2-zA*1e2, 'linewidth', 2)
    hold on
    plot(real(eps_diel_soil), L_depth*1e2, 'o')
    grid
    axis([0 30 0 z(end)])
    set(gca,'YDir','reverse')
    set(gca,'YLim',[0-zA*1e2 z(end, :) * 1e2-zA*1e2])
    % set(gca,'YTick',[0 zz z(end)] * 1e2)
    % aa = sort([zz, z(end), zA + L_depth]) - zA ;
    aa = [0, L_depth] ;
    set(gca,'YTick',aa * 1e2, 'FontSize', 6)
    bb = real([eps_diel_air, eps_diel_soil]) ;
    bb = sort(unique(bb)) ;
    set(gca,'XTick',bb)
    xlabel('\epsilon\prime - real part')
    ylabel('z [cm]')
    hold off
    title('Logistic function fit')
%     saveas(gcf, strcat(pwd, '\', 'Re0'), 'jpg')
%     close(gcf)
end

end

