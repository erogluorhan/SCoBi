
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


function MSU_Multi_Layer_Model

%% Clear
clear ; clc ;

%% global
global th0d polT polR
% th0d = 40 ;

%% System Parameters

% Operating Frequency [Hz]
FHZ = 1575.42e6 ;
c = 3e8 ;
lambda0 = c / FHZ ;

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
vsm = [0.21 0.34 0.38 0.44]  ; % cm^3/cm^3
vsm = [0.21 0.31 0.35 0.44]  ; % cm^3/cm^3
% vsm = [0.35 0.33 0.36 0.44]  ; % cm^3/cm^3

% number of layers
N_layer = length(rho_b) ;

% Soil Dielectric Constant
eps_diel_soil = round(10 * dielg(vsm, FHZ, sand, clay, rho_b)) / 10 ;

% Air Dielectric Constant
eps_diel_air = 1.0 - 0.0 * 1i ;


% Roughness
lambda = lambda0 * 1e2 ;       % in cm
ko = 2 * pi / lambda ;
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
if FigON == 1
    figure
    plot(real(eps_diel_z2), (z-zA)*1e2, 'linewidth', 2)
    hold
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
    saveas(gcf, strcat(pwd, '\', 'Re2'), 'jpg')
    close(gcf)
end
%% 3rd order polyfit
zz = zA + L_depth ;
p3 = polyfit(zz, eps_diel_soil, 3) ;
eps_diel_z3 = polyval(p3, z) ;
eps_diel_z3(z <= zA) = eps_diel_air ;
eps_diel_z3(real(eps_diel_z3) < 1) = 1 ;

FigON = 1 ;
if FigON == 1
    figure
    plot(real(eps_diel_z3), (z-zA)*1e2, 'linewidth', 2)
    hold
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
    saveas(gcf, strcat(pwd, '\', 'Re3'), 'jpg')
    close(gcf)
end


%% smoother transition
[z, eps_diel_z0] ...
    = DielPrf0(eps_diel_air, eps_diel_soil, L_thickness, L_depth, sL_depth, zA, z) ;


%% Discrete Slab

zzb = (zA + L_bottom)' ;
Lzb = diff(zzb)' ; % / lambda0 ; % complex optical length in units of lambda0
nA = sqrte(eps_diel_air) ;
nS = sqrte(eps_diel_soil) ;

na = [nA; nA; nA] ;
ns = [nS; nS; nS] ;

% input to multidiel
n = [na, ns] ;

THD0 = 10 : 10 : 80 ;
Rp1 = zeros(length(THD0), 1) ;
Rp2 = Rp1 ;

for ii = 1 : length(THD0)
    th0d = THD0(ii) ;
    % Reflection Coefficient
    rh = multidiel(n, Lzb, 1, th0d, 'te') ;
    rv = multidiel(n, Lzb, 1, th0d, 'th') ;
    r_s = [rv 0; 0 rh] ;
    [r0_coh1b, r0_coh2b] = AntennaPolVeg(r_s) ;
    
    if polT == 'X' && polR == 'X'
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh2b(2)) .^ 2 ; % at ?/?0 = 1
    else
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh1b(2)) .^ 2 ; % at ?/?0 = 1
    end
    
end

FigON = 1 ;
if FigON == 1
    figure
    plot(THD0, Rp1, '-or', THD0, Rp2, '-sb')
    axis([0 90 0 1])
    xlabel('Angle of Observation')
    ylabel('Reflectivity')
    title('Discrete Slab')
    
    saveas(gcf, strcat(pwd, '\', 'R_slab_', polT, polR), 'jpg')
    close(gcf)
end

%% calculate reflection coefficients


THD0 = 10 : 10 : 80 ;
Rp1 = zeros(length(THD0), 1) ;
Rp2 = Rp1 ;

for ii = 1 : length(THD0)
    
    [rh, rv] = Calc_RelectionCoef(lambda0, THD0(ii), z, eps_diel_z0, eps_diel_air) ;
    r_s = [rv 0; 0 rh] ;
    [r0_coh1b, r0_coh2b] = AntennaPolVeg(r_s) ;
    
    if polT == 'X' && polR == 'X'
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh2b(2)) .^ 2 ; % at ?/?0 = 1
    else
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh1b(2)) .^ 2 ; % at ?/?0 = 1
    end
    
    
end

FigON = 1 ;
if FigON == 1
    
    figure
    plot(THD0, Rp1, '-or', THD0, Rp2, '-sb')
    axis([0 90 0 1])
    xlabel('Angle of Observation')
    ylabel('Reflectivity')
    title('Logistic function')
    saveas(gcf, strcat(pwd, '\', 'R_smooth_', polT, polR), 'jpg')
    close(gcf)
end

%%
THD0 = 10 : 10 : 80 ;
Rp1 = zeros(length(THD0), 1) ;
Rp2 = Rp1 ;

for ii = 1 : length(THD0)
    
    [rh, rv] = Calc_RelectionCoef(lambda0, THD0(ii), z, eps_diel_z2, eps_diel_air) ;
    r_s = [rv 0; 0 rh] ;
    [r0_coh1b, r0_coh2b] = AntennaPolVeg(r_s) ;
    
    if polT == 'X' && polR == 'X'
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh2b(2)) .^ 2 ; % at ?/?0 = 1
    else
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh1b(2)) .^ 2 ; % at ?/?0 = 1
    end
    
    
end

FigON = 1 ;
if FigON == 1
    
    figure
    plot(THD0, Rp1, '-or', THD0, Rp2, '-sb')
    axis([0 90 0 1])
    xlabel('Angle of Observation')
    ylabel('Reflectivity')
    title('polynomial 2nd order')
    saveas(gcf, strcat(pwd, '\', 'R_2ndorder_', polT, polR), 'jpg')
    close(gcf)
end

%%
THD0 = 10 : 10 : 80 ;
Rp1 = zeros(length(THD0), 1) ;
Rp2 = Rp1 ;

for ii = 1 : length(THD0)
    
    [rh, rv] = Calc_RelectionCoef(lambda0, THD0(ii), z, eps_diel_z3, eps_diel_air) ;
    r_s = [rv 0; 0 rh] ;
    [r0_coh1b, r0_coh2b] = AntennaPolVeg(r_s) ;
    
    if polT == 'X' && polR == 'X'
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh2b(2)) .^ 2 ; % at ?/?0 = 1
    else
        % Reflectivity
        Rp1(ii) = abs(r0_coh1b(1)) .^ 2 ; % at ?/?0 = 1
        Rp2(ii) = abs(r0_coh1b(2)) .^ 2 ; % at ?/?0 = 1
    end
    
    
end

FigON = 1 ;
if FigON == 1
    
    figure
    plot(THD0, Rp1, '-or', THD0, Rp2, '-sb')
    axis([0 90 0 1])
    xlabel('Angle of Observation')
    ylabel('Reflectivity')
    title('polynomial 3rd order')
    saveas(gcf, strcat(pwd, '\', 'R_3rdorder_', polT, polR), 'jpg')
    close(gcf)
end


end

function [rh, rv] = Calc_RelectionCoef(lambda0, THD0, z, eps_diel_z, eps_diel_air)


Lz = diff(z)' / lambda0 ; % complex optical length in units of lambda0

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
rh = multidiel(n, Lz, 1, THD0, 'te') ;
rv = multidiel(n, Lz, 1, THD0, 'th') ;

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

FigON = 1 ;
if FigON == 1
    figure
    plot(real(eps_diel_z), z*1e2-zA*1e2, 'linewidth', 2)
    hold
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
    saveas(gcf, strcat(pwd, '\', 'Re0'), 'jpg')
    close(gcf)
end

end

