% Mehmet Kurm
% March 2, 2017

function calcRxAntPatMatrixGG

%% GET GLOBAL DIRECTORIES
dir_ant_lookup = SimulationFolders.getInstance.ant_lookup;


%% GET GLOBAL PARAMETERS
% Receiver Parameters (for Generalized-Gaussian pattern)
SLL_dB = RxGGParams.getInstance.SLL_dB;
XPL_dB = RxGGParams.getInstance.XPL_dB;
hpbw_deg = RxGGParams.getInstance.hpbw_deg;


%% CALCULATIONS
% TO-DO: Make a continous range
if SLL_dB == 15 ;
    a = 0.18 ; % 0.35 ;  % 15 dB    
elseif SLL_dB == 20 ;
    a = 0.11 ; % 0.25 ;  % 20 dB    
elseif SLL_dB == 30 ;
    a = 0.04 ; % 0.12 ;  % 30 dB    
elseif SLL_dB == 40 ;
    a = 0.01 ; % 0.055 ; % 40 dB    
end


% X-pol level
if XPL_dB == 15
    V = sqrt(0.0316) ; % -15 dB
elseif XPL_dB == 25
    V = sqrt(0.0100);  % -20 dB
elseif XPL_dB == 30
    V = sqrt(0.0010) ;  % -30 dB
elseif XPL_dB == 40
    V = sqrt(0.0001) ; % -40 dB
end

[th, ph, gg] = GGpattern( hpbw_deg, a) ;
gXX = gg ; gYY = gg ;

[~, ~, gg] = GGpattern( 2 * hpbw_deg, a ) ;
gXY = V * gg ; gYX = V * gg ;
    


%% complex (voltage) co- and x- patterns : X-PORT - v-pol
magg = sqrt(abs(gXX) .^2 + abs(gXY) .^2) ;
maxg = max(max(magg)) ;

% complex normalized (voltage) pattern
gnXX = gXX / maxg ;  
gnXY = gXY / maxg ;

%% complex (voltage) co- and x- patterns : Y-PORT - h-pol
magg = sqrt(abs(gYX) .^2 + abs(gYY) .^2) ;
maxg = max(max(magg)) ;

% complex normalized (voltage) pattern
gnYX = gYX / maxg ;  
gnYY = gYY / maxg ;

%% Antenna Pattern Matrix Elements

g11 = abs(gnXX) .^ 2 ;
g12 = abs(gnXY) .^ 2 ;
g13 = real(gnXX .* conj(gnXY)) ;
g14 = -imag(gnXX .* conj(gnXY)) ;

g21 = abs(gnYX) .^ 2 ;
g22 = abs(gnYY) .^ 2 ;
g23 = real(gnYX .* conj(gnYY)) ;
g24 = -imag(gnYX .* conj(gnYY)) ;

g31 = 2 * real(gnXX .* conj(gnYX)) ; 
g32 = 2 * real(gnXY .* conj(gnYY)) ; 
g33 = real(gnXX .* conj(gnYY) + gnXY .* conj(gnYX)) ;
g34 = -imag(gnXX .* conj(gnYY) - gnXY .* conj(gnYX)) ;

g41 = 2 * imag(gnXX .* conj(gnYX)) ; 
g42 = 2 * imag(gnXY .* conj(gnYY)) ; 
g43 = imag(gnXX .* conj(gnYY) + gnXY .* conj(gnYX)) ;
g44 = real(gnXX .* conj(gnYY) - gnXY .* conj(gnYX)) ;

%% Antenna Pattern Matrix

G = cell(4) ;

G{1,1} = g11 ; G{1,2} = g12 ; G{1,3} = g13 ; G{1,4} = g14 ;
G{2,1} = g21 ; G{2,2} = g22 ; G{2,3} = g23 ; G{2,4} = g24 ;
G{3,1} = g31 ; G{3,2} = g32 ; G{3,3} = g33 ; G{3,4} = g34 ;
G{4,1} = g41 ; G{4,2} = g42 ; G{4,3} = g43 ; G{4,4} = g44 ;

%%

g = cell(2) ;

g{1,1} = gnXX ; g{1,2} = gnXY ;
g{2,1} = gnYX ; g{2,2} = gnYY ; 

excelFileName = "default_ant_pat";

excelFileName = strcat( Directories.getInstance.input_sys, '\', excelFileName);

%% Half-Power Beanwidth

% X-port
Gn_co = G{1, 1} ;
Gn_x = G{1, 2} ;

Gn = Gn_co + Gn_x ; % co- + x- pols

indhpbw = (Gn >= 0.49) & (Gn <= 0.51) ;
hpbwX = mean(mean(th(indhpbw))) ;

del = 0 ;
while isnan(hpbwX)
    indhpbw = (Gn >= (0.50 - del)) & (Gn <= (0.50 + del)) ;
    hpbwX = mean(mean(th(indhpbw))) ;
    del = del + 0.01 ;
end

% Y-port
Gn_co = G{2, 2} ;
Gn_x = G{2, 1} ;

Gn = Gn_co + Gn_x ; % co- + x- pols

indhpbw = (Gn >= 0.49) & (Gn <= 0.51) ;
hpbwY = mean(mean(th(indhpbw))) ;

del = 0 ;
while isnan(hpbwY)
    indhpbw = (Gn >= (0.50 - del)) & (Gn <= (0.50 + del)) ;
    hpbwY = mean(mean(th(indhpbw))) ;
    del = del + 0.01 ;
end

hpbw2 = hpbwX + hpbwY ;  % Main beam

hpbw2 = hpbw2 * 180 / pi ;


%% SAVE
save([dir_ant_lookup '\AntPat.mat'], 'G', 'g', 'th', 'ph')

       
end


function [th, ph, gg] = GGpattern(ths_deg, a)

%% GET GLOBAL PARAMETERS
ant_pat_res_deg = RxParams.getInstance.ant_pat_res_deg;

ant_pat_th_range_deg = Constants.ant_pat_th_range_deg;
ant_pat_th_range_rad = ant_pat_th_range_deg * Constants.deg2rad;
ant_pat_ph_range_deg = Constants.ant_pat_ph_range_deg;
ant_pat_ph_range_rad = ant_pat_ph_range_deg * Constants.deg2rad;

% Beamwidth
% ths_deg = 12, 6, 3

% Sidelobe levels
% a = 0.35 (15 dB), 0.25 (20 dB), 0.12 (30 dB), 0.055 (40 db)

% TO-DO: Test for bad values (e.g. non-integer result)?
Nth = floor( ant_pat_th_range_deg / ant_pat_res_deg ) + 1;
Nph = floor( ant_pat_ph_range_deg / ant_pat_res_deg ) + 1;

th_rad = linspace(0, ant_pat_th_range_rad, Nth) ;
ph_rad = linspace(0, ant_pat_ph_range_rad, Nph) ;

[th, ph] =  meshgrid(th_rad, ph_rad);

% Angle span (Beamwidth)
% ths_deg = 12 ;                     
ths_rad = degtorad(ths_deg) ; 

% Generalized Gaussian Pattern parameters
alpha = 0.2 ; % 0.5 ;  % sidelobe width                 
% a = 0.35 ;  % sidelobe level

% Generelized Gaussian
gg = abs(1 / (1 - a) * exp( -(tan(th) / tan(ths_rad)) .^ 2)...
        - a / (1 - a) * exp(-(alpha * tan(th) / tan(ths_rad)) .^ 2)) ;

% To eliminate the angles that are over 0.4pi degrees from main beam
XX = min( min( gg( :, (th_rad > 0.4*pi) & (th_rad < 0.45*pi) ) ) ) ;
gg(:, th_rad > pi/2) = XX ;

gg(gg < XX) = XX ;

end