% Orhan Eroglu
% July 6, 2018

function calcRxAntPatMatrixUserDefined

%% GET GLOBAL DIRECTORIES
dir_ant_lookup = SimulationFolders.getInstance.ant_lookup;


%% GET GLOBAL PARAMETERS
% Receiver Parameters (for User-defined (Custom) pattern)
ant_pat_fullfilename = RxUserDefinedParams.getInstance.ant_pat_fullfilename;

ant_pat_th_range_deg = Constants.ant_pat_th_range_deg;
ant_pat_th_range_rad = ant_pat_th_range_deg * Constants.deg2rad;
ant_pat_ph_range_deg = Constants.ant_pat_ph_range_deg;
ant_pat_ph_range_rad = ant_pat_ph_range_deg * Constants.deg2rad;


%% READ ANTENNA PATTERN FILE
[gnXX, ~, ~] = xlsread( ant_pat_fullfilename, 1 );

[gnXY, ~, ~] = xlsread( ant_pat_fullfilename, 2 );

[gnYX, ~, ~] = xlsread( ant_pat_fullfilename, 3 );

[gnYY, ~, ~] = xlsread( ant_pat_fullfilename, 4 );


% Get the numbers of azimuth (phi) and incidence (theta) angles of the antenna pattern
[Nph, Nth] = size(gnXX);

% Calculate the antenna pattern resolution in degrees
ant_pat_res_deg = ant_pat_th_range_deg / (Nth - 1);
% Set the Receiver Parameters' antenna pattern resolution property
RxParams.getInstance.set_ant_pat_res_deg( ant_pat_res_deg );

% Calculate theta and phi angle values for the antenna pattern 
th_rad = linspace(0, ant_pat_th_range_rad, Nth) ;
ph_rad = linspace(0, ant_pat_ph_range_rad, Nph) ;
% Construct meshgrids for theta and phi angles
[th, ph] =  meshgrid(th_rad, ph_rad);


%% ANTENNA PATTERN MATRIX ELEMENTS
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

% Antenna Pattern Matrix
G = cell(4) ;

G{1,1} = g11 ; G{1,2} = g12 ; G{1,3} = g13 ; G{1,4} = g14 ;
G{2,1} = g21 ; G{2,2} = g22 ; G{2,3} = g23 ; G{2,4} = g24 ;
G{3,1} = g31 ; G{3,2} = g32 ; G{3,3} = g33 ; G{3,4} = g34 ;
G{4,1} = g41 ; G{4,2} = g42 ; G{4,3} = g43 ; G{4,4} = g44 ;

%%
g = cell(2) ;

g{1,1} = gnXX ; g{1,2} = gnXY ;
g{2,1} = gnYX ; g{2,2} = gnYY ; 


%% SAVE
save([dir_ant_lookup '\AntPat.mat'], 'G', 'g', 'th', 'ph')

       
end