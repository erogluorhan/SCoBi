function plotStalkPos


%% ANALYSIS INPUTS
inputFile_sys_tag = 'sysInput-CORN_PAPER-row_';
inputFile_veg_tag = 'vegVirRowInput-Corn-row';
% inputFile_sys_tag = 'sysInput-CORN_PAPER-HOMOGENOUS';
% inputFile_veg_tag = 'vegHomInput-CORN_PAPER';

cornStages = { 'V1_V9', 'V10_VT', 'R1_R4', 'R5', 'R6' };

% Corn field row azimuth angles
rowPhis = [0 30 45 60 90];


%% Choices
% TO-DO: Check here
rowPhi_choice = 5;
TH_choice = 7;
PH_choice = 1;
SM_choice = 1 ;
RMSH_choice = 1 ;
stage_choice = 3;


%% GET INPUT
% Just for getting the initial values of some parameters
inputFile_sys_ext = strcat( num2str( rowPhis(rowPhi_choice) ), '.xml' );
inputFile_veg_ext = strcat( num2str( rowPhis(rowPhi_choice) ), '-', cornStages{stage_choice}, '.xml' );
% inputFile_sys_ext = '-Tricky-TH70-PH0_45_90.xml';
% inputFile_veg_ext = strcat( '-', cornStages{stage_choice}, '.xml' );

inputFile_sys = strcat( inputFile_sys_tag, inputFile_sys_ext );
inputFile_veg = strcat( inputFile_veg_tag, inputFile_veg_ext );

getInput( inputFile_sys, inputFile_veg );

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();
% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();   
    
% Set theta index
ParamsManager.index_Th( TH_choice );
ParamsManager.index_VSM( SM_choice );
ParamsManager.index_Ph( PH_choice );
ParamsManager.index_RMSH( RMSH_choice );

% Initialize the directories depending on theta phi, VSM, and RMSH
SimulationFolders.getInstance.initializeDynamicDirs();


%% GET GLOBAL DIRECTORIES
dir_position = SimulationFolders.getInstance.position;
dir_fzones = SimulationFolders.getInstance.fzones;
dir_config = SimulationFolders.getInstance.config;


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nfz = SimParams.getInstance.Nfz;
% Bistatic Parameters
Tgs = BistaticParams.getInstance.Tgs;   % Transformation Gnd -> Specular
AntRotZ_Tx = BistaticParams.getInstance.AntRotZ_Tx;   % Transmitter Rotation along Z-Axis


%% READ META-DATA
% Scattering Particles (Only stalk)
fileName = 'R1_L1_T3_K1';
pP = readVar( dir_position, fileName );
Npart = readVar( dir_fzones, fileName );
% Fresnel zones
fileName = 'ellipses_FZ_m' ;
ellipses_FZ_m = readVar(dir_config, fileName) ;
fileName = 'ellipse_centers_FZ_m' ;
ellipse_centers_FZ_m = readVar(dir_config, fileName) ;


%% FIGURE
figure
axis equal
hold on

%Draw Stalk Positions
N2 = 1 ;
for fz = 1 : Nfz

    color = 'r.';
    if mod(fz,2) == 0
        color = 'b.';
    end
    N1 = N2 ;
    N2 = Npart(fz) ;
    if N2 > N1
        scatter(pP( N1:N2, 1 ), pP( N1:N2, 2 ), color);
    elseif N2 == 0
        N2 = 1;
    end

end


% Draw fresnel zones
majorAxis = 2 * ellipses_FZ_m(:, 1);
minorAxis = 2 * ellipses_FZ_m(:, 2); 

fz_thetas = -pi : 0.01 : pi;

lenR = length(fz_thetas);
tempR = zeros(3,lenR);

for ii = 1 : Nfz

    tempR(1,:) = majorAxis(ii) / 2 * cos(fz_thetas);
    tempR(2,:) = minorAxis(ii) / 2 * sin(fz_thetas);
    tempR = AntRotZ_Tx * tempR;
    fz_x = ellipse_centers_FZ_m(1,ii) + tempR(1,:);
    fz_y = ellipse_centers_FZ_m(2,ii) + tempR(2,:);

    fz_sf = Tgs * [fz_x; fz_y; zeros(1, length(fz_x))];

    plot(fz_x, fz_y, 'r--')
    %plot(fz_sf(1,:), fz_sf(2,:), 'r--')
end
xlabel('East ->')
ylabel('North ->')
% 
% %% OE: End of figures for testing


end

