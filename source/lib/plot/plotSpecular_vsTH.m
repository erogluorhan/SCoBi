% Mehmet Kurum
% modified 11/13/2017
% Reflectivity vs. Observation Angle

function plotSpecular_vsTH


%% INPUT FILES
% inputFile_sys = 'sysInput-Paulownia-PAPER_PBAND-hr_20.xml';
% inputFile_veg = 'vegHomInput-Paulownia.xml';
inputFile_sys = 'sysInput-CORN_PAPER-row_0.xml';
inputFile_veg = 'vegVirRowInput-Corn.xml';


%% GET INPUT FOR INITIAL START
getInput( inputFile_sys, inputFile_veg );

% Initialize but not create simulations' directories
SimulationFolders.getInstance.initializeStaticDirs();

% [isInputValid, isTerminate, terminateMsg] = ParamsManager.isInputValid();


%% GET GLOBAL PARAMETERS
% Satellite Parameters
th0_list_deg = SatParams.getInstance.th0_list_deg;
PH0_list_deg = SatParams.getInstance.PH0_list_deg;
polT = SatParams.getInstance.polT;
% Receiver Parameters
polR = RecParams.getInstance.polR;
% Ground Parameters
VSM_list_cm3cm3 = GndParams.getInstance.VSM_list_cm3cm3;
RMSH_list_cm = GndParams.getInstance.RMSH_list_cm;

num_Th = length( th0_list_deg );
num_Ph = length( PH0_list_deg );
num_VSM = length( VSM_list_cm3cm3 );
num_RMSH = length( RMSH_list_cm );

% over vegetation
P_coh1veg = zeros(2, num_Th) ;    
P_coh2veg = zeros(2, num_Th) ; 
% over bare soil
P_coh1bare = zeros(2, num_Th) ; 
P_coh2bare = zeros(2, num_Th) ; 


%% Reading ...
if SimSettings.getInstance.sim_mode == Constants.sim_mode.SNAPSHOT

    % For each phi (azimuth angle)
    for pp = 1 : num_Ph

        % Set phi index
        ParamsManager.index_Ph( pp );

        % For each VSM (volumetric soil moisture)
        for vv = 1 : num_VSM

            % Set VSM index
            ParamsManager.index_VSM( vv );

            % For each RMSH (root mean square height - roughness)
            for rr = 1 : num_RMSH

                % Set RMSH index
                ParamsManager.index_RMSH( rr );

                
                KKc_dB = [] ;

                % For each theta (looking angle)
                for tt = 1 : num_Th 

                    % Set theta index
                    ParamsManager.index_Th( tt );

                    % Initialize the directories depending on theta,
                    % phi, VSM, and RMSH
                    SimulationFolders.getInstance.initializeDynamicDirs();


                    %% GET GLOBAL DIRECTORIES
                    dir_out_specular = SimulationFolders.getInstance.out_specular;
                    
                    dir_fig_specular_P_vsTH = SimulationFolders.getInstance.fig_specular_P_vsTH;
                    if ~exist(dir_fig_specular_P_vsTH, 'dir')
                        mkdir(dir_fig_specular_P_vsTH)
                    end
                    

                    % read Kc
                    filename1 = 'Kc';
                    Kc = readComplexVar(dir_out_specular, filename1) ;
                    KKc_dB = [KKc_dB, 20 * log10(abs(Kc))] ;

                    [~, ~, ~, ~, ~, ~, ~, ~, ...
                    P_coh1vegx, P_coh2vegx, ~, ~, ...
                    P_coh1barex, P_coh2barex, ~, ~] = readSpecular;

                    P_coh1veg(:, tt) = P_coh1vegx(1 : 2) ; %#ok<*AGROW>
                    P_coh2veg(:, tt) = P_coh2vegx(1 : 2) ;
                    P_coh1bare(:, tt) = P_coh1barex(1 : 2);
                    P_coh2bare(:, tt) = P_coh2barex(1 : 2) ;

                end
                
                
                %%
                FigON = 1 ;
                if FigON == 1
                %% Received Power in dB

                    figure

                    if ( polT == 'X') && (polT == 'X' )

                        bareVV_dB = KKc_dB + 10 * log10(squeeze(P_coh1bare(1, vv, :))) ;
                        bareHH_dB = KKc_dB + 10 * log10(squeeze(P_coh2bare(2, vv, :))) ;
                        vegVV_dB = KKc_dB + 10 * log10(squeeze(P_coh1veg(1, vv, :))) ;
                        vegHH_dB = KKc_dB + 10 * log10(squeeze(P_coh2veg(2, vv, :))) ;

                        plot(th0_list_deg, bareVV_dB, ':or') % co-pol VV
                        hold
                        plot(th0_list_deg, bareHH_dB, '-or') % co-pol HH
                        plot(th0_list_deg, vegVV_dB, ':sg') % co-pol VV
                        plot(th0_list_deg, vegHH_dB, '-sg') % co-pol HH
                    
                    else
                        
                        BareCO_dB = KKc_dB + 10 * log10( squeeze(P_coh1bare(1, :)) ) ;
                        BareX_dB = KKc_dB + 10 * log10( squeeze(P_coh1bare(2, :)) ) ;
                        VegCO_dB = KKc_dB + 10 * log10( squeeze(P_coh1veg(1, :)) ) ;
                        VegX_dB = KKc_dB + 10 * log10( squeeze(P_coh1veg(2, :)) ) ;

                        plot(th0_list_deg, BareCO_dB, ':or') % co-pol
                        hold
                        plot(th0_list_deg, BareX_dB, '-or') % x-pol
                        plot(th0_list_deg, VegCO_dB, ':sg') % co-pol
                        plot(th0_list_deg, VegX_dB, '-sg') % x-pol
                    end

                    grid
                    xlabel('Incidence Angle [degrees]')
                    ylabel('Received Power [dB]')

                    if (polT == 'R') && (polR == 'R')
                        legend('bare: RR-pol', 'bare: RL-pol', 'veg: RR-pol', 'veg: RL-pol', 'Location', 'southeast')
                        title(strcat('Tx: RHCP, Rx: RHCP', '- VSM:', num2str( VSM_list_cm3cm3(vv) ), '%'))
                    elseif (polT == 'R') && (polR == 'X')
                        legend('bare: RX-pol', 'bare: RY-pol', 'veg: RX-pol', 'veg: RY-pol', 'Location', 'southeast')
                        title(strcat('Tx: RHCP, Rx: LINEAR', '- VSM:', num2str( VSM_list_cm3cm3(vv) ), '%'))
                    elseif (polT == 'X') && (polT == 'X')
                        legend('bare: XX-pol', 'bare: YY-pol', 'veg: XX-pol', 'veg: YY-pol', 'Location', 'southeast')
                        title(strcat('Tx: LINEAR, Rx: LINEAR', '- VSM:', num2str( VSM_list_cm3cm3(vv) ), '%'))
                    end

                    axis([0 80 -250 -180])

                    fname = strcat('ReceivedPower_dB', '-PH_', num2str( PH0_list_deg(pp) ), '-VSM_', num2str( VSM_list_cm3cm3(vv) ), '-RMSH_', num2str( RMSH_list_cm(rr) ) );
                    fname = strrep( fname, '.', 'dot' );
                    saveas( gcf, strcat(dir_fig_specular_P_vsTH, '\', fname ), 'tiff')
                    close
                    
                end

                % Reflectivity in dB 
                figure

                if (polT == 'X') && (polT == 'X')

                    bareVV_dB = 10 * log10(squeeze(P_coh1bare(1, :))) ;
                    bareHH_dB = 10 * log10(squeeze(P_coh2bare(2, :))) ;
                    vegVV_dB = 10 * log10(squeeze(P_coh1veg(1, :))) ;
                    vegHH_dB = 10 * log10(squeeze(P_coh2veg(2, :))) ;

                    plot(th0_list_deg, bareVV_dB, ':or') % co-pol VV
                    hold
                    plot(th0_list_deg, bareHH_dB, '-or') % co-pol HH
                    plot(th0_list_deg, vegVV_dB, ':sg') % co-pol VV
                    plot(th0_list_deg, vegHH_dB, '-sg') % co-pol HH
                else
                    BareCO_dB = 10 * log10(squeeze(P_coh1bare(1, :))) ;
                    BareX_dB = 10 * log10(squeeze(P_coh1bare(2, :))) ;
                    VegCO_dB = 10 * log10(squeeze(P_coh1veg(1, :))) ;
                    VegX_dB = 10 * log10(squeeze(P_coh1veg(2, :))) ;

                    plot(th0_list_deg, BareCO_dB, ':or') % co-pol
                    hold
                    plot(th0_list_deg, BareX_dB, '-or') % x-pol
                    plot(th0_list_deg, VegCO_dB, ':sg') % co-pol
                    plot(th0_list_deg, VegX_dB, '-sg') % x-pol
                end

                grid
                xlabel('Incidence Angle [degrees]')
                ylabel('Reflectivity [dB]')

                if (polT == 'R') && (polR == 'R')
                    legend('bare: RR-pol', 'bare: RL-pol', 'veg: RR-pol', 'veg: RL-pol', 'Location', 'southeast')
                    title(strcat('No Ks - Tx: RHCP, Rx: RHCP', '- VSM:', num2str( VSM_list_cm3cm3(vv) ), '%'))
                elseif (polT == 'R') && (polR == 'X')
                    legend('bare: RX-pol', 'bare: RY-pol', 'veg: RX-pol', 'veg: RY-pol', 'Location', 'southeast')
                    title(strcat('No Ks - Tx: RHCP, Rx: LINEAR', '- VSM:', num2str( VSM_list_cm3cm3(vv) ), '%'))
                elseif (polT == 'X') && (polT == 'X')
                    legend('bare: XX-pol', 'bare: YY-pol', 'veg: XX-pol', 'veg: YY-pol', 'Location', 'southeast')
                    title(strcat('No Ks - Tx: LINEAR, Rx: LINEAR', '- VSM:', num2str( VSM_list_cm3cm3(vv) ), '%'))
                end

                axis([0 80 -50 0])

                fname = strcat('Reflectivity_dB', '-PH_', num2str( PH0_list_deg(pp) ), '-VSM_', num2str( VSM_list_cm3cm3(vv) ), '-RMSH_', num2str( RMSH_list_cm(rr) ) );
                fname = strrep( fname, '.', 'dot' );
                saveas(gcf, strcat(dir_fig_specular_P_vsTH, '\',fname), 'tiff')
                close

            end

        end

    end
    
% Time-series simulation
else
    
    % For each corresponding tuple of theta (looking angle), 
    % phi (azimuth angle), VSM (volumetric soil moisture), and 
    % RMSH (root mean square height - roughness)
    for ii = 1 : num_Th  % The length of each is equal

        % Set theta, phi, VSM, ad RMSH index the same
        ParamsManager.index_Th( ii );
        ParamsManager.index_Ph( ii );
        ParamsManager.index_VSM( ii );
        ParamsManager.index_RMSH( ii );

        % Initialize the directories depending on theta, phi, VSM, and 
        % RMSH
        SimulationFolders.getInstance.initializeDynamicDirs();

        %% TO-DO: PLOT

    end
    
end

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