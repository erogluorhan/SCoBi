classdef SimulationFolders < handle
    % SIMULATIONFOLDERS Class to keep track of all output directories
    %   This class has one attribute for each source code or input
    %   directory. Every attribute can be reached by a static getter method.
    
    
properties (SetAccess = private, GetAccess = public)

simulations

analysis

hr
gnd
veg
sim_mode

th0_deg
afsa
fscat
geo
position
fzones
distance
incidence
scattering
observation

pol
config
ant
ant_real

rot
rot_lookup
rot_real

freq
freqdiff
freqdiff_b1
freqdiff_b1_tuple
freqdiff_b2
freqdiff_b2_tuple
freqdiff_b3
freqdiff_b3_tuple
freqdiff_b4
freqdiff_b4_tuple
freqdiff_P1
freqdiff_P1_tuple
freqdiff_P2
freqdiff_P2_tuple
freqdiff_P3
freqdiff_P3_tuple
freqdiff_P4
freqdiff_P4_tuple

out
out_direct
out_specular
out_specular_tuple
out_diffuse
out_diffuse_b1
out_diffuse_b2
out_diffuse_b3
out_diffuse_b4
out_diffuse_P1
out_diffuse_P1_tuple
out_diffuse_P2
out_diffuse_P2_tuple
out_diffuse_P3
out_diffuse_P3_tuple
out_diffuse_P4
out_diffuse_P4_tuple
out_diffuse_NBRCS
out_diffuse_NBRCS_tuple

out_ml_ref

fig
fig_direct
fig_specular
fig_specular_P
fig_specular_P_vsTH
fig_diffuse
fig_diffuse_P1_vsTH
fig_diffuse_NBRCS_vsTH
fig_diffuse_P1_vsFZ
fig_diffuse_NBRCS_vsFZ

end
    
    
methods (Access = private)

function obj = SimulationFolders
end

end
    
    
methods (Static)
    
function singleObj = getInstance
persistent localObj

if isempty(localObj) || ~isvalid(localObj)
   localObj = SimulationFolders;
end

 singleObj = localObj;
end
end
        
    
methods
        
function initializeStaticDirs(obj)


%% GET GLOBAL DIRECTORIES
main_dir = Directories.getInstance.main_dir;


%% GET GLOBAL PARAMETERS
% Simulation Settings
simulator_id = SimSettings.getInstance.simulator_id;
simulator = Constants.simulators{1, simulator_id};
sim_mode_id = SimSettings.getInstance.sim_mode_id;
sim_mode_str = Constants.sim_modes{1, sim_mode_id};
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;
gnd_cover = Constants.gnd_covers{1, gnd_cover_id};
% Simulation Parameters
sim_name = SimParams.getInstance.sim_name;
campaign = SimParams.getInstance.campaign;
campaign_date = SimParams.getInstance.campaign_date;
plot = SimParams.getInstance.plot;
% Receiver Parameters
hr_m = RxParams.getInstance.hr_m;


%% INITIALIZE DIRECTORIES
% Initialize static simulation directories
simulation_dir = strcat( main_dir, '\sims\', simulator);
simulation_dir = strcat( simulation_dir, '\', gnd_cover);

% If gnd_cover is Vegetation, then add the veg_method and
% vegetation_plant
if gnd_cover_id == Constants.id_veg_cover


    %% GET GLOBAL PARAMETERS
    %Simulation Parameters
    veg_method_id = SimParams.getInstance.veg_method_id;
    veg_method = Constants.veg_methods{1, veg_method_id};
    vegetation_plant = SimParams.getInstance.vegetation_plant;
    % Vegetation Parameters
    if veg_method_id == Constants.id_veg_vir

        vegetation_stage = VegVirRowParams.getInstance.vegetation_stage;

        veg_vir_orientation_id = SimParams.getInstance.veg_vir_orientation_id;
        veg_vir_orientation = Constants.veg_vir_orientations{1, veg_vir_orientation_id};

    end

    simulation_dir = strcat( simulation_dir, '\', vegetation_plant);
    simulation_dir = strcat( simulation_dir, '\', campaign, ...
                                '-', campaign_date, ...
                                '-plot_', plot );

    % If the simulator is SCoBi-Veg, then add veg_method and
    % related information to the folder system
    if simulator_id == Constants.id_veg_agr ...
            || simulator_id == Constants.id_veg_for

        simulation_dir = strcat( simulation_dir, '\', veg_method);

        % If veg_method is Virtual
        if veg_method_id == Constants.id_veg_vir

            simulation_dir = strcat( simulation_dir, '\', veg_vir_orientation);

            simulation_dir = strcat( simulation_dir, '\', vegetation_stage);

        end

    end

    simulation_dir = strcat( simulation_dir, '\', sim_name);

end

obj.simulations = simulation_dir;

%% Analysis
obj.analysis = strcat( obj.simulations, '\', 'Analysis') ;  

%% Receiver Height
hr_folder = strcat( '\hr_', num2str(hr_m) );
obj.hr = strcat( obj.simulations, hr_folder) ;  

%% Average Forward Scattering Amplitude
obj.afsa = strcat(obj.hr, '\', 'AFSA') ;

if gnd_cover_id == Constants.id_veg_cover
    %% Vegetation
    obj.veg = strcat(obj.hr, '\', 'VEG') ;
end

%% Simulation Mode
obj.sim_mode = strcat( obj.hr, '\', sim_mode_str );

end
        
        
        
function initializeDynamicDirs(obj)
% Initialize dynamically changing simulation directories


%% GET GLOBAL PARAMETERS
sim_counter = ParamsManager.sim_counter;
% Simulation Settings
simulator_id = SimSettings.getInstance.simulator_id;
% Dynamic Parameters
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;
VSM_cm3cm3 = VSM_list_cm3cm3( sim_counter );
RMSH_list_cm = DynParams.getInstance.RMSH_list_cm;
RMSH_cm = RMSH_list_cm( sim_counter );
th0_Tx_list_deg = DynParams.getInstance.th0_Tx_list_deg;
ph0_Tx_list_deg = DynParams.getInstance.ph0_Tx_list_deg;
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
% Receiver Parameters
pol_Rx = RxParams.getInstance.pol_Rx;
ant_pat_Rx_id = RxParams.getInstance.ant_pat_Rx_id;

if ant_pat_Rx_id == Constants.id_Rx_GG

    hpbw_deg = RxGGParams.getInstance.hpbw_deg;
    SLL_dB = RxGGParams.getInstance.SLL_dB;
    XPL_dB = RxGGParams.getInstance.XPL_dB;

elseif ant_pat_Rx_id == Constants.id_Rx_cos_pow_n

    % TO-DO: Should be implemented when cos-to-the-power-n
    % pattern is added

elseif ant_pat_Rx_id == Constants.id_Rx_user_defined

    ant_pat_fullfilename = RxUserDefinedParams.getInstance.ant_pat_fullfilename;

end



%% Ground
obj.gnd = strcat(obj.hr, '\', 'GND', '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ; 


%% Angle of incidence
th0d_folder = strcat('th0d_', num2str( th0_Tx_list_deg( sim_counter ) ), ...
                     '-ph0d_', num2str( ph0_Tx_list_deg( sim_counter ) ) ) ;
obj.th0_deg = strcat(obj.sim_mode, '\', th0d_folder) ;


% TO-DO: Determine for different simulators and settings
%% Bistatic Scattering Amplitudes
obj.fscat = strcat(obj.th0_deg, '\', 'FSCAT') ;


%% Geometry
obj.geo = strcat(obj.th0_deg, '\', 'GEO') ;
obj.position = strcat(obj.geo, '\', 'POSITION') ;                        
obj.fzones = strcat(obj.geo, '\', 'FZONES') ;            
obj.distance = strcat(obj.geo, '\', 'DISTANCE') ;            
obj.incidence = strcat(obj.geo, '\', 'INCIDENCE') ;            
obj.scattering = strcat(obj.geo, '\', 'SCATTERING') ;            
obj.observation = strcat(obj.geo, '\', 'OBSERVATION') ;


%% Polarization
obj.pol = strcat(obj.th0_deg, '\', pol_Tx, pol_Rx) ;

%Configuration
obj.config = strcat(obj.th0_deg, '\', 'CONFIG') ;


%% Antenna
FolderNameAnt = [];
if ant_pat_Rx_id == Constants.id_Rx_GG

    FolderNameAnt = strcat('AntPat-GG-thsd_', num2str( hpbw_deg ),...
        '-SLL_', num2str( SLL_dB ), '-XPL_', num2str( XPL_dB )) ;

elseif ant_pat_Rx_id == Constants.id_Rx_cos_pow_n

    % TO-DO: Should be implemented when cos-to-the-power-n
    % pattern is added

elseif ant_pat_Rx_id == Constants.id_Rx_user_defined

%                 FolderNameAnt = strcat('User_defined-', ant_pat_fullfilename ) ;
    FolderNameAnt = 'AntPat-User_defined';

end

obj.ant = strcat(obj.th0_deg, '\', 'ANT\', FolderNameAnt) ;     
obj.ant_real = strcat(obj.ant, '\', 'REALIZATION') ;


%% Rotation
obj.rot = strcat(obj.pol, '\', 'ROT') ;            
obj.rot_lookup = strcat(obj.rot, '\', 'LOOK-UP') ;            
obj.rot_real = strcat(obj.rot, '\', 'REALIZATION') ; 


%% Output
obj.out = strcat(obj.pol, '\OUTPUT\', FolderNameAnt) ;


%% Figure
obj.fig = strcat(obj.sim_mode, '\', 'FIGURE\', pol_Tx, pol_Rx, '\', FolderNameAnt ) ;

% SCoBi-Veg specific folders
if simulator_id == Constants.id_veg_agr ...
    || simulator_id == Constants.id_veg_for


    %% Frequency Response
    obj.freq = strcat(obj.pol, '\FREQ\', FolderNameAnt) ;            
    obj.freqdiff = strcat(obj.freq, '\', 'DIFFUSE') ;            
    obj.freqdiff_b1 = strcat(obj.freqdiff, '\', 'b1') ;          
    obj.freqdiff_b1_tuple = strcat(obj.freqdiff_b1, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;            
    obj.freqdiff_b2 = strcat(obj.freqdiff, '\', 'b2') ; 
    obj.freqdiff_b2_tuple = strcat(obj.freqdiff_b2, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;             
    obj.freqdiff_b3 = strcat(obj.freqdiff, '\', 'b3') ; 
    obj.freqdiff_b3_tuple = strcat(obj.freqdiff_b3, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;             
    obj.freqdiff_b4 = strcat(obj.freqdiff, '\', 'b4') ;  
    obj.freqdiff_b4_tuple = strcat(obj.freqdiff_b4, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;            
    obj.freqdiff_P1 = strcat(obj.freqdiff, '\', 'P1') ;  
    obj.freqdiff_P1_tuple = strcat(obj.freqdiff_P1, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;            
    obj.freqdiff_P2 = strcat(obj.freqdiff, '\', 'P2') ; 
    obj.freqdiff_P2_tuple = strcat(obj.freqdiff_P2, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;             
    obj.freqdiff_P3 = strcat(obj.freqdiff, '\', 'P3') ; 
    obj.freqdiff_P3_tuple = strcat(obj.freqdiff_P3, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ;             
    obj.freqdiff_P4 = strcat(obj.freqdiff, '\', 'P4') ;
    obj.freqdiff_P4_tuple = strcat(obj.freqdiff_P4, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ; 

    obj.out_direct = strcat(obj.out, '\', 'DIRECT') ;
    obj.out_specular = strcat(obj.out, '\', 'SPECULAR') ;
    obj.out_specular_tuple = strcat(obj.out_specular, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ; 
    obj.out_diffuse = strcat(obj.out, '\', 'DIFFUSE') ;
    obj.out_diffuse_b1 = strcat(obj.out_diffuse, '\', 'b1') ;
    obj.out_diffuse_b2 = strcat(obj.out_diffuse, '\', 'b2') ;
    obj.out_diffuse_b3 = strcat(obj.out_diffuse, '\', 'b3') ;
    obj.out_diffuse_b4 = strcat(obj.out_diffuse, '\', 'b4') ;
    obj.out_diffuse_P1 = strcat(obj.out_diffuse, '\', 'P1') ;
    obj.out_diffuse_P1_tuple = strcat(obj.out_diffuse_P1, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ; 
    obj.out_diffuse_P2 = strcat(obj.out_diffuse, '\', 'P2') ;
    obj.out_diffuse_P2_tuple = strcat(obj.out_diffuse_P2, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ; 
    obj.out_diffuse_P3 = strcat(obj.out_diffuse, '\', 'P3') ;
    obj.out_diffuse_P3_tuple = strcat(obj.out_diffuse_P3, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ; 
    obj.out_diffuse_P4 = strcat(obj.out_diffuse, '\', 'P4') ;
    obj.out_diffuse_P4_tuple = strcat(obj.out_diffuse_P4, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ; 
    obj.out_diffuse_NBRCS = strcat(obj.out_diffuse, '\', 'NBRCS') ;
    obj.out_diffuse_NBRCS_tuple = strcat(obj.out_diffuse_NBRCS, '\VSM_', num2str( VSM_cm3cm3 ), '-RMSH_', num2str(RMSH_cm)) ; 

    obj.fig_direct = strcat(obj.fig, '\', 'DIRECT') ;
    obj.fig_specular = strcat(obj.fig, '\', 'SPECULAR') ;
    obj.fig_specular_P = strcat(obj.fig_specular, '\', 'P') ;
    obj.fig_specular_P_vsTH = strcat(obj.fig_specular_P, '\', 'vs_TH') ;
    obj.fig_diffuse = strcat(obj.fig, '\', 'DIFFUSE') ;
    obj.fig_diffuse_P1_vsTH = strcat(obj.fig_diffuse, '\', 'P1', '\', 'vs_TH');
    obj.fig_diffuse_P1_vsFZ = strcat(obj.fig_diffuse, '\', 'P1', '\', 'vs_FZ');
    obj.fig_diffuse_NBRCS_vsTH = strcat(obj.fig_diffuse, '\', 'NBRCS', '\', 'vs_TH');
    obj.fig_diffuse_NBRCS_vsFZ = strcat(obj.fig_diffuse, '\', 'NBRCS', '\', 'vs_FZ');

% SCoBi-ML specific folders
elseif simulator_id == Constants.id_multi_layer

    obj.out_ml_ref = strcat(obj.out, '\', 'Reflectivity') ;

end

obj.makeDynamicDirs();

end
        
        
        function makeStaticDirs(obj)
            
            %% GET GLOBAL PARAMETERS
            % Simulation Parameters
            veg_method_id = SimParams.getInstance.veg_method_id;
        
            
            %% Analysis
            if ~exist(obj.analysis, 'dir')
                mkdir(obj.analysis);
            end
        
            %% Receiver Height
            if ~exist(obj.hr, 'dir')
                mkdir(obj.hr);
            end

            %% Average Forward Scattering Amplitude
            if ~exist(obj.afsa, 'dir')
                mkdir(obj.afsa)
            end
            
            %% Vegetation
            if ~exist(obj.veg, 'dir')
                mkdir(obj.veg)
            end   
            
            %% Simulation Mode
            if ~exist(obj.sim_mode, 'dir')
                mkdir(obj.sim_mode)
            end
            
            % TO-DO
            if veg_method_id == Constants.id_veg_hom
                VegParams.getInstance.write;
            end
        
        end
        
        
        function makeDynamicDirs(obj)
            

            %% Ground
            if ~exist(obj.gnd, 'dir')
                mkdir(obj.gnd)
            end
            
            %% Angle of incidence
            if ~exist(obj.th0_deg, 'dir')
                mkdir(obj.th0_deg)
            end

            %% Bistatic Scattering Amplitudes
            if ~exist(obj.fscat, 'dir')
                mkdir(obj.fscat)
            end

            %% Geometry
            if ~exist(obj.geo, 'dir')
                mkdir(obj.geo)
            end

            if ~exist(obj.position, 'dir')
                mkdir(obj.position)
            end

            if ~exist(obj.fzones, 'dir')
                mkdir(obj.fzones)
            end

            if ~exist(obj.distance, 'dir')
                mkdir(obj.distance)
            end

            if ~exist(obj.incidence, 'dir')
                mkdir(obj.incidence)
            end

            if ~exist(obj.scattering, 'dir')
                mkdir(obj.scattering)
            end

            if ~exist(obj.observation, 'dir')
                mkdir(obj.observation)
            end
            
            %% Polarization
            if ~exist(obj.pol, 'dir')
                mkdir(obj.pol)
            end

            %% Configuration
            if ~exist(obj.config, 'dir')
                mkdir(obj.config)
            end

            %% Antenna
            if ~exist(obj.ant, 'dir')
                mkdir(obj.ant)
            end

            if ~exist(obj.ant_real, 'dir')
                mkdir(obj.ant_real)
            end

            %% Rotation
            if ~exist(obj.rot, 'dir')
                mkdir(obj.rot)
            end

            if ~exist(obj.rot_lookup, 'dir')
                mkdir(obj.rot_lookup)
            end

            if ~exist(obj.rot_real, 'dir')
                mkdir(obj.rot_real)
            end

            %% Frequency Response
            if ~exist(obj.freq, 'dir')
                mkdir(obj.freq)
            end

            if ~exist(obj.freqdiff, 'dir')
                mkdir(obj.freqdiff)
            end

            if ~exist(obj.freqdiff_b1, 'dir')
                mkdir(obj.freqdiff_b1)
            end

            if ~exist(obj.freqdiff_b2, 'dir')
                mkdir(obj.freqdiff_b2)
            end

            if ~exist(obj.freqdiff_b3, 'dir')
                mkdir(obj.freqdiff_b3)
            end

            if ~exist(obj.freqdiff_b4, 'dir')
                mkdir(obj.freqdiff_b4)
            end

            if ~exist(obj.freqdiff_P1, 'dir')
                mkdir(obj.freqdiff_P1)
            end

            if ~exist(obj.freqdiff_P2, 'dir')
                mkdir(obj.freqdiff_P2)
            end

            if ~exist(obj.freqdiff_P3, 'dir')
                mkdir(obj.freqdiff_P3)
            end

            if ~exist(obj.freqdiff_P4, 'dir')
                mkdir(obj.freqdiff_P4)
            end

            %% Output
            if ~exist(obj.out, 'dir')
                mkdir(obj.out)
            end

            if ~exist(obj.out_direct, 'dir')
                mkdir(obj.out_direct)
            end

            if ~exist(obj.out_specular, 'dir')
                mkdir(obj.out_specular)
            end

            if ~exist(obj.out_diffuse, 'dir')
                mkdir(obj.out_diffuse)
            end

            if ~exist(obj.out_diffuse_b1, 'dir')
                mkdir(obj.out_diffuse_b1)
            end

            if ~exist(obj.out_diffuse_b2, 'dir')
                mkdir(obj.out_diffuse_b2)
            end

            if ~exist(obj.out_diffuse_b3, 'dir')
                mkdir(obj.out_diffuse_b3)
            end

            if ~exist(obj.out_diffuse_b4, 'dir')
                mkdir(obj.out_diffuse_b4)
            end

            if ~exist(obj.out_diffuse_P1, 'dir')
                mkdir(obj.out_diffuse_P1)
            end

            if ~exist(obj.out_diffuse_P2, 'dir')
                mkdir(obj.out_diffuse_P2)
            end
            if ~exist(obj.out_diffuse_P3, 'dir')
                mkdir(obj.out_diffuse_P3)
            end

            if ~exist(obj.out_diffuse_P4, 'dir')
                mkdir(obj.out_diffuse_P4)
            end

            if ~exist(obj.out_diffuse_NBRCS, 'dir')
                mkdir(obj.out_diffuse_NBRCS)
            end

            % Multilayer Reflectivity Output
            if ~exist(obj.out_ml_ref, 'dir')
                mkdir(obj.out_ml_ref)
            end

            % Figure
            if ~exist(obj.fig, 'dir')
                mkdir(obj.fig)
            end

            if ~exist(obj.fig_direct, 'dir')
                mkdir(obj.fig_direct)
            end

            if ~exist(obj.fig_specular, 'dir')
                mkdir(obj.fig_specular)
            end

            if ~exist(obj.fig_specular_P, 'dir')
                mkdir(obj.fig_specular_P)
            end

            if ~exist(obj.fig_specular_P_vsTH, 'dir')
                mkdir(obj.fig_specular_P_vsTH)
            end

            if ~exist(obj.fig_diffuse, 'dir')
                mkdir(obj.fig_diffuse)
            end

            if ~exist(obj.fig_diffuse_P1_vsTH, 'dir')
                mkdir(obj.fig_diffuse_P1_vsTH)
            end

            if ~exist(obj.fig_diffuse_P1_vsFZ, 'dir')
                mkdir(obj.fig_diffuse_P1_vsFZ)
            end

            if ~exist(obj.fig_diffuse_NBRCS_vsTH, 'dir')
                mkdir(obj.fig_diffuse_NBRCS_vsTH)
            end

            if ~exist(obj.fig_diffuse_NBRCS_vsFZ, 'dir')
                mkdir(obj.fig_diffuse_NBRCS_vsFZ)
            end
        
        end
        
        function out = get.simulations(obj)
            out = obj.simulations;        
        end
        
        function out = get.analysis(obj)
            out = obj.analysis;
        end 
        
        function out = get.hr(obj)
            out = obj.hr;
        end 
        
        function out = get.veg(obj)
            out = obj.veg;
        end 
        
        function out = get.th0_deg(obj)
            out = obj.th0_deg;
        end 
        
        function out = get.afsa(obj)
            out = obj.afsa;
        end 
        
        function out = get.fscat(obj)
            out = obj.fscat;
        end 
        
        function out = get.position(obj)
            out = obj.position;
        end 
        
        function out = get.fzones(obj)
            out = obj.fzones;
        end 
        
        function out = get.distance(obj)
            out = obj.distance;
        end 
        
        function out = get.incidence(obj)
            out = obj.incidence;
        end 
        
        function out = get.scattering(obj)
            out = obj.scattering;
        end 
        
        function out = get.observation(obj)
            out = obj.observation;
        end 
        
        function out = get.gnd(obj)
            out = obj.gnd;
        end 
        
        function out = get.pol(obj)
            out = obj.pol;
        end 
        
        function out = get.config(obj)
            out = obj.config;
        end 
        
        function out = get.ant(obj)
            out = obj.ant;
        end 
        
        function out = get.ant_real(obj)
            out = obj.ant_real;
        end
        
        function out = get.rot(obj)
            out = obj.rot;
        end
        
        function out = get.rot_lookup(obj)
            out = obj.rot_lookup;
        end
        
        function out = get.rot_real(obj)
            out = obj.rot_real;
        end
        
        function out = get.freq(obj)
            out = obj.freq;
        end
        
        function out = get.freqdiff(obj)
            out = obj.freqdiff;
        end
        
        function out = get.freqdiff_b1(obj)
            out = obj.freqdiff_b1;
        end
        
        function out = get.freqdiff_b1_tuple(obj)
            out = obj.freqdiff_b1_tuple;
        end
        
        function out = get.freqdiff_b2(obj)
            out = obj.freqdiff_b2;
        end
        
        function out = get.freqdiff_b2_tuple(obj)
            out = obj.freqdiff_b2_tuple;
        end
        
        function out = get.freqdiff_b3(obj)
            out = obj.freqdiff_b3;
        end
        
        function out = get.freqdiff_b3_tuple(obj)
            out = obj.freqdiff_b3_tuple;
        end
        
        function out = get.freqdiff_b4(obj)
            out = obj.freqdiff_b4;
        end
        
        function out = get.freqdiff_b4_tuple(obj)
            out = obj.freqdiff_b4_tuple;
        end
        
        function out = get.freqdiff_P1(obj)
            out = obj.freqdiff_P1;
        end
        
        function out = get.freqdiff_P1_tuple(obj)
            out = obj.freqdiff_P1_tuple;
        end
        
        function out = get.freqdiff_P2(obj)
            out = obj.freqdiff_P2;
        end
        
        function out = get.freqdiff_P2_tuple(obj)
            out = obj.freqdiff_P2_tuple;
        end
        
        function out = get.freqdiff_P3(obj)
            out = obj.freqdiff_P3;
        end
        
        function out = get.freqdiff_P3_tuple(obj)
            out = obj.freqdiff_P3_tuple;
        end
        
        function out = get.freqdiff_P4(obj)
            out = obj.freqdiff_P4;
        end
        
        function out = get.freqdiff_P4_tuple(obj)
            out = obj.freqdiff_P4_tuple;
        end
        
        function out = get.out(obj)
            out = obj.out;
        end
        
        function out = get.out_direct(obj)
            out = obj.out_direct;
        end
        
        function out = get.out_specular(obj)
            out = obj.out_specular;
        end
        
        function out = get.out_diffuse(obj)
            out = obj.out_diffuse;
        end
        
        function out = get.out_diffuse_b1(obj)
            out = obj.out_diffuse_b1;
        end
        
        function out = get.out_diffuse_b2(obj)
            out = obj.out_diffuse_b2;
        end
        
        function out = get.out_diffuse_b3(obj)
            out = obj.out_diffuse_b3;
        end
        
        function out = get.out_diffuse_b4(obj)
            out = obj.out_diffuse_b4;
        end
        
        function out = get.out_diffuse_P1(obj)
            out = obj.out_diffuse_P1;
        end
        
        function out = get.out_diffuse_P1_tuple(obj)
            out = obj.out_diffuse_P1_tuple;
        end
        
        function out = get.out_diffuse_P2(obj)
            out = obj.out_diffuse_P2;
        end
        
        function out = get.out_diffuse_P2_tuple(obj)
            out = obj.out_diffuse_P2_tuple;
        end
        
        function out = get.out_diffuse_P3(obj)
            out = obj.out_diffuse_P3;
        end
        
        function out = get.out_diffuse_P3_tuple(obj)
            out = obj.out_diffuse_P3_tuple;
        end
        
        function out = get.out_diffuse_P4(obj)
            out = obj.out_diffuse_P4;
        end
        
        function out = get.out_diffuse_P4_tuple(obj)
            out = obj.out_diffuse_P4_tuple;
        end
        
        function out = get.out_diffuse_NBRCS(obj)
            out = obj.out_diffuse_NBRCS;
        end
        
        function out = get.out_diffuse_NBRCS_tuple(obj)
            out = obj.out_diffuse_NBRCS_tuple;
        end
        
        function out = get.out_ml_ref(obj)
            out = obj.out_ml_ref;
        end
        
        function out = get.fig(obj)
            out = obj.fig;
        end
        
        function out = get.fig_direct(obj)
            out = obj.fig_direct;
        end
        
        function out = get.fig_specular(obj)
            out = obj.fig_specular;
        end
        
        function out = get.fig_specular_P(obj)
            out = obj.fig_specular_P;
        end
        
        function out = get.fig_specular_P_vsTH(obj)
            out = obj.fig_specular_P_vsTH;
        end
        
        function out = get.fig_diffuse(obj)
            out = obj.fig_diffuse;
        end
        
        function out = get.fig_diffuse_P1_vsTH(obj)
            out = obj.fig_diffuse_P1_vsTH;
        end
        
        function out = get.fig_diffuse_P1_vsFZ(obj)
            out = obj.fig_diffuse_P1_vsFZ;
        end
        
        function out = get.fig_diffuse_NBRCS_vsTH(obj)
            out = obj.fig_diffuse_NBRCS_vsTH;
        end
        
        function out = get.fig_diffuse_NBRCS_vsFZ(obj)
            out = obj.fig_diffuse_NBRCS_vsFZ;
        end
        
    end
    
end

