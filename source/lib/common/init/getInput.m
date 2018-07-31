function getInput(simulator_id, inputStruct)


%% SIMULATION SETTINGS 
setSimSettings( simulator_id, inputStruct );


% TO-DO: Fix the bugs with the Random-spread sim run
% TO-DO: Handle the Timestamps for the Time-series sim mode
% TO-DO: Handle the Bare-soil case
%% SIMULATION INPUTS 
setSimParams( inputStruct );


%% TRANSMITTER (Tx) INPUTS
setTxParams( inputStruct );


%% RECEIVER (Rx) ANTENNA INPUTS
setRxParams( inputStruct );


% TO-DO: Diel_model added
%% GROUND INPUTS
setGndParams( simulator_id, inputStruct );


%% DYNAMIC INPUTS
setDynParams(inputStruct);


%% VEGETATION INPUTS
% It checks if ground cover is Vegetation
setVegParams( inputStruct );

end


function setDynParams( inputStruct )


%% GET GLOBAL PARAMETERS
% Simulation Settings
simulator_id = SimSettings.getInstance.simulator_id;
sim_mode_id = SimSettings.getInstance.sim_mode_id;
% Ground Parameters
num_gnd_layers = GndParams.getInstance.num_layers;


% Get the dynamic inputs full file path and name
dynInputFullFile = inputStruct.dyn_inputs_file;

% Read the dynamic inputs file
[num, ~, ~] = xlsread( dynInputFullFile, 1 );

ind = 0;

% Timestamp - Day-of-year
DOYs = [];

% If sim_mode is Time-series, then input contains timestamps
if simulator_id == Constants.id_multi_layer ...
    || sim_mode_id == Constants.id_time_series
    
    ind = ind + 1;
    DOYs = num(:, ind);
    
end


% TO-DO: This is for satGeometryManual. There should be an option for
% satGeometry
% Incidence angle
ind = ind + 1;
th0_Tx_list_deg = num(:, ind);    
th0_Tx_list_deg( isnan(th0_Tx_list_deg) ) = [];

% Azimuth angle
ind = ind + 1;
ph0_Tx_list_deg = num(:, ind);    
ph0_Tx_list_deg( isnan(ph0_Tx_list_deg) ) = [];

% Surface rms height (cm)
ind = ind + 1;
RMSH_list_cm = num(:, ind);       
RMSH_list_cm( isnan(RMSH_list_cm) ) = [];

% Volumetric soil moisture
ind = ind + 1;
VSM_list_cm3cm3 = num(:, ind : ind + (num_gnd_layers - 1) );  
VSM_list_cm3cm3( isnan(VSM_list_cm3cm3) ) = [];


% TO-DO: Think about SCoB-ML date interval input
if simulator_id == Constants.id_multi_layer   
    
    % Custom date Interval
    % d1 = '05-01-2017' ; d2 = '08-31-2017' ;
    d1 = '05-25-2017' ; d2 = '06-07-2017' ;
    % d1 = '05-25-2017' ; d2 = '05-27-2017' ;
    startDate = datenum(d1) ;
    endDate = datenum(d2) ;
    DoY1 = date2doy(startDate) ;
    DoY2 = date2doy(endDate) ;

    % DOYs and SM of interval of interest
    VSM_list_cm3cm3 = VSM_list_cm3cm3(DoY1<DOYs & DOYs<DoY2, :);
    DOYs = DOYs(DoY1<DOYs & DOYs<DoY2);
    
end


% Initialize Dynamic Parameters
DynParams.getInstance.initialize( DOYs, th0_Tx_list_deg, ph0_Tx_list_deg, VSM_list_cm3cm3, RMSH_list_cm );

end


function setGndParams( simulator_id, inputStruct ) 

% TO-DO: Handle multi-laer ground for SCoBi-Veg

% If the simulator is SCoBi-Veg (Agriculture OR Forest), set the single
% ground layer parameters (For now only! It can be multi-layer as well in
% the future)
if simulator_id == Constants.id_veg_agr ...
        || simulator_id == Constants.id_veg_for
    
    setGndSingleParams( inputStruct );
    
% Else if simulator is SCoBi-ML (Multilayer), then set the multi-layer
% ground inputs
elseif simulator_id == Constants.id_multi_layer
    
    setGndMLParams( inputStruct );
    
end

end


function setGndSingleParams( inputStruct ) 


% Get the string diel_model and convert it to an integer index
diel_model_id = findElementIdInCell( Constants.diel_models, inputStruct.diel_model );

% Sand ratio of the soil texture
sand_ratio = inputStruct.sand_ratio;

% Clay ratio of the soil texture 
clay_ratio = inputStruct.clay_ratio;  

% Soil bulk density   
rhob_gcm3 = inputStruct.rhob_gcm3;   

% Default ground is a single-layered infinite soil. No need for layer depth
layer_depth_m = [];

% Initialize Ground Parameters
GndParams.getInstance.initialize( layer_depth_m, sand_ratio, clay_ratio, rhob_gcm3, diel_model_id );

end


function setGndMLParams( inputStruct ) 


%% GENERAL GROUND PARAMETERS

% Get the string diel_model and convert it to an integer index
diel_model_id = findElementIdInCell( Constants.diel_models, inputStruct.diel_model );

% Get the dynamic inputs file
dynInputFullFile = inputStruct.dyn_inputs_file;

% Read file
[num, ~, ~] = xlsread( dynInputFullFile, 2);

% Ground layer depth (meters)
ind = 1;
layer_depth_m = num(:, ind);
layer_depth_m( isnan(layer_depth_m) ) = [];

% Sand ratio of the soil texture
ind = ind + 1;
sand_ratio = num(:, ind);
sand_ratio( isnan(sand_ratio) ) = [];

% Clay ratio of the soil texture 
ind = ind + 1;
clay_ratio = num(:, ind);
clay_ratio( isnan(clay_ratio) ) = [];

% Soil bulk density
ind = ind + 1;
rhob_gcm3 = num(:, ind);
rhob_gcm3( isnan(rhob_gcm3) ) = []; 

% Initialize Ground Parameters
GndParams.getInstance.initialize( layer_depth_m, sand_ratio, clay_ratio, rhob_gcm3, diel_model_id );

 
%% SPECIFIC MULTI-LAYERED GROUND PARAMETERS
% Layer discritization
ind = ind + 1;
delZ_m = num(1, ind);
delZ_m( isnan(delZ_m) ) = []; 

% Air layer
ind = ind + 1;
zA_m = num(1, ind);
zA_m( isnan(zA_m) ) = []; 

% Bottom-most layer
ind = ind + 1;
zB_m = num(1, ind);
zB_m( isnan(zB_m) ) = []; 

% Initialize Ground Parameters
GndMLParams.getInstance.initialize( delZ_m, zA_m, zB_m );

end


function setRxParams( inputStruct )

% Antenna Height (m)
hr_m = inputStruct.hr_m;

% Receive Antenna Gain (dB) 
G0r_dB = inputStruct.G0r_dB;

% Receiver antenna polarization
pol_Rx = inputStruct.pol_Rx;

% Get the string orientation_Rx and convert it to an integer index
orientation_Rx_id = findElementIdInCell( Constants.Rx_orientations, inputStruct.orientation_Rx );


%% DETERMINE RECEIVER ANGLES
% If receiver has a fixed orientation (observation and azimuh angles)
if orientation_Rx_id == Constants.id_Rx_fixed

    th0_Rx_deg = inputStruct.th0_Rx_deg;   % Receiver observation (theta) angle
    ph0_Rx_deg = inputStruct.ph0_Rx_deg;   % Receiver azimuth (phi) angle

% Else if receiver always faces the specular point. Then it will always 
% have equal theta and azimuth angles with the transmitter. Ignore here!
elseif orientation_Rx_id == Constants.id_Rx_specular_facing

    th0_Rx_deg = [];
    ph0_Rx_deg = [];
    
end


%% DETERMINE ANTENNA PATTERN
% Get the string ant_pat_Rx and convert it to an integer index    
ant_pat_Rx_id = findElementIdInCell( Constants.Rx_ant_pats, inputStruct.ant_pat_Rx );

% If receiver antenna pattern is Generalized-Gaussian
if ant_pat_Rx_id == Constants.id_Rx_GG
    
    ant_pat_res_deg_Rx = inputStruct.ant_pat_res_deg_Rx;   % Antenna pattern resolution in degrees
    
    % Set GG pattern parameters
    setRxGGParams( inputStruct );
    
elseif ant_pat_Rx_id == Constants.id_Rx_user_defined
    
    % TO-DO: May be the user defined pattern and resolution is calculated
    % here
    
    % Antenna pattern resolution will be calculated later by using the
    % pattern file input
    ant_pat_res_deg_Rx = [];
    
    % Set GG pattern parameters
    setRxUserDefinedParams( inputStruct );
  
% Else if Antenna pattern is Cosine to the power n
elseif ant_pat_Rx_id == Constants.id_Rx_cos_pow_n 
    
    % Should be implemented when this pattern added
end


%% INITIALIZE
% Initialize Receiver Parameters
RxParams.getInstance.initialize( hr_m, G0r_dB, pol_Rx, ant_pat_Rx_id, ant_pat_res_deg_Rx, orientation_Rx_id, th0_Rx_deg, ph0_Rx_deg );

end


function setRxGGParams( inputStruct )
    
% Beamwidth (degrees)
hpbw_deg = inputStruct.hpbw_deg;      

% Sidelobe Level (dB)
SLL_dB = inputStruct.SLL_dB;      

% X-pol level (dB)
XPL_dB = inputStruct.XPL_dB;   

% Initialize Generalized-Gaussian Receiver Parameters
RxGGParams.getInstance.initialize( hpbw_deg, SLL_dB, XPL_dB );

end


function setRxUserDefinedParams( inputStruct )

% Filename with the full path
ant_pat_Rx_file = inputStruct.ant_pat_Rx_file;    

% Initialize Generalized-Gaussian Receiver Parameters
RxUserDefinedParams.getInstance.initialize( ant_pat_Rx_file );

end


function setSimParams( inputStruct )


%% GET GLOBAL PARAMETERS
% Simulation Settings
sim_mode_id = SimSettings.getInstance.sim_mode_id;
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;


% Simulation name
sim_name = inputStruct.sim_name;

% Simulation campaign
campaign = inputStruct.campaign; 

% Simulation campaign date
campaign_date = inputStruct.campaign_date; 

% Simulation campaign plot
plot = inputStruct.plot;

% Vegetation method
veg_method_id = [];
% If gnd_cover is Vegetation, then veg_method is taken from inputs
if gnd_cover_id == Constants.id_veg_cover
    
    % If sim_mode is Snapshot, then veg_method depends on input
    if sim_mode_id == Constants.id_snapshot

        veg_method_id = findElementIdInCell( Constants.veg_methods, inputStruct.veg_method );
        
    % Else if sim_mode is Time-series, veg_method can only be Homogenous
    elseif sim_mode_id == Constants.id_time_series

        veg_method_id = Constants.id_veg_hom;

    end
    
end    

% Virtual Vegetation Orientation
veg_vir_orientation_id = [];
% If veg_method is 'Virtual' only
if veg_method_id == Constants.id_veg_vir
    veg_vir_orientation_id = findElementIdInCell( Constants.veg_vir_orientations, inputStruct.veg_vir_orientation );
end

% Vegetation plant name
vegetation_plant = [];

% If gnd_cover is Vegetation, then veg_plant is used and taken from inputs
if gnd_cover_id == Constants.id_veg_cover
    
    vegetation_plant = inputStruct.veg_plant;    
    
end

% Number of Realizations
Nr = inputStruct.Nr;

% Number of Fresnel Zones
Nfz = inputStruct.Nfz;

% Initialize Simulation Parameters
SimParams.getInstance.initialize(sim_name, campaign, campaign_date, ...
                plot, veg_method_id, veg_vir_orientation_id, vegetation_plant, Nr, Nfz );

end


function setSimSettings( simulator_id, inputStruct )

sim_mode_id = findElementIdInCell( Constants.sim_modes, inputStruct.sim_mode );

gnd_cover_id = findElementIdInCell( Constants.gnd_covers, inputStruct.gnd_cover );

if simulator_id == Constants.id_veg_agr ...
        || simulator_id == Constants.id_veg_for
    
    % Flag to write Attenuation to Excel file 
    write_attenuation = inputStruct.write_attenuation;      

    % Flag to calculate direct term  
    calc_direct_term = inputStruct.calc_direct_term;

    % Flag to calculate specular term    
    calc_specular_term = inputStruct.calc_specular_term; 

    % Flag to calculate diffuse term
    % If sim_mode is Snapshot AND gnd_cover is Vegetation, then get the flag from inputs
    if sim_mode_id == Constants.id_snapshot ...
            && gnd_cover_id == Constants.id_veg_cover

        calc_diffuse_term = inputStruct.calc_diffuse_term;

    % Else, no way for diffuse calculations
    else

        calc_diffuse_term = 0;

    end
    
    draw_live_plots = 0;

elseif simulator_id == Constants.id_multi_layer
    
    calc_direct_term = 0;
    calc_specular_term = 1;
    calc_diffuse_term = 0;
    write_attenuation = 0;
    
    draw_live_plots = 1;
    
end

% Initialize Simulation Settings
SimSettings.getInstance.initialize(simulator_id, sim_mode_id, ...
                gnd_cover_id, write_attenuation, calc_direct_term, ...
                calc_specular_term, calc_diffuse_term, draw_live_plots );

end


function setTxParams( inputStruct )

f_MHz = inputStruct.f_MHz;      % Operating frequncy (MHz)

r_Tx_m = inputStruct.r_Tx_km * Constants.km2m;    % Transmitter range from Earth"s center (km - > m) 

EIRP_dB = inputStruct.EIRP_dB;    % Equivalent Isotropic Radiated Power

pol_Tx = inputStruct.pol_Tx;  % Transmitter polarization                        

% Initialize Transmitter Parameters
TxParams.getInstance.initialize( f_MHz, r_Tx_m, EIRP_dB, pol_Tx );

end


function setVegParams( inputStruct )

%% GET GLOBAL PARAMETERS
% Simulation Parameters
veg_method_id = SimParams.getInstance.veg_method_id;
veg_vir_orientation_id = SimParams.getInstance.veg_vir_orientation_id;
% Simulation Settings
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;


% If ground cover is Vegetation, then Vegetation Parameters are set
if gnd_cover_id == Constants.id_veg_cover

    % Virtual vegetation
    if veg_method_id == Constants.id_veg_vir

        % Virtual, Row-crop vegetation (e.g. crop fields like corn, soybean etc.)
        if veg_vir_orientation_id == Constants.id_veg_vir_row_crop

            setVegVirRowParams( inputStruct );

        % Virtual, Random-spread vegetation
        elseif veg_vir_orientation_id == Constants.id_veg_vir_random_spread

            setVegVirRndParams( inputStruct );

        end

    % Homogenous vegetation
    elseif SimParams.getInstance.veg_method_id == Constants.id_veg_hom

        setVegHomParams( inputStruct );

    end

end

end


function setVegHomParams( inputStruct )


%% EXCEL FILE
vegInputFullFile = inputStruct.veg_inputs_file;

% Read the vegetation layer structure first 
[numLayers, txtLayers, rawLayers] = xlsread( vegInputFullFile, 1 );
[numKinds, txtKinds, rawKinds] = xlsread( vegInputFullFile, 2 );

[num_veg_layers, ~] = size(numLayers);

%  Layer dimensions vector (in meters)
dim_layers_m = numLayers;

% Get the total number of particles of all kinds
[~, num_particles] = size(rawKinds);


%% PARTICLES
% Initialize particles and particleIDs cells, and particle index
particlesCell = cell( 1, num_particles );
particleIDs = cell( num_particles );
particle_ind = 1;

for ii = 1 : num_particles
        
    % Get the particle id string
    particleID = char( rawKinds(1, ii) );

    % Get the flag that shows if the particle is a scaterrer
    is_scatterer = numKinds(1, ii);

    % Get the density of the particle
    density = numKinds(2, ii);
    
    % Get the average dimension 1 of the particle
    dim1_m = numKinds(3, ii);

    % Get the average dimension 2 of the particle
    dim2_m = numKinds(4, ii);

    % Get the average dimension 3 of the particle
    dim3_m = numKinds(5, ii);

    % Get the average dielectric permittivity real part of the particle
    epsr_Re = numKinds(6, ii);

    % Get the average dielectric permittivity imaginary part of the particle
    epsr_Im = numKinds(7, ii);

    % Get the average beginning angle of the particle
    prob1_deg = numKinds(8, ii);

    % Get the average ending angle of the particle
    prob2_deg = numKinds(9, ii);

    % Generate a particle struct
    particle = generateParticle( particleID, is_scatterer, density, dim1_m, dim2_m, dim3_m, epsr_Re, epsr_Im, prob1_deg, prob2_deg);

    % Add the particle to the particles cell
    particlesCell{1, particle_ind} = particle;
    
    % Add the particle's ID to the particleIDs list 
    particleIDs{1, particle_ind} = particleID;
    
    % Increment the particle index to get the next particle
    particle_ind = particle_ind + 1;
end


%% LAYERS
% Initialize layers cell
layersCell = cell(num_veg_layers, 1);

[~, max_num] = size(txtLayers);

for ii = 1 : num_veg_layers

    kk = 1;
    
    partsCell = [];

    while ( (kk) <= max_num ) && ( ~isempty( char( txtLayers(ii, kk ) ) ) )

        % Get the particle
        part = char( txtLayers(ii, kk ) );

        % Find the particle's index in the particl IDs
        part_index = strfind( particleIDs(1,:), part );
        part_index = find(not(cellfun('isempty', part_index)));

        % Update the particle IDs list regarding the layer info
        index_down = 1;
        while ~isempty( particleIDs{index_down, part_index} )
            index_down = index_down + 1;
        end
        particleIDs{index_down, part_index} = ii;

        partsCell{1,kk} = part;
        
        kk = kk + 1;

    end

    layersCell{ii, 1} = partsCell;
end

% Initialize Vegetation Parameters
VegParams.getInstance.initialize( dim_layers_m, particleIDs, particlesCell, layersCell );

end


function setVegVirRndParams( inputFile )

% TO-DO: Should be implemented when Virtual Random-spraad Vegetation added

end


function setVegVirRowParams( inputStruct )
    

%% EXCEL FILE
vegInputFullFile = inputStruct.veg_inputs_file;

% Read the vegetation layer structure first 
[num, ~, raw] = xlsread( vegInputFullFile, 1 );

% Get vegetation stage for feeding the plugin
vegetation_stage = char( raw(1,1) );

% Plugin name to be run for virtual vegetation generation    
plugin = char( raw(2,1) );

%  Distance without vegetation between two rows (m)
row_space_m = num(1,1);

% Distance without vegetation within a row (m)
col_space_m = num(2,1);

% Azimuth angle of field rows from local North (degrees)
phi_row_deg = num(3,1);

% Max scattering dist. of a plant pos within a row (m)
seed_fluctuation_m = num(4,1);

% Initialize Virtual Row-Structured Vegetation Parameters
VegVirRowParams.getInstance.initialize( row_space_m, col_space_m, ...
    phi_row_deg, seed_fluctuation_m, plugin, vegetation_stage);

end


%% Generate particle structs for types (disk or cylinder)
function particle = generateParticle( particleID, is_scatterer, dnsty, dim1_m, dim2_m, dim3_m, epsr_Re, epsr_Im, prob1_deg, prob2_deg  )

particle = struct( 'PARTICLE_ID', particleID, ...
        'IS_SCATTERER', is_scatterer, 'DENSITY', dnsty, ...
        'DIM1', dim1_m, 'DIM2', dim2_m, 'DIM3', dim3_m, ...
        'EPSILON', epsr_Re+1i*epsr_Im, ...
        'PARM1', prob1_deg, 'PARM2', prob2_deg ) ;

end




