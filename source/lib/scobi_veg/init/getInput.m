function getInput(inputFile_sys, inputFile_veg, inputFile_dyn)

%% SIMULATION SETTINGS 
setSimSettings( inputFile_sys );


% TO-DO: Fix the bugs with the Random-spread sim run
% TO-DO: Handle the Timestamps for the Time-series sim mode
% TO-DO: Handle the Bare-soil case
%% SIMULATION INPUTS 
setSimParams( inputFile_sys );


%% TO-DO: Dynamic inputs added
setDynParams(inputFile_dyn);


%% TRANSMITTER (Tx) INPUTS
setTxParams( inputFile_sys );


%% RECEIVER (Rx) ANTENNA INPUTS
setRxParams( inputFile_sys );


% TO-DO: Diel_model added
%% OE: GROUND INPUTS
setGndParams( inputFile_sys );


% TO-DO: Moved to Excel file
%% VEGETATION INPUTS
setVegParams( inputFile_veg );

end


function setDynParams( inputFile )


%% GET GLOBAL DIRECTORIES
dir_input_dyn = Directories.getInstance.input_dyn;


%% GET GLOBAL PARAMETERS
sim_mode_id = SimSettings.getInstance.sim_mode_id;


[num, txt, raw] = xlsread( strcat( dir_input_dyn, '\', inputFile ),1);

ind = 1;

% TO-DO: This is for satGeometryManual. There should be an option for
% satGeometry

th0_Tx_list_deg = num(:, ind);    % Incidence angle
th0_Tx_list_deg( isnan(th0_Tx_list_deg) ) = [];
ind = ind + 1;

ph0_Tx_list_deg = num(:, ind);    % Azimuth angle
ph0_Tx_list_deg( isnan(ph0_Tx_list_deg) ) = [];
ind = ind + 1;

% TO-DO: VSM and RMSH moved to Excel file
VSM_list_cm3cm3 = num(:, ind);    % Theta probe 
VSM_list_cm3cm3( isnan(VSM_list_cm3cm3) ) = [];
ind = ind + 1;

RMSH_list_cm = num(:, ind);       % Surface rms height (cm)
RMSH_list_cm( isnan(RMSH_list_cm) ) = [];
ind = ind + 1;

if sim_mode_id == Constants.id_time_series
   
    time_stamps = num(:, ind);
    ind = ind + 1;
    
end

% Initialize Dynamic Parameters
DynParams.getInstance.initialize( th0_Tx_list_deg, ph0_Tx_list_deg, VSM_list_cm3cm3, RMSH_list_cm );

end


function setGndParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

sand_ratio = getDoubleFromXML(xDoc, ConstantNames.gnd_sand_ratio);  % Sand ratio of the soil texture

clay_ratio = getDoubleFromXML(xDoc, ConstantNames.gnd_clay_ratio);  % Clay ratio of the soil texture 

rhob_gcm3 = getDoubleFromXML(xDoc, ConstantNames.gnd_rhob_gcm3);  % Soil bulk density    

% Get the string diel_model and convert it to an integer index
diel_model = getStringFromXML(xDoc, ConstantNames.gnd_diel_model);              % The dielectric permittivity model
is_diel_model = cellfun(@(x)isequal(x, diel_model), Constants.diel_models );
[~, diel_model_id] = find(is_diel_model);  

% Initialize Ground Parameters
GndParams.getInstance.initialize(sand_ratio, clay_ratio, rhob_gcm3, diel_model_id );

end


function setRxParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

hr_m = getDoubleFromXML(xDoc, ConstantNames.Rx_hr_m);      % Antenna Height (m)

G0r_dB = getDoubleFromXML(xDoc, ConstantNames.Rx_G0r_dB);    % Receive Antenna Gain (dB) 

pol_Rx = getStringFromXML(xDoc, ConstantNames.Rx_pol_Rx);  % Receiver antenna polarization

% Get the string ant_pat_Rx and convert it to an integer index
ant_pat_Rx = getStringFromXML(xDoc, ConstantNames.Rx_ant_pat_Rx);        % The virtual vegetation method to be employed
is_ant_pat_Rx = cellfun(@(x)isequal(x, ant_pat_Rx), Constants.Rx_ant_pats );
[~, ant_pat_Rx_id] = find(is_ant_pat_Rx);     

ant_pat_res_deg = getDoubleFromXML(xDoc, ConstantNames.Rx_ant_pat_res_deg);   % Antenna pattern resolution in degrees

% Get the string orientation_Rx and convert it to an integer index
orientation_Rx = getStringFromXML(xDoc, ConstantNames.Rx_orientation_Rx);        % The virtual vegetation method to be employed
is_orientation_Rx = cellfun(@(x)isequal(x, orientation_Rx), Constants.Rx_orientations );
[~, orientation_Rx_id] = find(is_orientation_Rx);


%% DETERMINE RECEIVER ANGLES
% If receiver has a fixed orientation (observation and azimuh angles)
if orientation_Rx_id == Constants.id_Rx_fixed

    th0_Rx_deg = getDoubleFromXML(xDoc, ConstantNames.Rx_th0_Rx_deg);   % Receiver observation (theta) angle
    ph0_Rx_deg = getDoubleFromXML(xDoc, ConstantNames.Rx_ph0_Rx_deg);   % Receiver azimuth (phi) angle

% Else if receiver always faces the specular point. Then it will always 
% have equal theta and azimuth angles with the transmitter. Ignore here!
elseif orientation_Rx_id == Constants.id_Rx_specular_facing

    th0_Rx_deg = [];
    ph0_Rx_deg = [];
    
end

% Initialize Receiver Parameters
RxParams.getInstance.initialize( hr_m, G0r_dB, pol_Rx, ant_pat_Rx_id, ant_pat_res_deg, orientation_Rx_id, th0_Rx_deg, ph0_Rx_deg );


%% SPECIFIC RECEIVER ANTENNA PATTERN PARAMETERS
% Antenna pattern: Generalized-Gaussian
if ant_pat_Rx_id == Constants.id_Rx_GG 

    setRxGGParams(inputFile);

% Antenna pattern: Cosine to the power n
elseif ant_pat_Rx_id == Constants.id_Rx_cos_pow_n 
    
    % TO-DO: Should be implemented for this pattern
    
elseif ant_pat_Rx_id == Constants.id_Rx_user_defined
    
    setRxUserDefinedParams( inputFile )
    
end

end


function setRxGGParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

hpbw_deg = getDoubleFromXML(xDoc, ConstantNames.Rx_GG_hpbw_deg);   % Beamwidth      

SLL_dB = getDoubleFromXML(xDoc, ConstantNames.Rx_GG_SLL_dB);    % Sidelobe Level

XPL_dB = getDoubleFromXML(xDoc, ConstantNames.Rx_GG_XPL_dB);    % X-pol level 

% Initialize Generalized-Gaussian Receiver Parameters
RxGGParams.getInstance.initialize( hpbw_deg, SLL_dB, XPL_dB );

end


function setRxUserDefinedParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

user_def_ant_pat_file = getStringFromXML(xDoc, ConstantNames.Rx_user_def_ant_pat_file);    % Filename with the full path

% Initialize Generalized-Gaussian Receiver Parameters
RxUserDefinedParams.getInstance.initialize( user_def_ant_pat_file );

end


function setSimParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

sim_name = getStringFromXML(xDoc, ConstantNames.sim_simName);  % Plot of the campaign  

campaign = getStringFromXML(xDoc, ConstantNames.sim_campaign);      % The campaign of the simulation 

campaign_date = getStringFromXML(xDoc, ConstantNames.sim_campaignDate);  % Date of the campaign  

plot = getStringFromXML(xDoc, ConstantNames.sim_plot);  % Plot of the campaign  

% Get the string veg_method and convert it to an integer index
veg_method = getStringFromXML(xDoc, ConstantNames.sim_vegMethod);        % The virtual vegetation method to be employed
is_veg_method = cellfun(@(x)isequal(x, veg_method), Constants.veg_methods );
[~, veg_method_id] = find(is_veg_method);      

% Get the string veg_vir_orientation and convert it to an integer index
veg_vir_orientation = getStringFromXML(xDoc, ConstantNames.sim_vegVirOrientation);        % If veg_method is virtual, then is it Row-crop or Random-spread?
is_veg_vir_orientation = cellfun(@(x)isequal(x, veg_vir_orientation), Constants.veg_vir_orientations );
[~, veg_vir_orientation_id] = find(is_veg_vir_orientation);

vegetation_plant = getStringFromXML(xDoc, ConstantNames.sim_vegetationPlant);  % The vegetation_plant that simulations to be run 

Nr = getDoubleFromXML(xDoc, ConstantNames.sim_Nr);       % Number of Realizations  

Nfz = getDoubleFromXML(xDoc, ConstantNames.sim_Nfz);    % Number of Fresnel Zones

% Initialize Simulation Parameters
SimParams.getInstance.initialize(sim_name, campaign, campaign_date, ...
                plot, veg_method_id, veg_vir_orientation_id, vegetation_plant, Nr, Nfz );

end


function setSimSettings( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

simulator = getStringFromXML(xDoc, ConstantNames.set_simMode);        % 
is_simulator = cellfun(@(x)isequal(x, simulator), Constants.simulators );
[~, simulator_id] = find(is_simulator);

sim_mode = getStringFromXML(xDoc, ConstantNames.set_simMode);        % 
is_sim_mode = cellfun(@(x)isequal(x, sim_mode), Constants.sim_modes );
[~, sim_mode_id] = find(is_sim_mode);

gnd_cover = getStringFromXML(xDoc, ConstantNames.set_gndCover);             % Bare-soil OR Vegetation 
is_gnd_cover = cellfun(@(x)isequal(x, gnd_cover), Constants.gnd_covers );
[~, gnd_cover_id] = find(is_gnd_cover);

write_attenuation = getDoubleFromXML(xDoc, ConstantNames.set_writeAttenuation);      % Flag to write Attenuation to Excel file 

calc_direct_term = getDoubleFromXML(xDoc, ConstantNames.set_calcDirectTerm);  % Flag to calculate direct term  

calc_specular_term = getDoubleFromXML(xDoc, ConstantNames.set_calcSpecularTerm);  % Flag to calculate specular term     

calc_diffuse_term = getDoubleFromXML(xDoc, ConstantNames.set_calcDiffuseTerm);       % Flag to calculate diffuse term 

% Initialize Simulation Settings
SimSettings.getInstance.initialize(simulator_id, sim_mode_id, ...
                gnd_cover_id, write_attenuation, calc_direct_term, ...
                calc_specular_term, calc_diffuse_term );

end


function setTxParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

f_MHz = getDoubleFromXML(xDoc, ConstantNames.Tx_f_MHz);      % Operating frequncy (MHz)

r_Tx_m = getDoubleFromXML(xDoc, ConstantNames.Tx_rsat_km) * Constants.km2m;    % Transmitter range from Earth"s center (km - > m) 

EIRP_dB = getDoubleFromXML(xDoc, ConstantNames.Tx_EIRP_dB);    % Equivalent Isotropic Radiated Power

pol_Tx = getStringFromXML(xDoc, ConstantNames.Tx_pol_Tx);  % Transmitter polarization                        

% Initialize Transmitter Parameters
TxParams.getInstance.initialize( f_MHz, r_Tx_m, EIRP_dB, pol_Tx );

end


function setVegParams( inputFile )

%% GET GLOBAL PARAMETERS
% Simulation Parameters
veg_method_id = SimParams.getInstance.veg_method_id;
veg_vir_orientation_id = SimParams.getInstance.veg_vir_orientation_id;

% Virtual vegetation
if veg_method_id == Constants.id_veg_vir
    
    % Virtual, Row-crop vegetation (e.g. crop fields like corn, soybean etc.)
    if veg_vir_orientation_id == Constants.id_veg_vir_row_crop
    
        setVegVirRowParams( inputFile );
    
    % Virtual, Random-spread vegetation
    elseif veg_vir_orientation_id == Constants.id_veg_vir_random_spread
        
        setVegVirRndParams( inputFile );
        
    end
    
% Homogenous vegetation
elseif SimParams.getInstance.veg_method_id == Constants.id_veg_hom
    
    setVegHomParams( inputFile );
    
end

end


function setVegHomParams( inputFile )
    
xDoc = xmlread( strcat( Directories.getInstance.input_veg_hom, '\', inputFile ) );

vegetation_stage = getStringFromXML(xDoc, ConstantNames.veg_vegetationStage);

dim_layers_m = getDoubleArrayFromXML(xDoc, ConstantNames.veg_hom_dimLayers_m);    %  Layer dimensions vector

types = xDoc.getElementsByTagName('types').item(0);

num_particles = 0;

leaf = types.getElementsByTagName('leaf').item(0);
leaf_kinds = leaf.getElementsByTagName('kinds').item(0);
if~isempty( leaf_kinds )
    leaf_kinds = leaf_kinds.getElementsByTagName('kind');
    num_particles = num_particles + leaf_kinds.getLength;
end

branch = types.getElementsByTagName('branch').item(0);
branch_kinds = branch.getElementsByTagName('kinds').item(0);
if~isempty( branch_kinds )
    branch_kinds = branch_kinds.getElementsByTagName('kind');
    num_particles = num_particles + branch_kinds.getLength;
end

trunk = types.getElementsByTagName('trunk').item(0);
trunk_kinds = trunk.getElementsByTagName('kinds').item(0);
if~isempty( trunk_kinds )
    trunk_kinds = trunk_kinds.getElementsByTagName('kind');
    num_particles = num_particles + trunk_kinds.getLength;
end

needle = types.getElementsByTagName('needle').item(0);
needle_kinds = needle.getElementsByTagName('kinds').item(0);
if~isempty( needle_kinds )
    needle_kinds = needle_kinds.getElementsByTagName('kind');
    num_particles = num_particles + needle_kinds.getLength;
end

all_kinds = [leaf_kinds, branch_kinds, trunk_kinds, needle_kinds];

particlesCell = cell( 1, num_particles );
particleIDs = cell( num_particles );
particle_ind = 1;

for kk = 1 : length( all_kinds )
    
    if all_kinds(kk).getLength >= 0

        for ii = 1 : all_kinds(kk).getLength
            
            particleID = char( all_kinds(kk).item(ii-1).getElementsByTagName('id').item(0).getFirstChild.getData );

            is_scatterer = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('is_scatterer').item(0).getFirstChild.getData );

            density = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('density').item(0).getFirstChild.getData );

            dim1_m = str2double( char(all_kinds(kk).item(ii-1).getElementsByTagName('dim1_m').item(0).getFirstChild.getData) );

            dim2_m = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('dim2_m').item(0).getFirstChild.getData );

            dim3_m = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('dim3_m').item(0).getFirstChild.getData );

            epsr = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('epsr').item(0).getFirstChild.getData );

            prob1_deg = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('prob1_deg').item(0).getFirstChild.getData );

            prob2_deg = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('prob2_deg').item(0).getFirstChild.getData );

            particle = generateParticle( particleID, is_scatterer, density, dim1_m, dim2_m, dim3_m, epsr, prob1_deg, prob2_deg);

            particlesCell{1, particle_ind} = particle;
            particleIDs{1, particle_ind} = particleID;
            particle_ind = particle_ind + 1;

        end

    else
        % TO-DO: display error!!!
    end
end

layers = xDoc.getElementsByTagName('layers').item(0);
layerList = layers.getElementsByTagName('layer');
if layerList.getLength >= 0
    
    layersCell = cell(layerList.getLength, 1);
    
    for jj = 1 : layerList.getLength
        layerInd = str2double( layerList.item(jj-1).getElementsByTagName('id').item(0).getFirstChild.getData );
        particleList = layerList.item(jj-1).getElementsByTagName('particle');
        
        partsCell = cell(1,particleList.getLength);
        
        for kk = 1 : particleList.getLength
            
            part = char( particleList.item(kk-1).getFirstChild.getData );
                   
            part_index = strfind( particleIDs(1,:), part );
            part_index = find(not(cellfun('isempty', part_index)));
            
            index_down = 1;
            while ~isempty( particleIDs{index_down, part_index} )
                index_down = index_down + 1;
            end
            particleIDs{index_down, part_index} = layerInd;
        
            partsCell{1,kk} = part;
                    
        end
        
        layersCell{jj, 1} = partsCell;
    end
end

% Initialize Vegetation Parameters
VegParams.getInstance.initialize( vegetation_stage, dim_layers_m, particleIDs, particlesCell, layersCell );

end


function setVegVirRndParams( inputFile )
    
xDoc = xmlread( strcat( Directories.getInstance.input_veg_vir_rnd, '\', inputFile ) );

% TO-DO: Should be implemented for Virtual Random Vegetation

end


function setVegVirRowParams( inputFile )
    
xDoc = xmlread( strcat( Directories.getInstance.input_veg_vir_row, '\', inputFile ) );

row_space_m = getDoubleFromXML(xDoc, ConstantNames.veg_vir_row_rowSpace_m);    %  Distance without vegetation between two rows (m)

col_space_m = getDoubleFromXML(xDoc, ConstantNames.veg_vir_row_colSpace_m);    %  Distance without vegetation within a row (m)

phi_row_deg = getDoubleFromXML(xDoc, ConstantNames.veg_vir_row_phiRow_deg);   % Azimuth angle of field rows from local North (degrees)

seed_fluctuation_m = getDoubleFromXML(xDoc, ConstantNames.veg_vir_row_seedFluctuation_m);  % Max scattering dist. of a plant pos within a row (m)

plugin = getStringFromXML(xDoc, ConstantNames.veg_vir_row_plugin);  % Plugin name to be run for virtual vegetation generation    

vegetation_stage = getStringFromXML( xDoc, ConstantNames.veg_vegetationStage );

% Initialize Virtual Row-Structured Vegetation Parameters
VegVirRowParams.getInstance.initialize( row_space_m, col_space_m, ...
    phi_row_deg, seed_fluctuation_m, plugin, vegetation_stage);

end


%% Generate particle structs for types (disk or cylinder)
function particle = generateParticle( particleID, is_scatterer, dnsty, dim1_m, dim2_m, dim3_m, epsr, prob1_deg, prob2_deg  )

particle = struct( 'PARTICLE_ID', particleID, ...
        'IS_SCATTERER', is_scatterer, 'DENSITY', dnsty, ...
        'DIM1', dim1_m, 'DIM2', dim2_m, 'DIM3', dim3_m, ...
        'EPSILON', epsr, 'PARM1', prob1_deg, 'PARM2', prob2_deg ) ;

end




