


function initVegParams( inputStruct )

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

            initVegVirRowParams( inputStruct );

        % Virtual, Random-spread vegetation
        elseif veg_vir_orientation_id == Constants.id_veg_vir_random_spread

            initVegVirRndParams( inputStruct );

        end

    % Homogenous vegetation
    elseif SimParams.getInstance.veg_method_id == Constants.id_veg_hom

        initVegHomParams( inputStruct );

    end

end

end


function initVegHomParams( inputStruct )


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


function initVegVirRndParams( inputStruct )

% TO-DO: Should be implemented when Virtual Random-spraad Vegetation added

end


function initVegVirRowParams( inputStruct )
    

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