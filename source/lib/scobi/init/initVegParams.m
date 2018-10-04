


function vegInputFullFile = initVegParams( inputStruct )

%% GET GLOBAL PARAMETERS
% Simulation Settings
gnd_cover_id = SimSettings.getInstance.gnd_cover_id;

vegInputFullFile = [];

% If ground cover is Vegetation, then Vegetation Parameters are set
if gnd_cover_id == Constants.id_veg_cover
    
    % Initialize workspace for vegetation
    initWSVegetation();

    vegInputFullFile = initVegHomParams( inputStruct );

end

end


function vegInputFullFile = initVegHomParams( inputStruct )


%% EXCEL FILE
vegInputFullFile = inputStruct.veg_inputs_file;

% Read the vegetation layer structure first 
[dim_layers_m, txtLayers, ~] = xlsread( vegInputFullFile, 1 );
[numKinds, ~, rawKinds] = xlsread( vegInputFullFile, 2 );

% Extract the required information from the file
[num_veg_layers, ~] = size(dim_layers_m);
kindsData = rawKinds( :, 2 : end );
layersData = txtLayers( 2 : end, 3 : end );

% Get the total number of particles of all kinds
[~, num_particles] = size(kindsData);


%% PARTICLES
% Initialize particles and particleIDs cells, and particle index
particlesCell = cell( 1, num_particles );
particleIDs = cell( num_particles );
particle_ind = 1;

for ii = 1 : num_particles
        
    % Get the particle id string
    particleID = char( kindsData(1, ii) );

    % Get the density of the particle
    density = numKinds(1, ii);
    
    % Get the average dimension 1 of the particle
    dim1_m = numKinds(2, ii);

    % Get the average dimension 2 of the particle
    dim2_m = numKinds(3, ii);

    % Get the average dimension 3 of the particle
    dim3_m = numKinds(4, ii);

    % Get the average dielectric permittivity real part of the particle
    epsr_Re = numKinds(5, ii);

    % Get the average dielectric permittivity imaginary part of the particle
    epsr_Im = numKinds(6, ii);

    % Get the average beginning angle of the particle
    prob1_deg = numKinds(7, ii);

    % Get the average ending angle of the particle
    prob2_deg = numKinds(8, ii);

    % Generate a particle struct
    particle = generateParticle( particleID, density, dim1_m, dim2_m, dim3_m, epsr_Re, epsr_Im, prob1_deg, prob2_deg);

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

[~, max_num] = size(layersData);

for ii = 1 : num_veg_layers

    kk = 1;
    
    partsCell = [];

    while ( (kk) <= max_num ) && ( ~isempty( char( layersData(ii, kk ) ) ) )

        % Get the particle
        part = char( layersData(ii, kk ) );

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
VegParams.getInstance.setup( dim_layers_m, particleIDs, particlesCell, layersCell );

end


%% Generate particle structs for types (disk or cylinder)
function particle = generateParticle( particleID, dnsty, dim1_m, dim2_m, dim3_m, epsr_Re, epsr_Im, prob1_deg, prob2_deg  )

particle = struct( 'PARTICLE_ID', particleID, ...
        'DENSITY', dnsty, ...
        'DIM1', dim1_m, 'DIM2', dim2_m, 'DIM3', dim3_m, ...
        'EPSILON', epsr_Re+1i*epsr_Im, ...
        'PARM1', prob1_deg, 'PARM2', prob2_deg ) ;

end