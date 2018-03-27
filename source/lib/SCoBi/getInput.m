function getInput

input_files_xml = xmlread( strcat( Directories.getInstance.input, '\', 'inputFiles-Corn.xml' ) );
% input_files_xml = xmlread( strcat( Directories.getInstance.input, '\', 'inputFiles-Paulownia.xml' ) );

inputFile_sys = getStringFromXML(input_files_xml, ConstantNames.sys_input);
inputFile_veg = getStringFromXML(input_files_xml, ConstantNames.veg_input);

%% SIMULATION SETTINGS 
setSimSettings( inputFile_sys );


%% SIMULATION INPUTS 
setSimParams( inputFile_sys );


%% RECEIVER ANTENNA INPUTS
setRecParams( inputFile_sys );


%% TRANSMITTER SATELLITE INPUTS
setSatParams( inputFile_sys );


%% OE: GROUND INPUTS
setGndParams( inputFile_sys );


%% VEGETATION INPUTS
setVegParams( inputFile_veg );

end


function setGndParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

VSM_cm3cm3 = getDoubleArrayFromXML(xDoc, ConstantNames.gnd_VSM_cm3cm3);        % Theta probe 

sand_ratio = getDoubleFromXML(xDoc, ConstantNames.gnd_sand_ratio);  % Sand ratio of the soil texture

clay_ratio = getDoubleFromXML(xDoc, ConstantNames.gnd_clay_ratio);  % Clay ratio of the soil texture 

rhob_gcm3 = getDoubleFromXML(xDoc, ConstantNames.gnd_rhob_gcm3);  % Soil bulk density    

RMSH_cm = getDoubleArrayFromXML(xDoc, ConstantNames.gnd_RMSH_cm);  % Surface rms height (cm)

% Initialize Ground Parameters
GndParams.getInstance.initialize(VSM_cm3cm3, sand_ratio, clay_ratio, rhob_gcm3, RMSH_cm );

end


function setRecParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

hr_m = getDoubleFromXML(xDoc, ConstantNames.rec_hr_m);      % Antenna Height (m)

G0r_dB = getDoubleFromXML(xDoc, ConstantNames.rec_G0r_dB);    % Receive Antenna Gain (dB) 

hpbw_deg = getDoubleFromXML(xDoc, ConstantNames.rec_hpbw_deg);   % Beamwidth      

SLL_dB = getDoubleFromXML(xDoc, ConstantNames.rec_SLL_dB);    % Sidelobe Level

XPL_dB = getDoubleFromXML(xDoc, ConstantNames.rec_XPL_dB);    % X-pol level 

polR = getStringFromXML(xDoc, ConstantNames.rec_polR);  % Antenna polarization

% Initialize Receiver Parameters
RecParams.getInstance.initialize(hr_m, G0r_dB, hpbw_deg, SLL_dB, XPL_dB, polR);

end


function setSatParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

f_MHz = getDoubleFromXML(xDoc, ConstantNames.sat_f_MHz);      % Operating frequncy (MHz)

rsat_m = getDoubleFromXML(xDoc, ConstantNames.sat_rsat_km) * Constants.km2m;    % Satellite radius (km - > m) 

% TO-DO: This is for satGeometryManual. There should be an option for
% satGeometry
th0_deg = getDoubleArrayFromXML(xDoc, ConstantNames.sat_th0_deg);   % Incidence angle      

PH0_deg = getDoubleArrayFromXML(xDoc, ConstantNames.sat_PH0_deg);    % Azimuth angle

EIRP_dB = getDoubleFromXML(xDoc, ConstantNames.sat_EIRP_dB);    % Equivalent Isotropic Radiated Power

polT = getStringFromXML(xDoc, ConstantNames.sat_polT);  % Satellite polarization                        

% Initialize Satellite Parameters
SatParams.getInstance.initialize( f_MHz, rsat_m, th0_deg, PH0_deg, ...
                                          EIRP_dB, polT)

end


function setSimParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

sim_name = getStringFromXML(xDoc, ConstantNames.sim_simName);  % Plot of the campaign  

campaign = getStringFromXML(xDoc, ConstantNames.sim_campaign);      % The campaign of the simulation 

campaign_date = getStringFromXML(xDoc, ConstantNames.sim_campaignDate);  % Date of the campaign  

plot = getStringFromXML(xDoc, ConstantNames.sim_plot);  % Plot of the campaign  

vegetation_method = getStringFromXML(xDoc, ConstantNames.sim_vegetationMethod);  % The virtual vegetation method to be employed    

vegetation_isRow = getDoubleFromXML(xDoc, ConstantNames.sim_vegetationIsRow);  % If method is virtual, then is it row-structured or not?      

vegetation_plant = getStringFromXML(xDoc, ConstantNames.sim_vegetationPlant);  % The vegetation_plant that simulations to be run 

Nr = getDoubleFromXML(xDoc, ConstantNames.sim_Nr);       % Number of Realizations  

Nfz = getDoubleFromXML(xDoc, ConstantNames.sim_Nfz);    % Number of Fresnel Zones

% Initialize Simulation Parameters
SimParams.getInstance.initialize(sim_name, campaign, campaign_date, ...
                plot, vegetation_method, vegetation_isRow, vegetation_plant, Nr, Nfz );

end


function setSimSettings( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

sim_mode = getDoubleFromXML(xDoc, ConstantNames.set_simMode);        % 

ground_cover = getDoubleFromXML(xDoc, ConstantNames.set_groundCover);        % Vegetation mode, e.g. Bare/Veg/Both 

write_attenuation = getDoubleFromXML(xDoc, ConstantNames.set_writeAttenuation);      % Flag to write Attenuation to Excel file 

calc_direct_term = getDoubleFromXML(xDoc, ConstantNames.set_calcDirectTerm);  % Flag to calculate direct term  

calc_specular_term = getDoubleFromXML(xDoc, ConstantNames.set_calcSpecularTerm);  % Flag to calculate specular term     

calc_diffuse_term = getDoubleFromXML(xDoc, ConstantNames.set_calcDiffuseTerm);       % Flag to calculate diffuse term 

% Initialize Simulation Settings
SimSettings.getInstance.initialize(sim_mode, ground_cover, write_attenuation, ...
                calc_direct_term, calc_specular_term, ...
                calc_diffuse_term );

end


function setVegParams( inputFile )

% Virtual vegetation
if strcmp( SimParams.getInstance.vegetation_method, Constants.veg_methods.VIRTUAL )
    
    % Virtual, row-structured vegetation (e.g. crop fields like corn, soybean etc.)
    if SimParams.getInstance.vegetation_isRow
    
        setVegVirRowParams( inputFile );
    
    % Virtual, random-positioned vegetation (e.g. forests)
    else
        
        setVegVirRndParams( inputFile );
        
    end
    
% Homogenous vegetation
elseif strcmp( SimParams.getInstance.vegetation_method, Constants.veg_methods.HOMOGENOUS)
    
    setVegHomParams( inputFile );
    
end

end


function setVegHomParams( inputFile )
    
xDoc = xmlread( strcat( Directories.getInstance.input_veg_hom, '\', inputFile ) );

vegetation_stage = getStringFromXML(xDoc, ConstantNames.veg_hom_vegetationStage);

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

plant_row_spread_m = getDoubleFromXML(xDoc, ConstantNames.veg_vir_row_plantRowSpread_m);  % Max scattering dist. of a plant pos between rows (m)                       

plant_col_spread_m = getDoubleFromXML(xDoc, ConstantNames.veg_vir_row_plantColSpread_m);  % Max scattering dist. of a plant pos within a row (m)

plugin = getStringFromXML(xDoc, ConstantNames.veg_vir_row_plugin);  % Plugin name to be run for virtual vegetation generation    

% Initialize Virtual Row-Structured Vegetation Parameters
VegVirRowParams.getInstance.initialize( row_space_m, col_space_m, ...
    phi_row_deg, plant_row_spread_m, plant_col_spread_m, plugin);

end


%% Generate particle structs for types (disk or cylinder)
function particle = generateParticle( particleID, is_scatterer, dnsty, dim1_m, dim2_m, dim3_m, epsr, prob1_deg, prob2_deg  )

particle = struct( 'PARTICLE_ID', particleID, ...
        'IS_SCATTERER', is_scatterer, 'DENSITY', dnsty, ...
        'DIM1', dim1_m, 'DIM2', dim2_m, 'DIM3', dim3_m, ...
        'EPSILON', epsr, 'PARM1', prob1_deg, 'PARM2', prob2_deg ) ;

end


function output = getDoubleFromXML( xmlFile, varName )

output = str2double( xmlFile.getElementsByTagName(varName).item(0).getFirstChild.getData );

end

function output = getDoubleArrayFromXML( xmlFile, varName )

list = xmlFile.getElementsByTagName(varName).item(0);
rows = str2double( list.getElementsByTagName('rows').item(0).getFirstChild.getData );
cols = str2double( list.getElementsByTagName('cols').item(0).getFirstChild.getData );
values = list.getElementsByTagName('val');

output = zeros(rows, cols);

if values.getLength == rows * cols

    for ii = 1 : rows
        for jj = 1 : cols
            ind = (ii-1) * cols + jj - 1;
            output(ii, jj) = str2double( values.item(ind).getFirstChild.getData );
        end
    end
   
else
    % TO-DO: display error!!!
end

end

function output = getStringFromXML( xmlFile, varName )

output = char( xmlFile.getElementsByTagName(varName).item(0).getFirstChild.getData );

end