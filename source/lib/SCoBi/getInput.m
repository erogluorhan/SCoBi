function getInput

% input_files_xml = xmlread( strcat( Directories.getInstance.input, '\', 'inputFiles-Corn.xml' ) );
input_files_xml = xmlread( strcat( Directories.getInstance.input, '\', 'inputFiles-Paulownia.xml' ) );

%% SIMULATION SETTINGS 
setSimSettings( getStringFromXML(input_files_xml, ConstantNames.sys_input) );


%% SIMULATION INPUTS 
setSimParams( getStringFromXML(input_files_xml, ConstantNames.sys_input) );


%% RECEIVER ANTENNA INPUTS
setRecParams( getStringFromXML(input_files_xml, ConstantNames.sys_input) );


%% TRANSMITTER SATELLITE INPUTS
setSatParams( getStringFromXML(input_files_xml, ConstantNames.sys_input) );


%% OE: GROUND INPUTS
setGndParams( getStringFromXML(input_files_xml, ConstantNames.sys_input) );


%% VEGETATION INPUTS
setVegParams( getStringFromXML(input_files_xml, ConstantNames.veg_input) );

end


function setGndParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

VSM = getDoubleArrayFromXML(xDoc, ConstantNames.gnd_VSM);        % Theta probe 

sand = getDoubleFromXML(xDoc, ConstantNames.gnd_sand);  % Sand ratio of the soil texture

clay = getDoubleFromXML(xDoc, ConstantNames.gnd_clay);  % Clay ratio of the soil texture 

rho_b = getDoubleFromXML(xDoc, ConstantNames.gnd_rhoB);  % Soil bulk density    

RMSH = getDoubleArrayFromXML(xDoc, ConstantNames.gnd_RMSH);  % Surface rms height (cm)

GndParams.getInstance.initialize(VSM, sand, clay, rho_b, RMSH );

end


function setRecParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

hr = getDoubleFromXML(xDoc, ConstantNames.rec_hr);      % Antenna Height (m)

G0r_dB = getDoubleFromXML(xDoc, ConstantNames.rec_G0rdB);    % Receive Antenna Gain (dB) 

hpbw_deg = getDoubleFromXML(xDoc, ConstantNames.rec_hpbwDeg);   % Beamwidth      

SLL_dB = getDoubleFromXML(xDoc, ConstantNames.rec_SLLdB);    % Sidelobe Level

XPL_dB = getDoubleFromXML(xDoc, ConstantNames.rec_XPLdB);    % X-pol level 

polR = getStringFromXML(xDoc, ConstantNames.rec_polR);  % Antenna polarization

RecParams.getInstance.initialize(hr, G0r_dB, hpbw_deg, SLL_dB, XPL_dB, polR);

end


function setSatParams( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

fMHz = getDoubleFromXML(xDoc, ConstantNames.sat_fMHz);      % Operating frequncy (MHz)

rsat = getDoubleFromXML(xDoc, ConstantNames.sat_rsatKm) * Constants.km2m;    % Satellite radius (km - > m) 

% TO-DO: This is for satGeometryManual. There should be an option for
% satGeometry
th0_deg = getDoubleArrayFromXML(xDoc, ConstantNames.sat_th0Deg);   % Incidence angle      

PH0_deg = getDoubleArrayFromXML(xDoc, ConstantNames.sat_PH0Deg);    % Azimuth angle

EIRP_dB = getDoubleFromXML(xDoc, ConstantNames.sat_EIRPdB);    % Equivalent Isotropic Radiated Power

polT = getStringFromXML(xDoc, ConstantNames.sat_polT);  % Satellite polarization                        

SatParams.getInstance.initialize( fMHz, rsat, th0_deg, PH0_deg, ...
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

SimParams.getInstance.initialize(sim_name, campaign, campaign_date, ...
                plot, vegetation_method, vegetation_isRow, vegetation_plant, Nr, Nfz );

end


function setSimSettings( inputFile )

xDoc = xmlread( strcat( Directories.getInstance.input_sys, '\', inputFile ) );

sim_mode = getDoubleFromXML(xDoc, ConstantNames.set_simMode);        % 

ground_cover = getDoubleFromXML(xDoc, ConstantNames.set_groundCover);        % Vegetation mode, e.g. Bare/Veg/Both 

calc_meta_data = getDoubleFromXML(xDoc, ConstantNames.set_calcMetaData);      % Flag to calculate meta-data 

calc_direct_term = getDoubleFromXML(xDoc, ConstantNames.set_calcDirectTerm);  % Flag to calculate direct term  

calc_specular_term = getDoubleFromXML(xDoc, ConstantNames.set_calcSpecularTerm);  % Flag to calculate specular term     

calc_diffuse_term = getDoubleFromXML(xDoc, ConstantNames.set_calcDiffuseTerm);       % Flag to calculate diffuse term 

SimSettings.getInstance.initialize(sim_mode, ground_cover, calc_meta_data, ...
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

dim_layers = getDoubleArrayFromXML(xDoc, ConstantNames.veg_hom_dimLayers);    %  Layer dimensions vector

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

            dim1 = str2double( char(all_kinds(kk).item(ii-1).getElementsByTagName('dim1').item(0).getFirstChild.getData) );

            dim2 = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('dim2').item(0).getFirstChild.getData );

            dim3 = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('dim3').item(0).getFirstChild.getData );

            epsr = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('epsr').item(0).getFirstChild.getData );

            prob1 = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('prob1').item(0).getFirstChild.getData );

            prob2 = str2double( all_kinds(kk).item(ii-1).getElementsByTagName('prob2').item(0).getFirstChild.getData );

            particle = generateParticle( particleID, is_scatterer, density, dim1, dim2, dim3, epsr, prob1, prob2);

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

VegParams.getInstance.initialize( vegetation_stage, dim_layers, particleIDs, particlesCell, layersCell );

end


function setVegVirRndParams( inputFile )
    
xDoc = xmlread( strcat( Directories.getInstance.input_veg_vir_rnd, '\', inputFile ) );

end


function setVegVirRowParams( inputFile )
    
xDoc = xmlread( strcat( Directories.getInstance.input_veg_vir_row, '\', inputFile ) );

row_space = getDoubleFromXML(xDoc, ConstantNames.veg_vir_row_rowSpace);    %  Distance without vegetation between two rows (m)

col_space = getDoubleFromXML(xDoc, ConstantNames.veg_vir_row_colSpace);    %  Distance without vegetation within a row (m)

phi_row = getDoubleFromXML(xDoc, ConstantNames.veg_vir_row_phiRow);   % Azimuth angle of field rows from local North (degrees)

plant_row_spread = getDoubleFromXML(xDoc, ConstantNames.veg_vir_row_plantRowSpread);  % Max scattering dist. of a plant pos between rows (m)                       

plant_col_spread = getDoubleFromXML(xDoc, ConstantNames.veg_vir_row_plantColSpread);  % Max scattering dist. of a plant pos within a row (m)

plugin = getStringFromXML(xDoc, ConstantNames.veg_vir_row_plugin);  % Plugin name to be run for virtual vegetation generation    

VegVirRowParams.getInstance.initialize( row_space, col_space, ...
    phi_row, plant_row_spread, plant_col_spread, plugin);

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


%% Generate particle structs for types (disk or cylinder)
function particle = generateParticle( particleID, is_scatterer, dnsty, dim1, dim2, dim3, epsr, prob1, prob2  )

particle = struct( 'PARTICLE_ID', particleID, ...
        'IS_SCATTERER', is_scatterer, 'DENSITY', dnsty, ...
        'DIM1', dim1, 'DIM2', dim2, 'DIM3', dim3, ...
        'EPSILON', epsr, 'PARM1', prob1, 'PARM2', prob2 ) ;

end