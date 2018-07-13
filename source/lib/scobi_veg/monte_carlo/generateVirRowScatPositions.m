

function generateVirRowScatPositions( ind_realization, Nr_current )
% GENERATEVIRROWSCATPOSITIONS: Creates a realization of virtual,
% row-structured vegetation field.
% ind_realization: Current index of the realization that is performed
% Nr_current: Existing number of realizations in the simulation directory


    %% GET GLOBAL DIRECTORIES
    dir_config = SimulationFolders.getInstance.config;
    dir_veg = SimulationFolders.getInstance.veg;


    %% GET GLOBAL PARAMETERS
    % Simulation Parameters
    Nr = SimParams.getInstance.Nr;
    Nfz = SimParams.getInstance.Nfz;
    % Vegetation Parameters
    TYPES = VegParams.getInstance.TYPES;
    % Virtual Row-Structured Vegetation Parameters
    plugin = VegVirRowParams.getInstance.plugin;
    row_space_m = VegVirRowParams.getInstance.row_space_m;
    col_space_m = VegVirRowParams.getInstance.col_space_m;
    seed_fluctuation_m = VegVirRowParams.getInstance.seed_fluctuation_m;
    phi_row_deg = VegVirRowParams.getInstance.phi_row_deg;
    % Use the crop row azimuth angle from local East instead of local North
    phi_row_rad = degtorad( 90 - phi_row_deg );
    
    % Constant Structs
    particleDataStruct = Constants.particleDataStruct;
    
    
    %% LOAD META-DATA
    % Fresnel Zone Ellipses on the ground    
    filenamex = 'ellipses_FZ_m' ;
    ellipses_FZ_m = readVar(dir_config, filenamex) ;

    filenamex = 'ellipse_centers_FZ_m' ;
    ellipse_centers_FZ_m = readVar(dir_config, filenamex) ;
    
    
    %% TO-DO:Change the names of the below variables that are equal to global params
    %% DEFINE PERSISTENT VARIABLES
    persistent highestPlantHeight 
    persistent dim_layers_m scat_cal_veg LTK TYPKND ...
        dsty dim1_m dim2_m dim3_m epsr parm1_deg ...
        parm2_deg
    persistent num_layers num_types num_kinds
  

    %% Slant Range

    majorAxis = 2 * ellipses_FZ_m(:, 1);
    minorAxis = 2 * ellipses_FZ_m(:, 2);

    if majorAxis > minorAxis
        fieldSide = majorAxis;
    else
        fieldSide = minorAxis;
    end

    numRows = ceil(fieldSide / row_space_m) + 1;
    numCols = ceil(fieldSide / col_space_m) + 1;

    % Rotation in Azimuth due to the angle between vegetation field rows and
    % the local east
    fieldRotZ = [cos(phi_row_rad) -sin(phi_row_rad) 0;
                     sin(phi_row_rad) cos(phi_row_rad)  0;
                     0 0 1] ;

    particleData = [];   % [Layer Type Kind posX posY posZ downAngle azimuthAngle ...
                    %  dim1 dim2 dim3 epsrRe epsrIm fresnelZoneIndex]
                    % dim1: bottom node radius for stem, length/2 for leaf
                    % dim2: top node radius for stem, width/2 for leaf
                    % dim3: length for stem, thickness for leaf
                    % epsrRe: dielectric constant real part for any type/kind
                    % epsrIm: dielectric constant imaginary part for any type/kind
    
    % Generate Vegetation Field in the selected Fresnel Zone
    
    % TO-DO: Below position actually is a position relative to the
        % center of Fresnel zone, instead of pure specular or ground frame. i.e.
        % there is a small difference between the specular point and FZ
        % centers
    startPos_fz = [-fieldSide(Nfz)/2;      % The bottom-left corner of the field. 
                -fieldSide(Nfz)/2;      % pSc2_m is assumed to be shifted to (0,0), it will be added later
                0];                                     % i.e. It will be shifted back to its original position.

    % Generate plants in the field if they are in any of the fresnel zones
    for ii = 1 : numRows(end)

        for jj = 1 : numCols(end)

            currentPos_fz = startPos_fz + [(jj-1) * col_space_m; (ii-1) * row_space_m; 0];
            
            % Fluctuate plant position
            fluctuationAngle = deg2rad(360 * rand);
            fluctuationRadius = seed_fluctuation_m * rand;
            currentPos_fz(1) = currentPos_fz(1) + fluctuationRadius * cos(fluctuationAngle);
            currentPos_fz(2) = currentPos_fz(2) + fluctuationRadius * sin(fluctuationAngle);

            currentPosRotated_fz = fieldRotZ * (currentPos_fz);
            currentPosFinal = currentPosRotated_fz + ellipse_centers_FZ_m(:,end); % ground frame

            % Determine which fresnel zone the position is in, if any
            ff = Nfz;
            while ( ff > 0 )

                fzAngle = atan( ellipse_centers_FZ_m(2,ff) / ellipse_centers_FZ_m(1,ff) );

                if isPointInEllipse(currentPosFinal, majorAxis(ff)/2, minorAxis(ff)/2, ellipse_centers_FZ_m(:,ff), fzAngle)
%                 if isPointInFZRect( currentPosRotated_fz, majorAxis(ff)/2, minorAxis(ff)/2, fzAngle )
                    ff = ff - 1;
                else
                    break;
                end

            end

            if ff < Nfz
                
                ff = ff + 1;

                % Call the virtual vegetation plugin's generatePlant function
                [plantHeight, plantParticles] = plugin.generatePlant( TYPES, currentPosFinal );
                
                plantParticles(:,particleDataStruct.fzIndex) = ff;
                
                particleData = [particleData; plantParticles];
                
                if plantHeight > highestPlantHeight
                    highestPlantHeight = plantHeight; 
                end

            end

        end % numCols

    end % numRows
    
    num_layers = 1;
    num_types = max( particleData(:, particleDataStruct.Type) );
    num_kinds = max( particleData(:, particleDataStruct.Kind) );
        
    if ( ind_realization == Nr_current + 1 )
                    
        TYPES_names = fieldnames( TYPES );
        num_types_all = length(TYPES_names);
    
        highestPlantHeight = 0;

        LTK = cell(num_kinds, num_types_all, num_layers);
        
        if Nr_current == 0
            % Initialize new statistics
            dim_layers_m = 0;

            scat_cal_veg = zeros(num_kinds, num_types_all, num_layers);

            TYPKND = zeros(num_layers, num_types_all );

            dsty = zeros(num_kinds, num_types_all, num_layers);

            dim1_m = zeros(num_kinds, num_types_all, num_layers);

            dim2_m = zeros(num_kinds, num_types_all, num_layers);

            dim3_m = zeros(num_kinds, num_types_all, num_layers);

            epsr = zeros(num_kinds, num_types_all, num_layers);

            parm1_deg = zeros(num_kinds, num_types_all, num_layers);

            parm2_deg = zeros(num_kinds, num_types_all, num_layers);  
        
        else
            % Use statistics from the past 
            
            dim_layers_m = VegParams.getInstance.dim_layers_m;

            scat_cal_veg = VegParams.getInstance.scat_cal_veg;

            TYPKND = VegParams.getInstance.TYPKND;

            
            
            VOL = zeros(num_layers,1);
            fresnelZoneAreas = pi * majorAxis/2 .* minorAxis/2;

            for ii = 1 : num_layers  % Layers
                % Calculate layer volumes to use in calculating avg parameters
                VOL(ii) = dim_layers_m(ii) * fresnelZoneAreas(end);
                for jj = 1 : num_types  % Layers
                    for kk = 1 : num_kinds  % Layers
                        dsty(kk, jj, ii) = VegParams.getInstance.dsty(kk, jj, ii) * Nr_current * VOL(ii);
                    end
                end
            end

            dim1_m = VegParams.getInstance.dim1_m * Nr_current;

            dim2_m = VegParams.getInstance.dim2_m * Nr_current;

            dim3_m = VegParams.getInstance.dim3_m * Nr_current;

            epsr = VegParams.getInstance.epsr * Nr_current;

            parm1_deg = VegParams.getInstance.parm1_deg * Nr_current;

            parm2_deg = VegParams.getInstance.parm2_deg * Nr_current;  
        end

    end

    % Adjust the vegetation top layer height
    if highestPlantHeight < sum( dim_layers_m )
        highestPlantHeight = sum( dim_layers_m ); 
    end
    
                
    % Calculate Fresnel Zone Areas and corrsponding volumes 
    fresnelZoneAreas = pi * majorAxis/2 .* minorAxis/2;

    for ii = 1 : num_layers  % Layers

        %Find indices to elements in the Layer, Type and Kind columns that satisfy the equality
        %Then use the logical indices to return required sub-matrices
%         indLayer = particleData( :, particleDataStruct.Layer ) == ii;
%         particleDataPerLayer = particleData(indLayer,:);
        particleDataPerLayer = particleData; % Because there is only one layer in virtual row vegetation

        for jj = 1 : num_types      % Types

            indType = particleDataPerLayer( :, particleDataStruct.Type ) == jj;
            particleDataPerType = particleDataPerLayer(indType,:);

            for kk = 1 : num_kinds  % Kinds

                indKind = particleDataPerType( :, particleDataStruct.Kind ) == kk;
                particleDataPerKind = particleDataPerType(indKind,:);  

                if ~isempty(particleDataPerKind)

                    % If the particle Type is not a Leaf
                    if jj ~= TYPES.L

                        % Sort the scatterers w.r.t. fresnel zones
                        particleDataPerKind = sortrows(particleDataPerKind,particleDataStruct.fzIndex);
                        freq = unique(particleDataPerKind(:,particleDataStruct.fzIndex));
                        numParticlesPerFZ = [freq,histc(particleDataPerKind(:,particleDataStruct.fzIndex),freq)];   
                        numParticlesPerFZ = cumsum(numParticlesPerFZ(:,2));

                        filename = strcat('R', num2str(ind_realization), '_L', num2str(ii), '_T', num2str(jj), '_K', num2str(kk)) ;
                        
                        writeVar( SimulationFolders.getInstance.position, filename, particleDataPerKind(:, particleDataStruct.posX : particleDataStruct.epsrIm) ) ;

                        writeVar(SimulationFolders.getInstance.fzones, filename, numParticlesPerFZ) ;

                        scat_cal_veg(kk, jj, ii) = 1; % Update the VegParams.scat_cal_veg (Main scatterers' flags)

                        
                        %% Incident, scattering angles
                        calcScatAngles(filename, particleDataPerKind(:, particleDataStruct.posX : particleDataStruct.posZ)) ;
%                         calcScatAngles( ind_realization, filename, particleDataPerKind );
                    end

                    % Determine TYPKND and LTK for this Kind
                    if ind_realization == 1
                        % Update VegParams.LTK and VegParams.TYPKND
                        LTK{kk, jj, ii} = strcat( TYPES_names{jj,1}, num2str(kk) ) ;
                        TYPKND(ii, jj) = TYPKND(ii, jj) + 1;

                    end

                    dsty(kk, jj, ii) = dsty(kk, jj, ii) + length( particleDataPerKind );
                    dim1_m(kk, jj, ii) = dim1_m(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.dim1) );
                    dim2_m(kk, jj, ii) = dim2_m(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.dim2) );
                    dim3_m(kk, jj, ii) = dim3_m(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.dim3) );
                    epsr(kk, jj, ii) = epsr(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.epsrRe) ) + 1i * mean( particleDataPerKind(:, particleDataStruct.epsrIm) );
                    parm1_deg(kk, jj, ii) = parm1_deg(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.downAngle) );
                    parm2_deg(kk, jj, ii) = parm2_deg(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.downAngle) );
                    
                end
            end
        end
    end

    if( ind_realization == Nr )
        
        VOL = zeros(num_layers,1);
        
        dim_layers_m = highestPlantHeight;
        
        for ii = 1 : num_layers  % Layers
            
            % Calculate layer volumes to use in calculating avg parameters
            VOL(ii) = dim_layers_m(ii) * fresnelZoneAreas(end);

            for jj = 1 : num_types  % Layers
                for kk = 1 : num_kinds  % Layers
                    
                    % TO-DO: Density is used per homogenous attenuation
                    % calculations for now; thus, stem density can be calculated per volume. 
                    % If another use is needed in the future, it might be calculated per area
                    dsty(kk, jj, ii) = dsty(kk, jj, ii) / Nr / VOL(ii);
                    dim1_m(kk, jj, ii) = dim1_m(kk, jj, ii) / Nr;
                    dim2_m(kk, jj, ii) = dim2_m(kk, jj, ii) / Nr;
                    dim3_m(kk, jj, ii) = dim3_m(kk, jj, ii) / Nr;
                    epsr(kk, jj, ii) = epsr(kk, jj, ii) / Nr;
                    parm1_deg(kk, jj, ii) = parm1_deg(kk, jj, ii) / Nr;
                    parm2_deg(kk, jj, ii) = parm2_deg(kk, jj, ii) / Nr;
                end
            end
        end
        
        VegParams.getInstance.initialize3( dim_layers_m, scat_cal_veg,...
            TYPKND, LTK, dsty, dim1_m, dim2_m, dim3_m, epsr, parm1_deg, parm2_deg );

        %% SAVE ALL
        writeVar( dir_veg, ConstantNames.veg_hom_scatCalVeg, scat_cal_veg) ;

        writeVar( dir_veg, ConstantNames.veg_hom_dimLayers_m, dim_layers_m ) ;

        writeVar( dir_veg, ConstantNames.veg_hom_TYPKND, TYPKND) ;
                        
        currentDir = pwd;
        cd( dir_veg )
        save LTK LTK
        cd( currentDir )

        writeVar( dir_veg, ConstantNames.veg_hom_dsty, dsty) ;

        writeVar( dir_veg, ConstantNames.veg_hom_dim1_m, dim1_m) ;

        writeVar( dir_veg, ConstantNames.veg_hom_dim2_m, dim2_m) ;

        writeVar( dir_veg, ConstantNames.veg_hom_dim3_m, dim3_m) ;

        writeComplexVar( dir_veg, ConstantNames.veg_hom_epsr, epsr) ;

        writeVar( dir_veg, ConstantNames.veg_hom_parm1_deg, parm1_deg) ;
        
        writeVar( dir_veg, ConstantNames.veg_hom_parm2_deg, parm2_deg) ;
        
    end

end


function result = isPointInEllipse( position, semiMajorAxis, semiMinorAxis, ellipseCenter, ellipseAngle) 

    x = position(1);
    y = position(2);
    h = ellipseCenter(1);
    k = ellipseCenter(2);

    result = ( cos(ellipseAngle)*(x - h) + sin(ellipseAngle)*(y - k) )^2 / semiMajorAxis^2 + ( sin(ellipseAngle)*(x - h) - cos(ellipseAngle)*(y - k) )^2 / semiMinorAxis^2 <= 1.0;

end


function result = isPointInFZRect( position_sf, fzSemiMajorAxis, fzSemiMinorAxis, fzAngle_rad) 

    fzRotZ = [cos(-fzAngle_rad) -sin(-fzAngle_rad) 0;
                 sin(-fzAngle_rad) cos(-fzAngle_rad)  0;
                 0 0 1] ;
    
    position_sf_fzRot = fzRotZ * position_sf;    

    x = position_sf_fzRot(1);
    y = position_sf_fzRot(2);
    
    rectSemiMajor = sqrt(pi) / 2 * fzSemiMajorAxis ;  % semi-major side
    rectSemiMinor = sqrt(pi) / 2 * fzSemiMinorAxis ;  % semi-minor side

    result = ( x >= -rectSemiMajor && x <= rectSemiMajor && y >= -rectSemiMinor && y <= rectSemiMinor );

end