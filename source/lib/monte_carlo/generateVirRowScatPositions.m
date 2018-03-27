

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
    % Vegetation Parameters
    TYPES = VegParams.getInstance.TYPES;
    % Virtual Row-Structured Vegetation Parameters
    plugin = VegVirRowParams.getInstance.plugin;
    row_space_m = VegVirRowParams.getInstance.row_space_m;
    col_space_m = VegVirRowParams.getInstance.col_space_m;
    plant_row_spread_m = VegVirRowParams.getInstance.plant_row_spread_m;
    plant_col_spread_m = VegVirRowParams.getInstance.plant_col_spread_m;
    % Use the row azimuth angle from local East (x axis) instead of local North
    phi_row_rad = degtorad( 90 - VegVirRowParams.getInstance.phi_row_deg );
    
    % Constant Structs
    particleDataStruct = Constants.particleDataStruct;
    
    
    %% LOAD META-DATA
    filenamex = 'AllPoints_m' ;
    AllPoints_m = readVar(dir_config, filenamex) ;

    % Fresnel Zone Ellipses on the ground    
    filenamex = 'ellipse_s_m' ;
    ellipse_s_m = readVar(dir_config, filenamex) ;

    filenamex = 'ellipse_s_centers_m' ;
    ellipse_s_centers_m = readVar(dir_config, filenamex) ;

    % Transmitter Satellite Rotation along Z-Axis
    filenamex = 'AntRotZt' ;
    AntRotZt = readVar(dir_config, filenamex) ;

    % Ground to Specular Frame Transformation
    filenamex = 'Tgs' ;
    Tgs = readVar(dir_config, filenamex) ;
    
    
    %% DEFINE PERSISTENT VARIABLES
    persistent highestPlantHeight 
    persistent dim_layers_m scat_cal_veg LTK TYPKND ...
        dsty dim1_m dim2_m dim3_m epsr parm1_deg ...
        parm2_deg
    persistent num_layers num_types num_kinds
  

    %% Slant Range
    
    pSc2_m = AllPoints_m(:, 9);    % Center of first Fresnel Zone
    
    pT_m = AllPoints_m(:, 1) ;          % Transmitter

    ht = pT_m(3) ;

    majorAxis = 2 * ellipse_s_m(:, 1);
    minorAxis = 2 * ellipse_s_m(:, 2);

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
    [numFresnelZones, ~] = size(ellipse_s_m);
    
    % TO-DO: Below position actually is a position relative to the
        % center of Fresnel zone, instead of pure specular or ground frame. i.e.
        % there is a small difference between the specular point and FZ
        % centers
    startPos_fz = [-fieldSide(numFresnelZones)/2;      % The bottom-left corner of the field. 
                -fieldSide(numFresnelZones)/2;      % pSc2_m is assumed to be shifted to (0,0), it will be added later
                0];                                     % i.e. It will be shifted back to its original position.

    % Generate plants in the field if they are in any of the fresnel zones
    for ii = 1 : numRows(end)

        for jj = 1 : numCols(end)

            currentPos_fz = startPos_fz + [(jj-1) * col_space_m; (ii-1) * row_space_m; 0];
            currentPos_fz(1) = currentPos_fz(1) + ( -plant_row_spread_m ) + 2 * plant_row_spread_m * rand;
            currentPos_fz(2) = currentPos_fz(2) + ( -plant_col_spread_m) + 2 * plant_col_spread_m * rand;

            currentPosRotated_fz = fieldRotZ * (currentPos_fz);
            currentPosRotated = currentPosRotated_fz + ellipse_s_centers_m(:,end); % ground frame

            % Determine which fresnel zone the position is in, if any
            kk = numFresnelZones;
            while ( kk > 0 )

                fzAngle = atan( ellipse_s_centers_m(2,kk) / ellipse_s_centers_m(1,kk) );

                if isPointInEllipse(currentPosRotated, majorAxis(kk)/2, minorAxis(kk)/2, ellipse_s_centers_m(:,kk), fzAngle)
                    kk = kk - 1;
                else
                    break;
                end

            end

            if kk < numFresnelZones

                % Call the virtual vegetation plugin's generatePlant function
                [plantHeight, plantParticles] = plugin.generatePlant( TYPES, currentPosRotated );
                
                plantParticles(:,particleDataStruct.fzIndex) = kk+1;
                
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
        
        particleData_realizations = cell(Nr,1);
        
        particleDataPerKind_realizations = cell(Nr,1);
        
        %numParticlesPerFZ_realizations = cell(Nr,1);
        numParticlesPerFZ_realizations = 0;

    end

    % Adjust the vegetation top layer height
    if highestPlantHeight < sum( dim_layers_m )
        highestPlantHeight = sum( dim_layers_m ); 
    end

%     % TO-DO: Should be deprecated since virtual vegetation stuies are planned to have only one layer !
%     % Adjust the corresponding layers to the types because it was not finalized
%     % in GenerateCorn() method
%     adjustLayersToTypes(ind_realization);

%     % TO-DO: Move from here to analysis part
%     % %% OE: Figures for testing
%     figure
%     axis equal
%     hold on
%     
%     % Draw scatterer positions
%     indexParticleDataExceptLeaf = particleData(:, particleDataStruct.Type) ~= TYPES.L; % Choose types other than leaf 
%     particleDataExceptLeaf = particleData(indexParticleDataExceptLeaf,:);
%     temp = [];
%     for ii = 1 : numFresnelZones
%         indFresnelZone = particleDataExceptLeaf(:,particleDataStruct.fzIndex) == ii;
%         temp = [temp; particleDataExceptLeaf(indFresnelZone,:)];
%         numParticlesPerFZExceptLeaf(ii,1) = length(temp);
%     end
%     particleDataExceptLeaf = temp;
%     
%     %particleDataExceptLeaf_sf = Tgs * particleDataExceptLeaf(:, particleDataStruct.posX:particleDataStruct.posZ)';
%     
%     N2 = 1 ;
%     for fz = 1 : numFresnelZones
%         
%         color = 'r.';
%         if mod(fz,2) == 0
%             color = 'b.';
%         end
%         N1 = N2 ;
%         N2 = numParticlesPerFZExceptLeaf(fz) ;
%         if N2 > N1
%             scatter(particleDataExceptLeaf( N1:N2, particleDataStruct.posX ), particleDataExceptLeaf( N1:N2, particleDataStruct.posY ), color);
%             %scatter(particleDataExceptLeaf_sf(1, N1:N2), particleDataExceptLeaf_sf(2, N1:N2), color);
%         elseif N2 == 0
%             N2 = 1;
%         end
%     
%     end
%     
%     
%     % Draw fresnel zones
%     fz_thetas = -pi : 0.01 : pi;
%     
%     lenR = length(fz_thetas);
%     tempR = zeros(3,lenR);
%     
%     for ii = 1 : numFresnelZones
%         
%         tempR(1,:) = majorAxis(ii) / 2 * cos(fz_thetas);
%         tempR(2,:) = minorAxis(ii) / 2 * sin(fz_thetas);
%         tempR = AntRotZt * tempR;
%         fz_x = ellipse_s_centers_m(1,ii) + tempR(1,:);
%         fz_y = ellipse_s_centers_m(2,ii) + tempR(2,:);
%         
%         fz_sf = Tgs * [fz_x; fz_y; zeros(1, length(fz_x))];
%         
%         plot(fz_x, fz_y, 'r--')
%         %plot(fz_sf(1,:), fz_sf(2,:), 'r--')
%     end
%     xlabel('East ->')
%     ylabel('North ->')
%     % 
%     % %% OE: End of figures for testing

   
    particleData_realizations{ind_realization} = particleData;
                
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