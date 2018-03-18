%% Mehmet Kurum
% 03/16/2017

function generateScatPos(indRealization, Nr_current)
% Defined over Fresnel Zones

if strcmp( SimParams.getInstance.vegetation_method, Constants.veg_methods.HOMOGENOUS )
    
    generateHomScatPositions(indRealization);

elseif strcmp( SimParams.getInstance.vegetation_method, Constants.veg_methods.VIRTUAL )
    
    if SimParams.getInstance.vegetation_isRow

        generateVirRowScatPositions( indRealization, Nr_current );
    
    else
        
        generateVirRndScatPositions();
    
    end
    
end

end

function generateHomScatPositions( indRealization )

 
AllPoints = readVar(SimulationFolders.getInstance.config, 'AllPoints') ;

pT = AllPoints(:, 1) ;          % Transmitter position

ht = pT(3) ;                    % Transmitter height

%% Reading vegetation parameters...
dim_layers = VegParams.getInstance.dim_layers;

scat_cal_veg = VegParams.getInstance.scat_cal_veg;

TYPKND = VegParams.getInstance.TYPKND;

dsty = VegParams.getInstance.dsty;

dim3 = VegParams.getInstance.dim3;

%% Layer parameters
[Nlayer, Ntype] = size(TYPKND) ;

%% Calculations...
D2 = 0 ;

for ii = 1 : Nlayer
    
    D1 = D2 ;               % Layer Top
    D2 = D1 + dim_layers(ii, 1) ;    % Layer Bottom
    
    for jj = 1 : Ntype
        
        Nkind = TYPKND(ii, jj) ;
        
        for kk = 1 : Nkind
            
            filename = strcat('R', num2str(indRealization), '_L', num2str(ii), '_T',...
                num2str(jj), '_K', num2str(kk)) ;
            disp(filename)
                      
            % If the particle is a scatterer
            if scat_cal_veg(kk, jj, ii) == 1 
                
                tic ;
            
                RHO = dsty(kk, jj, ii) ;
                LEN = dim3(kk, jj, ii) ;
                
                disp('Calculating scatterer position...')
                
                scatPosition(ht, D1, D2, RHO, LEN, filename, jj)
                
                disp('done...')
                
                toc
                
            else
                
                disp('Skipped calculating scatterer position for the particle...')
                
            end                      
            
        end % Nkind
        
    end % Ntype
    
end % Nlayer
    
end


function generateVirRowScatPositions( indRealization, Nr_current )

    particleDataStruct = Constants.particleDataStruct;
    
    persistent highestPlantHeight 
    persistent dim_layers scat_cal_veg LTK TYPKND ...
        dsty dim1 dim2 dim3 epsr parm1 ...
        parm2
    
    persistent num_layers num_types num_kinds

    % Use the row azimuth angle from local East (x axis) instead of local North
    phi_row = degtorad( 90 - VegVirRowParams.getInstance.phi_row );  

    %% Slant Range
    pathname = SimulationFolders.getInstance.config;
    
    filenamex = 'AllPoints' ;
    AllPoints = readVar(pathname, filenamex) ;

    % Fresnel Zone Ellipses on the ground    
    filenamex = 'ellipse_s' ;
    ellipse_s = readVar(pathname, filenamex) ;

    filenamex = 'ellipse_s_centers' ;
    ellipse_s_centers = readVar(pathname, filenamex) ;

    % Transmitter Satellite Rotation along Z-Axis
    filenamex = 'AntRotZt' ;
    AntRotZt = readVar(pathname, filenamex) ;

    % Ground to Specular Frame Transformation
    filenamex = 'Tgs' ;
    Tgs = readVar(pathname, filenamex) ;

    pSc2 = AllPoints(:, 9);    % Center of first Fresnel Zone
    
    pT = AllPoints(:, 1) ;          % Transmitter

    ht = pT(3) ;

    majorAxis = 2 * ellipse_s(:, 1);
    minorAxis = 2 * ellipse_s(:, 2);

    if majorAxis > minorAxis
        fieldSide = majorAxis;
    else
        fieldSide = minorAxis;
    end

    numRows = ceil(fieldSide / VegVirRowParams.getInstance.row_space) + 1;
    numCols = ceil(fieldSide / VegVirRowParams.getInstance.col_space) + 1;

    % Rotation in Azimuth due to the angle between vegetation field rows and
    % the local east
    fieldRotZ = [cos(phi_row) -sin(phi_row) 0;
                     sin(phi_row) cos(phi_row)  0;
                     0 0 1] ;

    particleData = [];   % [Layer Type Kind posX posY posZ downAngle azimuthAngle ...
                    %  dim1 dim2 dim3 epsrRe epsrIm fresnelZoneIndex]
                    % dim1: bottom node radius for stem, length/2 for leaf
                    % dim2: top node radius for stem, width/2 for leaf
                    % dim3: length for stem, thickness for leaf
                    % epsrRe: dielectric constant real part for any type/kind
                    % epsrIm: dielectric constant imaginary part for any type/kind
    
    % Generate Vegetation Field in the selected Fresnel Zone
    [numFresnelZones, ~] = size(ellipse_s);
    
    % TO-DO: Below position actually is a position relative to the
        % center of Fresnel zone, instead of pure specular or ground frame. i.e.
        % there is a small difference between the specular point and FZ
        % centers
    startPos_fz = [-fieldSide(numFresnelZones)/2;      % The bottom-left corner of the field. 
                -fieldSide(numFresnelZones)/2;      % pSc2 is assumed to be shifted to (0,0), it will be added later
                0];                                     % i.e. It will be shifted back to its original position.

    % Generate plants in the field if they are in any of the fresnel zones
    for ii = 1 : numRows(end)

        for jj = 1 : numCols(end)

            currentPos_fz = startPos_fz + [(jj-1) * VegVirRowParams.getInstance.col_space; (ii-1) * VegVirRowParams.getInstance.row_space; 0];
            currentPos_fz(1) = currentPos_fz(1) + ( -VegVirRowParams.getInstance.plant_row_spread ) + 2 * VegVirRowParams.getInstance.plant_row_spread * rand;
            currentPos_fz(2) = currentPos_fz(2) + ( -VegVirRowParams.getInstance.plant_col_spread) + 2 * VegVirRowParams.getInstance.plant_col_spread * rand;

            currentPosRotated_fz = fieldRotZ * (currentPos_fz);
            currentPosRotated = currentPosRotated_fz + ellipse_s_centers(:,end); % ground frame

            % Determine which fresnel zone the position is in, if any
            kk = numFresnelZones;
            while ( kk > 0 )

                fzAngle = atan( ellipse_s_centers(2,kk) / ellipse_s_centers(1,kk) );

                if isPointInEllipse(currentPosRotated, majorAxis(kk)/2, minorAxis(kk)/2, ellipse_s_centers(:,kk), fzAngle)
                    kk = kk - 1;
                else
                    break;
                end

            end

            if kk < numFresnelZones

                % Call the virtual vegetation plugin's generatePlant function
                [plantHeight, plantParticles] = VegVirRowParams.getInstance.plugin.generatePlant( VegParams.getInstance.TYPES, currentPosRotated );
                
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
        
    if ( indRealization == Nr_current + 1 )
                    
        TYPES_struct = VegParams.getInstance.TYPES;
        TYPES_names = fieldnames( TYPES_struct );
        num_types_all = length(TYPES_names);
    
        highestPlantHeight = 0;

        LTK = cell(num_kinds, num_types_all, num_layers);
        
        if Nr_current == 0
            % Initialize new statistics
            dim_layers = 0;

            scat_cal_veg = zeros(num_kinds, num_types_all, num_layers);

            TYPKND = zeros(num_layers, num_types_all );

            dsty = zeros(num_kinds, num_types_all, num_layers);

            dim1 = zeros(num_kinds, num_types_all, num_layers);

            dim2 = zeros(num_kinds, num_types_all, num_layers);

            dim3 = zeros(num_kinds, num_types_all, num_layers);

            epsr = zeros(num_kinds, num_types_all, num_layers);

            parm1 = zeros(num_kinds, num_types_all, num_layers);

            parm2 = zeros(num_kinds, num_types_all, num_layers);  
        
        else
            % Use statistics from the past 
            
            dim_layers = VegParams.getInstance.dim_layers;

            scat_cal_veg = VegParams.getInstance.scat_cal_veg;

            TYPKND = VegParams.getInstance.TYPKND;

            
            
            VOL = zeros(num_layers,1);
            fresnelZoneAreas = pi * majorAxis/2 .* minorAxis/2;

            for ii = 1 : num_layers  % Layers
                % Calculate layer volumes to use in calculating avg parameters
                VOL(ii) = dim_layers(ii) * fresnelZoneAreas(end);
                for jj = 1 : num_types  % Layers
                    for kk = 1 : num_kinds  % Layers
                        dsty(kk, jj, ii) = VegParams.getInstance.dsty(kk, jj, ii) * Nr_current * VOL(ii);
                    end
                end
            end

            dim1 = VegParams.getInstance.dim1 * Nr_current;

            dim2 = VegParams.getInstance.dim2 * Nr_current;

            dim3 = VegParams.getInstance.dim3 * Nr_current;

            epsr = VegParams.getInstance.epsr * Nr_current;

            parm1 = VegParams.getInstance.parm1 * Nr_current;

            parm2 = VegParams.getInstance.parm2 * Nr_current;  
        end
        
        particleData_realizations = cell(SimParams.getInstance.Nr,1);
        
        particleDataPerKind_realizations = cell(SimParams.getInstance.Nr,1);
        
        %numParticlesPerFZ_realizations = cell(SimParams.getInstance.Nr,1);
        numParticlesPerFZ_realizations = 0;

    end

    % Adjust the vegetation top layer height
    if highestPlantHeight < sum( VegParams.getInstance.dim_layers )
        highestPlantHeight = sum( VegParams.getInstance.dim_layers ); 
    end

%     % TO-DO: Should be deprecated since virtual vegetation stuies are planned to have only one layer !
%     % Adjust the corresponding layers to the types because it was not finalized
%     % in GenerateCorn() method
%     adjustLayersToTypes(indRealization);

%     % TO-DO: Move from here to analysis part
%     % %% OE: Figures for testing
%     figure
%     axis equal
%     hold on
%     
%     % Draw scatterer positions
%     indexParticleDataExceptLeaf = particleData(:, particleDataStruct.Type) ~= VegParams.getInstance.TYPES.L; % Choose types other than leaf 
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
%         fz_x = ellipse_s_centers(1,ii) + tempR(1,:);
%         fz_y = ellipse_s_centers(2,ii) + tempR(2,:);
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

   
    particleData_realizations{indRealization} = particleData;
                
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
                    if jj ~= VegParams.getInstance.TYPES.L

                        % Sort the scatterers w.r.t. fresnel zones
                        particleDataPerKind = sortrows(particleDataPerKind,particleDataStruct.fzIndex);
                        freq = unique(particleDataPerKind(:,particleDataStruct.fzIndex));
                        numParticlesPerFZ = [freq,histc(particleDataPerKind(:,particleDataStruct.fzIndex),freq)];   
                        numParticlesPerFZ = cumsum(numParticlesPerFZ(:,2));

                        filename = strcat('R', num2str(indRealization), '_L', num2str(ii), '_T', num2str(jj), '_K', num2str(kk)) ;
                        
                        writeVar( SimulationFolders.getInstance.position, filename, particleDataPerKind(:, particleDataStruct.posX : particleDataStruct.epsrIm) ) ;

                        writeVar(SimulationFolders.getInstance.fzones, filename, numParticlesPerFZ) ;

                        scat_cal_veg(kk, jj, ii) = 1; % Update the VegParams.scat_cal_veg (Main scatterers' flags)

                        
                        %% Incident, scattering angles
                        calcScatAngles(filename, particleDataPerKind(:, particleDataStruct.posX : particleDataStruct.posZ)) ;
%                         calcScatAngles( indRealization, filename, particleDataPerKind );
                    end

                    % Determine TYPKND and LTK for this Kind
                    if indRealization == 1
                        % Update VegParams.LTK and VegParams.TYPKND
                        LTK{kk, jj, ii} = strcat( TYPES_names{jj,1}, num2str(kk) ) ;
                        TYPKND(ii, jj) = TYPKND(ii, jj) + 1;

                    end

                    dsty(kk, jj, ii) = dsty(kk, jj, ii) + length( particleDataPerKind );
                    dim1(kk, jj, ii) = dim1(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.dim1) );
                    dim2(kk, jj, ii) = dim2(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.dim2) );
                    dim3(kk, jj, ii) = dim3(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.dim3) );
                    epsr(kk, jj, ii) = epsr(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.epsrRe) ) + 1i * mean( particleDataPerKind(:, particleDataStruct.epsrIm) );
                    parm1(kk, jj, ii) = parm1(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.downAngle) );
                    parm2(kk, jj, ii) = parm2(kk, jj, ii) + mean( particleDataPerKind(:, particleDataStruct.downAngle) );
                    
                end
            end
        end
    end

    if( indRealization == SimParams.getInstance.Nr )
        
        VOL = zeros(num_layers,1);
        
        dim_layers = highestPlantHeight;
        
        for ii = 1 : num_layers  % Layers
            
            % Calculate layer volumes to use in calculating avg parameters
            VOL(ii) = dim_layers(ii) * fresnelZoneAreas(end);

            for jj = 1 : num_types  % Layers
                for kk = 1 : num_kinds  % Layers
                    
                    % TO-DO: Density is used per homogenous attenuation
                    % calculations for now; thus, stem density can be calculated per volume. 
                    % If another use is needed in the future, it might be calculated per area
                    dsty(kk, jj, ii) = dsty(kk, jj, ii) / SimParams.getInstance.Nr / VOL(ii);
                    dim1(kk, jj, ii) = dim1(kk, jj, ii) / SimParams.getInstance.Nr;
                    dim2(kk, jj, ii) = dim2(kk, jj, ii) / SimParams.getInstance.Nr;
                    dim3(kk, jj, ii) = dim3(kk, jj, ii) / SimParams.getInstance.Nr;
                    epsr(kk, jj, ii) = epsr(kk, jj, ii) / SimParams.getInstance.Nr;
                    parm1(kk, jj, ii) = parm1(kk, jj, ii) / SimParams.getInstance.Nr;
                    parm2(kk, jj, ii) = parm2(kk, jj, ii) / SimParams.getInstance.Nr;
                end
            end
        end
        
        VegParams.getInstance.initialize3( dim_layers, scat_cal_veg,...
            TYPKND, LTK, dsty, dim1, dim2, dim3, epsr, parm1, parm2 );

        %% Saving...
        filename = 'scat_cal_veg' ;
        writeVar(SimulationFolders.getInstance.veg, filename, VegParams.getInstance.scat_cal_veg) ;

        filename = strcat( 'dim_layers' );
        writeVar(SimulationFolders.getInstance.veg, filename, VegParams.getInstance.dim_layers ) ;

        filename = 'TYPKND' ;
        writeVar(SimulationFolders.getInstance.veg, filename, VegParams.getInstance.TYPKND) ;
                        
        currentDir = pwd;
        cd( SimulationFolders.getInstance.veg )
        save LTK LTK
        cd(currentDir)

        filename = 'dsty' ;
        writeVar(SimulationFolders.getInstance.veg, filename, VegParams.getInstance.dsty) ;

        filename = 'dim1' ;
        writeVar(SimulationFolders.getInstance.veg, filename, VegParams.getInstance.dim1) ;

        filename = 'dim2' ;
        writeVar(SimulationFolders.getInstance.veg, filename, VegParams.getInstance.dim2) ;

        filename = 'dim3' ;
        writeVar(SimulationFolders.getInstance.veg, filename, VegParams.getInstance.dim3) ;

        filename = 'epsr' ;
        writeComplexVar(SimulationFolders.getInstance.veg, filename, VegParams.getInstance.epsr) ;

        filename = 'parm1' ;
        writeVar(SimulationFolders.getInstance.veg, filename, VegParams.getInstance.parm1) ;
        
        filename = 'parm2' ;
        writeVar(SimulationFolders.getInstance.veg, filename, VegParams.getInstance.parm2) ;
        
    end

end



%% Generate Scatterer positions
function scatPosition(ht, ds, de, rho, LEN, filename, ifT)

vegetation_plant = SimParams.getInstance.vegetation_plant;

%% Reading
dim_layers = VegParams.getInstance.dim_layers;

%% Layer Depth and Height
d = sum(dim_layers) ;        % thickness of layer in meter
hr = RecParams.getInstance.hr ;            % Receiver height
hr = hr - d ;         % heigth of antenna from tree tops in meter
ht = ht - d ;         % heigth of transmitter from tree tops in meter

%% Ellipse on the layer where z = -ds (layer top)
[~, ~, ax1, by1] = calcFresnelZones(ht + ds, hr + ds) ;
aas = ax1 ; % major axis
bbs = by1 ; % minor axis

%% Ellipse on the layer where z = -de (layer bottom)
[~, ~, ax1, by1] = calcFresnelZones(ht + de, hr + de) ;
aae = ax1 ; % major axis
bbe = by1 ; % minor axis

%% Illuminated Volume
As = pi * aas .* bbs ; % area on top layer
Ae = pi * aae .* bbe ; % area on bottom layer
VOL = (de - ds) * (As + Ae) / 2 ;

% Above VOL produces the same volume results with the following 
% aas_r = sqrt(pi) / 2 * aas ; % major side
% bbs_r = sqrt(pi) / 2 * bbs ; % minor side
% aae_r = sqrt(pi) / 2 * aae ; % major side
% bbe_r = sqrt(pi) / 2 * bbe ; % minor side
% As_r = 4 * aas_r .* bbs_r ; % area on top layer
% Ae_r = 4 * aae_r .* bbe_r ; % area on bottom layer
% VOL_r = (de - ds) * (As_r + Ae_r) / 2 ;


%% Generate random positions

% Number of particles
Npart = ceil(VOL * rho) ;
disp(strcat('# of kind in Fresnel Zones :', num2str(Npart)))

num_fresnel_zones = SimParams.getInstance.Nfz

% TO-DO: Decide what to do with this
if strcmp(vegetation_plant, 'OneTrunkCenter')
    Npart(:) = 1 ;
    SimParams.getInstance.Nfz = 1 ;
end

% positions of particles (xp_sf, yp_sf, zp)
xp_sf = zeros(Npart(end), 1) ;
yp_sf = xp_sf ;
zp = xp_sf ;

N2 = 1 ;
for fz = 1 : num_fresnel_zones
    
    N1 = N2 ;
    N2 = Npart(fz) ;
    
    for ii = N1 : N2
        
        if ifT == VegParams.getInstance.TYPES.T
            % The trunk is attached to the ground
            % Calculate the position for center of gravity
            zp(ii, 1) = -de + LEN / 2 ;            
        else            
            % Generating zp between ds (layer top) and de (layer bottom)
            % Calculate the position for center of gravity
            s = -ds ; e = -de ;
            zp(ii, 1) = s + (e - s) * rand(1) ;            
        end
        
        [S0x, x1, ax1, by1] = calcFresnelZones(ht + abs(zp(ii, 1)), hr + abs(zp(ii, 1))) ;        

        % outer ellipse
        a2 = ax1(fz) ;  % major axis
        b2 = by1(fz) ;  % minor axis
        
        % inner ellipse
        if fz == 1   % first Fresnel zone         
            a1 = 0 ;  
            b1 = 0 ; 
        else
            a1 = ax1(fz - 1) ;
            b1 = by1(fz - 1) ;            
        end 
        
        a1_r = sqrt(pi) / 2 * a1 ;
        b1_r = sqrt(pi) / 2 * b1 ;
        a2_r = sqrt(pi) / 2 * a2 ;  % major side
        b2_r = sqrt(pi) / 2 * b2 ;  % minor side
  
        outRect = 0 ;    % outside inner rectangle
        while ~outRect
            s = -a2_r ;  e = a2_r ;
            xx_r = s + (e - s) * rand(1) ;
            s = -b2_r ; e = b2_r ;
            yy_r = s + (e - s) * rand(1) ;
            outRect = (abs(xx_r)>= a1_r) || (abs(yy_r)>= b1_r) ;
        end
        
        % TO-DO: Below summation actually gives a position relative to the
        % center of Fresnel zone, instead of pure specular frame. i.e.
        % there is a small difference between the specular point and FZ
        % centers
        xp_sf(ii, 1) = x1(fz) + xx_r ;
        yp_sf(ii, 1) = yy_r ;
        
        % TO-DO: Decide what to do with this
        if strcmp(vegetation_plant, 'OneTrunkCenter')
            xp_sf(ii, 1) = S0x ;
            yp_sf(ii, 1) = 2 * b2_r ;
        end
        
    end
    
end

% Convert zp to absolute z-values (that increases upwards starting from the zero-ground) from relative values within the
% vegetation layer (that regards the vegetation top as zero, and goes negative downwards)
zp = zp + d;

% The position of particle center of gravity in specular (local) frame
pP_sf = [xp_sf, yp_sf, zp] ;

% Convert specular frame positions to ground frame
Tgs = readVar(SimulationFolders.getInstance.config, 'Tgs') ;   % G -> S
pP = (Tgs' * pP_sf')';


%% Saving . . .

disp('saving...')

% Save particle center  of gravity positions in the ground frame with 
% absolute z-values starting from zero-ground
writeVar(SimulationFolders.getInstance.position, filename, pP) ;

% Number of particles per each Fresnel zone with a cumulative array
writeVar(SimulationFolders.getInstance.fzones, filename, Npart) ;

%% Incident, scattering angles
calcScatAngles(filename, pP) ;


end

%% calculates  incident and scattering angles

function calcScatAngles(filename, pP)
% pP: Particle center of gravity positions in ground frame with absolute z-values

%% Recall transformations and configuration
% incident unit vectors (only in one direction - transmitter is far away)
pathname = SimulationFolders.getInstance.config;

filenamex = 'isn' ;
isn = readVar(pathname, filenamex) ;   % propagation vector (i_s^-)
filenamex = 'osp' ;
isp = readVar(pathname, filenamex) ;   % propagation vector (i_s^+=o_s^+)

% transformations
filenamex = 'Tgs' ;
Tgs = readVar(pathname, filenamex) ;   % G -> S
filenamex = 'Tgr' ;
Tgr = readVar(pathname, filenamex) ;   % G -> R
filenamex = 'TgrI' ;
TgrI = readVar(pathname, filenamex) ;  % G -> RI

filenamex = 'AllPoints' ;
AllPoints = readVar(pathname, filenamex) ;

pT = AllPoints(:, 1);        % Transmitter
pTI = AllPoints(:, 2);        % Transmitter Image

pR = AllPoints(:, 4);        % Receiver
pRI = AllPoints(:, 5);       % Receiver Image

%% Number of Particles
[Npart, ~] = size(pP) ; 

%% repmat
pR = repmat(pR', Npart, 1) ;         % Receiver
pRI = repmat(pRI', Npart, 1) ;       % Image Receiver

pT = repmat(pT', Npart, 1) ;         % Transmitter
pTI = repmat(pTI', Npart, 1) ;       % Image Transmitter


%% Scatteting directions and distances to Transmitter / Image Transmitter
ki = pT - pP ;                  % vector from particle to antenna
ri = sqrt(sum(ki' .^ 2))' ;     % distance from particle to antenna

kiI = pTI - pP ;                    % vector from particle to image antenna
riI = sqrt(sum(kiI' .^ 2))' ;       % distance from particle to image antenna

%% Scatteting directions and distances to receiver/image reciever
ko = pR - pP ;                  % vector from particle to antenna
ro = sqrt(sum(ko' .^ 2))' ;     % distance from particle to antenna
ko = ko ./ repmat(ro, 1, 3) ;   % unit vector

koI = pRI - pP ;                    % vector from particle to image antenna
roI = sqrt(sum(koI' .^ 2))' ;       % distance from particle to image antenna
koI = koI ./ repmat(roI, 1, 3) ;    % unit vector

% TO-DO: Decide what to do with this
% figure
% riroI = ri + roI ;
%  riIro = riI + ro ;
% plot(riIro, 'or')
% hold
% plot(riroI, 'ob')

% figure
% subplot(2,2,1)
% plot(ro, ':o')
% title('r_o')
% subplot(2,2,2)
% plot(roI, ':o')
% title('r_o_I')
% subplot(2,2,3)
% plot(ri, ':o')
% title('r_i')
% subplot(2,2,4)
% plot(riI, ':o')
% title('r_i_I')
% 
% figure
% subplot(2,2,1)
% plot(ro+ri, ':o')
% title('dd : r_o + r_i')
% subplot(2,2,2)
% plot(roI+ri, ':o')
% title('rd : r_o_I + r_i')
% subplot(2,2,3)
% plot(ro+riI, ':o')
% title('dr : r_o + r_i_I')
% subplot(2,2,4)
% plot(riI+roI, ':o')
% title('rr : r_o_I + r_i_I')

%% Scattering angles from particles 

% propagation vector to receiver in local (specular) frame
ko_sf = (Tgs * ko')' ;

% FROM PARTICLE TO REC
% incidence angle
thsd = acos(ko_sf(:, 3)) * 180 / pi ;
% azimuth angle
phsd = atan2(ko_sf(:, 2), ko_sf(:, 1)) * 180 / pi ;

% propagation vector to image of the receiver in local (specular) frame
koI_sf = (Tgs * koI')' ;

% FROM PARTICLE TO IM of REC
% incidence angle
thsdI = acos(koI_sf(:, 3)) * 180 / pi ;
% azimuth angle
phsdI = atan2(koI_sf(:, 2), koI_sf(:, 1)) * 180 / pi ;


%% Incident angles on particles

% propagation vector TRANS TO PARTICLE in local (specular) frame
ki_sf = (Tgs * isn)' ;

% incidence angle
thid = acos(-ki_sf(:, 3)) * 180 / pi ;
% azimuth angle
phid = atan2(-ki_sf(:, 2), -ki_sf(:, 1)) * 180 / pi ;

% propagation vector IM OF TRANS TO PARTICLE in local (specular) frame
kiI_sf = (Tgs * isp)' ;

% incidence angle
thidI = acos(-kiI_sf(:, 3)) * 180 / pi ;
% azimuth angle
phidI = atan2(-kiI_sf(:, 2), -kiI_sf(:, 1)) * 180 / pi ;

%% incident angles to receiver and image of receiver
% propagation vector in receiver antenna system
ko_rf = (Tgr * ko')' ;

% FROM PARTICLE TO REC
% off-axis angle of zr towards particle
thrd = mod(acos(-ko_rf(:, 3)) * 180 / pi + 360, 360) ;
% Reciever orientation - azimuth
phrd = mod(atan2(-ko_rf(:, 2), -ko_rf(:, 1)) * 180 / pi + 360, 360) ;
% we take mod to make sure the angle is between 0 and 360

% propagation vector in image receiver antenna system
koI_rIf = (TgrI * koI')' ;

% FROM PARTICLE TO IM of REC
% off-axis angle of zr towards particle
thrdI = mod(acos(-koI_rIf(:, 3)) * 180 / pi + 360, 360) ;
% Reciever orientation - azimuth
phrdI = mod(atan2(-koI_rIf(:, 2), -koI_rIf(:, 1)) * 180 / pi + 360, 360) ; 
% we take mod to make sure the angle is between 0 and 360

%% Saving Incidence Angles

disp('Saving incidence angles...')
tic;
pathname = SimulationFolders.getInstance.incidence;
filenamex = strcat('thid_', filename) ;
writeVar(pathname, filenamex, thid)
filenamex = strcat('thidI_', filename) ;
writeVar(pathname, filenamex, thidI)
filenamex = strcat('phid_', filename) ;
writeVar(pathname, filenamex, phid)
filenamex = strcat('phidI_', filename) ;
writeVar(pathname, filenamex, phidI)
toc


%% Saving Scattering Angles

disp('Saving scattering angles...')
tic;
pathname = SimulationFolders.getInstance.scattering ;
filenamex = strcat('thsd_', filename) ;
writeVar(pathname, filenamex, thsd)
filenamex = strcat('thsdI_', filename) ;
writeVar(pathname, filenamex, thsdI)
filenamex = strcat('phsd_', filename) ;
writeVar(pathname, filenamex, phsd)
filenamex = strcat('phsdI_', filename) ;
writeVar(pathname, filenamex, phsdI)
toc

%% Saving Receiver Observation Angles

disp('Saving observation angles...')
tic;
pathname = SimulationFolders.getInstance.observation ;
filenamex = strcat('thrd_', filename) ;
writeVar(pathname, filenamex, thrd)
filenamex = strcat('thrdI_', filename) ;
writeVar(pathname, filenamex, thrdI)
filenamex = strcat('phrd_', filename) ;
writeVar(pathname, filenamex, phrd)
filenamex = strcat('phrdI_', filename) ;
writeVar(pathname, filenamex, phrdI)
toc


%%  Saving particle distance to antennas and image antennas
disp('Saving distances...')
tic;
pathname = SimulationFolders.getInstance.distance ;
filenamex = strcat('ri_', filename) ;
writeVar(pathname, filenamex, ri)
filenamex = strcat('riI_', filename) ;
writeVar(pathname, filenamex, riI)
filenamex = strcat('ro_', filename) ;
writeVar(pathname, filenamex, ro)
filenamex = strcat('roI_', filename) ;
writeVar(pathname, filenamex, roI)
toc


end


function result = isPointInEllipse( position, semiMajorAxis, semiMinorAxis, ellipseCenter, ellipseAngle) 

    x = position(1);
    y = position(2);
    h = ellipseCenter(1);
    k = ellipseCenter(2);

    result = ( cos(ellipseAngle)*(x - h) + sin(ellipseAngle)*(y - k) )^2 / semiMajorAxis^2 + ( sin(ellipseAngle)*(x - h) - cos(ellipseAngle)*(y - k) )^2 / semiMinorAxis^2 <= 1.0;

end