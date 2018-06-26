function generateHomScatPositions( ind_realization )
 
%% GET GLOBAL DIRECTORIES
dir_config = SimulationFolders.getInstance.config;


%% GET GLOBAL PARAMETERS
% Vegetation Parameters
dim_layers_m = VegParams.getInstance.dim_layers_m;
scat_cal_veg = VegParams.getInstance.scat_cal_veg;
TYPKND = VegParams.getInstance.TYPKND;
dsty = VegParams.getInstance.dsty;
dim3_m = VegParams.getInstance.dim3_m;


%% READ META-DATA
% All Positions
AllPoints_m = readVar(dir_config, 'AllPoints_m') ;
pT_m = AllPoints_m(:, 1) ;          % Transmitter position
ht = pT_m(3) ;                    % Transmitter height


%% CALCULATIONS

% Layer parameters
[Nlayer, Ntype] = size(TYPKND) ;

%% Calculations...
D2 = 0 ;

for ii = 1 : Nlayer
    
    D1 = D2 ;               % Layer Top
    D2 = D1 + dim_layers_m(ii, 1) ;    % Layer Bottom
    
    for jj = 1 : Ntype
        
        Nkind = TYPKND(ii, jj) ;
        
        for kk = 1 : Nkind
            
            filename = strcat('R', num2str(ind_realization), '_L', num2str(ii), '_T',...
                num2str(jj), '_K', num2str(kk)) ;
            disp(filename)
                      
            % If the particle is a scatterer
            if scat_cal_veg(kk, jj, ii) == 1 
                
                tic ;
            
                RHO = dsty(kk, jj, ii) ;
                LEN = dim3_m(kk, jj, ii) ;
                
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


%% Generate Scatterer positions
function scatPosition(ht, ds, de, rho, LEN, filename, ifT)

%% GET GLOBAL DIRECTORIES
dir_config = SimulationFolders.getInstance.config;
dir_position = SimulationFolders.getInstance.position;
dir_fzones = SimulationFolders.getInstance.fzones;


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nfz = SimParams.getInstance.Nfz;
% Receiver Parameters
hr_m = RecParams.getInstance.hr_m ;     % Receiver height
% Vegetation Parameters
TYPES = VegParams.getInstance.TYPES;
dim_layers_m = VegParams.getInstance.dim_layers_m;


%% READ META-DATA
Tgs = readVar( dir_config, 'Tgs' ) ;   % G -> SP


%% CALCULATIONS
% Layer Depth and Height
d = sum(dim_layers_m) ;        % thickness of layer in meter

hr_m = hr_m - d ;         % heigth of antenna from veg. layer top in meter
ht = ht - d ;         % heigth of transmitter from veg. layer top in meter

% Ellipse on the layer where z = -ds (layer top)
[~, ~, ax1, by1] = calcFresnelZones(ht + ds, hr_m + ds) ;
aas = ax1 ; % major axis
bbs = by1 ; % minor axis

% Ellipse on the layer where z = -de (layer bottom)
[~, ~, ax1, by1] = calcFresnelZones(ht + de, hr_m + de) ;
aae = ax1 ; % major axis
bbe = by1 ; % minor axis

% Illuminated Volume
As = pi * aas .* bbs ; % area on top layer
Ae = pi * aae .* bbe ; % area on bottom layer

% TO-DO: Finalize Homogenous volume filtering method
VOL = (de - ds) * (As + Ae) / 2 ;

% Above VOL produces the same volume results with the following 
% aas_r = sqrt(pi) / 2 * aas ; % major side
% bbs_r = sqrt(pi) / 2 * bbs ; % minor side
% aae_r = sqrt(pi) / 2 * aae ; % major side
% bbe_r = sqrt(pi) / 2 * bbe ; % minor side
% As_r = 4 * aas_r .* bbs_r ; % area on top layer
% Ae_r = 4 * aae_r .* bbe_r ; % area on bottom layer
% VOL_r = (de - ds) * (As_r + Ae_r) / 2 ;


%% GENERATE RANDOM POSITIONS
% Number of particles
Npart = ceil(VOL * rho) ;
disp(strcat('# of kind in Fresnel Zones :', num2str(Npart)))

% Positions of particles (xp_sf, yp_sf, zp)
xp_sf = zeros(Npart(end), 1) ;
yp_sf = xp_sf ;
zp = xp_sf ;

% Generate particles w.r.t. each Fresnel zone
N2 = 1 ;
for fz = 1 : Nfz
    
    N1 = N2 ;
    N2 = Npart(fz) ;
    
    for ii = N1 : N2
        
        if ifT == TYPES.T
            % The trunk is attached to the ground
            % Calculate the position for center of gravity
            zp(ii, 1) = -de + LEN / 2 ;            
        else            
            % Generating zp between ds (layer top) and de (layer bottom)
            % Calculate the position for center of gravity
            s = -ds ; e = -de ;
            zp(ii, 1) = s + (e - s) * rand(1) ;
        end
        
        % TO-DO: Finalize Homogenous volume filtering method
        [S0x, x1, ax1, by1] = calcFresnelZones(ht + abs(zp(ii, 1)), hr_m + abs(zp(ii, 1))) ;

        % outer ellipse
        a2 = ax1(fz) ;  % semi-major axis
        b2 = by1(fz) ;  % semi-minor axis
        
        % inner ellipse
        if fz == 1   % first Fresnel zone         
            a1 = 0 ;  
            b1 = 0 ; 
        else
            a1 = ax1(fz - 1) ;
            b1 = by1(fz - 1) ;            
        end 
        
        % TO-DO: Finalize Homogenous Fresnel Zone representation method:
        % Rectanbular OR Elliptical?
        a1_r = sqrt(pi) / 2 * a1 ;
        b1_r = sqrt(pi) / 2 * b1 ;
        a2_r = sqrt(pi) / 2 * a2 ;  % semi-major side
        b2_r = sqrt(pi) / 2 * b2 ;  % semi-minor side
  
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
        
    end
    
end

% Convert zp to absolute z-values (that increases upwards starting from the zero-ground) from relative values within the
% vegetation layer (that regards the vegetation top as zero, and goes negative downwards)
zp = zp + d;

% The position of particle center of gravity in specular (local) frame
pP_sf = [xp_sf, yp_sf, zp] ;

% Convert specular frame positions to ground frame
pP = (Tgs' * pP_sf')';  % Tsg = Tgs'  (SP -> G)


%% SAVE ALL
disp('saving...')

% Save particle center  of gravity positions in the ground frame with 
% absolute z-values starting from zero-ground
writeVar( dir_position, filename, pP) ;

% Number of particles per each Fresnel zone with a cumulative array
writeVar( dir_fzones, filename, Npart) ;

%% Incident, scattering angles
calcScatAngles(filename, pP) ;


end

