
function calcFScatAmp(rr)


%% GET GLOBAL PARAMETERS
% Vegetation Parameters
vegetation_method = SimParams.getInstance.vegetation_method;
TYPKND = VegParams.getInstance.TYPKND;
scat_cal_veg = VegParams.getInstance.scat_cal_veg;
epsr = VegParams.getInstance.epsr;
parm1_deg = VegParams.getInstance.parm1_deg;
parm2_deg = VegParams.getInstance.parm2_deg;
dim1_m = VegParams.getInstance.dim1_m;
dim3_m = VegParams.getInstance.dim3_m;


%% Read Particle Positions
% Option #1 - Uniform Layered Vegetation. - current implementation
% Option #2 - Virtual Vegetation


%% Option 2
[Nlayer, Ntype] = size(TYPKND) ;

for ii = 1 : Nlayer
    
    for jj = 1 : Ntype
        
        Nkind = TYPKND(ii, jj) ;
        
        for kk = 1 : Nkind
            
            filename = strcat('R', num2str(rr), '_L', num2str(ii), '_T',...
                num2str(jj), '_K', num2str(kk)) ;
            disp(filename)
            
            if scat_cal_veg(kk, jj, ii) == 1
                
                if strcmp( vegetation_method, Constants.veg_methods.HOMOGENOUS )
                    
                    EPS = epsr(kk, jj, ii) ;
                    PROB = [ parm1_deg(kk, jj, ii), parm2_deg(kk, jj, ii) ] ;
                    % TO-DO: Should be revised for the case leaves are
                    % scatterer
                    RAD = dim1_m(kk, jj, ii) ;
                    LEN = dim3_m(kk, jj, ii) ; % length  % for cylinders (trunks and branches)
                    
                    fScatAmp(filename, PROB, RAD, LEN, EPS ) ;
                        
                elseif strcmp( vegetation_method, Constants.veg_methods.VIRTUAL )
                    fScatAmp(filename);
                end
                
                disp('done...')
            else
                disp('skipped...')
            end

        end % Nkind

    end % Ntype

end % Nlayer
    
% end % realization

end

function fScatAmp(filename, prob, RAD, LEN, epsr)

%% GET GLOBAL DIRECTORIES
dir_position = SimulationFolders.getInstance.position;
dir_incidence = SimulationFolders.getInstance.incidence;
dir_scattering = SimulationFolders.getInstance.scattering;
dir_fscat = SimulationFolders.getInstance.fscat;


%% GET GLOBAL PARAMETERS
% Constant Structs
scattererParamsStruct = Constants.scattererParamsStruct;
% Satellite Parameters
f_MHz = SatParams.getInstance.f_MHz;


f_Hz = f_MHz * Constants.MHz2Hz ;

%% Reading the particle positions
pP = readVar( dir_position, filename ) ;
[N, ~] = size(pP) ;

%% Bistatic scattering amplitudes
BISTATIC1 = zeros(N, 2, 2) ;   % fpq(o_a^+, i_a^-)     - direct direct
BISTATIC2 = zeros(N, 2, 2) ;   % fpq(oI_a^-, i_a^-)    - direcr reflected
BISTATIC3 = zeros(N, 2, 2) ;   % fpq(o_a^+, iI_a^+)    - reflected direct
BISTATIC4 = zeros(N, 2, 2) ;   % fpq(oI_a^-, iI_a^+)   - reflected reflected

%% Reading Incidence Angles

disp('Reading incidence angles...')
tic;
filenamex = strcat('thid_', filename) ;
thid = readVar(dir_incidence, filenamex) ;
filenamex = strcat('thidI_', filename) ;
thidI = readVar(dir_incidence, filenamex) ;
filenamex = strcat('phid_', filename) ;
phid = readVar(dir_incidence, filenamex) ;
filenamex = strcat('phidI_', filename) ;
phidI = readVar(dir_incidence, filenamex) ;
toc

%% Reading Scattering Angles

disp('Reading scattering angles...')
tic;
filenamex = strcat('thsd_', filename) ;
thsd = readVar(dir_scattering, filenamex) ;
filenamex = strcat('thsdI_', filename) ;
thsdI = readVar(dir_scattering, filenamex) ;
filenamex = strcat('phsd_', filename) ;
phsd = readVar(dir_scattering, filenamex) ;
filenamex = strcat('phsdI_', filename) ;
phsdI = readVar(dir_scattering, filenamex) ;
toc


%% Scattering Amplitudes
tic ;
disp('Calculating scattering amplitudes. . .')
for ii = 1 : N
    
    % Virtual vegetation
    if (nargin == 1)
        % Get the down angle from zenith of the scatterer 
        thcyl = pP(ii, scattererParamsStruct.downAngle);  % Down angle 4

        % Get the azimuth angle of the scatterer 
        phcyl = pP(ii, scattererParamsStruct.azimuthAngle);  % Azimuth angle

        % TO-DO: If leaf is counted as a scatterer in the future, change it!
        % Get the radius and length of the scatterer
        RAD = (pP(ii, scattererParamsStruct.dim1) + ...
            pP(ii, scattererParamsStruct.dim2)) / 2 ;  % Radius
        LEN = pP(ii, scattererParamsStruct.dim3 );  % Length
        epsr = pP(ii, scattererParamsStruct.epsrRe) + ...
            1i * pP(ii, scattererParamsStruct.epsrIm);  % Dielectric Constant
   
    % Homogenous vegetation
    else
        % Generating an angle from the uniform disribution [th1, th2]
        s = prob(1, 1) ;
        e = prob(1, 2) ;
        thcyl = s + (e - s) * rand(1) ;

        % Generating an angle from the uniform disribution [0, 360]
        s = 0 ;
        e = 360 ;
        phcyl = s + (e - s) * rand(1) ;
    end
    
    tin = thid(1, 1) * Constants.deg2rad ;
    pin = phid(1, 1) * Constants.deg2rad ;
    ts = thsd(ii, 1) * Constants.deg2rad ;
    ps = phsd(ii, 1) * Constants.deg2rad ;

    % TO-DO: If leaf is counted as a scatterer in the future, add ellipse!
    % dd : direct - direct
    x = CYLINDER(tin ,pin, ts, ps, thcyl, phcyl, f_Hz, RAD, LEN, epsr) ;
    
    BISTATIC1(ii, 1, 1) = x(4, 1) ; % VV
    BISTATIC1(ii, 1, 2) = x(2, 1) ; % VH    
    BISTATIC1(ii, 2, 1) = x(3, 1) ; % HV
    BISTATIC1(ii, 2, 2) = x(1, 1) ; % HH
    
    tin = thid(1, 1) * Constants.deg2rad ;
    pin = phid(1, 1) * Constants.deg2rad ;
    ts = thsdI(ii, 1) * Constants.deg2rad ;
    ps = phsdI(ii, 1) * Constants.deg2rad ;

    % rd : reflected - direct
    x = CYLINDER(tin ,pin, ts, ps, thcyl, phcyl, f_Hz, RAD, LEN, epsr) ;
    
    BISTATIC2(ii, 1, 1) = x(4, 1) ; % VV
    BISTATIC2(ii, 1, 2) = x(2, 1) ; % VH    
    BISTATIC2(ii, 2, 1) = x(3, 1) ; % HV
    BISTATIC2(ii, 2, 2) = x(1, 1) ; % HH
  
    tin = thidI(1, 1) * Constants.deg2rad ;
    pin = phidI(1, 1) * Constants.deg2rad ;
    ts = thsd(ii, 1) * Constants.deg2rad ;
    ps = phsd(ii, 1) * Constants.deg2rad ;

    % dr : direct - reflected
    x = CYLINDER(tin ,pin, ts, ps, thcyl, phcyl, f_Hz, RAD, LEN, epsr) ;
    
    BISTATIC3(ii, 1, 1) = x(4, 1) ; % VV
    BISTATIC3(ii, 1, 2) = x(2, 1) ; % VH    
    BISTATIC3(ii, 2, 1) = x(3, 1) ; % HV
    BISTATIC3(ii, 2, 2) = x(1, 1) ; % HH
   
    tin = thidI(1, 1) * Constants.deg2rad ;
    pin = phidI(1, 1) * Constants.deg2rad ;
    ts = thsdI(ii, 1) * Constants.deg2rad ;
    ps = phsdI(ii, 1) * Constants.deg2rad ;
    
    % rr : reflected - reflected    
    x = CYLINDER(tin ,pin, ts, ps, thcyl, phcyl, f_Hz, RAD, LEN, epsr) ;
    
    BISTATIC4(ii, 1, 1) = x(4, 1) ; % VV
    BISTATIC4(ii, 1, 2) = x(2, 1) ; % VH    
    BISTATIC4(ii, 2, 1) = x(3, 1) ; % HV
    BISTATIC4(ii, 2, 2) = x(1, 1) ; % HH
        
end
toc

%% Saving scattering amplitudes..
disp('Saving scattering amplitudes...')
tic ;
filenamex = strcat('BISTATIC1_', filename) ;
writeComplexVar(dir_fscat, filenamex, BISTATIC1)
filenamex = strcat('BISTATIC2_', filename) ;
writeComplexVar(dir_fscat, filenamex, BISTATIC2)
filenamex = strcat('BISTATIC3_', filename) ;
writeComplexVar(dir_fscat, filenamex, BISTATIC3)
filenamex = strcat('BISTATIC4_', filename) ;
writeComplexVar(dir_fscat, filenamex, BISTATIC4)
toc

end

