%% Mehmet Kurum
% 03/13/2017
% modified - 11/09/2017

% Diffuse (Incoherent) Term
function diffuseTerm(Nr)

% Get global parameters
Nfz = SimParams.getInstance.Nfz;
fMHz = SatParams.getInstance.fMHz;
EIRP = convertDecibelToNatural( SatParams.getInstance.EIRP_dB );
G0r = convertDecibelToNatural( RecParams.getInstance.G0r_dB );
VSM = GndParams.getInstance.VSM( ParamsManager.index_VSM );
RMSH = GndParams.getInstance.RMSH( ParamsManager.index_RMSH );
scat_cal_veg = VegParams.getInstance.scat_cal_veg ;
TYPKND = VegParams.getInstance.TYPKND;
d = sum( VegParams.getInstance.dim_layers ) ;           % thickness of layer in meter

 
%% Local
% Nr : Nth Realization (different than global Nr)


%% Shifting coordinate system to the top of the vegetation
d_layer = [0; 0; d] ;

%% Positions relative to ground
% pT: Transmitter, pS2: specular point, pR2: receiver, pG2: ground (reference), pBr2:boresight,
% pCr2: center of footprint, pSc2: center of fresnel zone
% AllPoints = [pT, pTI, pS2, pR, pRI, pG2, pBr2, pCr2, pSc2] ;
filenamex = 'AllPoints' ;
AllPoints = readVar(SimulationFolders.getInstance.config, filenamex) ;

pT = AllPoints(:, 1) ; % - d_layer ;        % Transmitter with respect to veg top
pS2 = AllPoints(:, 3) ;       % Specular point  (at top of vegetation)
pR = AllPoints(:, 4) ; % - d_layer ;        % Receiver with respect to veg top

%% Slant Range - relative to ground
ST = pS2 - pT ;         % Transmitter to Specular
r_st = vectorMagnitude(ST) ;  % slant range

RS = pR - pS2 ;         % Specular to Receiver
r_sr = vectorMagnitude(RS) ;  % slant range

f0hz = fMHz * Constants.MHz2Hz ;
lambda = Constants.c / f0hz ;     % Wavelength
k0 = 2 * pi * f0hz / Constants.c ;    % Wave number


%% Factor Ki
K = 1i * sqrt(EIRP) * sqrt(G0r) * lambda / (4 * pi) ;

Ki = K * exp(1i * k0 * (r_st + r_sr)) / (r_st * r_sr) ; 


%% Layer Parameters
TYPKND %#ok<NOPRT>
[Nlayer, Ntype] = size(TYPKND) ;
NkindMax = max(max(TYPKND)) ;
% sTYPKND = sum(TYPKND) ;
% Ntype = length(sTYPKND(sTYPKND ~= 0)) ; % L, B, T

%% Initializing...
% per kind
b4_inc1_t1 = zeros(2, Nfz, NkindMax, Ntype, Nlayer) ;  % dd
b4_inc2_t1 = b4_inc1_t1 ;                                   % rd
b4_inc3_t1 = b4_inc1_t1 ;                                   % dr
b4_inc4_t1 = b4_inc1_t1 ;                                   % rr
b4_inc1_t2 = zeros(2, Nfz, NkindMax, Ntype, Nlayer) ;
b4_inc2_t2 = b4_inc1_t2 ;
b4_inc3_t2 = b4_inc1_t2 ;
b4_inc4_t2 = b4_inc1_t2 ;

% per type
b3_inc1_t1 = zeros(2, Nfz, Ntype, Nlayer) ;
b3_inc2_t1 = b3_inc1_t1 ;
b3_inc3_t1 = b3_inc1_t1 ;
b3_inc4_t1 = b3_inc1_t1 ;
b3_inc1_t2 = zeros(2, Nfz, Ntype, Nlayer) ;
b3_inc2_t2 = b3_inc1_t2 ;
b3_inc3_t2 = b3_inc1_t2 ;
b3_inc4_t2 = b3_inc1_t2 ;

% per layer
b2_inc1_t1 = zeros(2, Nfz, Nlayer) ;
b2_inc2_t1 = b2_inc1_t1 ;
b2_inc3_t1 = b2_inc1_t1 ;
b2_inc4_t1 = b2_inc1_t1 ;
b2_inc1_t2 = zeros(2, Nfz, Nlayer) ;
b2_inc2_t2 = b2_inc1_t2 ;
b2_inc3_t2 = b2_inc1_t2 ;
b2_inc4_t2 = b2_inc1_t2 ;

% medium
b1_inc1_t1 = zeros(2, Nfz) ;
b1_inc2_t1 = b1_inc1_t1 ;
b1_inc3_t1 = b1_inc1_t1 ;
b1_inc4_t1 = b1_inc1_t1 ;
b1_inc1_t2 = zeros(2, Nfz) ;
b1_inc2_t2 = b1_inc1_t2 ;
b1_inc3_t2 = b1_inc1_t2 ;
b1_inc4_t2 = b1_inc1_t2 ;

%%
% per kind
P4_inc1_t1 = zeros(4, Nfz, NkindMax, Ntype, Nlayer) ;  % dd
P4_inc2_t1 = P4_inc1_t1 ;                                   % rd
P4_inc3_t1 = P4_inc1_t1 ;                                   % dr
P4_inc4_t1 = P4_inc1_t1 ;                                   % rr
P4_inc1_t2 = zeros(4, Nfz, NkindMax, Ntype, Nlayer) ;
P4_inc2_t2 = P4_inc1_t2 ;
P4_inc3_t2 = P4_inc1_t2 ;
P4_inc4_t2 = P4_inc1_t2 ;

% per type
P3_inc1_t1 = zeros(4, Nfz, Ntype, Nlayer) ;
P3_inc2_t1 = P3_inc1_t1 ;
P3_inc3_t1 = P3_inc1_t1 ;
P3_inc4_t1 = P3_inc1_t1 ;
P3_inc1_t2 = zeros(4, Nfz, Ntype, Nlayer) ;
P3_inc2_t2 = P3_inc1_t2 ;
P3_inc3_t2 = P3_inc1_t2 ;
P3_inc4_t2 = P3_inc1_t2 ;

% per layer
P2_inc1_t1 = zeros(4, Nfz, Nlayer) ;
P2_inc2_t1 = P2_inc1_t1 ;
P2_inc3_t1 = P2_inc1_t1 ;
P2_inc4_t1 = P2_inc1_t1 ;
P2_inc1_t2 = zeros(4, Nfz, Nlayer) ;
P2_inc2_t2 = P2_inc1_t2 ;
P2_inc3_t2 = P2_inc1_t2 ;
P2_inc4_t2 = P2_inc1_t2 ;

% medium
P1_inc1_t1 = zeros(4, Nfz) ;
P1_inc2_t1 = P1_inc1_t1 ;
P1_inc3_t1 = P1_inc1_t1 ;
P1_inc4_t1 = P1_inc1_t1 ;
P1_inc1_t2 = zeros(4, Nfz) ;
P1_inc2_t2 = P1_inc1_t2 ;
P1_inc3_t2 = P1_inc1_t2 ;
P1_inc4_t2 = P1_inc1_t2 ;

%% Calculations...

for ii = 1 : Nlayer
    
    for jj = 1 : Ntype
        
        Nkind = TYPKND(ii, jj) ;
        
        if Nkind ~= 0
            
            for kk = 1 : Nkind
                
                % +++++++++++++++++++++++++
                disp('reading...')
                filename = strcat('R', num2str(Nr), '_L', num2str(ii), '_T',...
                    num2str(jj), '_K', num2str(kk)) ;
                disp(filename)                
                
                % If the particle is a scatterer
                if scat_cal_veg(kk, jj, ii) == 1
                    
                    tic ;
                    
                    [b4_inc1_t1(:, :, kk, jj, ii), b4_inc2_t1(:, :, kk, jj, ii), ...
                        b4_inc3_t1(:, :, kk, jj, ii), b4_inc4_t1(:, :, kk, jj, ii), ...
                        b4_inc1_t2(:, :, kk, jj, ii), b4_inc2_t2(:, :, kk, jj, ii), ...
                        b4_inc3_t2(:, :, kk, jj, ii), b4_inc4_t2(:, :, kk, jj, ii), ...
                        P4_inc1_t1(:, :, kk, jj, ii), P4_inc2_t1(:, :, kk, jj, ii), ...
                        P4_inc3_t1(:, :, kk, jj, ii), P4_inc4_t1(:, :, kk, jj, ii), ...
                        P4_inc1_t2(:, :, kk, jj, ii), P4_inc2_t2(:, :, kk, jj, ii), ...
                        P4_inc3_t2(:, :, kk, jj, ii), P4_inc4_t2(:, :, kk, jj, ii)] ...
                        = ScatMech(ii, filename, r_st, r_sr) ;
                    
                    disp('done...')
                    
                    toc
                    
                else
                    disp('Skipped to calculate diffuse term for the particle...')
                end

                
                % +++++++++++++++++++++++++
                % sum over all kinds of each type in each layer
                % 2 by ....
                b3_inc1_t1(:, :, jj, ii) = b3_inc1_t1(:, :, jj, ii) + b4_inc1_t1(:, :, kk, jj, ii) ;
                b3_inc2_t1(:, :, jj, ii) = b3_inc2_t1(:, :, jj, ii) + b4_inc2_t1(:, :, kk, jj, ii) ;
                b3_inc3_t1(:, :, jj, ii) = b3_inc3_t1(:, :, jj, ii) + b4_inc3_t1(:, :, kk, jj, ii) ;
                b3_inc4_t1(:, :, jj, ii) = b3_inc4_t1(:, :, jj, ii) + b4_inc4_t1(:, :, kk, jj, ii) ;
             
                b3_inc1_t2(:, :, jj, ii) = b3_inc1_t2(:, :, jj, ii) + b4_inc1_t2(:, :, kk, jj, ii) ;
                b3_inc2_t2(:, :, jj, ii) = b3_inc2_t2(:, :, jj, ii) + b4_inc2_t2(:, :, kk, jj, ii) ;
                b3_inc3_t2(:, :, jj, ii) = b3_inc3_t2(:, :, jj, ii) + b4_inc3_t2(:, :, kk, jj, ii) ;
                b3_inc4_t2(:, :, jj, ii) = b3_inc4_t2(:, :, jj, ii) + b4_inc4_t2(:, :, kk, jj, ii) ;
                
                % 4 by ....
                P3_inc1_t1(:, :, jj, ii) = P3_inc1_t1(:, :, jj, ii) + P4_inc1_t1(:, :, kk, jj, ii) ;
                P3_inc2_t1(:, :, jj, ii) = P3_inc2_t1(:, :, jj, ii) + P4_inc2_t1(:, :, kk, jj, ii) ;
                P3_inc3_t1(:, :, jj, ii) = P3_inc3_t1(:, :, jj, ii) + P4_inc3_t1(:, :, kk, jj, ii) ;
                P3_inc4_t1(:, :, jj, ii) = P3_inc4_t1(:, :, jj, ii) + P4_inc4_t1(:, :, kk, jj, ii) ;
             
                P3_inc1_t2(:, :, jj, ii) = P3_inc1_t2(:, :, jj, ii) + P4_inc1_t2(:, :, kk, jj, ii) ;
                P3_inc2_t2(:, :, jj, ii) = P3_inc2_t2(:, :, jj, ii) + P4_inc2_t2(:, :, kk, jj, ii) ;
                P3_inc3_t2(:, :, jj, ii) = P3_inc3_t2(:, :, jj, ii) + P4_inc3_t2(:, :, kk, jj, ii) ;
                P3_inc4_t2(:, :, jj, ii) = P3_inc4_t2(:, :, jj, ii) + P4_inc4_t2(:, :, kk, jj, ii) ;
                % +++++++++++++++++++++++++
                
            end % Nkind -- kk
            
            % +++++++++++++++++++++++++
            % sum over all types in each layer
            % 2 by ....
            b2_inc1_t1(:, :, ii) = b2_inc1_t1(:, :, ii) + b3_inc1_t1(:, :, jj, ii) ;
            b2_inc2_t1(:, :, ii) = b2_inc2_t1(:, :, ii) + b3_inc2_t1(:, :, jj, ii) ;
            b2_inc3_t1(:, :, ii) = b2_inc3_t1(:, :, ii) + b3_inc3_t1(:, :, jj, ii) ;
            b2_inc4_t1(:, :, ii) = b2_inc4_t1(:, :, ii) + b3_inc4_t1(:, :, jj, ii) ;
            
            b2_inc1_t2(:, :, ii) = b2_inc1_t2(:, :, ii) + b3_inc1_t2(:, :, jj, ii) ;
            b2_inc2_t2(:, :, ii) = b2_inc2_t2(:, :, ii) + b3_inc2_t2(:, :, jj, ii) ;
            b2_inc3_t2(:, :, ii) = b2_inc3_t2(:, :, ii) + b3_inc3_t2(:, :, jj, ii) ;
            b2_inc4_t2(:, :, ii) = b2_inc4_t2(:, :, ii) + b3_inc4_t2(:, :, jj, ii) ;
            
            % 4 by ....
            P2_inc1_t1(:, :, ii) = P2_inc1_t1(:, :, ii) + P3_inc1_t1(:, :, jj, ii) ;
            P2_inc2_t1(:, :, ii) = P2_inc2_t1(:, :, ii) + P3_inc2_t1(:, :, jj, ii) ;
            P2_inc3_t1(:, :, ii) = P2_inc3_t1(:, :, ii) + P3_inc3_t1(:, :, jj, ii) ;
            P2_inc4_t1(:, :, ii) = P2_inc4_t1(:, :, ii) + P3_inc4_t1(:, :, jj, ii) ;
            
            P2_inc1_t2(:, :, ii) = P2_inc1_t2(:, :, ii) + P3_inc1_t2(:, :, jj, ii) ;
            P2_inc2_t2(:, :, ii) = P2_inc2_t2(:, :, ii) + P3_inc2_t2(:, :, jj, ii) ;
            P2_inc3_t2(:, :, ii) = P2_inc3_t2(:, :, ii) + P3_inc3_t2(:, :, jj, ii) ;
            P2_inc4_t2(:, :, ii) = P2_inc4_t2(:, :, ii) + P3_inc4_t2(:, :, jj, ii) ;
            % +++++++++++++++++++++++++
            
        end % if Nkind ~= 0
        
    end % Ntype -- jj
    
    % +++++++++++++++++++++++++
    % sum over all layers
    % 2 by ....
    b1_inc1_t1(:, :) = b1_inc1_t1(:, :) + b2_inc1_t1(:, :, ii) ;
    b1_inc2_t1(:, :) = b1_inc2_t1(:, :) + b2_inc2_t1(:, :, ii) ;
    b1_inc3_t1(:, :) = b1_inc3_t1(:, :) + b2_inc3_t1(:, :, ii) ;
    b1_inc4_t1(:, :) = b1_inc4_t1(:, :) + b2_inc4_t1(:, :, ii) ;
    
    b1_inc1_t2(:, :) = b1_inc1_t2(:, :) + b2_inc1_t2(:, :, ii) ;
    b1_inc2_t2(:, :) = b1_inc2_t2(:, :) + b2_inc2_t2(:, :, ii) ;
    b1_inc3_t2(:, :) = b1_inc3_t2(:, :) + b2_inc3_t2(:, :, ii) ;
    b1_inc4_t2(:, :) = b1_inc4_t2(:, :) + b2_inc4_t2(:, :, ii) ;
    
    % 4 by ....
    P1_inc1_t1(:, :) = P1_inc1_t1(:, :) + P2_inc1_t1(:, :, ii) ;
    P1_inc2_t1(:, :) = P1_inc2_t1(:, :) + P2_inc2_t1(:, :, ii) ;
    P1_inc3_t1(:, :) = P1_inc3_t1(:, :) + P2_inc3_t1(:, :, ii) ;
    P1_inc4_t1(:, :) = P1_inc4_t1(:, :) + P2_inc4_t1(:, :, ii) ;
    
    P1_inc1_t2(:, :) = P1_inc1_t2(:, :) + P2_inc1_t2(:, :, ii) ;
    P1_inc2_t2(:, :) = P1_inc2_t2(:, :) + P2_inc2_t2(:, :, ii) ;
    P1_inc3_t2(:, :) = P1_inc3_t2(:, :) + P2_inc3_t2(:, :, ii) ;
    P1_inc4_t2(:, :) = P1_inc4_t2(:, :) + P2_inc4_t2(:, :, ii) ;
    % +++++++++++++++++++++++++
    
end % Nlayer -- ii

%% Saving...

% Ki
filename3 = strcat('Ki') ;
writeComplexVar(SimulationFolders.getInstance.freqdiff, filename3, Ki)

% ====================================================
% b1
pathname = strcat(SimulationFolders.getInstance.freqdiff_b1, '\VSM_', num2str( VSM ), '-RMSH_', num2str(RMSH)) ;
filename1 = strcat('b1_inc1_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b1_inc1_t1))
filename1 = strcat('b1_inc2_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b1_inc2_t1))
filename1 = strcat('b1_inc3_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b1_inc3_t1))
filename1 = strcat('b1_inc4_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b1_inc4_t1))

filename1 = strcat('b1_inc1_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b1_inc1_t2))
filename1 = strcat('b1_inc2_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b1_inc2_t2))
filename1 = strcat('b1_inc3_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b1_inc3_t2))
filename1 = strcat('b1_inc4_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b1_inc4_t2))

% b2
pathname = strcat(SimulationFolders.getInstance.freqdiff_b2, '\VSM_', num2str( VSM ), '-RMSH_', num2str(RMSH)) ;
filename1 = strcat('b2_inc1_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b2_inc1_t1))
filename1 = strcat('b2_inc2_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b2_inc2_t1))
filename1 = strcat('b2_inc3_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b2_inc3_t1))
filename1 = strcat('b2_inc4_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b2_inc4_t1))

filename1 = strcat('b2_inc1_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b2_inc1_t2))
filename1 = strcat('b2_inc2_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b2_inc2_t2))
filename1 = strcat('b2_inc3_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b2_inc3_t2))
filename1 = strcat('b2_inc4_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b2_inc4_t2))

% b3
pathname = strcat(SimulationFolders.getInstance.freqdiff_b3, '\VSM_', num2str( VSM ), '-RMSH_', num2str(RMSH)) ;
filename1 = strcat('b3_inc1_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b3_inc1_t1))
filename1 = strcat('b3_inc2_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b3_inc2_t1))
filename1 = strcat('b3_inc3_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b3_inc3_t1))
filename1 = strcat('b3_inc4_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b3_inc4_t1))

filename1 = strcat('b3_inc1_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b3_inc1_t2))
filename1 = strcat('b3_inc2_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b3_inc2_t2))
filename1 = strcat('b3_inc3_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b3_inc3_t2))
filename1 = strcat('b3_inc4_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b3_inc4_t2))

% b4
pathname = strcat(SimulationFolders.getInstance.freqdiff_b4, '\VSM_', num2str( VSM ), '-RMSH_', num2str(RMSH)) ;
filename1 = strcat('b4_inc1_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b4_inc1_t1))
filename1 = strcat('b4_inc2_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b4_inc2_t1))
filename1 = strcat('b4_inc3_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b4_inc3_t1))
filename1 = strcat('b4_inc4_t1', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b4_inc4_t1))

filename1 = strcat('b4_inc1_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b4_inc1_t2))
filename1 = strcat('b4_inc2_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b4_inc2_t2))
filename1 = strcat('b4_inc3_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b4_inc3_t2))
filename1 = strcat('b4_inc4_t2', '_R', num2str(Nr)) ;
writeComplexVar(pathname, filename1, (b4_inc4_t2))
% =====================================================

% P1
pathname = strcat(SimulationFolders.getInstance.freqdiff_P1, '\VSM_', num2str( VSM ), '-RMSH_', num2str(RMSH)) ;
filename1 = strcat('P1_inc1_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P1_inc1_t1))
filename1 = strcat('P1_inc2_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P1_inc2_t1))
filename1 = strcat('P1_inc3_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P1_inc3_t1))
filename1 = strcat('P1_inc4_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P1_inc4_t1))

filename1 = strcat('P1_inc1_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P1_inc1_t2))
filename1 = strcat('P1_inc2_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P1_inc2_t2))
filename1 = strcat('P1_inc3_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P1_inc3_t2))
filename1 = strcat('P1_inc4_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P1_inc4_t2))

% P2
pathname = strcat(SimulationFolders.getInstance.freqdiff_P2, '\VSM_', num2str( VSM ), '-RMSH_', num2str(RMSH)) ;
filename1 = strcat('P2_inc1_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P2_inc1_t1))
filename1 = strcat('P2_inc2_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P2_inc2_t1))
filename1 = strcat('P2_inc3_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P2_inc3_t1))
filename1 = strcat('P2_inc4_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P2_inc4_t1))

filename1 = strcat('P2_inc1_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P2_inc1_t2))
filename1 = strcat('P2_inc2_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P2_inc2_t2))
filename1 = strcat('P2_inc3_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P2_inc3_t2))
filename1 = strcat('P2_inc4_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P2_inc4_t2))

% P3
pathname = strcat(SimulationFolders.getInstance.freqdiff_P3, '\VSM_', num2str( VSM ), '-RMSH_', num2str(RMSH)) ;
filename1 = strcat('P3_inc1_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P3_inc1_t1))
filename1 = strcat('P3_inc2_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P3_inc2_t1))
filename1 = strcat('P3_inc3_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P3_inc3_t1))
filename1 = strcat('P3_inc4_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P3_inc4_t1))

filename1 = strcat('P3_inc1_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P3_inc1_t2))
filename1 = strcat('P3_inc2_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P3_inc2_t2))
filename1 = strcat('P3_inc3_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P3_inc3_t2))
filename1 = strcat('P3_inc4_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P3_inc4_t2))

% P4
pathname = strcat(SimulationFolders.getInstance.freqdiff_P4, '\VSM_', num2str( VSM ), '-RMSH_', num2str(RMSH)) ;
filename1 = strcat('P4_inc1_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P4_inc1_t1))
filename1 = strcat('P4_inc2_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P4_inc2_t1))
filename1 = strcat('P4_inc3_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P4_inc3_t1))
filename1 = strcat('P4_inc4_t1', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P4_inc4_t1))

filename1 = strcat('P4_inc1_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P4_inc1_t2))
filename1 = strcat('P4_inc2_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P4_inc2_t2))
filename1 = strcat('P4_inc3_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P4_inc3_t2))
filename1 = strcat('P4_inc4_t2', '_R', num2str(Nr)) ;
writeVar(pathname, filename1, (P4_inc4_t2))


end

%% Scattering Matrices
function [b_dd_t1, b_rd_t1, b_dr_t1, b_rr_t1,...
    b_dd_t2, b_rd_t2, b_dr_t2, b_rr_t2, ...
    P_dd_t1, P_rd_t1, P_dr_t1, P_rr_t1,...
    P_dd_t2, P_rd_t2, P_dr_t2, P_rr_t2] = ScatMech(layerIndex, filename, r_st, r_sr)


fMHz = SatParams.getInstance.fMHz
Nfz = SimParams.getInstance.Nfz;

dim_layers = VegParams.getInstance.dim_layers;

%% Ground Parameters
disp('Reading ground parameters...')
filenamex = 'G' ;
grnd_par = readVar(SimulationFolders.getInstance.gnd, filenamex) ;
h = grnd_par(1, 1) ;

epsg = grnd_par(1, 2) + ( 1i * grnd_par(1, 3) ) ;

%% Incremental Propagation Constant
disp('Reading Incremental propagation constants...')
filenamex = 'dKz' ;
dKz = readComplexVar(SimulationFolders.getInstance.afsa, filenamex) ;
filenamex = 'ANGDEG' ;
ANGDEG = readVar(SimulationFolders.getInstance.afsa, filenamex) ;

%% Particle Positions
disp('Reading the particle positions...')

pP = readVar(SimulationFolders.getInstance.position, filename) ;

Npart = readVar(SimulationFolders.getInstance.fzones, filename) ;

%% Incidence Angles
disp('Reading incidence angles...')
pathname = SimulationFolders.getInstance.incidence ;
filenamex = strcat('thid_', filename) ;
thid = readVar(pathname, filenamex) ;
filenamex = strcat('thidI_', filename) ;
thidI = readVar(pathname, filenamex) ;

%% Reading Scattering Angles
disp('Reading scattering angles...')
pathname = SimulationFolders.getInstance.scattering ;
filenamex = strcat('thsd_', filename) ;
thsd = readVar(pathname, filenamex) ;
filenamex = strcat('thsdI_', filename) ;
thsdI = readVar(pathname, filenamex) ;
% no need for phsd and phsdI

%% Reading Antenna gain values
disp('Reading antenna values...')
pathname = SimulationFolders.getInstance.ant_real ;
filenamex = strcat('gr_', filename) ;
gr = readComplexVar(pathname, filenamex) ;
filenamex = strcat('grI_', filename) ;
grI = readComplexVar(pathname, filenamex) ;

%% Reading Rotation Matrix values
disp('Reading rotation matrices...')
pathname = SimulationFolders.getInstance.rot_real ;
filenamex = strcat('u_gar_', filename) ;
u_gar2 = readComplexVar(pathname, filenamex) ;
filenamex = strcat('u_garI_', filename) ;
u_garI2 = readComplexVar(pathname, filenamex) ;
pathname = SimulationFolders.getInstance.rot_lookup ;
load([pathname '\u_ts.mat'], 'u_ts') 
load([pathname '\u_tIs.mat'], 'u_tIs')

U_ts = calc_Muller(u_ts) ;
U_tIs = calc_Muller(u_tIs) ;


%% Reading scattering amplitudes..
disp('Reading scattering amplitudes...')
pathname = SimulationFolders.getInstance.fscat;
filenamex = strcat('BISTATIC1_', filename) ;
BISTATIC1 = readComplexVar(pathname, filenamex) ;
filenamex = strcat('BISTATIC2_', filename) ;
BISTATIC2 = readComplexVar(pathname, filenamex) ;
filenamex = strcat('BISTATIC3_', filename) ;
BISTATIC3 = readComplexVar(pathname, filenamex) ;
filenamex = strcat('BISTATIC4_', filename) ;
BISTATIC4 = readComplexVar(pathname, filenamex) ;

%%  Reading particle distances to receiver and image receiver
disp('Reading distances...')
pathname = SimulationFolders.getInstance.distance ;
filenamex = strcat('ro_', filename) ;
ro = readVar(pathname, filenamex) ;
filenamex = strcat('roI_', filename) ;
roI = readVar(pathname, filenamex) ;
filenamex = strcat('ri_', filename) ;
ri = readVar(pathname, filenamex) ;
filenamex = strcat('riI_', filename) ;
riI = readVar(pathname, filenamex) ;

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

%% Wave Number
f = fMHz * Constants.MHz2Hz;  % Frequency points (Hz)
k0 = 2 * pi * f / Constants.c ;    % Wave number

%% Scattering Matrix Calculations
disp(strcat('Number of Scatters:', num2str(Npart)))

%% Transmitter Antenna pattern
g_t = SatParams.getInstance.g_t ; % ideal
G_t = calc_Muller(g_t) ;

%% Transmitter Pol State
e_t1 = SatParams.getInstance.e_t1 ;
e_t2 = SatParams.getInstance.e_t2 ;
%TO-DO: Calculate E_t1 from e_t1?
E_t1 = [1; 0; 0; 0] ;
E_t2 = [0; 0; 0; 1] ;
Q = [1 0 0 0; 0 0 0 1; 0 1 1 0; 0 -1i -1i 0] ;
E_t1 = Q * E_t1 ;
E_t2 = Q * E_t2 ;
%% Intilizing. . . 
b_dd_t1 = zeros(2, Nfz) ;
b_dr_t1 = b_dd_t1 ;
b_rd_t1 = b_dd_t1 ;
b_rr_t1 = b_dd_t1 ;
b_dd_t2 = zeros(2, Nfz) ;
b_dr_t2 = b_dd_t2 ;
b_rd_t2 = b_dd_t2 ;
b_rr_t2 = b_dd_t2 ;

bfz_dd_t1 = zeros(2, Nfz) ;
bfz_dr_t1 = bfz_dd_t1 ;
bfz_rd_t1 = bfz_dd_t1 ;
bfz_rr_t1 = bfz_dd_t1 ;
bfz_dd_t2 = zeros(2, Nfz) ;
bfz_dr_t2 = bfz_dd_t2 ;
bfz_rd_t2 = bfz_dd_t2 ;
bfz_rr_t2 = bfz_dd_t2 ;

%
P_dd_t1 = zeros(4, Nfz) ;
P_dr_t1 = P_dd_t1 ;
P_rd_t1 = P_dd_t1 ;
P_rr_t1 = P_dd_t1 ;
P_dd_t2 = zeros(4, Nfz) ;
P_dr_t2 = P_dd_t2 ;
P_rd_t2 = P_dd_t2 ;
P_rr_t2 = P_dd_t2 ;

Pfz_dd_t1 = zeros(4, Nfz) ;
Pfz_dr_t1 = Pfz_dd_t1 ;
Pfz_rd_t1 = Pfz_dd_t1 ;
Pfz_rr_t1 = Pfz_dd_t1 ;
Pfz_dd_t2 = zeros(4, Nfz) ;
Pfz_dr_t2 = Pfz_dd_t2 ;
Pfz_rd_t2 = Pfz_dd_t2 ;
Pfz_rr_t2 = Pfz_dd_t2 ;

% +++++++++++++++++++++++++
disp('calculating...')
tic ;

N2 = 1 ;
for fz = 1 : Nfz    % Fresnel Zones
    
    N1 = N2 ;
    N2 = Npart(fz) ;
    
    for pp = N1 : N2
        
        %% Calculate Transmissivity Matrices
        dKz_i = squeeze(dKz(:, ANGDEG == round(thid), :)) ;
        dKz_iI = squeeze(dKz(:, ANGDEG == round(180 - thidI), :)) ;
        dKz_o = squeeze(dKz(:, ANGDEG == round(thsd(pp)), :)) ;
        dKz_oI = squeeze(dKz(:, ANGDEG == round(180 - thsdI(pp)), :)) ;
        
        zp = pP(pp, 3) ;
        
        [tp_i, tp_iI, tp_o, tp_oI] = CalcTransMat(layerIndex, zp, dim_layers, dKz_i, dKz_iI, dKz_o, dKz_oI) ;

%         tp_i = [1 0; 0 1] ;
%         tp_iI = [1 0; 0 1] ;
%         tp_o = [1 0; 0 1] ;
%         tp_oI = [1 0; 0 1] ;
        %% Calculate Reflection Coefficient
        thiI = degtorad(180 - thidI) ;        % incident angle
        thsI = degtorad(180 - thsdI(pp)) ;    % scattered angle
        [RGH_iI, RGV_iI, RGH_oI, RGV_oI] = reflectionCoeff(thiI, thsI, epsg, h) ;

        %% Antenna Pattern Matrix
        % 2 X 2
        g_ro = squeeze(gr(pp, :, :)) ;
        g_roI = squeeze(grI(pp, :, :)) ;
% %         g_ro = [1 0; 0 1] ;
% %         g_roI = [1 0; 0 1] ;
        % 4 X 4
        G_ro = calc_Muller(g_ro) ;
        G_roI = calc_Muller(g_roI) ;                  

        %% Rotation Matrix
        % 2 X 2
        u_gar = squeeze(u_gar2(pp, :, :)) ;
        u_garI = squeeze(u_garI2(pp, :, :)) ;
        % 4 X 4
        U_gar = calc_Muller(u_gar) ;
        U_garI = calc_Muller(u_garI) ;
        %% Bistatic Scattering Amplitude
        f_oi = squeeze(BISTATIC1(pp, :, :)) ;
        f_oIi = squeeze(BISTATIC2(pp, :, :)) ;
        f_oiI = squeeze(BISTATIC3(pp, :, :)) ;
        f_oIiI = squeeze(BISTATIC4(pp, :, :)) ;
%         f_oi = [1 0; 0 1] ;
%         f_oIi = [1 0; 0 1] ;
%         f_oiI = [1 0; 0 1] ;
%         f_oIiI = [1 0; 0 1] ;

        %% B-factor        
        B_dd = (exp(+1i * k0 * (ro(pp) + ri(pp))) / (ro(pp) * ri(pp))) ...
            / (exp(+1i * k0 * (r_st + r_sr)) / (r_st * r_sr)) ;
        B_rd = (exp(+1i * k0 * (roI(pp) + ri(pp))) / (roI(pp) * ri(pp))) ...
            / (exp(+1i * k0 * (r_st + r_sr)) / (r_st * r_sr)) ;
        B_dr = (exp(+1i * k0 * (ro(pp) + riI(pp))) / (ro(pp) * riI(pp))) ...
            / (exp(+1i * k0 * (r_st + r_sr)) / (r_st * r_sr)) ;
        B_rr = (exp(+1i * k0 * (roI(pp) + riI(pp))) / (roI(pp) * riI(pp))) ...
            / (exp(+1i * k0 * (r_st + r_sr)) / (r_st * r_sr)) ;      
%         B_dd = 1 ; %  
%         B_rd = 1 ; % 
%         B_dr = 1 ; % 
%         B_rr = 1 ; %
            
        %% Ground Reflection Matrices ===============
        rg_iI = [RGV_iI 0; 0 RGH_iI] ;
        rg_oI = [RGV_oI 0; 0 RGH_oI] ;            
% %             rg_iI = [1 0; 0 1] ;
% %             rg_oI = [1 0; 0 1] ;

        %% Calculate Scattering Matrices =========================
        % 2 X 2
        s_dd = tp_o * f_oi * tp_i ;
        s_rd = tp_oI * rg_oI * f_oIi * tp_i ;
        s_dr = tp_o * f_oiI * rg_iI * tp_iI ;
        s_rr = tp_oI * rg_oI * f_oIiI * rg_iI * tp_iI ;
        % 4 X 4
        sig_dd = 4 * pi * calc_Muller(s_dd) ;
        sig_rd = 4 * pi * calc_Muller(s_rd) ;
        sig_dr = 4 * pi * calc_Muller(s_dr) ;
        sig_rr = 4 * pi * calc_Muller(s_rr) ;

        %% ======== summing over indivdual contributions=================
        % 2 X 1
        bfz_dd_t1(:, fz) = bfz_dd_t1(:, fz) + B_dd * g_ro * u_gar * s_dd * u_ts * g_t * e_t1 ;
        bfz_rd_t1(:, fz) = bfz_rd_t1(:, fz) + B_rd * g_roI * u_garI * s_rd * u_ts * g_t * e_t1 ;
        bfz_dr_t1(:, fz) = bfz_dr_t1(:, fz) + B_dr * g_ro * u_gar * s_dr * u_tIs * g_t * e_t1 ;
        bfz_rr_t1(:, fz) = bfz_rr_t1(:, fz) + B_rr * g_roI * u_garI * s_rr * u_tIs * g_t * e_t1 ;

        bfz_dd_t2(:, fz) = bfz_dd_t2(:, fz) + B_dd * g_ro * u_gar * s_dd * u_ts * g_t * e_t2 ;
        bfz_rd_t2(:, fz) = bfz_rd_t2(:, fz) + B_rd * g_roI * u_garI * s_rd * u_ts * g_t * e_t2 ;
        bfz_dr_t2(:, fz) = bfz_dr_t2(:, fz) + B_dr * g_ro * u_gar * s_dr * u_tIs * g_t * e_t2 ;
        bfz_rr_t2(:, fz) = bfz_rr_t2(:, fz) + B_rr * g_roI * u_garI * s_rr * u_tIs * g_t * e_t2 ;

        % 4 X 1
        Pfz_dd_t1(:, fz) = Pfz_dd_t1(:, fz) + abs(B_dd) ^ 2 * G_ro * U_gar * sig_dd * U_ts * G_t * E_t1 ;
        Pfz_rd_t1(:, fz) = Pfz_rd_t1(:, fz) + abs(B_rd) ^ 2 * G_roI * U_garI * sig_rd * U_ts * G_t * E_t1 ;
        Pfz_dr_t1(:, fz) = Pfz_dr_t1(:, fz) + abs(B_dr) ^ 2 * G_ro * U_gar * sig_dr * U_tIs * G_t * E_t1 ;
        Pfz_rr_t1(:, fz) = Pfz_rr_t1(:, fz) + abs(B_rr) ^ 2 * G_roI * U_garI * sig_rr * U_tIs * G_t * E_t1 ;

        Pfz_dd_t2(:, fz) = Pfz_dd_t2(:, fz) + abs(B_dd) ^ 2 * G_ro * U_gar * sig_dd * U_ts * G_t * E_t2 ;
        Pfz_rd_t2(:, fz) = Pfz_rd_t2(:, fz) + abs(B_rd) ^ 2 * G_roI * U_garI * sig_rd * U_ts * G_t * E_t2 ;
        Pfz_dr_t2(:, fz) = Pfz_dr_t2(:, fz) + abs(B_dr) ^ 2 * G_ro * U_gar * sig_dr * U_tIs * G_t * E_t2 ;
        Pfz_rr_t2(:, fz) = Pfz_rr_t2(:, fz) + abs(B_rr) ^ 2 * G_roI * U_garI * sig_rr * U_tIs * G_t * E_t2 ;
        
    end % Npart -- pp
    
    %% === summing over fresnel zones =================    
    for ii = 1 : fz
        % 2 X 1
        b_dd_t1(:, fz) = b_dd_t1(:, fz) + bfz_dd_t1(:, ii) ;
        b_rd_t1(:, fz) = b_rd_t1(:, fz) + bfz_rd_t1(:, ii) ;
        b_dr_t1(:, fz) = b_dr_t1(:, fz) + bfz_dr_t1(:, ii) ;
        b_rr_t1(:, fz) = b_rr_t1(:, fz) + bfz_rr_t1(:, ii) ;
        
        b_dd_t2(:, fz) = b_dd_t2(:, fz) + bfz_dd_t2(:, ii) ;
        b_rd_t2(:, fz) = b_rd_t2(:, fz) + bfz_rd_t2(:, ii) ;
        b_dr_t2(:, fz) = b_dr_t2(:, fz) + bfz_dr_t2(:, ii) ;
        b_rr_t2(:, fz) = b_rr_t2(:, fz) + bfz_rr_t2(:, ii) ;
        
        % 4 X 1
        P_dd_t1(:, fz) = P_dd_t1(:, fz) + Pfz_dd_t1(:, ii) ;
        P_rd_t1(:, fz) = P_rd_t1(:, fz) + Pfz_rd_t1(:, ii) ;
        P_dr_t1(:, fz) = P_dr_t1(:, fz) + Pfz_dr_t1(:, ii) ;
        P_rr_t1(:, fz) = P_rr_t1(:, fz) + Pfz_rr_t1(:, ii) ;
        
        P_dd_t2(:, fz) = P_dd_t2(:, fz) + Pfz_dd_t2(:, ii) ;
        P_rd_t2(:, fz) = P_rd_t2(:, fz) + Pfz_rd_t2(:, ii) ;
        P_dr_t2(:, fz) = P_dr_t2(:, fz) + Pfz_dr_t2(:, ii) ;
        P_rr_t2(:, fz) = P_rr_t2(:, fz) + Pfz_rr_t2(:, ii) ;
        
    end
    %%=================================================
    
end % Nfz -- fz
toc
       
end

%% Calculate Transmissivity Matrices

function [tp_i, tp_iI, tp_o, tp_oI] = CalcTransMat(layerIndex, zp, D, dKz_i, dKz_iI, dKz_o, dKz_oI)
% zp: particle center of gravity position with absolute z-values starting 
% from zero-ground  


% Convert absoulte z-alues to relative ones within the vegetation layer
% where vegetation top is zero.
zp = zp - sum(D);

ArgH_i = 0 ;        ArgV_i = 0 ;
ArgH_iI = 0 ;       ArgV_iI = 0 ;
ArgH_o = 0 ;        ArgV_o = 0 ;
ArgH_oI = 0 ;       ArgV_oI = 0 ;

Nlayer = length(D) ;
D2 = 0 ;
for ii = 1 : Nlayer
    
    D1 = D2 ;               % Layer Top
    D2 = D1 + D(ii, 1) ;    % Layer Bottom
    
    if ii < layerIndex % attenuation above layerIndex-th layer till layer top
        % ==========================================
        ArgH_i = ArgH_i + dKz_i(1, ii) * D(ii, 1) ;
        ArgV_i = ArgV_i + dKz_i(2, ii) * D(ii, 1) ;
        
        ArgH_o = ArgH_o + dKz_o(1, ii) * D(ii, 1) ;
        ArgV_o = ArgV_o + dKz_o(2, ii) * D(ii, 1) ;
        % =========================================
    elseif ii == layerIndex % partial attenuation within layerIndex-th layer
        % ==================================================
        ArgH_i = ArgH_i + dKz_i(1, ii) * (abs(zp) - D1) ;
        ArgV_i = ArgV_i + dKz_i(2, ii) * (abs(zp) - D1) ;
        
        ArgH_o = ArgH_o + dKz_o(1, ii) * (abs(zp) - D1) ;
        ArgV_o = ArgV_o + dKz_o(2, ii) * (abs(zp) - D1) ;
        
        ArgH_iI = ArgH_iI + dKz_iI(1, ii) * (D2 - abs(zp)) ;
        ArgV_iI = ArgV_iI + dKz_iI(2, ii) * (D2 - abs(zp)) ;
        
        ArgH_oI = ArgH_oI + dKz_oI(1, ii) * (D2 - abs(zp)) ;
        ArgV_oI = ArgV_oI + dKz_oI(2, ii) * (D2 - abs(zp)) ;
        % =================================================
    elseif ii > layerIndex % attenuation below layerIndex-th layer till ground
        % =============================================
        ArgH_iI = ArgH_iI + dKz_iI(1, ii) * D(ii, 1) ;
        ArgV_iI = ArgV_iI + dKz_iI(2, ii) * D(ii, 1) ;
        
        ArgH_oI = ArgH_oI + dKz_oI(1, ii) * D(ii, 1) ;
        ArgV_oI = ArgV_oI + dKz_oI(2, ii) * D(ii, 1) ;
        % =============================================
    end
    % all-layer [attenuation over the image layer]
    % ============================================
    ArgH_iI = ArgH_iI + dKz_iI(1, ii) * D(ii, 1) ;
    ArgV_iI = ArgV_iI + dKz_iI(2, ii) * D(ii, 1) ;
    
    ArgH_oI = ArgH_oI + dKz_oI(1, ii) * D(ii, 1) ;
    ArgV_oI = ArgV_oI + dKz_oI(2, ii) * D(ii, 1) ;
    % ===========================================
end

%% Forming Transmissivity Matrices
tp_i = [exp(+1i * ArgV_i) 0; 0 exp(+1i * ArgH_i)] ;
tp_iI = [exp(+1i * ArgV_iI) 0; 0 exp(+1i * ArgH_iI)] ;
tp_o = [exp(+1i * ArgV_o) 0; 0 exp(+1i * ArgH_o)] ;
tp_oI = [exp(+1i * ArgV_oI) 0; 0 exp(+1i * ArgH_oI)] ;

end