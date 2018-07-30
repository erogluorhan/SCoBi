%% Mehmet Kurum
% 03/13/2017
% modified - 11/09/2017

% Diffuse (Incoherent) Term
function diffuseTerm_B1_G1(ind_realization)

%% GET GLOBAL DIRECTORIES
dir_config = SimulationFolders.getInstance.config;
dir_freqdiff = SimulationFolders.getInstance.freqdiff;
dir_freqdiff_b1_tuple = SimulationFolders.getInstance.freqdiff_b1_tuple;
dir_freqdiff_b2_tuple = SimulationFolders.getInstance.freqdiff_b2_tuple;
dir_freqdiff_b3_tuple = SimulationFolders.getInstance.freqdiff_b3_tuple;
dir_freqdiff_b4_tuple = SimulationFolders.getInstance.freqdiff_b4_tuple;
dir_freqdiff_P1_tuple = SimulationFolders.getInstance.freqdiff_P1_tuple;
dir_freqdiff_P2_tuple = SimulationFolders.getInstance.freqdiff_P2_tuple;
dir_freqdiff_P3_tuple = SimulationFolders.getInstance.freqdiff_P3_tuple;
dir_freqdiff_P4_tuple = SimulationFolders.getInstance.freqdiff_P4_tuple;


%% GET GLOBAL PARAMETERS
% Simulation Parameters
Nfz = SimParams.getInstance.Nfz;
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
EIRP_dB = TxParams.getInstance.EIRP_dB;
EIRP = convertDecibelToNatural( EIRP_dB );
% Receiver Parameters
G0r_dB = RxParams.getInstance.G0r_dB;
G0r = convertDecibelToNatural( G0r_dB );
% Vegetation Parameters
scat_cal_veg = VegParams.getInstance.scat_cal_veg ;
TYPKND = VegParams.getInstance.TYPKND;


%% READ OR LOAD META-DATA
% Positions relative to ground
% pos_FP_Rx_m: center of footprint, pos_FZ_m: center of fresnel zone
% AllPoints_m = [pos_Tx_m, pos_TxI_m, pos_SP_m, pos_Rx_m, pos_RxI_m, pos_Gnd_m, pos_B_Rx_m, pos_FP_Rx_m, pos_FZ_m] ;
filenamex = 'AllPoints_m' ;
AllPoints_m = readVar(dir_config, filenamex) ;
pos_Tx_m = AllPoints_m(:, 1) ; % - d_layer_m ;        % Transmitter with respect to veg top
pos_SP_m = AllPoints_m(:, 3) ;       % Specular point  (at top of vegetation)
pos_Rx_m = AllPoints_m(:, 4) ; % - d_layer_m ;        % Receiver with respect to veg top


%% CALCULATIONS
% Slant Range - relative to ground
ST = pos_SP_m - pos_Tx_m ;         % Transmitter to Specular
r_st = vectorMagnitude(ST) ;  % slant range

RS = pos_Rx_m - pos_SP_m ;         % Specular to Receiver
r_sr = vectorMagnitude(RS) ;  % slant range

f_Hz = f_MHz * Constants.MHz2Hz ;
lambda_m = Constants.c / f_Hz ;     % Wavelength
k0 = 2 * pi * f_Hz / Constants.c ;    % Wave number


% Factor Ki
K = 1i * sqrt(EIRP) * sqrt(G0r) * lambda_m / (4 * pi) ;

Ki = K * exp(1i * k0 * (r_st + r_sr)) / (r_st * r_sr) ; 


%% Layer Parameters
[Nlayer, Ntype] = size(TYPKND) ;
NkindMax = max(max(TYPKND)) ;


%% INITIALIZE REQUIRED VARIABLES
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
                filename = strcat('R', num2str(ind_realization), '_L', num2str(ii), '_T',...
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

%% SAVE ALL
% Ki
filename3 = strcat('Ki') ;
writeComplexVar(dir_freqdiff, filename3, Ki)

% ====================================================
% b1
filename1 = strcat('b1_inc1_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b1_tuple, filename1, (b1_inc1_t1))
filename1 = strcat('b1_inc2_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b1_tuple, filename1, (b1_inc2_t1))
filename1 = strcat('b1_inc3_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b1_tuple, filename1, (b1_inc3_t1))
filename1 = strcat('b1_inc4_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b1_tuple, filename1, (b1_inc4_t1))

filename1 = strcat('b1_inc1_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b1_tuple, filename1, (b1_inc1_t2))
filename1 = strcat('b1_inc2_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b1_tuple, filename1, (b1_inc2_t2))
filename1 = strcat('b1_inc3_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b1_tuple, filename1, (b1_inc3_t2))
filename1 = strcat('b1_inc4_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b1_tuple, filename1, (b1_inc4_t2))

% b2
filename1 = strcat('b2_inc1_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b2_tuple, filename1, (b2_inc1_t1))
filename1 = strcat('b2_inc2_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b2_tuple, filename1, (b2_inc2_t1))
filename1 = strcat('b2_inc3_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b2_tuple, filename1, (b2_inc3_t1))
filename1 = strcat('b2_inc4_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b2_tuple, filename1, (b2_inc4_t1))

filename1 = strcat('b2_inc1_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b2_tuple, filename1, (b2_inc1_t2))
filename1 = strcat('b2_inc2_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b2_tuple, filename1, (b2_inc2_t2))
filename1 = strcat('b2_inc3_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b2_tuple, filename1, (b2_inc3_t2))
filename1 = strcat('b2_inc4_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b2_tuple, filename1, (b2_inc4_t2))

% b3
filename1 = strcat('b3_inc1_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b3_tuple, filename1, (b3_inc1_t1))
filename1 = strcat('b3_inc2_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b3_tuple, filename1, (b3_inc2_t1))
filename1 = strcat('b3_inc3_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b3_tuple, filename1, (b3_inc3_t1))
filename1 = strcat('b3_inc4_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b3_tuple, filename1, (b3_inc4_t1))

filename1 = strcat('b3_inc1_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b3_tuple, filename1, (b3_inc1_t2))
filename1 = strcat('b3_inc2_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b3_tuple, filename1, (b3_inc2_t2))
filename1 = strcat('b3_inc3_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b3_tuple, filename1, (b3_inc3_t2))
filename1 = strcat('b3_inc4_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b3_tuple, filename1, (b3_inc4_t2))

% b4
filename1 = strcat('b4_inc1_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b4_tuple, filename1, (b4_inc1_t1))
filename1 = strcat('b4_inc2_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b4_tuple, filename1, (b4_inc2_t1))
filename1 = strcat('b4_inc3_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b4_tuple, filename1, (b4_inc3_t1))
filename1 = strcat('b4_inc4_t1', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b4_tuple, filename1, (b4_inc4_t1))

filename1 = strcat('b4_inc1_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b4_tuple, filename1, (b4_inc1_t2))
filename1 = strcat('b4_inc2_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b4_tuple, filename1, (b4_inc2_t2))
filename1 = strcat('b4_inc3_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b4_tuple, filename1, (b4_inc3_t2))
filename1 = strcat('b4_inc4_t2', '_R', num2str(ind_realization)) ;
writeComplexVar(dir_freqdiff_b4_tuple, filename1, (b4_inc4_t2))

% =====================================================
% P1
filename1 = strcat('P1_inc1_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P1_tuple, filename1, (P1_inc1_t1))
filename1 = strcat('P1_inc2_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P1_tuple, filename1, (P1_inc2_t1))
filename1 = strcat('P1_inc3_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P1_tuple, filename1, (P1_inc3_t1))
filename1 = strcat('P1_inc4_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P1_tuple, filename1, (P1_inc4_t1))

filename1 = strcat('P1_inc1_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P1_tuple, filename1, (P1_inc1_t2))
filename1 = strcat('P1_inc2_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P1_tuple, filename1, (P1_inc2_t2))
filename1 = strcat('P1_inc3_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P1_tuple, filename1, (P1_inc3_t2))
filename1 = strcat('P1_inc4_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P1_tuple, filename1, (P1_inc4_t2))

% P2
filename1 = strcat('P2_inc1_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P2_tuple, filename1, (P2_inc1_t1))
filename1 = strcat('P2_inc2_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P2_tuple, filename1, (P2_inc2_t1))
filename1 = strcat('P2_inc3_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P2_tuple, filename1, (P2_inc3_t1))
filename1 = strcat('P2_inc4_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P2_tuple, filename1, (P2_inc4_t1))

filename1 = strcat('P2_inc1_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P2_tuple, filename1, (P2_inc1_t2))
filename1 = strcat('P2_inc2_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P2_tuple, filename1, (P2_inc2_t2))
filename1 = strcat('P2_inc3_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P2_tuple, filename1, (P2_inc3_t2))
filename1 = strcat('P2_inc4_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P2_tuple, filename1, (P2_inc4_t2))

% P3
filename1 = strcat('P3_inc1_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P3_tuple, filename1, (P3_inc1_t1))
filename1 = strcat('P3_inc2_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P3_tuple, filename1, (P3_inc2_t1))
filename1 = strcat('P3_inc3_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P3_tuple, filename1, (P3_inc3_t1))
filename1 = strcat('P3_inc4_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P3_tuple, filename1, (P3_inc4_t1))

filename1 = strcat('P3_inc1_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P3_tuple, filename1, (P3_inc1_t2))
filename1 = strcat('P3_inc2_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P3_tuple, filename1, (P3_inc2_t2))
filename1 = strcat('P3_inc3_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P3_tuple, filename1, (P3_inc3_t2))
filename1 = strcat('P3_inc4_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P3_tuple, filename1, (P3_inc4_t2))

% P4
filename1 = strcat('P4_inc1_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P4_tuple, filename1, (P4_inc1_t1))
filename1 = strcat('P4_inc2_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P4_tuple, filename1, (P4_inc2_t1))
filename1 = strcat('P4_inc3_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P4_tuple, filename1, (P4_inc3_t1))
filename1 = strcat('P4_inc4_t1', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P4_tuple, filename1, (P4_inc4_t1))

filename1 = strcat('P4_inc1_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P4_tuple, filename1, (P4_inc1_t2))
filename1 = strcat('P4_inc2_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P4_tuple, filename1, (P4_inc2_t2))
filename1 = strcat('P4_inc3_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P4_tuple, filename1, (P4_inc3_t2))
filename1 = strcat('P4_inc4_t2', '_R', num2str(ind_realization)) ;
writeVar(dir_freqdiff_P4_tuple, filename1, (P4_inc4_t2))


end

%% Scattering Matrices
function [b_dd_t1, b_rd_t1, b_dr_t1, b_rr_t1,...
    b_dd_t2, b_rd_t2, b_dr_t2, b_rr_t2, ...
    P_dd_t1, P_rd_t1, P_dr_t1, P_rr_t1,...
    P_dd_t2, P_rd_t2, P_dr_t2, P_rr_t2] = ScatMech(layerIndex, filename, r_st, r_sr)


%% GET GLOBAL DIRECTORIES
dir_afsa = SimulationFolders.getInstance.afsa;
dir_gnd = SimulationFolders.getInstance.gnd;
dir_position = SimulationFolders.getInstance.position;
dir_fzones = SimulationFolders.getInstance.fzones;
dir_incidence = SimulationFolders.getInstance.incidence;
dir_scattering = SimulationFolders.getInstance.scattering;
dir_ant_real = SimulationFolders.getInstance.ant_real;
dir_rot_real = SimulationFolders.getInstance.rot_real;
dir_rot_lookup = SimulationFolders.getInstance.rot_lookup;
dir_fscat = SimulationFolders.getInstance.fscat;
dir_distance = SimulationFolders.getInstance.distance ;


%% GET GLOBAL PARAMETERS
%Simulation Parameters
Nfz = SimParams.getInstance.Nfz;
% Transmitter Parameters
f_MHz = TxParams.getInstance.f_MHz;
f_Hz = f_MHz * Constants.MHz2Hz;  % Frequency points (Hz)
g_t = TxParams.getInstance.g_t ; % ideal
e_t1 = TxParams.getInstance.e_t1 ;
e_t2 = TxParams.getInstance.e_t2 ;
% Vegetation Parameters
dim_layers_m = VegParams.getInstance.dim_layers_m;


%% READ META-DATA
% Ground Parameters
disp('Reading ground parameters...')
filenamex = 'G' ;
grnd_par = readComplexVar( dir_gnd, filenamex );
h = real(grnd_par(1, 1)) ;
% TO-DO: Solve for multilayered ground
% epsg = grnd_par(1, 2) + 1i * grnd_par(1, 3) ;
epsg = grnd_par(2, 1);

% Incremental Propagation Constant
disp('Reading Incremental propagation constants...')
filenamex = 'dKz' ;
dKz = readComplexVar( dir_afsa, filenamex) ;
filenamex = 'ANGDEG' ;
ANGDEG = readVar( dir_afsa, filenamex) ;

% Particle Positions
disp('Reading the particle positions...')

pP = readVar( dir_position, filename );

Npart = readVar( dir_fzones, filename ) ;

% Incidence Angles
disp('Reading incidence angles...')
filenamex = strcat('thid_', filename) ;
thid = readVar(dir_incidence, filenamex) ;
filenamex = strcat('thidI_', filename) ;
thidI = readVar(dir_incidence, filenamex) ;

% Scattering Angles
disp('Reading scattering angles...')
filenamex = strcat('thsd_', filename) ;
thsd = readVar(dir_scattering, filenamex) ;
filenamex = strcat('thsdI_', filename) ;
thsdI = readVar(dir_scattering, filenamex) ;
% no need for phsd and phsdI

% Antenna gain values
disp('Reading antenna values...')
filenamex = strcat('gr_', filename) ;
gr = readComplexVar(dir_ant_real, filenamex) ;
filenamex = strcat('grI_', filename) ;
grI = readComplexVar(dir_ant_real, filenamex) ;

% Rotation Matrix values
disp('Reading rotation matrices...')
filenamex = strcat('u_gar_', filename) ;
u_gar2 = readComplexVar(dir_rot_real, filenamex) ;
filenamex = strcat('u_garI_', filename) ;
u_garI2 = readComplexVar(dir_rot_real, filenamex);
load([dir_rot_lookup '\u_ts.mat'], 'u_ts') 
load([dir_rot_lookup '\u_tIs.mat'], 'u_tIs')

U_ts = calc_Muller(u_ts) ;
U_tIs = calc_Muller(u_tIs) ;

% Scattering amplitudes..
disp('Reading scattering amplitudes...')
filenamex = strcat('BISTATIC1_', filename) ;
BISTATIC1 = readComplexVar(dir_fscat, filenamex) ;
filenamex = strcat('BISTATIC2_', filename) ;
BISTATIC2 = readComplexVar(dir_fscat, filenamex) ;
filenamex = strcat('BISTATIC3_', filename) ;
BISTATIC3 = readComplexVar(dir_fscat, filenamex) ;
filenamex = strcat('BISTATIC4_', filename) ;
BISTATIC4 = readComplexVar(dir_fscat, filenamex) ;

%  Particle distances to receiver and image receiver
disp('Reading distances...')
filenamex = strcat('ro_', filename) ;
ro = readVar(dir_distance, filenamex) ;
filenamex = strcat('roI_', filename) ;
roI = readVar(dir_distance, filenamex) ;
filenamex = strcat('ri_', filename) ;
ri = readVar(dir_distance, filenamex) ;
filenamex = strcat('riI_', filename) ;
riI = readVar(dir_distance, filenamex) ;

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


%% CALCULATIONS
% Wave Number
k0 = 2 * pi * f_Hz / Constants.c ;    % Wave number

% Scattering Matrix Calculations
disp(strcat('Number of Scatters:', num2str(Npart)))

%% Transmitter Antenna pattern
G_t = calc_Muller(g_t) ;


%% INITIALIZE REQUIRED VARIABLES
% Transmitter Pol State
E_t1 = [1; 0; 0; 0] ;
E_t2 = [0; 0; 0; 1] ;
Q = [1 0 0 0; 0 0 0 1; 0 1 1 0; 0 -1i -1i 0] ;
E_t1 = Q * E_t1 ;
E_t2 = Q * E_t2 ;
 
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
        
        [tp_i, tp_iI, tp_o, tp_oI] = CalcTransMat(layerIndex, zp, dim_layers_m, dKz_i, dKz_iI, dKz_o, dKz_oI) ;

%         tp_i = [1 0; 0 1] ;
%         tp_iI = [1 0; 0 1] ;
%         tp_o = [1 0; 0 1] ;
%         tp_oI = [1 0; 0 1] ;
        % Calculate Reflection Coefficient
        thiI = degtorad(180 - thidI) ;        % incident angle
        thsI = degtorad(180 - thsdI(pp)) ;    % scattered angle
        [RGH_iI, RGV_iI, RGH_oI, RGV_oI] = reflectionCoeff(thiI, thsI, epsg, h) ;

        % Antenna Pattern Matrix
        % 2 X 2
%         g_ro = squeeze(gr(pp, :, :)) ;
%         g_roI = squeeze(grI(pp, :, :)) ;
        g_ro = [1 0; 0 1] ;
        g_roI = [1 0; 0 1] ;
        % 4 X 4
        G_ro = calc_Muller(g_ro) ;
        G_roI = calc_Muller(g_roI) ;                  

        % Rotation Matrix
        % 2 X 2
        u_gar = squeeze(u_gar2(pp, :, :)) ;
        u_garI = squeeze(u_garI2(pp, :, :)) ;
        % 4 X 4
        U_gar = calc_Muller(u_gar) ;
        U_garI = calc_Muller(u_garI) ;
        
        % Bistatic Scattering Amplitude
        f_oi = squeeze(BISTATIC1(pp, :, :)) ;
        f_oIi = squeeze(BISTATIC2(pp, :, :)) ;
        f_oiI = squeeze(BISTATIC3(pp, :, :)) ;
        f_oIiI = squeeze(BISTATIC4(pp, :, :)) ;
%         f_oi = [1 0; 0 1] ;
%         f_oIi = [1 0; 0 1] ;
%         f_oiI = [1 0; 0 1] ;
%         f_oIiI = [1 0; 0 1] ;

        %% B-factor        
%         B_dd = (exp(+1i * k0 * (ro(pp) + ri(pp))) / (ro(pp) * ri(pp))) ...
%             / (exp(+1i * k0 * (r_st + r_sr)) / (r_st * r_sr)) ;
%         B_rd = (exp(+1i * k0 * (roI(pp) + ri(pp))) / (roI(pp) * ri(pp))) ...
%             / (exp(+1i * k0 * (r_st + r_sr)) / (r_st * r_sr)) ;
%         B_dr = (exp(+1i * k0 * (ro(pp) + riI(pp))) / (ro(pp) * riI(pp))) ...
%             / (exp(+1i * k0 * (r_st + r_sr)) / (r_st * r_sr)) ;
%         B_rr = (exp(+1i * k0 * (roI(pp) + riI(pp))) / (roI(pp) * riI(pp))) ...
%             / (exp(+1i * k0 * (r_st + r_sr)) / (r_st * r_sr)) ;      
        B_dd = 1 ; %  
        B_rd = 1 ; % 
        B_dr = 1 ; % 
        B_rr = 1 ; %
            
        % Ground Reflection Matrices
        rg_iI = [RGV_iI 0; 0 RGH_iI] ;
        rg_oI = [RGV_oI 0; 0 RGH_oI] ;            
% %             rg_iI = [1 0; 0 1] ;
% %             rg_oI = [1 0; 0 1] ;

        % Calculate Scattering Matrices
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

        
        %% SUM OVER INDIVIDUAL CONTRIBUTIONS
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
    
    
    %% SUM OVER FRESNEL ZONES   
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
    
end % Nfz -- fz
toc
       
end


%% Calculate Transmissivity Matrices
function [tp_i, tp_iI, tp_o, tp_oI] = CalcTransMat(layerIndex, zp, D, dKz_i, dKz_iI, dKz_o, dKz_oI)
% zp: particle center of gravity position with absolute z-values starting 
% from zero-ground  


% Convert absoulte z-values to relative ones within the vegetation layer
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