function plotSpecular3

%% Global
global th0d tp FolderPath_fig_spec polT polR Ks

%% Observation angles
th0ds = 10 : 10 : 70 ;

%% Intilization
% over vegetation
b_coh1veg = [] ; P_coh1veg = [] ;    % with non-ideal antenna pattern
b_coh2veg = [] ; P_coh2veg = [] ;
b0_coh1veg = [] ; P0_coh1veg = [] ;  % ideal
b0_coh2veg = [] ; P0_coh2veg = [] ;
% over bare soil
b_coh1bare = [] ; P_coh1bare = [] ;
b_coh2bare = [] ; P_coh2bare = [] ;
b0_coh1bare = [] ; P0_coh1bare = [] ;
b0_coh2bare = [] ; P0_coh2bare = [] ;
Kss = [] ;  % Constant Ks for each angle of observation

%%
FolderPath_write_b = strcat(FolderPath_fig_spec, '\b') ;
if ~exist(FolderPath_write_b, 'dir')
    mkdir(FolderPath_write_b)
end

FolderPath_write_P = strcat(FolderPath_fig_spec, '\P') ;
if ~exist(FolderPath_write_P, 'dir')
    mkdir(FolderPath_write_P)
end

%% Reading ...
for ii = 1 : length(th0ds) 
    
    th0d = th0ds(ii) 
    inputVals ;
    calc_Ks ;
    
    [b_coh1vegx, b_coh2vegx, b0_coh1vegx, b0_coh2vegx, ...
    b_coh1barex, b_coh2barex, b0_coh1barex, b0_coh2barex, ...
    P_coh1vegx, P_coh2vegx, P0_coh1vegx, P0_coh2vegx, ...
    P_coh1barex, P_coh2barex, P0_coh1barex, P0_coh2barex] = readSpecular3 ;

    b_coh1veg = [b_coh1veg b_coh1vegx] ; %#ok<*AGROW>
    b_coh2veg = [b_coh2veg b_coh2vegx] ;
    b0_coh1veg = [b0_coh1veg b0_coh1vegx] ;
    b0_coh2veg = [b0_coh2veg b0_coh2vegx] ;
    b_coh1bare = [b_coh1bare b_coh1barex] ;
    b_coh2bare = [b_coh2bare b_coh2barex] ;
    b0_coh1bare = [b0_coh1bare b0_coh1barex] ;
    b0_coh2bare = [b0_coh2bare b0_coh2barex] ;
    
    P_coh1veg = [P_coh1veg P_coh1vegx] ; %#ok<*AGROW>
    P_coh2veg = [P_coh2veg P_coh2vegx] ;
    P0_coh1veg = [P0_coh1veg P0_coh1vegx] ;
    P0_coh2veg = [P0_coh2veg P0_coh2vegx] ;
    P_coh1bare = [P_coh1bare P_coh1barex] ;
    P_coh2bare = [P_coh2bare P_coh2barex] ;
    P0_coh1bare = [P0_coh1bare P0_coh1barex] ;
    P0_coh2bare = [P0_coh2bare P0_coh2barex] ;
    
    Kss = [Kss Ks] ;
end

[M, ~] = size(b_coh1veg) ;
Nsm = M / 2 ;
% Nang = N ;

%%

FigON = 1 ;
if FigON == 1
    for ii = 1 : Nsm   % for each soil moisture
        
        figure
        
        if (polT == 'X') && (polT == 'X')
            
            bareVV = b_coh1bare(2 * ii - 1, :) ;
            bareHH = b_coh2bare(2 * ii, :) ;
            vegVV = b_coh1veg(2 * ii - 1, :) ;
            vegHH = b_coh2veg(2 * ii, :) ;
            % bare
            plot(th0ds, abs(bareVV) .^ 2, ':or')     % co-pol VV
            hold
            plot(th0ds, abs(bareHH) .^ 2, '-or')         % co-pol HH
            % veg
            plot(th0ds, abs(vegVV) .^ 2, ':og')      % co-pol VV
            plot(th0ds, abs(vegHH) .^ 2, '-og')          % co-pol HH
        else
            BareCO = b_coh1bare(2 * ii - 1, :) ;
            BareX = b_coh1bare(2 * ii, :) ;
            VegCO = b_coh1veg(2 * ii - 1, :) ;
            VegX = b_coh1veg(2 * ii, :) ;
            % bare
            plot(th0ds, abs(BareCO) .^ 2, ':or')     % co-pol
            hold
            plot(th0ds, abs(BareX) .^ 2, '-or')         % x-pol
            % veg
            plot(th0ds, abs(VegCO) .^ 2, ':og')      % co-pol
            plot(th0ds, abs(VegX) .^ 2, '-og')          % x-pol
        end
        
        axis([0 90 0 0.5])
        set(gca,'XTick', th0ds);
        set(gca,'XTickLabel',th0ds)
        set(gca,'YTick', 0.1 : 0.1 : 0.5);
        set(gca,'YTickLabel', 0.1 : 0.1 : 0.5)
        grid
        xlabel('Incidence Angle')
        ylabel('Reflectivity')
        
        if (polT == 'R') && (polR == 'R')
            legend('bare: RR-pol', 'bare: RL-pol', 'veg: RR-pol', 'veg: RL-pol')
            title(strcat('Transmit: RHCP, Receive: RHCP', '- VSM:', num2str(tp(ii)), '%'))
        elseif (polT == 'R') && (polR == 'X')
            legend('bare: RX-pol', 'bare: RY-pol', 'veg: RX-pol', 'veg: RY-pol')
            title(strcat('Transmit: RHCP, Receive: LINEAR', '- VSM:', num2str(tp(ii)), '%'))
        elseif (polT == 'X') && (polT == 'X')
            legend('bare: XX-pol', 'bare: YY-pol', 'veg: XX-pol', 'veg: YY-pol')
            title(strcat('Transmit: LINEAR, Receive: LINEAR', '- VSM:', num2str(tp(ii)), '%'))
        end
        
        fname = strcat(polT, polR, num2str(tp(ii))) ;
        saveas(gcf, strcat(FolderPath_write_b, '\',fname), 'jpg')
        close
        
    end
end

FigON = 1 ;
if FigON == 1
    % in dB and includes Ks
    for ii = 1 : Nsm
        
        figure
        
        if (polT == 'X') && (polT == 'X')
            
            bareVV = Kss .* b_coh1bare(2 * ii - 1, :) ;
            bareHH = Kss .* b_coh2bare(2 * ii, :) ;
            vegVV = Kss .* b_coh1veg(2 * ii - 1, :) ;
            vegHH = Kss .* b_coh2veg(2 * ii, :) ;
            
            plot(th0ds, 10 * log10(abs(bareVV) .^ 2), ':or') % co-pol VV
            hold
            plot(th0ds, 10 * log10(abs(bareHH) .^ 2), '-or') % co-pol HH
            plot(th0ds, 10 * log10(abs(vegVV) .^ 2), ':og') % co-pol VV
            plot(th0ds, 10 * log10(abs(vegHH) .^ 2), '-og') % co-pol HH
        else
            BareCO = Kss .* b_coh1bare(2 * ii - 1, :) ;
            BareX = Kss .* b_coh1bare(2 * ii, :) ;
            VegCO = Kss .* b_coh1veg(2 * ii - 1, :) ;
            VegX = Kss .* b_coh1veg(2 * ii, :) ;
            
            plot(th0ds, 10 * log10(abs(BareCO) .^ 2), ':or') % co-pol
            hold
            plot(th0ds, 10 * log10(abs(BareX) .^ 2), '-or') % x-pol
            plot(th0ds, 10 * log10(abs(VegCO) .^ 2), ':og') % co-pol
            plot(th0ds, 10 * log10(abs(VegX) .^ 2), '-og') % x-pol
        end
        
        
        grid
        xlabel('Incidence Angle')
        ylabel('Received Power [dB]')
        
        if (polT == 'R') && (polR == 'R')
            legend('bare: RR-pol', 'bare: RL-pol', 'veg: RR-pol', 'veg: RL-pol')
            title(strcat('Transmit: RHCP, Receive: RHCP', '- VSM:', num2str(tp(ii)), '%'))
        elseif (polT == 'R') && (polR == 'X')
            legend('bare: RX-pol', 'bare: RY-pol', 'veg: RX-pol', 'veg: RY-pol')
            title(strcat('Transmit: RHCP, Receive: LINEAR', '- VSM:', num2str(tp(ii)), '%'))
        elseif (polT == 'X') && (polT == 'X')
            legend('bare: XX-pol', 'bare: YY-pol', 'veg: XX-pol', 'veg: YY-pol')
            title(strcat('Transmit: LINEAR, Receive: LINEAR', '- VSM:', num2str(tp(ii)), '%'))
        end
        
        axis([0 80 -280 -210])
        
        fname = strcat(polT, polR, num2str(tp(ii)), '_dB') ;
        saveas(gcf, strcat(FolderPath_write_b, '\',fname), 'jpg')
        close
        
    end
end

FigON = 1 ;
if FigON == 1
    % in dB and no Ks
    for ii = 1 : Nsm
        
        figure
        
        if (polT == 'X') && (polT == 'X')
            
            bareVV = b_coh1bare(2 * ii - 1, :) ;
            bareHH = b_coh2bare(2 * ii, :) ;
            vegVV = b_coh1veg(2 * ii - 1, :) ;
            vegHH = b_coh2veg(2 * ii, :) ;
            
            plot(th0ds, 10 * log10(abs(bareVV) .^ 2), ':or') % co-pol VV
            hold
            plot(th0ds, 10 * log10(abs(bareHH) .^ 2), '-or') % co-pol HH
            plot(th0ds, 10 * log10(abs(vegVV) .^ 2), ':og') % co-pol VV
            plot(th0ds, 10 * log10(abs(vegHH) .^ 2), '-og') % co-pol HH
        else
            BareCO = b_coh1bare(2 * ii - 1, :) ;
            BareX = b_coh1bare(2 * ii, :) ;
            VegCO = b_coh1veg(2 * ii - 1, :) ;
            VegX = b_coh1veg(2 * ii, :) ;
            
            plot(th0ds, 10 * log10(abs(BareCO) .^ 2), ':or') % co-pol
            hold
            plot(th0ds, 10 * log10(abs(BareX) .^ 2), '-or') % x-pol
            plot(th0ds, 10 * log10(abs(VegCO) .^ 2), ':og') % co-pol
            plot(th0ds, 10 * log10(abs(VegX) .^ 2), '-og') % x-pol
        end
        
        
        grid
        xlabel('Incidence Angle')
        ylabel('Received Power [dB]')
        
        if (polT == 'R') && (polR == 'R')
            legend('bare: RR-pol', 'bare: RL-pol', 'veg: RR-pol', 'veg: RL-pol', 'Location', 'southwest')
            title(strcat('No Ks - Transmit: RHCP, Receive: RHCP', '- VSM:', num2str(tp(ii)), '%'))
        elseif (polT == 'R') && (polR == 'X')
            legend('bare: RX-pol', 'bare: RY-pol', 'veg: RX-pol', 'veg: RY-pol', 'Location', 'southwest')
            title(strcat('No Ks - Transmit: RHCP, Receive: LINEAR', '- VSM:', num2str(tp(ii)), '%'))
        elseif (polT == 'X') && (polT == 'X')
            legend('bare: XX-pol', 'bare: YY-pol', 'veg: XX-pol', 'veg: YY-pol', 'Location', 'southwest')
            title(strcat('No Ks - Transmit: LINEAR, Receive: LINEAR', '- VSM:', num2str(tp(ii)), '%'))
        end
        
        axis([0 80 -50 0])
        
        fname = strcat(polT, polR, num2str(tp(ii)), '_dB_noKs') ;
        saveas(gcf, strcat(FolderPath_write_b, '\',fname), 'jpg')
        close
        
    end
end

%%

FigON = 1 ;
if FigON == 1
    % in dB and includes Ks
    for ii = 1 : Nsm
        
        figure
        
        if (polT == 'X') && (polT == 'X')
            
            bareVV = abs(Kss) .^ 2 .* P_coh1bare(4 * ii - 3, :) / (4 * pi ) ;
            bareHH = abs(Kss) .^ 2 .* P_coh2bare(4 * ii - 2, :) / (4 * pi ) ;
            vegVV = abs(Kss) .^ 2 .* P_coh1veg(4 * ii - 3, :) / (4 * pi ) ;
            vegHH = abs(Kss) .^ 2 .* P_coh2veg(4 * ii - 2, :) / (4 * pi ) ;
            
            plot(th0ds, 10 * log10(bareVV), ':or') % co-pol VV
            hold
            plot(th0ds, 10 * log10(bareHH), '-or') % co-pol HH
            plot(th0ds, 10 * log10(vegVV), ':og') % co-pol VV
            plot(th0ds, 10 * log10(vegHH), '-og') % co-pol HH
        else
            BareCO = abs(Kss) .^ 2 .* P_coh1bare(4 * ii - 3, :) / (4 * pi ) ;
            BareX = abs(Kss) .^ 2 .* P_coh1bare(4 * ii - 2, :) / (4 * pi ) ;
            VegCO = abs(Kss) .^ 2 .* P_coh1veg(4 * ii - 3, :) / (4 * pi ) ;
            VegX = abs(Kss) .^ 2 .* P_coh1veg(4 * ii - 2, :) / (4 * pi ) ;
            
            plot(th0ds, 10 * log10(BareCO), ':or') % co-pol
            hold
            plot(th0ds, 10 * log10(BareX), '-or') % x-pol
            plot(th0ds, 10 * log10(VegCO), ':og') % co-pol
            plot(th0ds, 10 * log10(VegX), '-og') % x-pol
        end
        
        
        grid
        xlabel('Incidence Angle')
        ylabel('Received Power [dB]')
        
        if (polT == 'R') && (polR == 'R')
            legend('bare: RR-pol', 'bare: RL-pol', 'veg: RR-pol', 'veg: RL-pol')
            title(strcat('Transmit: RHCP, Receive: RHCP', '- VSM:', num2str(tp(ii)), '%'))
        elseif (polT == 'R') && (polR == 'X')
            legend('bare: RX-pol', 'bare: RY-pol', 'veg: RX-pol', 'veg: RY-pol')
            title(strcat('Transmit: RHCP, Receive: LINEAR', '- VSM:', num2str(tp(ii)), '%'))
        elseif (polT == 'X') && (polT == 'X')
            legend('bare: XX-pol', 'bare: YY-pol', 'veg: XX-pol', 'veg: YY-pol')
            title(strcat('Transmit: LINEAR, Receive: LINEAR', '- VSM:', num2str(tp(ii)), '%'))
        end
        
        axis([0 80 -280 -210])
        
        fname = strcat(polT, polR, num2str(tp(ii)), '_dB') ;
        saveas(gcf, strcat(FolderPath_write_P, '\',fname), 'jpg')
        close
        
    end
end

% in dB and no Ks
for ii = 1 : Nsm
    
    figure
    
    if (polT == 'X') && (polT == 'X')
        
        bareVV = P_coh1bare(4 * ii - 3, :) ; % 1st element
        bareHH = P_coh2bare(4 * ii - 2, :) ; % 2nd element
        vegVV = P_coh1veg(4 * ii - 3, :) ;
        vegHH = P_coh2veg(4 * ii - 2, :) ;
        
        plot(th0ds, 10 * log10(bareVV), ':or') % co-pol VV
        hold
        plot(th0ds, 10 * log10(bareHH), '-or') % co-pol HH
        plot(th0ds, 10 * log10(vegVV), ':og') % co-pol VV
        plot(th0ds, 10 * log10(vegHH), '-og') % co-pol HH
    else
        BareCO = P_coh1bare(4 * ii - 3, :) ; % 1st element
        BareX = P_coh1bare(4 * ii - 2, :) ;  % 2nd element
        VegCO = P_coh1veg(4 * ii - 3, :) ;
        VegX = P_coh1veg(4 * ii - 2, :) ;
        
        plot(th0ds, 10 * log10(BareCO), ':or') % co-pol
        hold
        plot(th0ds, 10 * log10(BareX), '-or') % x-pol
        plot(th0ds, 10 * log10(VegCO), ':og') % co-pol
        plot(th0ds, 10 * log10(VegX), '-og') % x-pol
    end
    
    
    grid
    xlabel('Incidence Angle')
    ylabel('Received Power [dB]')
    
    if (polT == 'R') && (polR == 'R')
        legend('bare: RR-pol', 'bare: RL-pol', 'veg: RR-pol', 'veg: RL-pol', 'Location', 'southwest')
        title(strcat('No Ks - Transmit: RHCP, Receive: RHCP', '- VSM:', num2str(tp(ii)), '%'))
    elseif (polT == 'R') && (polR == 'X')
        legend('bare: RX-pol', 'bare: RY-pol', 'veg: RX-pol', 'veg: RY-pol', 'Location', 'southwest')
        title(strcat('No Ks - Transmit: RHCP, Receive: LINEAR', '- VSM:', num2str(tp(ii)), '%'))
    elseif (polT == 'X') && (polT == 'X')
        legend('bare: XX-pol', 'bare: YY-pol', 'veg: XX-pol', 'veg: YY-pol', 'Location', 'southwest')
        title(strcat('No Ks - Transmit: LINEAR, Receive: LINEAR', '- VSM:', num2str(tp(ii)), '%'))
    end
    
    axis([0 80 -30 20])
    
    fname = strcat(polT, polR, num2str(tp(ii)), '_dB_noKs') ;
    saveas(gcf, strcat(FolderPath_write_P, '\',fname), 'jpg')
    close
    
end


end


function calc_Ks

%% global

global FolderPath_Config fMHz EIRP G0r Ks 

%% points
% pT: Transmitter, pS2: specular point, pR2: receiver, pG2: ground (reference), pBr2:boresight,
% pCr2: center of footprint, pSc2: center of fresnel zone
% AllPoints = [pT, pTI, pS2, pR, pRI, pG2, pBr2, pCr2, pSc2] ;
filenamex = 'AllPoints' ;
AllPoints = read_var(FolderPath_Config, filenamex) ;

pT = AllPoints(:, 1) ;        % Transmitter
pS2 = AllPoints(:, 3) ;       % Specular point
pR = AllPoints(:, 4) ;        % Receiver

%% Free space Loss
ST = pS2 - pT ;         % Satellite to Specular
r_st = vectorMag(ST) ;  % slant range

RS = pR - pS2 ;         % Specular to Receiver
r_sr = vectorMag(RS) ;  % slant range

f0hz = fMHz * 1e6 ;
c = 3e8 ;
lambda = c / f0hz ;     % Wavelength

Ls = (4 * pi * r_st / lambda) ^ 2 ;  % free space loss

%% Factor Ks
Ks = 1i * sqrt(EIRP) * sqrt(G0r) / (r_sr * sqrt(Ls)) ;

end

function inputVals

%% Input values
global hr th0d hpbw SLL XPL polT polR
global typeOfCanopy

inputmaster2;

%% Config veg/grnd
global FolderPath_Config FolderPath_Veg FolderPath_gnd FolderPath_afsa
global FolderPath_hr FolderPath_th0d FolderPath_pol 

% Receiver Height
Folder_hr = strcat('hr-', num2str(hr)) ;
FolderPath_hr = strcat('SIM\', typeOfCanopy, '\', Folder_hr) ;
% Angle of incidence
Folder_th0d = strcat('th0d-', num2str(th0d)) ;
FolderPath_th0d = strcat(FolderPath_hr, '\', Folder_th0d) ;
% Polarization
FolderPath_pol = strcat(FolderPath_th0d, '\', polT, polR) ;
% Average Forward Scattering Amplitude
FolderPath_afsa = strcat(FolderPath_hr, '\', 'AFSA') ;
% Vegetation
FolderPath_Veg = strcat(FolderPath_hr, '\', 'VEG') ;
% Ground
FolderPath_gnd = strcat(FolderPath_hr, '\', 'GND') ;
% Configuration
FolderPath_Config = strcat(FolderPath_th0d, '\', 'CONFIG') ;

%% Antenna
global FolderPath_Ant FolderPath_Ant_lookup  

%% Rotation
global FolderPath_rot FolderPath_rot_lookup 

%% Output
global FolderPath_out FolderPath_out_dir FolderPath_out_spec

%% Figure
global FolderPath_fig FolderPath_fig_dir FolderPath_fig_spec 

%% Antenna
FolderNameAnt = strcat('thsd-', num2str(hpbw),...
    '_SLL-', num2str(SLL), '_XPL-', num2str(XPL)) ;
FolderPath_Ant = strcat(FolderPath_th0d, '\', 'ANT\', FolderNameAnt) ;
FolderPath_Ant_lookup = strcat(FolderPath_Ant, '\', 'LOOK-UP') ;

%% Rotation
FolderPath_rot = strcat(FolderPath_pol, '\', 'ROT') ;
FolderPath_rot_lookup = strcat(FolderPath_rot, '\', 'LOOK-UP') ;

%% Output
FolderPath_out = strcat(FolderPath_pol, '\', 'OUTPUT\', FolderNameAnt) ;
FolderPath_out_dir = strcat(FolderPath_out, '\', 'DIRECT') ;
FolderPath_out_spec = strcat(FolderPath_out, '\', 'SPECULAR') ;

%% Figure
FolderPath_fig = strcat(FolderPath_hr, '\', 'FIGURE\', FolderNameAnt) ;
FolderPath_fig_dir = strcat(FolderPath_fig, '\', 'DIRECT') ;
FolderPath_fig_spec = strcat(FolderPath_fig, '\', 'SPECULAR') ;

end

%% read results
function [b_coh1veg, b_coh2veg, b0_coh1veg, b0_coh2veg, ...
    b_coh1bare, b_coh2bare, b0_coh1bare, b0_coh2bare, ...
    P_coh1veg, P_coh2veg, P0_coh1veg, P0_coh2veg, ...
    P_coh1bare, P_coh2bare, P0_coh1bare, P0_coh2bare] = readSpecular3


global FolderPath_out polT polR

% saved output
pathname = strcat(FolderPath_out, '\SPECULAR\') ;

% 2 X 2
filename1 = strcat('Veg1', polT, polR) ;
filename2 = strcat('Veg2', polT, polR) ;
filename01 = strcat('Veg01', polT, polR) ;
filename02 = strcat('Veg02', polT, polR) ;
b_coh1veg = read_cplxvar(pathname, filename1) ;
b_coh2veg = read_cplxvar(pathname, filename2) ;
b0_coh1veg = read_cplxvar(pathname, filename01) ;
b0_coh2veg = read_cplxvar(pathname, filename02) ;

filename1 = strcat('Bare1', polT, polR) ;
filename2 = strcat('Bare2', polT, polR) ;
filename01 = strcat('Bare01', polT, polR) ;
filename02 = strcat('Bare02', polT, polR) ;
b_coh1bare = read_cplxvar(pathname, filename1) ;
b_coh2bare = read_cplxvar(pathname, filename2) ;
b0_coh1bare = read_cplxvar(pathname, filename01) ;
b0_coh2bare = read_cplxvar(pathname, filename02) ;

% 4 X 4
filename1 = strcat('P_Veg1', polT, polR) ;
filename2 = strcat('P_Veg2', polT, polR) ;
filename01 = strcat('P_Veg01', polT, polR) ;
filename02 = strcat('P_Veg02', polT, polR) ;
P_coh1veg = read_var(pathname, filename1) ;
P_coh2veg = read_var(pathname, filename2) ;
P0_coh1veg = read_var(pathname, filename01) ;
P0_coh2veg = read_var(pathname, filename02) ;

filename1 = strcat('P_Bare1', polT, polR) ;
filename2 = strcat('P_Bare2', polT, polR) ;
filename01 = strcat('P_Bare01', polT, polR) ;
filename02 = strcat('P_Bare02', polT, polR) ;
P_coh1bare = read_var(pathname, filename1) ;
P_coh2bare = read_var(pathname, filename2) ;
P0_coh1bare = read_var(pathname, filename01) ;
P0_coh2bare = read_var(pathname, filename02) ;

end

%% Function that calculates magnitude of the given vector

function abs = vectorMag(vec)

[m,n]=size(vec);
if (m~=1)&&(n~=1)  % or unit colomn or unit row
    abs = 0;
    disp('Error - vector is not of proper dimensions');
else
    abs = sqrt(sum (vec.^2));
end;

end