
function writeAttenuation
% function writeAttenuation 
%
%   Writes the alculated vegetation propagation values into an Excel file, 
%   if this option is selected in the simulation preferences.  
%
%   See also mainSCoBi, calcPropagation.

%   Copyright © 2017-2018 Mehmet Kurum, Orhan Eroglu, Dylan R. Boyd

%   This program is free software: You can redistribute it and/or 
%   modify it under the terms of the GNU General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.

%   Version: 1.0.0



%% GET GLOBAL DIRECTORIES
dir_afsa = SimulationFolders.getInstance.afsa;


%% GET GLOBAL PARAMETERS
% Vegetation Parameters
LTK = VegParams.getInstance.LTK;
TYPKND = VegParams.getInstance.TYPKND;


%% READ META-DATA
disp('reading...')

Desired_Medium = 1 ;
Desired_Layer = 1 ; 
Desired_Type = 1 ;
Desired_Kind = 1 ;

% Layer Combined
if Desired_Medium == 1

    filename = 'ATTENH' ;
    ATTENH = readVar(dir_afsa, filename) ;
    filename = 'ATTENV' ;
    ATTENV = readVar(dir_afsa, filename) ;
    
end

% Type Combined
if Desired_Layer == 1
    
    filename = 'ATTENPIQS' ;
    ATTENPIQS = readVar(dir_afsa, filename) ;
    
end

% Kind Combined
if Desired_Type == 1
    
    filename = 'ATTPIQS' ;
    ATTPIQS = readVar(dir_afsa, filename) ;
    
end

% Scatterers Combined
if  Desired_Kind == 1
    
    filename = 'ATPIQS' ;
    ATPIQS = readVar(dir_afsa, filename) ;
    
end


% Save Sig0
outFileName = ConstantNames.ATTENUATION_OUT_FILENAME ;
outFileName = strcat(dir_afsa, outFileName) ;

% polarization
pols = {'HH'; 'VV'} ;

% Header
Mech = {'One-Way'} ;
TypeNames = {'Leaf', 'Branch', 'Trunk', 'Needle'} ;

% Layer parameters
[Nlayer, Ntype] = size(TYPKND);

for ii = 1 : Nlayer
    
    ofsetnum1 = 63 ; 
    ofsetnum2 = 63 ;
    Lyrnum = {strcat('Layer_ ', num2str(ii))} ;
    
    for jj = 1 : Ntype
        
        Nkind = TYPKND(ii, jj) ;
        
        if Nkind ~= 0
            
            ofsetnum2 = ofsetnum2 + 5 ;
            
            for kk = 1 : Nkind
                
                ofsetnum1 = ofsetnum1 + 5 ;
                
                
                %% SCATTERERS COMBINED
                if  Desired_Kind == 1
                    
                    disp(strcat('L', num2str(ii), 'T', num2str(jj), 'K', num2str(kk)))                    
                    
                    sheetName = 'Kind' ;
                    
                    cellnum = 4 + 6 * (ii - 1) ;
                    
                    % Write Headers
                    cellname = strcat(char(ofsetnum1 - 2), num2str(cellnum + 1)) ;
                    xlswrite(outFileName, strcat(Lyrnum), sheetName, cellname)
                    
                    cellname = strcat(char(ofsetnum1 - 1), num2str(cellnum - 1)) ;
                    xlswrite(outFileName, strcat(LTK(kk, jj, ii)),...
                        sheetName, cellname)
                    
                    cellname = strcat(char(ofsetnum1 - 1), num2str(cellnum)) ;
                    xlswrite(outFileName, strcat(pols), sheetName, cellname)
                    
                    cellname = strcat(char(ofsetnum1), num2str(cellnum - 1)) ;
                    xlswrite(outFileName, strcat(Mech), sheetName, cellname)
                    
                    % Write Data
                    cellname = strcat(char(ofsetnum1), num2str(cellnum)) ;
                    a = squeeze(ATPIQS(:, 35, kk, jj, ii)) ;
                    xlswrite(outFileName, num2cell(a), sheetName, cellname)
                  
                end % Desired Kind check
                
            end % End of Kinds
            
            
            %% KIND COMBINED
            if  Desired_Type == 1
                
                disp(strcat('L', num2str(ii), 'T', num2str(jj)))                
                
                sheetName = 'Type' ;
                
                cellnum = 4 + 6 * (ii - 1) ;
                
                %% Write Headers
                % Layers
                cellname = strcat(char(ofsetnum2 - 2), num2str(cellnum + 1)) ;
                xlswrite(outFileName, strcat(Lyrnum), sheetName, cellname)
                % Type name
                cellname = strcat(char(ofsetnum2 - 1), num2str(cellnum - 1)) ;
                xlswrite(outFileName, strcat(TypeNames(jj)),...
                    sheetName, cellname)
                % Polarization
                cellname = strcat(char(ofsetnum2 - 1), num2str(cellnum)) ;
                xlswrite(outFileName, strcat(pols), sheetName, cellname)
                % Mechanisms
                cellname = strcat(char(ofsetnum2), num2str(cellnum - 1)) ;
                xlswrite(outFileName, strcat(Mech), sheetName, cellname)
                
                %% Write Data
                cellname = strcat(char(ofsetnum2), num2str(cellnum)) ;
                a = squeeze(ATTPIQS(:, 35, jj, ii)) ;
                xlswrite(outFileName, num2cell(a), sheetName, cellname)
                
            end  % Desired Type check
            
        end % End of Kinds
        
    end % End of Types
    
    
    %% TYPE COMBINED
    if  Desired_Layer == 1
     
        disp(strcat('L', num2str(ii)))
        
        sheetName = strcat('Layer') ;
        ofsetnum = 68 ;
        cellnum = 4 + 6 * (ii - 1) ;
        
        cellname = strcat(char(ofsetnum), num2str(cellnum)) ;
        a = squeeze(ATTENPIQS(:, 35, ii)) ;
        xlswrite(outFileName, num2cell(a), sheetName, cellname)
                    
        % Layers
        cellname = strcat(char(ofsetnum - 2), num2str(cellnum + 1)) ;
        xlswrite(outFileName, strcat(Lyrnum), sheetName, cellname)
        % Mechanisms
        cellname = strcat(char(ofsetnum), num2str(cellnum - 1)) ;
        xlswrite(outFileName, strcat(Mech(1)), sheetName, cellname)
        % Polarization
        cellname = strcat(char(ofsetnum - 1), num2str(cellnum)) ;
        xlswrite(outFileName, strcat(pols), sheetName, cellname)
        
    end
    
end


%% LAYER COMBINED
if  Desired_Medium == 1

    sheetName = strcat('Layer') ;
    cellnum = 4 + 6 * Nlayer ;
    ofsetnum = 68 ;
    cellname = strcat(char(ofsetnum), num2str(cellnum)) ;
    a = [ATTENH(45, 1); ATTENV(45, 1)] ;
    xlswrite(outFileName, num2cell(a), sheetName, cellname)
    
    % Layers
    cellname = strcat(char(ofsetnum - 2), num2str(cellnum + 1)) ;
    xlswrite(outFileName, strcat({'Medium'}), sheetName, cellname)
    % Mechanisms
    cellname = strcat(char(ofsetnum), num2str(cellnum - 1)) ;
    xlswrite(outFileName, strcat(Mech(1)), sheetName, cellname)
    % Polarization
    cellname = strcat(char(ofsetnum - 1), num2str(cellnum)) ;
    xlswrite(outFileName, strcat(pols), sheetName, cellname)
            
end


end