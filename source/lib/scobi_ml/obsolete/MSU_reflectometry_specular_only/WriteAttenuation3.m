function WriteAttenuation3

%% Global
global FolderPath_afsa FolderPath_Veg FolderPath_hr
%%
Desired_Medium = 1 ;
Desired_Layer = 1 ; 
Desired_Type = 1 ;
Desired_Kind = 1 ;

% Save Sig0
outFileName = '\AFSA\Attenuation.xls' ;

%% polarization
pols = {'HH'; 'VV'} ;

%% Header
Mech = {'One-Way'} ;
TypeNames = {'Leaf', 'Branch', 'Trunk', 'Needle'} ;

%% Save Sig0
outFileName = strcat(FolderPath_hr, outFileName) ;

%% Reading vegetation parameters...

filename = 'D' ;
D = read_var(FolderPath_Veg, filename) ;

filename = 'TYPKND' ;
TYPKND = read_var(FolderPath_Veg, filename) ;

load LTK LTK

%% Layer parameters
D  %#ok<NOPRT>
TYPKND %#ok<NOPRT>
[Nlayer, Ntype] = size(TYPKND) ;
% sTYPKND = sum(TYPKND) ;
% Ntype = length(sTYPKND(sTYPKND ~= 0)) ; % L, B, T

%% Reading

disp('reading...')
tic ;

% Layer Combined
if Desired_Medium == 1

    filename = 'ATTENH' ;
    ATTENH = read_var(FolderPath_afsa, filename) ;
    filename = 'ATTENV' ;
    ATTENV = read_var(FolderPath_afsa, filename) ;
    
end

% Type Combined
if Desired_Layer == 1
    
    filename = 'ATTENPIQS' ;
    ATTENPIQS = read_var(FolderPath_afsa, filename) ;
    
end

% Kind Combined
if Desired_Type == 1
    
    filename = 'ATTPIQS' ;
    ATTPIQS = read_var(FolderPath_afsa, filename) ;
    
end

% Scatterers Combined
if  Desired_Kind == 1
    
    filename = 'ATPIQS' ;
    ATPIQS = read_var(FolderPath_afsa, filename) ;
    
end


%%

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
                
                % ++++++++++++++++++++++++++++++
                % Scatterers Combined
                if  Desired_Kind == 1
                    tic ;
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
                    
                                        
                   toc ;
                end
                % +++++++++++++++++++++++++++++
            end
            
            % ++++++++++++++++++++++++++++++
            % Kind Combined
            if  Desired_Type == 1
                tic ;
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
                
                toc ;
            end
            % +++++++++++++++++++++++++++++
            
        end
    end
    
    
    % ++++++++++++++++++++++++++++++
    % Type Combined
    if  Desired_Layer == 1
        tic ;
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
        
        
        toc ;
    end
    % +++++++++++++++++++++++++++++
    
end


% ++++++++++++++++++++++++++++++
% Layer Combined
if  Desired_Medium == 1
    tic ;
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
        
    toc ;
end
% +++++++++++++++++++++++++++++


end