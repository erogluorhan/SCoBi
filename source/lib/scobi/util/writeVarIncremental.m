
function writeVarIncremental(pathname, filename, index, var)
% var: Column vector    

% First read the existing variable, if any    
currentVar = readVar(pathname, filename);
    
% If no current variable, write var as initial
if isnan( currentVar )

    [N, M] = size(var);

% Else if current variable exists, append var to the end
else
    % First read the current variable's size
    [N, M] = size(currentVar);
    
    if index > M
        % Enlarge the current var w.r.t var's end index
        M = index; 
    end
    
    currentVar(:, index) = var;
    var = currentVar;

end

% Create the folder if not exists, set the full filename
if ~isempty(pathname)

    if exist(pathname, 'dir') ~= 7
        mkdir( pathname )
    end
    
    filename = strcat( pathname, '\', filename, '.dat' ) ;

end

% Open the file to write
fid = fopen(filename, 'w') ;


fprintf(fid, '%6.4f\n', ndims(var)) ;
fprintf(fid, '%6.4f\n', size(var)) ;


for nn = 1 : N

    for mm = 1 : M

        fprintf(fid, '%6.16f  ', var(nn, mm)) ;

    end

    fprintf(fid, '\n') ;

end

fclose(fid) ;
    
end
    