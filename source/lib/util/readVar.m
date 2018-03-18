function var = readVar(pathname, filename)

folderNotEmpty = any(size(dir([strcat( pathname ) '/*.dat' ] ), 1 ));

if folderNotEmpty
    
    if ~isempty(pathname), filename = strcat( pathname, '\', filename, '.dat' ) ; end

    fid = fopen(filename,'r') ;

    NofDims = str2num(fgetl(fid)) ; %#ok<ST2NM>

    SizeVar = zeros(1, NofDims) ;

    for ii = 1 : NofDims

        SizeVar(ii) = str2num(fgetl(fid)) ; %#ok<ST2NM>

    end
    var = zeros(SizeVar) ;

    N = SizeVar(1) ;

    for n = 1 : N

        var(n, :) = str2num(fgetl(fid)) ; %#ok<ST2NM>

    end

    fclose(fid);
    
else
    
    var = NaN;
    
end


end