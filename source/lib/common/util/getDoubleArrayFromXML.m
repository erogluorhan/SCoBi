function output = getDoubleArrayFromXML( xmlFile, varName )

list = xmlFile.getElementsByTagName(varName).item(0);
rows = str2double( list.getElementsByTagName('rows').item(0).getFirstChild.getData );
cols = str2double( list.getElementsByTagName('cols').item(0).getFirstChild.getData );
values = list.getElementsByTagName('val');

output = zeros(rows, cols);

if values.getLength == rows * cols

    for ii = 1 : rows
        for jj = 1 : cols
            ind = (ii-1) * cols + jj - 1;
            output(ii, jj) = str2double( values.item(ind).getFirstChild.getData );
        end
    end
   
else
    % TO-DO: display error!!!
end

end