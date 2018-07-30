function [ element_id ] = findElementIdInCell( cellArr, element )

is_element = cellfun(@(x)isequal(x, element), cellArr );

[~, element_id] = find(is_element);

end

