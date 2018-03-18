
function [result] = convertDecibelToNatural( valDB )

result = power( 10, valDB / 10 );

end

