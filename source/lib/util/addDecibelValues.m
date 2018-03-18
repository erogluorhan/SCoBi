
function [result] = addDecibelValues( val1, val2 )

pow1 = power( 10, val1 / 10 );
pow2 = power( 10, val2 / 10 );

result = 10 * log10( pow1 + pow2 );

end

