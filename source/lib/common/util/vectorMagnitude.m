%% Function that calculates magnitude of the given vector

function abs = vectorMagnitude(vec)

[m,n] = size(vec);
if (m ~= 1)&&(n ~= 1)  % or unit colomn or unit row
    abs = 0;
    disp('Error - vector is not of proper dimensions');
else
    abs = sqrt(sum (vec.^2));
end;

end