

function UU = calc_Muller(uu)

u11 = uu(1, 1) ; u12 = uu(1, 2) ;
u21 = uu(2, 1) ; u22 = uu(2, 2) ;

U11 = abs(u11) .^ 2 ;
U12 = abs(u12) .^ 2 ;
U13 = real(u11 .* conj(u12)) ;
U14 = -imag(u11 .* conj(u12)) ;

U21 = abs(u21) .^ 2 ;
U22 = abs(u22) .^ 2 ;
U23 = real(u21 .* conj(u22)) ;
U24 = -imag(u21 .* conj(u22)) ;

U31 = 2 * real(u11 .* conj(u21)) ; 
U32 = 2 * real(u12 .* conj(u22)) ; 
U33 = real(u11 .* conj(u22) + u12 .* conj(u21)) ;
U34 = -imag(u11 .* conj(u22) - u12 .* conj(u21)) ;

U41 = 2 * imag(u11 .* conj(u21)) ; 
U42 = 2 * imag(u12 .* conj(u22)) ; 
U43 = imag(u11 .* conj(u22) + u12 .* conj(u21)) ;
U44 = real(u11 .* conj(u22) - u12 .* conj(u21)) ;

UU = [U11, U12, U13, U14; ...
    U21, U22, U23, U24; ...
    U31, U32, U33, U34; ...
    U41, U42, U43, U44] ;

end

