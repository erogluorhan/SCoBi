function result = isPointInEllipse( position, semiMajorAxis, semiMinorAxis, ellipseCenter, ellipseAngle) 

    x = position(1);
    y = position(2);
    h = ellipseCenter(1);
    k = ellipseCenter(2);

    result = ( cos(ellipseAngle)*(x - h) + sin(ellipseAngle)*(y - k) )^2 / semiMajorAxis^2 + ( sin(ellipseAngle)*(x - h) - cos(ellipseAngle)*(y - k) )^2 / semiMinorAxis^2 <= 1.0;

end

