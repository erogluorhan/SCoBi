%% Orhan Eroglu
% 13 July, 2017

%% Function that converts any angle to the [0,360] range
function newAngle = convertAngleTo360Range( angle )

factor = floor( angle / 360 );

newAngle = angle - factor * 360;

end