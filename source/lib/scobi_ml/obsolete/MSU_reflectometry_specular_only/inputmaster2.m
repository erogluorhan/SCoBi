function inputmaster2


%% Input values
global fMHz EIRP G0r hr PH0d hpbw SLL XPL polT polR tp
global typeOfCanopy

typeOfCanopy = 'CornACRE2017' ;

% % Nfz = 10 ;                      % Number of Fresnel Zones
fMHz = 370 ;                % Operating frequncy in MHz
G0r = 1 ;                       % Receive Antenna Gain
EIRP = 1 ;                      % Pt * G0t - equivalent isotropic radiated power
hr = 20 ;                       % Antenna Height (m)
PH0d = 0 ;                      % Azimuth angle
hpbw = 30 ;                     % beamwidth
SLL = 15 ;                      % Sidelobe Level
XPL = 25 ;                      % X-pol level
% % Nr = 20 ;                        % Number of Realization
tp = [5 15 25] ;                % Theta Probe

polT = 'R' ;
polR = 'X' ;


end