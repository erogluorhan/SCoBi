function initSCoBiWS( simulator_id )

% Store current debug breakpoints before doing clear all
myBreakpoints = dbstatus;
save('myBreakpoints.mat', 'myBreakpoints');

% TO-DO: Decooment it out when simulator_id is read from GUI
% % Clear all the workspace
% clear all;
% clc ;

% Restore debug breakpoints
load('myBreakpoints.mat');
dbstop(myBreakpoints);
clear myBreakpoints;
if (exist('myBreakpoints.mat','file')) delete('myBreakpoints.mat'); end

% Add all subdirectories to the path
addpath( genpath( strcat(pwd, '/common') ) ); % This is required to first addpath
addpath( genpath( Directories.getInstance.input ) );

% If SCoBi-Veg
if simulator_id == Constants.id_veg

    initSCoBiVegWS();
    
% Else if SCoBi-ML
elseif simulator_id == Constants.id_multi_layer

    initSCoBiMLWS();
    
end

end


% Initialize SCoBi-Veg workspace
function initSCoBiVegWS

addpath( genpath( Directories.getInstance.scobi_veg ) );

% addpath( genpath( Directories.getInstance.scobi_veg_analysis ) );
% addpath( genpath( Directories.getInstance.scobi_veg_bistatic ) );
% addpath( genpath( Directories.getInstance.scobi_veg_gui ) );
% addpath( genpath( Directories.getInstance.scobi_veg_init ) );
% addpath( genpath( Directories.getInstance.scobi_veg_monte_carlo ) );
% addpath( genpath( Directories.getInstance.scobi_veg_param ) );
% addpath( genpath( Directories.getInstance.scobi_veg_products ) );
% addpath( genpath( Directories.getInstance.scobi_veg_SCoBi ) );

end


% Initialize SCoBi-ML workspace
function initSCoBiMLWS

addpath( genpath( Directories.getInstance.scobi_ml ) );

end