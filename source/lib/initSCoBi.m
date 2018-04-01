% Store current debug breakpoints before doing clear all
myBreakpoints = dbstatus;
save('myBreakpoints.mat', 'myBreakpoints');

% Clear all the workspace
clear all;
clc ;

% Restore debug breakpoints
load('myBreakpoints.mat');
dbstop(myBreakpoints);
clear myBreakpoints;
if (exist('myBreakpoints.mat','file')) delete('myBreakpoints.mat'); end

% Add all subdirectories to the path
addpath( genpath( strcat(pwd, '/constants') ) ); % This is required to first addpath
addpath( genpath( Directories.getInstance.analysis ) );
addpath( genpath( Directories.getInstance.bistatic ) );
addpath( genpath( Directories.getInstance.gui ) );
addpath( genpath( Directories.getInstance.init ) );
addpath( genpath( Directories.getInstance.input ) );
addpath( genpath( Directories.getInstance.monte_carlo ) );
addpath( genpath( Directories.getInstance.param ) );
addpath( genpath( Directories.getInstance.products ) );
addpath( genpath( Directories.getInstance.SCoBi ) );
addpath( genpath( Directories.getInstance.util ) );