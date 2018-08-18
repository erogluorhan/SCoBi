function resetWS

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

if ( exist('myBreakpoints.mat','file') ) 
    delete('myBreakpoints.mat'); 
end

end