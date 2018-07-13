% read data from csv file
function [doy, THD0, SM5, SM10, SM20, SM40] = read_CSV_MSU
file_dir = strcat(pwd, '\Dylan\', 'output', '\') ;
file = 'L1b_summary.dat';
fn = strcat(file_dir,'\',file);

% open data
x = importdata(fn);
data = x.data;

% separate data
% ----------------------------------------------------------------------- %
% I currently have the data separated by the following columns. This 
% information can be viewed in the "x" variable above, but it hurts 
% the eyes:
%
% doy = col 1, delay = col 2, reflectivity = col 3, THD0 = col 4
% SM5 = col 5, SM10  = col 6,         SM20 = col 7, SM40 = col 8
% ----------------------------------------------------------------------- %
doy    =  data(:,1);
THD0   =  data(:,4);
SM5    =  data(:,5);
SM10   =  data(:, 6);
SM20   =  data(:,7);
SM40   =  data(:,8);
end