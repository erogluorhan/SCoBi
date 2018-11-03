function plotMovie( Mdata )
 

hf = figure ;
movie( hf, Mdata );

% mplay(Mdata)

%% add Dylan file path (containing movie support functions) and all subdirs
addpath( genpath( strcat(pwd, '\Dylan') ) );
% Mdata = delEmptyFrames(Mdata);
% writeMovie('test_movie.mat', Mdata);
mplay(Mdata);


end

