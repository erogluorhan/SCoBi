function plotMovie( Mdata )


%% GET GLOBAL PARAMETER
% Simulation Settings
draw_live_plots = SimSettings.getInstance.draw_live_plots;


if draw_live_plots    

    hf = figure ;
    movie( hf, Mdata );

    % mplay(Mdata)

    %% add Dylan file path (containing movie support functions) and all subdirs
    addpath( genpath( strcat(pwd, '\Dylan') ) );
    % Mdata = delEmptyFrames(Mdata);
    % writeMovie('test_movie.mat', Mdata);
    mplay(Mdata);

end


end

