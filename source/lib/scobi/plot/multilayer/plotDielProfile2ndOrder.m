function plotDielProfile2ndOrder( fig, eps_diel_z2nd, z )


%% GET GLOBAL PARAMETER
% Simulation Settings
draw_live_plots = SimSettings.getInstance.draw_live_plots;
% Ground Parameters
gnd_layer_depth_m = GndParams.getInstance.layer_depth_m;
% Surface Dynamic Params
eps_g = SurfaceDynParams.getInstance.eps_g;
% Ground-ML Parameters
zA_m = GndMLParams.getInstance.zA_m;


if draw_live_plots
    
    FigON = 1 ;
    if FigON == 1
        figure(fig)
        subplot(3,4,1)
        plot(real(eps_diel_z2nd), (z - zA_m)*1e2, 'k', 'linewidth', 2)
        hold on
        plot(real(eps_g), gnd_layer_depth_m*1e2, 'o')
        grid
        axis([0 30 0 z(end)])
        set(gca,'YDir','reverse')
        set(gca,'YLim',[0-zA_m*1e2 z(end, :) * 1e2-zA_m*1e2])
        % set(gca,'YTick',[0 zz z(end)] * 1e2)
        % aa = sort([zz, z(end), zA_m + gnd_layer_depth_m]) - zA_m ;
        aa = [0; gnd_layer_depth_m] ;
        set(gca,'YTick',aa * 1e2)
        bb = real([Constants.eps_diel_air; eps_g]) ;
        bb = sort(unique(bb)) ;
        set(gca,'XTick',bb)
        xlabel('\epsilon\prime - real part')
        ylabel('z [cm]')
        title('2nd order polynomial')
        hold off
    %     saveas(gcf, strcat(pwd, '\', 'Re2'), 'jpg')
    %     close(gcf)
    end

    FigON = 0 ;
    if FigON == 1
        figure(10)
        subplot(2,2,1)
        plot(real(eps_diel_z2nd), (z-zA_m)*1e2, 'k', 'linewidth', 2)
        hold on
        plot(real(eps_g), gnd_layer_depth_m*1e2, 'o')
        grid
        axis([0 30 0 z(end)])
        set(gca,'YDir','reverse')
        set(gca,'YLim',[0-zA_m*1e2 z(end, :) * 1e2-zA_m*1e2])
        % set(gca,'YTick',[0 zz z(end)] * 1e2)
        % aa = sort([zz, z(end), zA_m + gnd_layer_depth_m]) - zA_m ;
        aa = [0, gnd_layer_depth_m] ;
        set(gca,'YTick',aa * 1e2, 'FontSize', 6)
        bb = real([Constants.eps_diel_air, eps_g]) ;
        bb = sort(unique(bb)) ;
        set(gca,'XTick',bb)
    %     xlabel('\epsilon\prime - real part')
        ylabel('z [cm]')
        title('2nd order polynomial')
        hold off
    %     saveas(gcf, strcat(pwd, '\', 'Re2'), 'jpg')
    %     close(gcf)
    end
    
end


end

