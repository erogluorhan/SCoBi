function plotDielProfileDiscreteSlab( fig, eps_diel_zs, z )


%% GET GLOBAL PARAMETER
% Simulation Settings
draw_live_plots = SimSettings.getInstance.draw_live_plots;
% Ground Parameters
gnd_layer_depth_m = GndParams.getInstance.layer_depth_m;
% Ground-ML Parameters
zA_m = GndMLParams.getInstance.zA_m;
% Ground Dynamic Params
eps_g = GndDynParams.getInstance.eps_g;


if draw_live_plots
    
    FigON = 1 ;
    
    if FigON == 1
        figure(fig);
        subplot(3,4,4)
        plot(real(eps_diel_zs), (z-zA_m)*1e2, 'r', 'linewidth', 2)
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
            title('Slab')
     hold off
     %     saveas(gcf, strcat(pwd, '\', 'Re3'), 'jpg')
    %     close(gcf)
    end

    FigON = 0 ;
    if FigON == 1
        figure(10)
        subplot(2,2,4)
        plot(real(eps_diel_zs), (z-zA_m)*1e2, 'r', 'linewidth', 2)
        hold on
        plot(real(eps_g), gnd_layer_depth_m*1e2, 'o')
        grid
        axis([0 30 0 z(end)])
        set(gca,'YDir','reverse')
        set(gca,'YLim',[0-zA_m*1e2 z(end, :) * 1e2-zA_m*1e2])
        % set(gca,'YTick',[0 zz z(end)] * 1e2)
        % aa = sort([zz, z(end), zA_m + gnd_layer_depth_m]) - zA_m ;
        aa = [0; gnd_layer_depth_m] ;
        set(gca,'YTick',aa * 1e2, 'FontSize', 6)
        bb = real([Constants.eps_diel_air; eps_g]) ;
        bb = sort(unique(bb)) ;
        set(gca,'XTick',bb)
        xlabel('\epsilon\prime - real part')
    %     ylabel('z [cm]')
            title('Slab')
     hold off
         saveas(gcf, strcat(pwd, '\', 'profiles'), 'jpg')
        close(gcf)
    end
    
end




end

