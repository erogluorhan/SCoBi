function plotDielProfile3rdOrder( fig, eps_diel_z3, eps_diel_soil, z )


%% GET GLOBAL PARAMETER
% Simulation Settings
draw_live_plots = SimSettings.getInstance.draw_live_plots;
% Ground Parameters
gnd_layer_depth_m = GndParams.getInstance.layer_depth_m;
% Ground-ML Parameters
zA_m = GndMLParams.getInstance.zA_m;


if draw_live_plots
    
    FigON = 1 ;
    if FigON == 1
        figure(fig);
        subplot(3,4,2)
        plot(real(eps_diel_z3), (z-zA_m)*Constants.m2cm, 'c', 'linewidth', 2)
        hold on
        plot(real(eps_diel_soil), gnd_layer_depth_m*Constants.m2cm, 'o')
        grid
        axis([0 30 0 z(end)])
        set(gca,'YDir','reverse')
        set(gca,'YLim',[0-zA_m*Constants.m2cm z(end, :)*Constants.m2cm - zA_m*Constants.m2cm])
        % set(gca,'YTick',[0 zz z(end)] * Constants.m2cm)
        % aa = sort([zz, z(end), zA_m + gnd_layer_depth_m]) - zA_m ;
        aa = [0; gnd_layer_depth_m] ;
        set(gca,'YTick',aa * Constants.m2cm)
        bb = real([Constants.eps_diel_air; eps_diel_soil]) ;
        bb = sort(unique(bb)) ;
        set(gca,'XTick',bb)
        xlabel('\epsilon\prime - real part')
        ylabel('z [cm]')
        title('3rd order polynomial')
        hold off

    %     saveas(gcf, strcat(pwd, '\', 'Re3'), 'jpg')
    %     close(gcf)
    end

    FigON = 0 ;
    if FigON == 1
        figure(10)
        subplot(2,2,2)
        plot(real(eps_diel_z3), (z-zA_m)*Constants.m2cm, 'c', 'linewidth', 2)
        hold on
        plot(real(eps_diel_soil), gnd_layer_depth_m*Constants.m2cm, 'o')
        grid
        axis([0 30 0 z(end)])
        set(gca,'YDir','reverse')
        set(gca,'YLim',[0-zA_m*Constants.m2cm z(end, :) * Constants.m2cm-zA_m*Constants.m2cm])
        % set(gca,'YTick',[0 zz z(end)] * Constants.m2cm)
        % aa = sort([zz, z(end), zA_m + gnd_layer_depth_m]) - zA_m ;
        aa = [0, gnd_layer_depth_m] ;
        set(gca,'YTick',aa * Constants.m2cm, 'FontSize', 6)
        bb = real([Constants.eps_diel_air, eps_diel_soil]) ;
        bb = sort(unique(bb)) ;
        set(gca,'XTick',bb)
    %     xlabel('\epsilon\prime - real part')
    %     ylabel('z [cm]')
        title('3rd order plynomial')
        hold off

    %     saveas(gcf, strcat(pwd, '\', 'Re3'), 'jpg')
    %     close(gcf)
    end
    
end


end

