function [M1, M2] = plotAddSMpoint( fig1, fig2, DOY, VSM_cm3cm3 )


%% GET GLOBAL PARAMETER
% Simulation Settings
draw_live_plots = SimSettings.getInstance.draw_live_plots;
% Dynamic Parameters
DOYs = DynParams.getInstance.DOYs;
VSM_list_cm3cm3 = DynParams.getInstance.VSM_list_cm3cm3;


M1 = [];
M2 = [];

if draw_live_plots

    DoY1 = DOYs(1);
    DoY2 = DOYs(end);
    
    figure(1)
    subplot(3,4,9:12)
    hold on
    plot(DOY, VSM_cm3cm3, 'o')
    axis([DoY1 DoY2 0 0.5])

    M1 = getframe(fig1);
    drawnow
    
    M2 = getframe(fig2);
    drawnow
    
end
    
end
