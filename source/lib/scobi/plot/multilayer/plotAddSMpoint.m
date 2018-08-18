function [M1, M2] = plotAddSMpoint


%% GET GLOBAL PARAMETER
sim_counter = ParamsManager.sim_counter;
% Simulation Settings
draw_live_plots = SimSettings.getInstance.draw_live_plots;
% Configuration Parameters
DoYs = ConfigParams.getInstance.DoYs;
VSM_list_cm3cm3 = ConfigParams.getInstance.VSM_list_cm3cm3;
VSM_cm3cm3 = VSM_list_cm3cm3(sim_counter,:)';
% Dielectric Parameters
fig1 = DielMLDynParams.getInstance.fig1;
fig2 = DielMLDynParams.getInstance.fig2;


M1 = [];
M2 = [];

if draw_live_plots
    
    DoY = DoYs( sim_counter );

    DoY1 = DoYs(1);
    DoY2 = DoYs(end);
    
    figure(1)
    subplot(3,4,9:12)
    hold on
    plot(DoY, VSM_cm3cm3, 'o')
    axis([DoY1 DoY2 0 0.5])

    M1 = getframe(fig1);
    drawnow
    
    M2 = getframe(fig2);
    drawnow
    
end
    
end
