function plotReflectivityForProfiles( DoY, Rp1, Rp2, Rp1_L, Rp2_L, Rp1_2, Rp2_2, Rp1_3, Rp2_3 )


%% GET GLOBAL PARAMETERS
% Simulation Settings
draw_live_plots = SimSettings.getInstance.draw_live_plots;
% Transmitter Parameters
pol_Tx = TxParams.getInstance.pol_Tx;
pol_Rx = RxParams.getInstance.pol_Rx;
% Dielectric Parameters
fig1 = DielParams.getInstance.fig1;
fig2 = DielParams.getInstance.fig2;


if draw_live_plots
    
    figure(fig1)
    subplot(3,4,5:8)
    plot(DoY, Rp1, '-or', 'MarkerFaceColor', 'r', 'markers', 3)
    plot(DoY, Rp2, '-or', 'MarkerFaceColor', 'w', 'markers', 3)
    plot(DoY, Rp1_L, '-ob', 'MarkerFaceColor', 'b', 'markers', 3)
    plot(DoY, Rp2_L, '-ob', 'MarkerFaceColor', 'w', 'markers', 3)
    plot(DoY, Rp1_2, '-ok', 'MarkerFaceColor', 'k', 'markers', 3)
    plot(DoY, Rp2_2, '-ok', 'MarkerFaceColor', 'w', 'markers', 3)
    plot(DoY, Rp1_3, '-oc', 'MarkerFaceColor', 'c', 'markers', 3)
    plot(DoY, Rp2_3, '-oc', 'MarkerFaceColor', 'w', 'markers', 3)
    ylabel('Reflectivity')
    title(strcat('Transmit/Receive Polarization:',  pol_Tx, pol_Rx))

    figure(fig2)
    subplot(2,1,1)
    plot(DoY, Rp1, '-or', 'MarkerFaceColor', 'r', 'markers', 1)
    plot(DoY, Rp2, '-or', 'MarkerFaceColor', 'w', 'markers', 1)
    plot(DoY, Rp1_L, '-ob', 'MarkerFaceColor', 'b', 'markers', 1)
    plot(DoY, Rp2_L, '-ob', 'MarkerFaceColor', 'w', 'markers', 1)
    plot(DoY, Rp1_2, '-ok', 'MarkerFaceColor', 'k', 'markers', 1)
    plot(DoY, Rp2_2, '-ok', 'MarkerFaceColor', 'w', 'markers', 1)
    plot(DoY, Rp1_3, '-oc', 'MarkerFaceColor', 'c', 'markers', 1)
    plot(DoY, Rp2_3, '-oc', 'MarkerFaceColor', 'w', 'markers', 1)
    ylabel('Reflectivity')
    % title(strcat('Transmit/Receive Polarization:',  pol_Tx, pol_Rx))

end

end

