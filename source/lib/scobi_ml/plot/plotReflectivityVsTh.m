function plotReflectivityVsTh( Rp1, Rp2, Rp1_L, Rp2_L, Rp1_2, Rp2_2, Rp1_3, Rp2_3 )


%% GET GLOBAL PARAMETER
% Simulation Settings
draw_live_plots = SimSettings.getInstance.draw_live_plots;


if draw_live_plots
    
    % FigON = 1 ;
    % if FigON == 1
    %     
    %     figure
    %     plot(th0_Tx_list_deg, Rp1_L, '-or', th0_Tx_list_deg, Rp2_L, '-sb')
    %     axis([0 90 0 1])
    %     xlabel('Angle of Observation')
    %     ylabel('Reflectivity')
    %     title('Logistic function')
    %     saveas(gcf, strcat(pwd, '\', 'R_smooth_', pol_Tx, pol_Rx), 'jpg')
    %     close(gcf)
    % end

    % FigON = 1 ;
    % if FigON == 1
    %     
    %     figure
    %     plot(th0_Tx_list_deg, Rp1_2, '-or', th0_Tx_list_deg, Rp2_2, '-sb')
    %     axis([0 90 0 1])
    %     xlabel('Angle of Observation')
    %     ylabel('Reflectivity')
    %     title('polynomial 2nd order')
    %     saveas(gcf, strcat(pwd, '\', 'R_2ndorder_', pol_Tx, pol_Rx), 'jpg')
    %     close(gcf)
    % end

    % FigON = 1 ;
    % if FigON == 1
    %     
    %     figure
    %     plot(th0_Tx_list_deg, Rp1_3, '-or', th0_Tx_list_deg, Rp2_3, '-sb')
    %     axis([0 90 0 1])
    %     xlabel('Angle of Observation')
    %     ylabel('Reflectivity')
    %     title('polynomial 3rd order')
    %     saveas(gcf, strcat(pwd, '\', 'R_3rdorder_', pol_Tx, pol_Rx), 'jpg')
    %     close(gcf)
    % end
    
end

end

