function plotSMdata(fig1, fig2)


%% GET GLOBAL PARAMETER
% Ground Parameters
layer_depth_m = GndMLParams.getInstance.layer_depth_m;
layer_depth_cm = layer_depth_m * Constants.M_TO_CM;
% Configuration Parameters
DoYs = ConfigParams.getInstance.DoYs;
VSM_list_cm3cm3 = ConfigParams.getInstance.VSM_list_cm3cm3;

    
DoY1 = DoYs(1);
DoY2 = DoYs(end);

[~, num_gnd_layers] = size(VSM_list_cm3cm3);

%%
% subplot(3,4,1)
% subplot(3,4,2)
% subplot(3,4,3)
% subplot(3,4,4)
% subplot(3,4,5:8)
% subplot(3,4,9:12)


%% PLOT SM DATA
figure(fig1) ;
subplot(3,4,9:12)
FigON = 1 ;
% TO-DO: Add enough number of colors
colors = {'.r', '.b', '.c', '.k'};
depths = [];

for ii = 1 : num_gnd_layers
    depths{1, ii} = strcat( num2str(layer_depth_cm(ii)), 'cm');
end

if FigON == 1

    figure(fig1)
    hold
    for ii = 1 : num_gnd_layers
        plot( DoYs, VSM_list_cm3cm3(:,ii), colors{1,ii}, 'MarkerSize', 6, 'MarkerFaceColor', 'blue', 'linewidth', 0.5);
    end

    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'YLim',[0 0.50], 'FontSize', 14)
    set(gca,'YTick', 0 : 0.10 : 0.50)
    set(gca,'XLim',[DoY1 - 1 DoY2 + 1], 'FontSize', 14)

    grid   
    legend( depths, 'location', 'SouthWest', ...
        'Orientation', 'horizontal', 'FontSize', 8, 'AutoUpdate', 'off');
    xlabel('Day of Year (DoY)', 'FontSize', 14)
    ylabel('VSM [cm^3/cm^3]', 'FontSize', 14)

    title_txt = strcat('ACRE 2017 - DoY from', {' '}, ...
        num2str(DoY1), {' '}, 'to', {' '}, num2str(DoY2)) ;
    title(title_txt)

end

figure(fig2) ;
subplot(2,1,2)
FigON = 1 ;
if FigON == 1

    hold
    for ii = 1 : num_gnd_layers
        plot( DoYs, VSM_list_cm3cm3(:,ii), '.k', 'MarkerSize', 3, 'MarkerFaceColor', 'black', 'linewidth', 0.25)
    end

    set(gca,'YLim',[0 0.50], 'FontSize', 12)
    set(gca,'YTick', 0 : 0.10 : 0.50)
    set(gca,'XLim',[DoY1 - 1 DoY2 + 1], 'FontSize', 12)

    grid    
    xlabel('Day of Year (DoY)', 'FontSize', 12)
    ylabel('VSM [cm^3/cm^3]', 'FontSize', 12)

end

%%
figure(fig1)
subplot(3,4,5:8)
hold on
% axis([DoY1 DoY2 0.45 0.55])
axis([DoY1 DoY2 0.0 0.2])
% axis([DoY1 DoY2 0.4 0.6])
grid

figure(fig2)
subplot(2,1,1)
hold on
axis([DoY1 DoY2 0.0 0.2])
set(gca,'XLim',[DoY1 - 1 DoY2 + 1], 'FontSize', 12)

grid


end

