

function Main_MSU_ML_process

%% GET GLOBAL DIRECTORIES
dirFig = SimulationFolders.getInstance.fig;

%% GET GLOBAL PARAMETERS
doy_SM = DynParams.getInstance.doy;
% TO-DO: Generalize the below hard code numbers
VSM5 = DynParams.getInstance.VSM_list_cm3cm3(:,1);
VSM10 = DynParams.getInstance.VSM_list_cm3cm3(:,2);
VSM20 = DynParams.getInstance.VSM_list_cm3cm3(:,3);
VSM40 = DynParams.getInstance.VSM_list_cm3cm3(:,4);

% %% Read Ground Data (SM at 5, 10, 20, and 40 cm)
% [doy_SM, VSM5, VSM10, VSM20, VSM40] = Read_Ground_Data ;

%% Custom date Interval
% d1 = '05-01-2017' ; d2 = '08-31-2017' ;
d1 = '05-25-2017' ; d2 = '06-07-2017' ;
% d1 = '05-25-2017' ; d2 = '05-27-2017' ;
startDate = datenum(d1) ;
endDate = datenum(d2) ;
DoY1 = date2doy(startDate) ;
DoY2 = date2doy(endDate) ;
ci_doy_SM = doy_SM(DoY1<doy_SM & doy_SM<DoY2) ;
ci_SM5 = VSM5(DoY1<doy_SM & doy_SM<DoY2) ;
ci_SM10 = VSM10(DoY1<doy_SM & doy_SM<DoY2) ;
ci_SM20 = VSM20(DoY1<doy_SM & doy_SM<DoY2) ;
ci_SM40 = VSM40(DoY1<doy_SM & doy_SM<DoY2) ;

%%
% subplot(3,4,1)
% subplot(3,4,2)
% subplot(3,4,3)
% subplot(3,4,4)
% subplot(3,4,5:8)
% subplot(3,4,9:12)

%% plot SM data
fig = figure(1) ;
subplot(3,4,9:12)
FigON = 1 ;
if FigON == 1
%     figname = strcat('\SM\SM-', num2str(DoY1), '-', num2str(DoY2)) ; 
    figure(1)
    plot(ci_doy_SM, ci_SM5, '.r', 'MarkerSize', 6, 'MarkerFaceColor', 'blue', 'linewidth', 0.5)
    hold
    plot(ci_doy_SM, ci_SM10, '.b', 'MarkerSize', 6, 'MarkerFaceColor', 'red', 'linewidth', 0.5)
    plot(ci_doy_SM, ci_SM20, '.c', 'MarkerSize', 6, 'MarkerFaceColor', 'red', 'linewidth', 0.5)
    plot(ci_doy_SM, ci_SM40, '.k', 'MarkerSize', 6, 'MarkerFaceColor', 'red', 'linewidth', 0.5)
    
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'YLim',[0 0.50], 'FontSize', 14)
    set(gca,'YTick', 0 : 0.10 : 0.50)
    set(gca,'XLim',[DoY1 - 1 DoY2 + 1], 'FontSize', 14)
%     set(gca,'XTick', [DoY1 DoY2])

    grid    
    legend({'5cm', '10cm','20cm', '40cm'}, 'location', 'SouthWest', ...
        'Orientation', 'horizontal', 'FontSize', 8, 'AutoUpdate', 'off');
    xlabel('Day of Year (DoY)', 'FontSize', 14)
    ylabel('VSM [cm^3/cm^3]', 'FontSize', 14)
    
    title_txt = strcat('ACRE 2017 - DoY from', {' '}, ...
        num2str(DoY1), {' '}, 'to', {' '}, num2str(DoY2)) ;
    title(title_txt)
    
%     saveas(gcf, strcat(dirFig, '\', figname), 'pdf')
%     close
end

fig2 = figure(2) ;
subplot(2,1,2)
FigON = 1 ;
if FigON == 1
%     figname = strcat('\SM\SM-', num2str(DoY1), '-', num2str(DoY2)) ; 
%     figure(1)
    plot(ci_doy_SM, ci_SM5, '.k', 'MarkerSize', 3, 'MarkerFaceColor', 'black', 'linewidth', 0.25)
    hold
    plot(ci_doy_SM, ci_SM10, '.k', 'MarkerSize', 3, 'MarkerFaceColor', 'black', 'linewidth', 0.25)
    plot(ci_doy_SM, ci_SM20, '.k', 'MarkerSize', 3, 'MarkerFaceColor', 'black', 'linewidth', 0.25)
    plot(ci_doy_SM, ci_SM40, '.k', 'MarkerSize', 3, 'MarkerFaceColor', 'black', 'linewidth', 0.25)
    
    set(gca,'YLim',[0 0.50], 'FontSize', 12)
    set(gca,'YTick', 0 : 0.10 : 0.50)
    set(gca,'XLim',[DoY1 - 1 DoY2 + 1], 'FontSize', 12)
%     set(gca,'XTick', [DoY1 DoY2])

    grid    
%     legend('5cm', '10cm','20cm', '40cm', 'location', 'SouthWest')
    xlabel('Day of Year (DoY)', 'FontSize', 12)
    ylabel('VSM [cm^3/cm^3]', 'FontSize', 12)
    
%     title_txt = strcat('ACRE 2017 - DoY from', num2str(DoY1), '-to-', num2str(DoY2)) ;
%     title(title_txt)
    
%     saveas(gcf, strcat(dirFig, '\', figname), 'pdf')
%     close
end

%%
figure(1)
subplot(3,4,5:8)
hold on
% axis([DoY1 DoY2 0.45 0.55])
axis([DoY1 DoY2 0.0 0.2])
% axis([DoY1 DoY2 0.4 0.6])
grid

figure(2)
subplot(2,1,1)
hold on
axis([DoY1 DoY2 0.0 0.2])
set(gca,'XLim',[DoY1 - 1 DoY2 + 1], 'FontSize', 12)

grid
    
Ninterval = length(ci_doy_SM) ;
for ii = 1 : 10 : Ninterval
    
    
    vsm = [ci_SM5(ii) ci_SM10(ii) ci_SM20(ii) ci_SM40(ii)] ;
    MSU_Multi_Layer_Model2(vsm, ci_doy_SM(ii)) ;
    
    figure(1)
    subplot(3,4,9:12)
    hold on
    plot(ci_doy_SM(ii), ci_SM5(ii), 'o')
        axis([DoY1 DoY2 0 0.5])

    M(ii) = getframe(fig);
    drawnow
    
    M2(ii) = getframe(fig2);
    drawnow
    
end

hf = figure ;
movie(hf,M);

% mplay(M)

%% add Dylan file path (containing movie support functions) and all subdirs
addpath( genpath( strcat(pwd, '\Dylan') ) );
M = delEmptyFrames(M);
writeMovie('test_movie.mat', M);
mplay(M);


end