

% load data_gasprices

data  = gasprices(6:end,8:11);
dates = datetime(datevec(labels_gasprices(6:end,2),'dd.mm.yyyy'),'Format','MM.yyyy');

t0 = datetime('01.2019','Format','MM.yyyy');
tT = datetime('04.2023','Format','MM.yyyy');

plot_periods = find(dates>=t0 & dates<=tT);
f = figure('Units','normalized','Position',[.2 .2 .6 .4]);
plot(dates,data(:,1),'k','LineWidth',2); hold on
plot(dates,data(:,4),'b','LineWidth',3); hold on
plot(dates,data(:,2),'r:','LineWidth',2); hold on
plot(dates,data(:,3),'Color',[255 153 51]./255,'LineWidth',2,'LineStyle','-.'); hold on
xlim([datetime('01.06.2020','InputFormat','dd.MM.yyyy')  datetime('31.12.2022','InputFormat','dd.MM.yyyy')]); ylim([0 1400]);
xtickformat('MM.yyyy'); grid on; box off; 
xticks(datetime('01.01.2020','InputFormat','dd.MM.yyyy'):calmonths(6):datetime('30.04.2023','InputFormat','dd.MM.yyyy'))
% xline(datetime('01.05.2020','InputFormat','dd.MM.yyyy'),'Label','Covid-recovery','LabelHorizontalAlignment','left','FontSize',12);
xline(datetime('24.02.2022','InputFormat','dd.MM.yyyy'),'Label','Russian Invasion of Ukraine','LabelHorizontalAlignment','left','FontSize',12);
xline(datetime('23.03.2022','InputFormat','dd.MM.yyyy'),'Label','Ramp-up gas of storages','LabelHorizontalAlignment','left','FontSize',12);
xline(datetime('14.07.2022','InputFormat','dd.MM.yyyy'),'Label',{'Stop of Russian','Gas supplies'},'LabelHorizontalAlignment','left','FontSize',12);
xline(datetime('22.09.2022','InputFormat','dd.MM.yyyy'),'Label',{'Gas storage','at 90%'},'LabelHorizontalAlignment','right','FontSize',12);
l = legend('TTF','Import price','Henry hub','JKM'); l.Location = 'Northwest'; l.Box = 'off'; l.FontSize = 14; l.NumColumns = 3;
set(gca,'FontSize',16,'FontName','Times')
rect = [0.005, 0.80, .5, .1];
set(l,'Position',rect)    
ylabel('€/MMBtu€ (2019=100)')
if save_figs == 1
    exportgraphics(f,'gasprices.eps','BackgroundColor','none'); close
end

