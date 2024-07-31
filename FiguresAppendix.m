
colors = [255 51 51;
          255 153 51;
          051 153 255
          000 204 000 
          160 160 160
          000 000 000]./255;      

if GLP == 1
    plot_hyperparameter   
end


%% A.1

data2022     = data_storages(:,1);
data2023     = data_storages(:,2);
data_min1821 = data_storages(:,3);
data_max1821 = data_storages(:,4);

ticks = {'Oct','Nov','Dez','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct'};


dates = datetime('01.10.2022','Format','dd.MM.yyyy') + caldays((0:364));
dates = datetime(dates,'Format','MMM');

f = figure('Units','normalized','Position',[.2 .2 .6 .4]);
plot(dates,data2022,'LineWidth',2,'Color',colors(1,:)); hold on
plot(dates,data2023,'LineWidth',2,'Color',colors(2,:));
plot(dates,data_min1821,'LineWidth',2,'Color',colors(3,:),'LineStyle','--');
plot(dates,data_max1821,'LineWidth',2,'Color',colors(4,:),'LineStyle','--');
ylabel('%'); grid on; box off; 
ax = gca; ax.FontSize = 12; ax.FontName = 'Times'; ax.XTickLabel = ticks; axis tight
l = legend('2022','2023','Max. 2018-2021','Min. 2018-2021');
l.Box = 'off'; l.FontSize = 12; l.Location = 'southwest';
if save_figs == 1
    exportgraphics(f,'speicherstand.pdf','BackgroundColor','none');
end


%% A.2

dates = datetime('01.01.1999','Format','dd-MMM-yyyy') + calmonths((0:286))';


f = figure('Units','normalized','Position',[.2 .2 .6 .8]);
tl = tiledlayout(4,1,'Padding','compact','TileSpacing','compact'); nexttile;
plot(dates,(data_model(:,1)),'-k','Linewidth',2); title('Gas Supply Growth','FontWeight','normal'); hold on;
yline(0,'k--'); ylabel('%'); grid on; box off; set(gca,'FontSize',12,'FontName','Times');
nexttile;
plot(dates,(data_model(:,2)),'-k','Linewidth',2);title('Industrial Production','FontWeight','normal'); hold on;
set(gca,'FontSize',12,'FontName','Times'); ylabel('log x 100, 2019=100'); box off; grid on; 
nexttile
plot(dates,(data_model(:,3)),'-k','Linewidth',2); title('Real Gas Price Growth','FontWeight','normal'); hold on;
yline(0,'k--'); ylabel('%'); grid on; box off; set(gca,'FontSize',12,'FontName','Times');
nexttile
plot(dates,(data_model(:,4)),'-k','Linewidth',2); title('Change in Gas Inventories','FontWeight','normal'); hold on;
yline(0,'k--'); ylabel('%'); grid on; box off; set(gca,'FontSize',12,'FontName','Times');
if save_figs == 1
    exportgraphics(f,'Data_for_Matlab_dp_new.pdf','BackgroundColor','none')
end


%% A.3
 
dates = datetime(datetime('07.01.2021','Format','dd.MM.yyyy') + caldays((0:783))','Format','MM.yyyy');

f = figure('Units','normalized','Position',[.2 .2 .6 .4]);
plot(dates,data_flows(:,1),'LineWidth',2,'Color',colors(1,:)); hold on
plot(dates,data_flows(:,2),'LineWidth',2,'Color',colors(2,:));
plot(dates,data_flows(:,3),'LineWidth',2,'Color',colors(3,:));
plot(dates,data_flows(:,4),'LineWidth',2,'Color',colors(4,:));
plot(dates,data_flows(:,5),'LineWidth',2,'Color',colors(5,:));
plot(dates,data_flows(:,6),'LineWidth',2,'Color',colors(6,:));
l = legend('Jamal','Transgas','Nordstream 1','Europipe 1','Europipe 2','Eynatten');
ylabel('GWh/day'); grid on; box off; axis tight
l.Box = 'off'; l.FontSize = 12; 
ax = gca; ax.FontSize = 12; ax.FontName = 'Times';
if save_figs == 1
    exportgraphics(f,'Data_gas_flows.pdf','BackgroundColor','none');
end

    
%% A.4

data_2021 = data_consumption(2:end,2);
data_2022 = data_consumption(2:end,3);
data_2023 = data_consumption(2:end,4);
data_avg  = data_consumption(2:end,5);

ticks = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dez'};

f = figure('Units','normalized','Position',[.2 .2 .6 .4]);
plot(data_2021,'LineWidth',2,'Color',colors(1,:)); hold on
plot(data_2022,'LineWidth',2,'Color',colors(2,:));
plot(data_2023,'LineWidth',2,'Color',colors(3,:));
plot(data_avg,'k--','LineWidth',2);
ax = gca; ax.FontSize = 12; ax.FontName = 'Times'; ax.XTickLabel = ticks(1:1:end);
l = legend('2021','2022','2023','Average 2018-2021');
l.Box = 'off'; l.FontSize = 12; 
ylabel('GWh/day'); grid on; box off; axis tight;
if save_figs == 1
    exportgraphics(f,'verbrauch_IND.pdf','BackgroundColor','none');
end