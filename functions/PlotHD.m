%% HD with Bands PointWise
% Compute HDs for sample MEAN temperature scenario (NO structural shocks)


Xtemp = exog;
Ta    = length(y);
exoga = ones(Ta,1);

time_trend = [];


%% baseline
% time_trend = [];
% 
% numSavedDraws = size(A0_save,3);
% 
% draws_HDs = nan(n,T-p,n,numSavedDraws);
% draws_DCs = nan(T-p,n,numSavedDraws);
% 
% for draw = 1:numSavedDraws
% 
%     A0    = A0_save(:,:,draw);
%     phi   = Beta_save(size(exog,2)+1:end,:,draw)';  
%     delta = Beta_save(1,:,draw)';
%     if trend == 1
%         time_trend = Beta_save(2,:,draw)'*Xtemp(:,2,1)';
%     end
%     temperature = Beta_save(2+trend,:,draw)'*Xtemp(:,2+trend,1)';
%     Dummies     = Beta_save(3+trend:13+trend,:,draw)'*Xtemp(:,3+trend:13+trend,1)';
% 
%     % get decomposition
%     [~,~,dc,~,sc]         = Get_SVAR_scenario(y,exoga,A0,delta,phi,p,0,[],time_trend,temperature,Dummies);
%     draws_HDs(:,:,:,draw) = sc(:,2:end,:);  % cumsum (for gassupply) n x T x shocks
%     draws_DCs(:,:,draw)   = dc(p+1:end,:);
%     if mod(draw,100) == 0
%         fprintf('Computing HD for traditional restrictions, draw %d of %d\n',draw,numSavedDraws)
%     end
% end
% 
% % output
% HD.dUnconditional_level = y(p+1:end,:) - draws_DCs;                         % total deviation from unconditional path
% HD.level                = permute(draws_HDs,[2,3,1,4]);

%% narrative restrictions

% point-wise
% numSavedDraws = size(A0_narrative,3);
% 
% draws_HDs_narrative = nan(n,T-p,n,numSavedDraws);
% draws_DCs_narrative = nan(T-p,n,numSavedDraws);
% 
% 
% for draw = 1:numSavedDraws
% 
%     A0    = A0_narrative(:,:,draw);
%     phi   = Beta_narrative(size(exog,2)+1:end,:,draw)';  % REPLACED CONSTANT WITH size(exog,2) WHICH INCLUDES THE CONSTANT
%     delta = Beta_narrative(1,:,draw)';
%     if trend == 1
%         time_trend = Beta_narrative(2,:,draw)'*Xtemp(:,2,1)';
%     end
%     temperature = Beta_narrative(2+trend,:,draw)'*Xtemp(:,2+trend,1)';
%     Dummies     = Beta_narrative(3+trend:13+trend,:,draw)'*Xtemp(:,3+trend:13+trend,1)';
% 
%     % get decomposition
%     [~,~,ic,dc,ec,tc,sc,shocks,~] = Get_SVAR_scenario(y,exoga,A0,delta,phi,p,0,[],time_trend,temperature,Dummies);
%     draws_HDs_narrative(:,:,:,draw)  = shocks;
%     draws_DCs_narrative(:,:,draw)    = ic + dc + ec + tc + sc;
%     if mod(draw,10) == 0
%         fprintf('Computing HD for narrative restrictions, draw %d of %d\n',draw,numSavedDraws)
%     end
% end

%----------------------------------------------------------------------------------------------------------------------------
% using Bayes-estimator of IRF
narrative = 1;
if narrative == 1
    bb  = Beta_narrative(:,:,In);
    aa  = A0_narrative(:,:,In);
else
    bb  = Beta_save(:,:,I);
    aa  = A0_save(:,:,I);
end
    
A0    = aa;
phi   = bb(size(exog,2)+1:end,:)'; 
delta = bb(1,:)';
if trend == 1
    time_trend = bb(2,:)'*Xtemp(:,2,1)';
end

clear temperature Dummies
if size(exog,2) > size(sa_dummies,2) + 1
    temperature = bb(1+trend+(1:p_temp+1),:)'*Xtemp(:,1+trend+(1:p_temp+1),1)';
    Dummies     = bb(2+trend+p_temp+(1:11),:)'*Xtemp(:,2+trend+p_temp+(1:11),1)';
else
    temperature = zeros(n,size(y,1));
    Dummies     = bb(2+trend+(1:11),:)'*Xtemp(:,2+trend+(1:11),1)';
end
  

% get decomposition
[~,~,ic,dc,ec,tc,sc,shocks,lrm] = Get_SVAR_scenario(y,exoga,A0,delta,phi,p,0,[],time_trend,temperature,Dummies);
draws_HDs_narrative(:,:,:,1)    = shocks(:,2:end,:);
draws_DCs_narrative(:,:,1)      = ic(2:end,:) + dc(2:end,:) + ec(2:end,:) + sc(2:end,:) + tc(2:end,:);

% output
HD.narrative.deterministic        = draws_DCs_narrative;
HD.narrative.dUnconditional_level = y(p+1:end,:) - draws_DCs_narrative;                         % total deviation from unconditional path
HD.narrative.stochastic           = permute(draws_HDs_narrative,[2,3,1,4]);

%% plotting
dates_plot     = datetime(datevec(dates(p+1:end)));
period_first   = '01.2020';
period_last    = '12.2022';
t0             = datetime(period_first,'Format','MM.yyyy');
tT             = datetime(period_last,'Format','MM.yyyy');
plot_periods   = find(dates_plot>=t0 & dates_plot<=tT);

colors = [255 51 51;
          255 153 51;
           51 51 255;
           160 160 160;
           32  32  32]./255;

cumHDs = 0;

f  = figure('Units','normalized','Position',[.1 .1 .5 .7]);   
tl = tiledlayout(floor(n/2),n-floor(n/2),'TileSpacing','compact','Padding','tight');
for i = 1:n
    nexttile;
    Hdtmp = squeeze(HD.narrative.stochastic(:,:,i,:));                         % bayes estimate
    unc   = mean(HD.narrative.dUnconditional_level(:,i,:),3);

    if ismember(i,cumulateWhich) && cumHDs == 1
        total      = cumsum(Hdtmp(plot_periods,:));
        unc        = cumsum(unc(plot_periods));
    else
        total      = (Hdtmp(plot_periods,:));   
        unc        = unc(plot_periods);
    end
    neg        = total;
    neg(neg>0) = 0;
    pos        = total;
    pos(pos<0) = 0;
    
    bP = bar(dates_plot(plot_periods),pos,'stack','Barwidth',.8,'EdgeColor','flat','FaceColor','flat','FaceAlpha',.7); hold on
    bN = bar(dates_plot(plot_periods),neg,'stack','Barwidth',.8,'EdgeColor','flat','FaceColor','flat','FaceAlpha',.7); hold on
    lU = plot(dates_plot(plot_periods),unc,'k-','LineWidth',2.5);
    pz  = yline(0,'k-','LineWidth',2);   grid on
    ax = gca;   ax.FontName = 'Times'; ax.FontSize = 14;
    ax.YLim(1) = round(ax.YLim(1) - sum(abs(ax.YLim))/2*.1,2);  ax.FontName = 'Times'; ax.XTickLabelRotation = 45; ax.GridLineStyle = ':'; grid on; 
    xtickformat('MM.yyyy');
    for j = 1:size(pos,2)
        bP(j).CData = colors(j,:);
        bN(j).CData = bP(j).CData;
    end
    title(sprintf(varNames{i}),'FontName','Times','FontSize',16,'FontWeight','normal');
    if ismember(i,[1,3])
        ylabel('m-o-m in %','FontName','Times','FontSize',14)
    end
    box off
    if i == n-1
        l = legend([bP,lU],[shockNames_for_HD 'Detrended Data'],'Location','SouthWest','Orientation','Horizontal','Box','off','NumColumns',1,'Interpreter','tex','FontSize',12);
    end
end
% title(tl,'Historical decomposition 2020-2022','FontName','Times','FontSize',20,'FontWeight','bold');
% l = legend([bP,lU],[shockNames_for_HD '\Delta Unconditional'],'Location','SouthWest','Orientation','Horizontal','Box','off','NumColumns',1,'Interpreter','tex','FontSize',12);
% rect = [0.0000000001, 0.08, .8, .1];
% set(l,'Position',rect)    

if save_figs == 1
    exportgraphics(f,'Output/Revision/HD_2019_2022_dlIP_tlag.pdf')
end

t0           = datetime('02.2000','Format','MM.yyyy');
tT           = datetime('12.2018','Format','MM.yyyy');
dates_plot   = datetime(dates_plot,'Format','MM.yyyy');
plot_periods = find(dates_plot>=t0 & dates_plot<=tT);


% cumulative HDs
% for i = 1:n
%     Hdtmp = squeeze(HD.narrative.stochastic(:,:,i,:));                         % bayes estimate    
%     data_demeaned = mean(HD.narrative.dUnconditional_level(:,i,:),3);    
%     fig   = figure('Units','normalized','Position',[.1 .1 .5 .5]);   
%     tl    = tiledlayout(n,1);
%     for j = 1:n
%         nexttile;
%         plot(dates_plot(plot_periods),Hdtmp(plot_periods,j),'r','Linewidth',2); grid on; hold on
%         plot(dates_plot(plot_periods),data_demeaned(plot_periods),'k--','Linewidth',2);
%         yline(0,'k-','LineWidth',2);
%         ax = gca; ax.FontSize = 14; ax.FontName = 'Times'; 
%         tt = title(sprintf('Cumulative Effect of %s Shocks on %s',shockNames{j},varNames{i}));
%         tt.Interpreter = 'latex'; tt.FontSize = 14; box off       
%     end
%     tl.TileSpacing    = 'compact';
%     tl.Padding        = 'compact';
%     tmp               = varNames{i};    % remove spaces from variable names for saving
%     tmp(isspace(tmp)) = '_';
%     if save_figs == 1
%         saveas(fig,sprintf('Output/New Baseline/HDcum_%s',tmp),'epsc'); close
%     end
% end

%% HDs as bar chart

dates_plot = datetime(dates_plot,'Format','MM.yyyy');
t0_1       = datetime('02.2000','Format','MM.yyyy');
tT_1       = datetime('12.2019','Format','MM.yyyy');
t0_2       = datetime('01.2021','Format','MM.yyyy');
tT_2       = datetime('12.2022','Format','MM.yyyy');

plot_periods_1 = find(dates_plot>=t0_1 & dates_plot<=tT_1);
plot_periods_2 = find(dates_plot>=t0_2 & dates_plot<=tT_2);

% estimates in levels (only for log-difference specification)
for i = 1:n
    data_demeaned = mean(squeeze(HD.narrative.dUnconditional_level(plot_periods_1(1):plot_periods_1(end),i,:)),2);                         % total cumulative change of non-deterministic components
    data_shocks   = HD.narrative.stochastic(plot_periods_1(1):plot_periods_1(end),:,i,:);
    if ismember(i,cumulateWhich) ~=0
        data_demeaned = cumsum(data_demeaned);
        data_shocks   = cumsum(data_shocks);
    end
    D_1(i,1)      = data_shocks(end,1) - data_shocks(1,1);                                                                                       
    D_1(i,2)      = data_shocks(end,2) - data_shocks(1,2);                                   
    D_1(i,3)      = data_shocks(end,3) - data_shocks(1,3);
    D_1(i,4)      = data_shocks(end,4) - data_shocks(1,4);
    D_1true(i)    = data_demeaned(end) - data_demeaned(1); 

    data_demeaned = mean(squeeze(HD.narrative.dUnconditional_level(plot_periods_2(1):plot_periods_2(end),i,:)),2);                         % total variance due to stochastic part of the model
    data_shocks   = HD.narrative.stochastic(plot_periods_2(1):plot_periods_2(end),:,i,:);
    if ismember(i,cumulateWhich) ~=0
        data_demeaned = cumsum(data_demeaned);
        data_shocks   = cumsum(data_shocks);
    end    
    D_2(i,1)      = data_shocks(end,1) - data_shocks(1,1);    
    D_2(i,2)      = data_shocks(end,2) - data_shocks(1,2);
    D_2(i,3)      = data_shocks(end,3) - data_shocks(1,3);
    D_2(i,4)      = data_shocks(end,4) - data_shocks(1,4);   
    D_2true(i)    = data_demeaned(end) - data_demeaned(1);     
end

if max(sum(D_1,2)' - D_1true) > 10E-3, warning('There seems to be an error in the bar charts'); end
if max(sum(D_2,2)' - D_2true) > 10E-3, warning('There seems to be an error in the bar charts'); end

% long sample
% figure('Units','normalized','Position',[.1 .1 .5 .7])
% tiledlayout(2,2,'TileSpacing','loose','Padding','compact');
% for i = 1:n 
%     nexttile;
%     bb = bar(D_1(i,:),'FaceColor','k','FaceAlpha',.5); grid on; box off;
%     ax = gca; 
% 
%     if ismember(i,[1,3]) == 1
%         ax.YLabel.String = 'Cumulative change';
%     end
%     ax.XTickLabel    = [shockNames(1),shockNames(2),shockNames(3),shockNames(4)];
%     ax.Title.String  = varNames{i};
% end
% 
% % short sample
% figure('Units','normalized','Position',[.1 .1 .5 .7])
% tiledlayout(2,2,'TileSpacing','loose','Padding','compact');
% for i = 1:n 
%     nexttile;
%     bb = bar(D_2(i,:),'FaceColor','k','FaceAlpha',.5); grid on; box off;
%     ax = gca; 
% 
%     if ismember(i,[1,3]) == 1
%         ax.YLabel.String = 'Cumulative change';
%     end
%     ax.XTickLabel    = [shockNames(1),shockNames(2),shockNames(3),shockNames(4)];
%     ax.Title.String  = varNames{i};
% end

% both in one graph
% short sample
f = figure('Units','normalized','Position',[.1 .1 .5 .7]);
tiledlayout(floor(n/2),n-floor(n/2),'TileSpacing','compact','Padding','tight');
for i = 1:n 
    nexttile;
    b1 = bar([D_1(i,:); D_2(i,:)]','FaceColor','k','FaceAlpha',.8,'BarLayout','grouped');
    b1(1).FaceColor = 'k'; b1(1).FaceAlpha = .7;
    b1(2).FaceColor = 'b'; b1(2).FaceAlpha = .3;
    grid on; box off;
    ax = gca; 
    if ismember(i,[1,3]) == 1
        ax.YLabel.String = 'Cumulative effect';
    end
    for j = 1:length(shockNames_for_HD)
        ax.XTickLabel(j) = shockNames_for_HD(j);
    end
    ax.Title.String  = varNames{i};
    if i == 1
        l = legend('2000-2019','2021-2022'); l.Box = 'off'; l.FontSize = 12; l.Location = 'Northeast'; l.Interpreter = 'tex';
    end
end

if save_figs == 1
    exportgraphics(f,'Output/Revision/HDbars_dlIP_tlag.pdf')
end
