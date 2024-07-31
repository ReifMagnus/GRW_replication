%% HD with Bands PointWise and Bayesian Estimate

%% housekeeping
clear Xtemp

temp_scen     = Scenario.temp;
dummies       = Scenario.dummies;
scenario      = Scenario.shocks;

baseline = zeros(size(scenario));
idxStart = find(ismember(dates,startDate));
idxEnd   = find(ismember(dates,endDate));
idxStand = find(ismember(dates_stand,startDate));

Ytemp = y(idxStart-p:idxEnd,:);

% construct matrix of lagged temperatures
clear temp_scen_lags
if lagged_temperature == 1
    for i = 1:size(temp_scen,2)
        temp_scen_lags(:,:,i) = lagmatrix(temp_scen(:,i),0:p_temp);
        for j = 2:p_temp+1
            temp_scen_lags(1:j-1,j,i) = exog(idxEnd-j+2:idxEnd,2);
        end
    end
else
    for i = 1:size(temp_scen,2)
        temp_scen_lags(:,:,i) = temp_scen(:,i);
    end
end

if trend == 1
    Xtemp = [exog(idxStart-p:idxEnd,:); ones(h,1), (exog(idxEnd,2)+1:exog(idxEnd,2)+h)', temp_scen_lags(:,:,1), dummies]; % Append MEAN temperature scenario and dummies to realized exogenous variables
else
    Xtemp = [exog(idxStart-p:idxEnd,:); ones(h,1), temp_scen_lags(:,:,1), dummies]; % Append MEAN temperature scenario and dummies to realized exogenous variables
end
Ta    = length(Ytemp);
exoga = ones(Ta,1);

%% Compute HDs
%% sign restrictions
numSavedDraws = size(A0_save,3);

[draws_HDs_base,draws_HDs_scen]                     = deal(nan(Ta-p+1+h,n,n,numSavedDraws)); % time, shock, draws
[draws_dc_base,draws_dc_scen,uncond_path,scen_path] = deal(nan(Ta-p+1+h,n,numSavedDraws)); % time, shock, draws

[supply_base,price_base,inv_base,supply_scen,ip_scen,price_scen,inv_scen] = deal(nan(h,1,numSavedDraws));
[gas_base,ip_base,gas_scen]                                               = deal(nan(h,1,numSavedDraws));

[gas_base_prc,gas_scen_prc,inv_base_prc,inv_scen_prc,price_base_prc,price_scen_prc,ip_base_prc,ip_scen_prc] = deal(nan(Ta,3));
[gas_base_BE,gas_scen_BE,inv_base_BE,inv_scen_BE,price_base_BE,price_scen_BE,ip_base_BE,ip_scen_BE,]        = deal(nan(Ta,1));

for draw = 1:numSavedDraws
    if mod(draw,100) == 0
        fprintf('Running (sign) draw %d\n',draw);
    end
    %% get coefficients (signs)
    A0    = A0_save(:,:,draw);
    phi   = Beta_save(size(exog,2)+1:end,:,draw)';
    delta = Beta_save(1,:,draw)';
    if trend == 1
        time_trend = Beta_save(2,:,draw)'*Xtemp(:,2)';
    else
        time_trend = [];
    end

    clear temperature Dummies
    for i = 1:size(Xtemp,3)
        if size(exog,2) > size(sa_dummies,2) + 1
            temperature(:,:,i) = Beta_save(1+trend+(1:p_temp+1),:,draw)'*Xtemp(:,1+trend+(1:p_temp+1),i)';
            Dummies(:,:,i)     = Beta_save(2+trend+p_temp+(1:11),:,draw)'*Xtemp(:,2+trend+p_temp+(1:11),i)';
        else
            temperature(:,:,i) = zeros(n,size(y,1));
            Dummies(:,:,i)     = Beta_save(2+trend+(1:11),:,draw)'*Xtemp(:,2+trend+(1:11),i)';
        end
    end

    % unconditional decomposition
    [~,~,ic,dc,ec,tc,sc,shocks,~] = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,baseline,time_trend,temperature(:,:,1),Dummies(:,:,1));
    draws_HDs_base(:,:,:,draw) = permute(shocks,[2,1,3]);
    draws_dc_base(:,:,draw)    = ic + dc + ec + tc + sc;

    % decomposition of scenario
    [~,~,ic,dc,ec,tc,sc,shocks,lrm] = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,scenario,time_trend,temperature(:,:,1),Dummies(:,:,1));
    draws_HDs_scen(:,:,:,draw)        = permute(shocks,[2,1,3]);  % cumsum (for gassupply)
    draws_dc_scen(:,:,draw)           = ic + dc + ec + tc + sc;

    %% Compute time paths of transformed variables

    
    stochastic            = squeeze(sum(draws_HDs_base(:,:,:,draw),3));    % stochastic component (w/o scenario)
    deterministic         = draws_dc_base(:,:,draw);                       % determinsitic component
    uncond_path(:,:,draw) = deterministic + stochastic;               % path without scenario

    stochastic_scen      = squeeze(sum(draws_HDs_scen(:,:,:,draw),3));    % stochastic component (with scenario)
    scen_path (:,:,draw) = deterministic + stochastic_scen;           % path with scenario

    %% Compute time paths of original variables
    supply_base(:,:,draw) = squeeze(uncond_path(2:end-1,1,draw))/100; % timing t-1 for h = 1,...,12
    ip_base(:,:,draw)     = squeeze(uncond_path(3:end,2,draw));       % timing t for h = 1,...,12
    price_base(:,:,draw)  = squeeze(uncond_path(3:end,3,draw));    % timing t for h = 1,...,12
    inv_base(:,:,draw)    = squeeze(uncond_path(3:end,4,draw));  % timing t for h = 1,...,12

    supply_scen(:,:,draw) = squeeze(scen_path(2:end-1,1,draw))/100;   % timing t-1 for h = 1,...,12
    ip_scen(:,:,draw)     = squeeze(scen_path(3:end,2,draw));         % timing t for h = 1,...,12
    price_scen(:,:,draw)  = squeeze(scen_path(3:end,3,draw));      % timing t for h = 1,...,12
    inv_scen(:,:,draw)    = squeeze(scen_path(3:end,4,draw));    % timing t for h = 1,...,12

    gas_base(1,:,draw) = ones(size(supply_base(1,:,draw)))*netgas_level(idxEnd,1);
    gas_scen(1,:,draw) = ones(size(supply_scen(1,:,draw)))*netgas_level(idxEnd,1);
    for i = 2:h
        gas_base(i,:,draw) = (1 + supply_base(i,:,draw)).*gas_base(i-1,:,draw);
        gas_scen(i,:,draw) = (1 + supply_scen(i,:,draw)).*gas_scen(i-1,:,draw);
    end
    gas_base(h+1,:,draw) = (1 + supply_base(h,:,draw)).*gas_base(h-1,:,draw);
    gas_scen(h+1,:,draw) = (1 + supply_scen(h,:,draw)).*gas_scen(h-1,:,draw);
    inv_base(:,:,draw)   = stand(idxStand) + cumsum(inv_base(:,:,draw),1)*277.778*1000*100/230000000;    % Transform to capacity (230 mio. MWh in 2020)
    inv_scen(:,:,draw)   = stand(idxStand) + cumsum(inv_scen(:,:,draw),1)*277.778*1000*100/230000000;
    price_base(:,:,draw) = cumsum(price_base(:,:,draw));
    price_scen(:,:,draw) = cumsum(price_scen(:,:,draw));
    ip_base(:,:,draw)    = cumsum(ip_base(:,:,draw));
    ip_scen(:,:,draw)    = cumsum(ip_scen(:,:,draw));

end

for i = 1:h
    gas_base_prc(i,1) = prctile(gas_base(i,1,:),50);
    gas_base_prc(i,2) = prctile(gas_base(i,1,:),84);
    gas_base_prc(i,3) = prctile(gas_base(i,1,:),16);
    gas_base_BE(i,1)  = gas_base(i,1,I);

    gas_scen_prc(i,1) = prctile(gas_scen(i,1,:),50);
    gas_scen_prc(i,2) = prctile(gas_scen(i,1,:),84);
    gas_scen_prc(i,3) = prctile(gas_scen(i,1,:),16);
    gas_scen_BE(i,1)  = gas_scen(i,1,I);

    inv_base_prc(i,1) = prctile(inv_base(i,1,:),50);
    inv_base_prc(i,2) = prctile(inv_base(i,1,:),84);
    inv_base_prc(i,3) = prctile(inv_base(i,1,:),16);
    inv_base_BE(i,1)  = inv_base(i,1,I);

    inv_scen_prc(i,1) = prctile(inv_scen(i,1,:),50);
    inv_scen_prc(i,2) = prctile(inv_scen(i,1,:),84);
    inv_scen_prc(i,3) = prctile(inv_scen(i,1,:),16);
    inv_scen_BE(i,1)  = inv_scen(i,1,I);

    price_base_prc(i,1) = prctile(price_base(i,1,:),50);
    price_base_prc(i,2) = prctile(price_base(i,1,:),84);
    price_base_prc(i,3) = prctile(price_base(i,1,:),16);
    price_base_BE(i,1)  = price_base(i,1,I);

    price_scen_prc(i,1) = prctile(price_scen(i,1,:),50);
    price_scen_prc(i,2) = prctile(price_scen(i,1,:),84);
    price_scen_prc(i,3) = prctile(price_scen(i,1,:),16);
    price_scen_BE(i,1)  = price_scen(i,1,I);

    ip_base_prc(i,1) = prctile(ip_base(i,1,:),50);
    ip_base_prc(i,2) = prctile(ip_base(i,1,:),84);
    ip_base_prc(i,3) = prctile(ip_base(i,1,:),16);
    ip_base_BE(i,1)  = ip_base(i,1,I);

    ip_scen_prc(i,1) = prctile(ip_scen(i,1,:),50);
    ip_scen_prc(i,2) = prctile(ip_scen(i,1,:),84);
    ip_scen_prc(i,3) = prctile(ip_scen(i,1,:),16);
    ip_scen_BE(i,1)  = ip_scen(i,1,I);
end
gas_base_prc(h+1,1) = prctile(gas_base(h+1,1,:),50);
gas_base_prc(h+1,2) = prctile(gas_base(h+1,1,:),84);
gas_base_prc(h+1,3) = prctile(gas_base(h+1,1,:),16);
gas_base_BE(h+1,1)  = gas_base(h+1,1,I);

gas_scen_prc(h+1,1) = prctile(gas_scen(h+1,1,:),50);
gas_scen_prc(h+1,2) = prctile(gas_scen(h+1,1,:),84);
gas_scen_prc(h+1,3) = prctile(gas_scen(h+1,1,:),16);
gas_scen_BE(h+1,1)  = gas_scen(h+1,1,I);

 
%% narrative restrictions
numSavedDraws = size(A0_narrative,3);

[narrative_draws_HDs_base,narrative_draws_HDs_scen]                                           = deal(nan(Ta-p+1+h,n,n,numSavedDraws)); % time, shock, draws
[narrative_draws_dc_base,narrative_draws_dc_scen,narrative_uncond_path,narrative_scen_path]   = deal(nan(Ta-p+1+h,n,numSavedDraws)); % time, shock, draws

[narrative_supply_base,narrative_price_base,narrative_inv_base,narrative_supply_scen,narrative_ip_scen,narrative_price_scen,narrative_inv_scen] = deal(nan(h,1,numSavedDraws));
[narrative_gas_base,narrative_ip_base,narrative_gas_scen]                                                       = deal(nan(h,1,numSavedDraws));

[narrative_gas_base_prc,narrative_gas_scen_prc,narrative_inv_base_prc,narrative_inv_scen_prc,narrative_price_base_prc,narrative_price_scen_prc,narrative_ip_base_prc,narrative_ip_scen_prc] = deal(nan(h,3));
[narrative_gas_base_BE,narrative_gas_scen_BE,narrative_inv_base_BE,narrative_inv_scen_BE,narrative_price_base_BE,narrative_price_scen_BE,narrative_ip_base_BE,narrative_ip_scen_BE,]        = deal(nan(h,1));

for draw = 1:numSavedDraws
    if mod(draw,100) == 0
        fprintf('Running (narrative) draw %d\n',draw);
    end

    %% get coefficients (narratives)
    A0    = A0_narrative(:,:,draw);
    phi   = Beta_narrative(size(exog,2)+1:end,:,draw)';
    delta = Beta_narrative(1,:,draw)';
    if trend == 1
        time_trend = Beta_narrative(2,:,draw)'*Xtemp(:,2)';
    else
        time_trend = [];
    end

    clear temperature Dummies
    for i = 1:size(Xtemp,3)
        if size(exog,2) > size(sa_dummies,2) + 1
            temperature(:,:,i) = Beta_narrative(1+trend+(1:p_temp+1),:,draw)'*Xtemp(:,1+trend+(1:p_temp+1),i)';
            Dummies(:,:,i)     = Beta_narrative(2+trend+p_temp+(1:11),:,draw)'*Xtemp(:,2+trend+p_temp+(1:11),i)';
        else
            temperature(:,:,i) = zeros(n,size(y,1));
            Dummies(:,:,i)     = Beta_narrative(2+trend+(1:11),:,draw)'*Xtemp(:,2+trend+(1:11),i)';
        end
    end

    % unconditional decomposition
    [~,~,ic,dc,ec,tc,sc,shocks,~]        = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,baseline,time_trend,temperature(:,:,1),Dummies(:,:,1));
    narrative_draws_HDs_base(:,:,:,draw) = permute(shocks,[2,1,3]);
    narrative_draws_dc_base(:,:,draw)    = ic + dc + ec + tc + sc;

    % decomposition of scenario
    [~,~,ic,dc,ec,tc,sc,shocks,~]               = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,scenario,time_trend,temperature(:,:,1),Dummies(:,:,1));
    narrative_draws_HDs_scen(:,:,:,draw)        = permute(shocks,[2,1,3]);
    narrative_draws_dc_scen(:,:,draw)           = ic + dc + ec + tc + sc;

    %% Compute time paths of transformed variables

    narrative_stochastic            = squeeze(sum(narrative_draws_HDs_base(:,:,:,draw),3));   % stochastic component (w/o scenario)
    narrative_deterministic         = narrative_draws_dc_base(:,:,draw);                      % deterministic component
    narrative_uncond_path(:,:,draw) = narrative_deterministic + narrative_stochastic;     % path without scenario

    narrative_stochastic_scen     = squeeze(sum(narrative_draws_HDs_scen(:,:,:,draw),3));   % stochastic component (with scenario)
    narrative_scen_path(:,:,draw) = narrative_deterministic + narrative_stochastic_scen;  % path with scenario

    %% Compute time paths of original variables   
    narrative_supply_base(:,:,draw) = squeeze(narrative_uncond_path(2:end-1,1,draw))/100;  % timing t-1 for h = 1,...,12
    narrative_ip_base(:,:,draw)     = squeeze(narrative_uncond_path(3:end,2,draw));        % timing t for h = 1,...,12
    narrative_price_base(:,:,draw)  = squeeze(narrative_uncond_path(3:end,3,draw));     % timing t for h = 1,...,12
    narrative_inv_base(:,:,draw)    = squeeze(narrative_uncond_path(3:end,4,draw));   % timing t for h = 1,...,12

    narrative_supply_scen(:,:,draw) = squeeze(narrative_scen_path(2:end-1,1,draw))/100;  % timing t-1 for h = 1,...,12
    narrative_ip_scen(:,:,draw)     = squeeze(narrative_scen_path(3:end,2,draw));        % timing t for h = 1,...,12
    narrative_price_scen(:,:,draw)  = squeeze(narrative_scen_path(3:end,3,draw));     % timing t for h = 1,...,12
    narrative_inv_scen(:,:,draw)    = squeeze(narrative_scen_path(3:end,4,draw));   % timing t for h = 1,...,12

    narrative_gas_base(1,:,draw) = ones(size(narrative_supply_base(1,:,draw)))*netgas_level(idxEnd,1);
    narrative_gas_scen(1,:,draw) = ones(size(narrative_supply_scen(1,:,draw)))*netgas_level(idxEnd,1);
    for i = 2:h
        narrative_gas_base(i,:,draw) = (1 + narrative_supply_base(i,:,draw)).*narrative_gas_base(i-1,:,draw);
        narrative_gas_scen(i,:,draw) = (1 + narrative_supply_scen(i,:,draw)).*narrative_gas_scen(i-1,:,draw);
    end
    narrative_gas_base(h+1,:,draw) = (1 + narrative_supply_base(h,:,draw)).*narrative_gas_base(h-1,:,draw);
    narrative_gas_scen(h+1,:,draw) = (1 + narrative_supply_scen(h,:,draw)).*narrative_gas_scen(h-1,:,draw);
    narrative_inv_base(:,:,draw)   = stand(idxStand) + cumsum(narrative_inv_base(:,:,draw),1)*277.778*1000*100/230000000;    % Transform to capacity (230 mio. MWh in 2020)
    narrative_inv_scen(:,:,draw)   = stand(idxStand) + cumsum(narrative_inv_scen(:,:,draw),1)*277.778*1000*100/230000000;
    narrative_price_base(:,:,draw) = cumsum(narrative_price_base(:,:,draw));
    narrative_price_scen(:,:,draw) = cumsum(narrative_price_scen(:,:,draw));
    narrative_ip_base(:,:,draw)    = cumsum(narrative_ip_base(:,:,draw));
    narrative_ip_scen(:,:,draw)    = cumsum(narrative_ip_scen(:,:,draw));

end

for i = 1:h
    narrative_gas_base_prc(i,1) = prctile(narrative_gas_base(i,1,:),50);
    narrative_gas_base_prc(i,2) = prctile(narrative_gas_base(i,1,:),84);
    narrative_gas_base_prc(i,3) = prctile(narrative_gas_base(i,1,:),16);
    narrative_gas_base_BE(i,1)  = narrative_gas_base(i,1,In);
    
    narrative_gas_scen_prc(i,1) = prctile(narrative_gas_scen(i,1,:),50);
    narrative_gas_scen_prc(i,2) = prctile(narrative_gas_scen(i,1,:),84);
    narrative_gas_scen_prc(i,3) = prctile(narrative_gas_scen(i,1,:),16);
    narrative_gas_scen_BE(i,1)  = narrative_gas_scen(i,1,In);
    
    narrative_inv_base_prc(i,1) = prctile(narrative_inv_base(i,1,:),50);
    narrative_inv_base_prc(i,2) = prctile(narrative_inv_base(i,1,:),84);
    narrative_inv_base_prc(i,3) = prctile(narrative_inv_base(i,1,:),16);
    narrative_inv_base_BE(i,1)  = narrative_inv_base(i,1,In);
    
    narrative_inv_scen_prc(i,1) = prctile(narrative_inv_scen(i,1,:),50);
    narrative_inv_scen_prc(i,2) = prctile(narrative_inv_scen(i,1,:),84);
    narrative_inv_scen_prc(i,3) = prctile(narrative_inv_scen(i,1,:),16);
    narrative_inv_scen_BE(i,1)  = narrative_inv_scen(i,1,In);
    
    narrative_price_base_prc(i,1) = prctile(narrative_price_base(i,1,:),50);
    narrative_price_base_prc(i,2) = prctile(narrative_price_base(i,1,:),84);
    narrative_price_base_prc(i,3) = prctile(narrative_price_base(i,1,:),16);
    narrative_price_base_BE(i,1)  = narrative_price_base(i,1,In);
    
    narrative_price_scen_prc(i,1) = prctile(narrative_price_scen(i,1,:),50);
    narrative_price_scen_prc(i,2) = prctile(narrative_price_scen(i,1,:),84);
    narrative_price_scen_prc(i,3) = prctile(narrative_price_scen(i,1,:),16);
    narrative_price_scen_BE(i,1)  = narrative_price_scen(i,1,In);
    
    narrative_ip_base_prc(i,1) = prctile(narrative_ip_base(i,1,:),50);
    narrative_ip_base_prc(i,2) = prctile(narrative_ip_base(i,1,:),84);
    narrative_ip_base_prc(i,3) = prctile(narrative_ip_base(i,1,:),16);
    narrative_ip_base_BE(i,1)  = narrative_ip_base(i,1,In);
    
    narrative_ip_scen_prc(i,1) = prctile(narrative_ip_scen(i,1,:),50);
    narrative_ip_scen_prc(i,2) = prctile(narrative_ip_scen(i,1,:),84);
    narrative_ip_scen_prc(i,3) = prctile(narrative_ip_scen(i,1,:),16);
    narrative_ip_scen_BE(i,1)  = narrative_ip_scen(i,1,In);
end
narrative_gas_base_prc(h+1,1) = prctile(narrative_gas_base(h+1,1,:),50);
narrative_gas_base_prc(h+1,2) = prctile(narrative_gas_base(h+1,1,:),84);
narrative_gas_base_prc(h+1,3) = prctile(narrative_gas_base(h+1,1,:),16);
narrative_gas_base_BE(h+1,1)  = narrative_gas_base(h+1,1,In);

narrative_gas_scen_prc(h+1,1) = prctile(narrative_gas_scen(h+1,1,:),50);
narrative_gas_scen_prc(h+1,2) = prctile(narrative_gas_scen(h+1,1,:),84);
narrative_gas_scen_prc(h+1,3) = prctile(narrative_gas_scen(h+1,1,:),16);
narrative_gas_scen_BE(h+1,1)  = narrative_gas_scen(h+1,1,In);

dates_plot = dates_f(idxStart+1:idxEnd+h);
dates_tick = dates_plot([1,7,12],1,1);


tt    = datetime(datevec(dates_plot));
ytemp = [y(idxEnd+1:ix_last,:); nan(idxEnd+h-ix_last,n)];
xband = [tt(1):calmonths(1):tt(end) tt(end):-calmonths(1):tt(1)];

f = figure('Units','normalized','Position',[.1 .1 .5 .3]);
tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile;
percentiles1 = narrative_gas_base_prc(2:end,2:3)/1000;
percentiles2 = narrative_gas_scen_prc(2:end,2:3)/1000;
yband1       = [percentiles1(:,2)' percentiles1(end:-1:1,1)'];
yband2       = [percentiles2(:,2)' percentiles2(end:-1:1,1)'];
p1 = fill(xband,yband1,LightBlue,'FaceAlpha',.3,'EdgeColor',LightBlue,'EdgeAlpha',.3); hold on
p2 = fill(xband,yband2,LightRed,'FaceAlpha',.3,'EdgeColor',LightRed,'EdgeAlpha',.3); 
p3 = plot(datetime(datevec(dates_plot)),narrative_gas_base_prc(2:end,1)/1000,'b-','LineWidth',2);
p4 = plot(datetime(datevec(dates_plot)),narrative_gas_scen_prc(2:end,1)/1000,'r:','LineWidth',2);
% add data
if add_data == 1
    gas_data(1,1) = (1+ytemp(1,1)/100)*netgas_level(idxEnd,1);
    for i = 2:h
        gas_data(i,1) = (1 + ytemp(i,1)/100).*gas_data(i-1,1);
    end
    p5 = plot(datetime(datevec(dates_plot)),gas_data/1000,'k--','linewidth',2,'DisplayName','Data');
end
ax = gca; ax.Box = 'off'; ax.FontSize = 10; ax.YLabel.String = '1000 TJ';
ax.Title.String = varNames(1); ax.Title.FontSize = 12; ax.Title.FontName = 'Times';
grid on;  xtickformat('MMM-yyyy');

nexttile;
percentiles1 = narrative_ip_base_prc(:,2:3);
percentiles2 = narrative_ip_scen_prc(:,2:3);
yband1       = [percentiles1(:,2)' percentiles1(end:-1:1,1)'];
yband2       = [percentiles2(:,2)' percentiles2(end:-1:1,1)'];
p1 = fill(xband,yband1,LightBlue,'FaceAlpha',.3,'EdgeColor',LightBlue,'EdgeAlpha',.3); hold on
p2 = fill(xband,yband2,LightRed,'FaceAlpha',.3,'EdgeColor',LightRed,'EdgeAlpha',.3); 
p3 = plot(datetime(datevec(dates_plot)),narrative_ip_base_prc(:,1),'b-','LineWidth',2);
p4 = plot(datetime(datevec(dates_plot)),narrative_ip_scen_prc(:,1),'r:','LineWidth',2);
% add data
if add_data == 1
    p5 = plot(datetime(datevec(dates_plot)),cumsum(ytemp(:,2)),'k--','linewidth',2,'DisplayName','Data');
end
ax = gca; ax.Box = 'off'; ax.FontSize = 10; ax.YLabel.String = 'cumulated %';
ax.Title.String = varNames(2); ax.Title.FontSize = 12; ax.Title.FontName = 'Times';
grid on;  xtickformat('MMM-yyyy');


nexttile;
percentiles1 = narrative_price_base_prc(:,2:3);
percentiles2 = narrative_price_scen_prc(:,2:3);
yband1       = [percentiles1(:,2)' percentiles1(end:-1:1,1)'];
yband2       = [percentiles2(:,2)' percentiles2(end:-1:1,1)'];
p1 = fill(xband,yband1,LightBlue,'FaceAlpha',.3,'EdgeColor',LightBlue,'EdgeAlpha',.3); hold on
p2 = fill(xband,yband2,LightRed,'FaceAlpha',.3,'EdgeColor',LightRed,'EdgeAlpha',.3); 
p3 = plot(datetime(datevec(dates_plot)),narrative_price_base_prc(:,1),'b-','LineWidth',2);
p4 = plot(datetime(datevec(dates_plot)),narrative_price_scen_prc(:,1),'r:','LineWidth',2);
% add data
if add_data == 1
    p5 = plot(datetime(datevec(dates_plot)),cumsum(ytemp(:,3)),'k--','linewidth',2,'DisplayName','Data');
end
ax = gca; ax.Box = 'off'; ax.FontSize = 10; ax.YLabel.String = 'cumulated %';
ax.Title.String = varNames(3); ax.Title.FontSize = 12; ax.Title.FontName = 'Times';
grid on; xtickformat('MMM yyyy');

nexttile;
percentiles1 = narrative_inv_base_prc(:,2:3);
percentiles2 = narrative_inv_scen_prc(:,2:3);
yband1       = [percentiles1(:,2)' percentiles1(end:-1:1,1)'];
yband2       = [percentiles2(:,2)' percentiles2(end:-1:1,1)'];
p1 = fill(xband,yband1,LightBlue,'FaceAlpha',.3,'EdgeColor',LightBlue,'EdgeAlpha',.3); hold on
p2 = fill(xband,yband2,LightRed,'FaceAlpha',.3,'EdgeColor',LightRed,'EdgeAlpha',.3); 
p3 = plot(datetime(datevec(dates_plot)),narrative_inv_base_prc(:,1),'b-','LineWidth',2);
p4 = plot(datetime(datevec(dates_plot)),narrative_inv_scen_prc(:,1),'r:','LineWidth',2);
if add_data == 1
    l = legend([p3,p4],'Baseline','Scenario','Location','southwest','Box','off','Orientation','horizontal','FontName','Times','FontSize',12);
    p5 = plot(datetime(datevec(dates_plot)),stand(idxStand+1:idxStand+h),'k--','linewidth',2,'DisplayName','Data');    
else
    l = legend([p3,p4],'Baseline','Scenario','Location','northwest','Box','off','Orientation','horizontal','FontName','Times','FontSize',12);
end
ax = gca; ax.Box = 'off'; ax.FontSize = 10; ax.YLabel.String = '% of capacity';
ax.Title.String = varNames(4); ax.Title.FontSize = 12; ax.Title.FontName = 'Times';
grid on;  xtickformat('MMM-yyyy');


    