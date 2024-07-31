%% HD with Bands PointWise and Bayesian Estimate

%% Compute HDs for sample MEAN temperature scenario (NO structural shocks)

clear Xtemp

% load data_speicher_Okt2021.mat

temp_scen     = Scenario.temp;
dummies       = Scenario.dummies;
scenario      = Scenario.shocks;

baseline = zeros(size(scenario));
idxStart = find(ismember(dates,startDate));
idxEnd   = find(ismember(dates,endDate));
idxStand = find(ismember(dates_stand,startDate));

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

% contruct matrix of lagged true temperature
clear temp_true_lags
temp_true      = [exog(idxEnd+1:end,2+trend); temp2023(1:3,:)];
temp_true_lags = lagmatrix(temp_true,0:p_temp);
for j = 2:p_temp+1
    temp_true_lags(1:j-1,j) =  exog(idxEnd-j+2:idxEnd,2);
end

Ytemp = y(idxStart-p:idxEnd,:);

if trend == 1
    Xtemp(:,:,1) = [exog(idxStart-p:idxEnd,:); ones(h,1), (exog(idxEnd,2)+1:exog(idxEnd,2)+h)', temp_scen_lags(:,:,1), dummies]; % Append MEAN temperature scenario and dummies to realized exogenous variables
    Xtemp(:,:,2) = [exog(idxStart-p:idxEnd,:); ones(h,1), (exog(idxEnd,2)+1:exog(idxEnd,2)+h)', temp_scen_lags(:,:,2), dummies]; % Append WARM temperature scenario and dummies to realized exogenous variables
    Xtemp(:,:,3) = [exog(idxStart-p:idxEnd,:); ones(h,1), (exog(idxEnd,2)+1:exog(idxEnd,2)+h)', temp_scen_lags(:,:,3), dummies]; % Append COLD temperature scenario and dummies to realized exogenous variables
    Xtemp(:,:,4) = [exog(idxStart-p:idxEnd,:); ones(h,1), (exog(idxEnd,2)+1:exog(idxEnd,2)+h)', [exog(idxEnd+1:end,2+trend); temperature(1:3,1)], dummies];    % Append TRUE temperature scenario and dummies to realized exogenous variables
else
    Xtemp(:,:,1) = [exog(idxStart-p:idxEnd,:); ones(h,1), temp_scen_lags(:,:,1), dummies]; % Append MEAN temperature scenario and dummies to realized exogenous variables
    Xtemp(:,:,2) = [exog(idxStart-p:idxEnd,:); ones(h,1), temp_scen_lags(:,:,2), dummies]; % Append WARM temperature scenario and dummies to realized exogenous variables
    Xtemp(:,:,3) = [exog(idxStart-p:idxEnd,:); ones(h,1), temp_scen_lags(:,:,3), dummies]; % Append COLD temperature scenario and dummies to realized exogenous variables
    Xtemp(:,:,4) = [exog(idxStart-p:idxEnd,:); ones(h,1), temp_true_lags, dummies];    % Append TRUE temperature scenario and dummies to realized exogenous variables
end
Ta    = length(Ytemp);
exoga = ones(Ta,1);

%% baseline and mean temperature (for conventional SRs)

numSavedDraws = size(A0_save,3);

%% pre-allocation

[draws_HDs_scen_avg,draws_HDs_scen_warm,draws_HDs_scen_cold,draws_HDs_scen_true]                                          = deal(nan(Ta-p+1+h,n,n,numSavedDraws)); % time, shock, draws
[draws_dc_scen_avg,draws_dc_scen_warm,draws_dc_scen_cold,draws_dc_scen_true]                                              = deal(nan(Ta-p+1+h,n,numSavedDraws)); % time, shock, draws
[temp_path_avg,temp_path_warm,temp_path_cold,temp_path_true]                                                              = deal(nan(Ta+1,n,numSavedDraws));
[supply_avg,price_avg,inv_avg,supply_cold,price_cold,inv_cold,supply_true,price_true,inv_true,gas_avg,gas_cold,gas_true]  = deal(nan(h,1,numSavedDraws));
[gas_avg_prc,gas_cold_prc,gas_true_prc,inv_avg_prc,inv_cold_prc,inv_true_prc,price_avg_prc,price_cold_prc,price_true_prc] = deal(nan(h,3));
[gas_avg_BE,gas_cold_BE,gas_true_BE,inv_avg_BE,inv_cold_BE,inv_true_BE,price_avg_BE,price_cold_BE,price_true_BE]          = deal(nan(Ta,1));

for draw = 1:numSavedDraws
    if mod(draw,100) == 0
        fprintf('Running draw %d\n',draw);
    end

    %% get coefficients (signs)
    A0    = A0_save(:,:,draw);
    phi   = Beta_save(size(exog,2)+1:end,:,draw)';
    delta = Beta_save(1,:,draw)';
    if trend == 1
        time_trend = Beta_save(2,:,draw)'*Xtemp(:,2,1)';
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

    % average temperature scenario
    [~,~,ic,dc,ec,tc,sc,shocks,~]   = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,scenario,time_trend,temperature(:,:,1),Dummies(:,:,1));
    draws_HDs_scen_avg(:,:,:,draw)  = permute(shocks,[2,1,3]);
    draws_dc_scen_avg(:,:,draw)     = dc + ic + ec + tc + sc;

    % warm temperature scenario
    [~,~,ic,dc,ec,tc,sc,shocks,~]   = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,scenario,time_trend,temperature(:,:,2),Dummies(:,:,2));
    draws_HDs_scen_warm(:,:,:,draw) = permute(shocks,[2,1,3]);
    draws_dc_scen_warm(:,:,draw)    = dc + ic + ec + tc + sc;

    % cold temperature scenario
    [~,~,ic,dc,ec,tc,sc,shocks,~]   = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,scenario,time_trend,temperature(:,:,3),Dummies(:,:,3));
    draws_HDs_scen_cold(:,:,:,draw) = permute(shocks,[2,1,3]);
    draws_dc_scen_cold(:,:,draw)    = dc + ic + ec + tc + sc;

    % true temperature scenario
    [~,~,ic,dc,ec,tc,sc,shocks,~]   = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,scenario,time_trend,temperature(:,:,4),Dummies(:,:,4));
    draws_HDs_scen_true(:,:,:,draw) = permute(shocks,[2,1,3]);
    draws_dc_scen_true(:,:,draw)    = dc + ic + ec + tc + sc;

    %% Compute time paths of transformed variables
    allshocks_avg               = squeeze(sum(draws_HDs_scen_avg(:,:,:,draw),3));   % stochastic component (unconditional)
    deterministic_avg           = draws_dc_scen_avg(:,:,draw);                      % deterministic component (including temperature!)
    temp_path_avg(:,:,draw)     = allshocks_avg + deterministic_avg;

    allshocks_warm              = squeeze(sum(draws_HDs_scen_warm(:,:,:,draw),3));  % stochastic component (warm scenario)
    deterministic_warm          = draws_dc_scen_warm(:,:,draw);                     % deterministic component (including temperature!)
    temp_path_warm(:,:,draw)    = allshocks_warm + deterministic_warm;

    allshocks_cold              = squeeze(sum(draws_HDs_scen_cold(:,:,:,draw),3));  % stochastic component (cold scenario)
    deterministic_cold          = draws_dc_scen_cold(:,:,draw);                     % deterministic component (including temperature!)
    temp_path_cold(:,:,draw)    = allshocks_cold + deterministic_cold;

    allshocks_true              = squeeze(sum(draws_HDs_scen_true(:,:,:,draw),3));  % stochastic component (true scenario)
    deterministic_true          = draws_dc_scen_true(:,:,draw);                     % deterministic component (including temperature!)
    temp_path_true(:,:,draw)    = allshocks_true + deterministic_true;

    %% Compute time paths of original variables
    supply_avg(:,:,draw)  = temp_path_avg(2:end-1,1,draw)/100;  % timing t-1 for h = 1,...,12
    price_avg(:,:,draw)   = temp_path_avg(3:end,3,draw);    % timing t for h = 1,...,12
    inv_avg(:,:,draw)     = temp_path_avg(3:end,4,draw);    % timing t for h = 1,...,12

    supply_cold(:,:,draw)  = temp_path_cold(2:end-1,1,draw)/100;  % timing t-1 for h = 1,...,12
    price_cold(:,:,draw)   = temp_path_cold(3:end,3,draw);    % timing t for h = 1,...,12
    inv_cold(:,:,draw)     = temp_path_cold(3:end,4,draw);    % timing t for h = 1,...,12

    supply_true(:,:,draw) = temp_path_true(2:end-1,1,draw)/100;  % timing t-1 for h = 1,...,12
    price_true(:,:,draw)  = temp_path_true(3:end,3,draw);    % timing t for h = 1,...,12
    inv_true(:,:,draw)    = temp_path_true(3:end,4,draw);  % timing t for h = 1,...,12

    gas_avg(1,:,draw)  = ones(size(supply_avg(1,:,draw)))*netgas_level(idxEnd,1);
    gas_cold(1,:,draw) = ones(size(supply_cold(1,:,draw)))*netgas_level(idxEnd,1);
    gas_true(1,:,draw) = ones(size(supply_true(1,:,draw)))*netgas_level(idxEnd,1);
    for i = 2:h
        gas_avg(i,:,draw)  = (1 + supply_avg(i,:,draw)).*gas_avg(i-1,:,draw);
        gas_cold(i,:,draw) = (1 + supply_cold(i,:,draw)).*gas_cold(i-1,:,draw);
        gas_true(i,:,draw) = (1 + supply_true(i,:,draw)).*gas_true(i-1,:,draw);
    end
    gas_avg(h+1,:,draw)  = (1 + supply_avg(h,:,draw)).*gas_avg(h-1,:,draw);
    gas_cold(h+1,:,draw) = (1 + supply_cold(h,:,draw)).*gas_cold(h-1,:,draw);
    gas_true(h+1,:,draw) = (1 + supply_true(h,:,draw)).*gas_true(h-1,:,draw);
    inv_avg(:,:,draw)    = stand(idxStand) + cumsum(inv_avg(:,:,draw),1)*277.778*1000*100/230000000;    % Transform to capacity (230 mio. MWh in 2020)
    inv_cold(:,:,draw)   = stand(idxStand) + cumsum(inv_cold(:,:,draw),1)*277.778*1000*100/230000000;
    inv_true(:,:,draw)   = stand(idxStand) + cumsum(inv_true(:,:,draw),1)*277.778*1000*100/230000000;
    price_avg(:,:,draw)  = cumsum(price_avg(:,:,draw));
    price_cold(:,:,draw) = cumsum(price_cold(:,:,draw));
    price_true(:,:,draw) = cumsum(price_true(:,:,draw));

end

for i = 1:h
    gas_avg_prc(i,1) = prctile(gas_avg(i,1,:),50);
    gas_avg_prc(i,2) = prctile(gas_avg(i,1,:),84);
    gas_avg_prc(i,3) = prctile(gas_avg(i,1,:),16);
    gas_avg_BE(i,1)  = gas_avg(i,1,I);

    gas_avg_prc(i,1) = prctile(gas_cold(i,1,:),50);
    gas_cold_prc(i,2) = prctile(gas_cold(i,1,:),84);
    gas_cold_prc(i,3) = prctile(gas_cold(i,1,:),16);
    gas_cold_BE(i,1)  = gas_cold(i,1,I);

    gas_true_prc(i,1) = prctile(gas_true(i,1,:),50);
    gas_true_prc(i,2) = prctile(gas_true(i,1,:),84);
    gas_true_prc(i,3) = prctile(gas_true(i,1,:),16);
    gas_true_BE(i,1)  = gas_true(i,1,I);

    inv_avg_prc(i,1) = prctile(inv_avg(i,1,:),50);
    inv_avg_prc(i,2) = prctile(inv_avg(i,1,:),84);
    inv_avg_prc(i,3) = prctile(inv_avg(i,1,:),16);
    inv_avg_BE(i,1)  = inv_avg(i,1,I);

    inv_cold_prc(i,1) = prctile(inv_cold(i,1,:),50);
    inv_cold_prc(i,2) = prctile(inv_cold(i,1,:),84);
    inv_cold_prc(i,3) = prctile(inv_cold(i,1,:),16);
    inv_cold_BE(i,1)  = inv_cold(i,1,I);

    inv_true_prc(i,1) = prctile(inv_true(i,1,:),50);
    inv_true_prc(i,2) = prctile(inv_true(i,1,:),84);
    inv_true_prc(i,3) = prctile(inv_true(i,1,:),16);
    inv_true_BE(i,1)  = inv_true(i,1,I);

    price_avg_prc(i,1) = prctile(price_avg(i,1,:),50);
    price_avg_prc(i,2) = prctile(price_avg(i,1,:),84);
    price_avg_prc(i,3) = prctile(price_avg(i,1,:),16);
    price_avg_BE(i,1)  = price_avg(i,1,I);

    price_cold_prc(i,1) = prctile(price_cold(i,1,:),50);
    price_cold_prc(i,2) = prctile(price_cold(i,1,:),84);
    price_cold_prc(i,3) = prctile(price_cold(i,1,:),16);
    price_cold_BE(i,1)  = price_cold(i,1,I);

    price_true_prc(i,1) = prctile(price_true(i,1,:),50);
    price_true_prc(i,2) = prctile(price_true(i,1,:),84);
    price_true_prc(i,3) = prctile(price_true(i,1,:),16);
    price_true_BE(i,1)  = price_true(i,1,I);
end
gas_avg_prc(h+1,1) = prctile(gas_avg(h+1,1,:),50);
gas_avg_prc(h+1,2) = prctile(gas_avg(h+1,1,:),84);
gas_avg_prc(h+1,3) = prctile(gas_avg(h+1,1,:),16);
gas_avg_BE(h+1,1)  = gas_avg(h+1,1,I);

gas_cold_prc(h+1,1) = prctile(gas_cold(h+1,1,:),50);
gas_cold_prc(h+1,2) = prctile(gas_cold(h+1,1,:),84);
gas_cold_prc(h+1,3) = prctile(gas_cold(h+1,1,:),16);
gas_cold_BE(h+1,1)  = gas_cold(h+1,1,I);

gas_true_prc(h+1,1) = prctile(gas_true(h+1,1,:),50);
gas_true_prc(h+1,2) = prctile(gas_true(h+1,1,:),84);
gas_true_prc(h+1,3) = prctile(gas_true(h+1,1,:),16);
gas_true_BE(h+1,1)  = gas_true(h+1,1,I);

%% plotting

% Effect of temperature on variable
dates_plot = dates_f(idxStart+1:idxEnd+h);

%% baseline and mean temperature (with narrative SRs)

numSavedDraws = size(A0_narrative,3);   

[draws_HDs_base,draws_HDs_scen_avg,draws_HDs_scen_warm,draws_HDs_scen_cold,draws_HDs_scen_true] = deal(nan(Ta-p+1+h,n,n,numSavedDraws)); % time, shock, draws
[draws_dc_base,draws_dc_scen_avg,draws_dc_scen_warm,draws_dc_scen_cold,draws_dc_scen_true]   = deal(nan(Ta-p+1+h,n,numSavedDraws)); % time, shock, draws

for draw = 1:numSavedDraws
    if mod(draw,100) == 0
        fprintf('Running draw %d\n',draw);
    end

    A0    = A0_narrative(:,:,draw);
    phi   = Beta_narrative(size(exog,2)+1:end,:,draw)';
    delta = Beta_narrative(1,:,draw)';
    if trend == 1
        time_trend = Beta_save(2,:,draw)'*Xtemp(:,2,1)';
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

    % average temperature scenario
    [~,~,ic,dc,ec,tc,sc,shocks,~]   = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,scenario,time_trend,temperature(:,:,1),Dummies(:,:,1));
    draws_HDs_scen_avg(:,:,:,draw)  = permute(shocks,[2,1,3]);
    draws_dc_scen_avg(:,:,draw)     = dc + ic + ec + tc + sc;

    % warm temperature scenario
    [~,~,ic,dc,ec,tc,sc,shocks,~]   = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,scenario,time_trend,temperature(:,:,2),Dummies(:,:,2));
    draws_HDs_scen_warm(:,:,:,draw) = permute(shocks,[2,1,3]);
    draws_dc_scen_warm(:,:,draw)    = dc + ic + ec + tc + sc;

    % cold temperature scenario
    [~,~,ic,dc,ec,tc,sc,shocks,~]   = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,scenario,time_trend,temperature(:,:,3),Dummies(:,:,3));
    draws_HDs_scen_cold(:,:,:,draw) = permute(shocks,[2,1,3]);
    draws_dc_scen_cold(:,:,draw)    = dc + ic + ec + tc + sc;

    % true temperature scenario
    [~,~,ic,dc,ec,tc,sc,shocks,~]   = Get_SVAR_scenario(Ytemp,exoga,A0,delta,phi,p,h,scenario,time_trend,temperature(:,:,4),Dummies(:,:,4));
    draws_HDs_scen_true(:,:,:,draw) = permute(shocks,[2,1,3]);
    draws_dc_scen_true(:,:,draw)    = dc + ic + ec + tc + sc;

    %% Compute time paths of transformed variables
    allshocks_avg               = squeeze(sum(draws_HDs_scen_avg(:,:,:,draw),3));   % stochastic component (unconditional)
    deterministic_avg           = draws_dc_scen_avg(:,:,draw);                      % deterministic component (including temperature!)
    temp_path_avg(:,:,draw)     = allshocks_avg + deterministic_avg;

    allshocks_warm              = squeeze(sum(draws_HDs_scen_warm(:,:,:,draw),3));  % stochastic component (warm scenario)
    deterministic_warm          = draws_dc_scen_warm(:,:,draw);                     % deterministic component (including temperature!)
    temp_path_warm(:,:,draw)    = allshocks_warm + deterministic_warm;

    allshocks_cold              = squeeze(sum(draws_HDs_scen_cold(:,:,:,draw),3));  % stochastic component (cold scenario)
    deterministic_cold          = draws_dc_scen_cold(:,:,draw);                     % deterministic component (including temperature!)
    temp_path_cold(:,:,draw)    = allshocks_cold + deterministic_cold;

    allshocks_true              = squeeze(sum(draws_HDs_scen_true(:,:,:,draw),3));  % stochastic component (true scenario)
    deterministic_true          = draws_dc_scen_true(:,:,draw);                     % deterministic component (including temperature!)
    temp_path_true(:,:,draw)    = allshocks_true + deterministic_true;

    %% Compute time paths of original variables

    supply_avg(:,:,draw)  = temp_path_avg(2:end-1,1,draw)/100;  % timing t-1 for h = 1,...,12
    price_avg(:,:,draw)   = temp_path_avg(3:end,3,draw);    % timing t for h = 1,...,12
    inv_avg(:,:,draw)     = temp_path_avg(3:end,4,draw);    % timing t for h = 1,...,12

    supply_cold(:,:,draw)  = temp_path_cold(2:end-1,1,draw)/100;  % timing t-1 for h = 1,...,12
    price_cold(:,:,draw)   = temp_path_cold(3:end,3,draw);    % timing t for h = 1,...,12
    inv_cold(:,:,draw)     = temp_path_cold(3:end,4,draw);    % timing t for h = 1,...,12

    supply_true(:,:,draw) = temp_path_true(2:end-1,1,draw)/100;  % timing t-1 for h = 1,...,12
    price_true(:,:,draw)  = temp_path_true(3:end,3,draw);    % timing t for h = 1,...,12
    inv_true(:,:,draw)    = temp_path_true(3:end,4,draw);  % timing t for h = 1,...,12

    gas_avg(1,:,draw)  = ones(size(supply_avg(1,:,draw)))*netgas_level(idxEnd,1);
    gas_cold(1,:,draw) = ones(size(supply_cold(1,:,draw)))*netgas_level(idxEnd,1);
    gas_true(1,:,draw) = ones(size(supply_true(1,:,draw)))*netgas_level(idxEnd,1);
    for i = 2:h
        gas_avg(i,:,draw)  = (1 + supply_avg(i,:,draw)).*gas_avg(i-1,:,draw);
        gas_cold(i,:,draw) = (1 + supply_cold(i,:,draw)).*gas_cold(i-1,:,draw);
        gas_true(i,:,draw) = (1 + supply_true(i,:,draw)).*gas_true(i-1,:,draw);
    end
    gas_avg(h+1,:,draw)  = (1 + supply_avg(h,:,draw)).*gas_avg(h-1,:,draw);
    gas_cold(h+1,:,draw) = (1 + supply_cold(h,:,draw)).*gas_cold(h-1,:,draw);
    gas_true(h+1,:,draw) = (1 + supply_true(h,:,draw)).*gas_true(h-1,:,draw);
    inv_avg(:,:,draw)    = stand(idxStand) + cumsum(inv_avg(:,:,draw),1)*277.778*1000*100/230000000;    % Transform to capacity (230 mio. MWh in 2020)
    inv_cold(:,:,draw)   = stand(idxStand) + cumsum(inv_cold(:,:,draw),1)*277.778*1000*100/230000000;
    inv_true(:,:,draw)   = stand(idxStand) + cumsum(inv_true(:,:,draw),1)*277.778*1000*100/230000000;
    price_avg(:,:,draw)  = cumsum(price_avg(:,:,draw));
    price_cold(:,:,draw) = cumsum(price_cold(:,:,draw));
    price_true(:,:,draw) = cumsum(price_true(:,:,draw));

end

for i = 1:h
    gas_avg_prc(i,1) = prctile(gas_avg(i,1,:),50);
    gas_avg_prc(i,2) = prctile(gas_avg(i,1,:),84);
    gas_avg_prc(i,3) = prctile(gas_avg(i,1,:),16);
    gas_avg_BE(i,1)  = gas_avg(i,1,In);

    gas_cold_prc(i,1) = prctile(gas_cold(i,1,:),50);
    gas_cold_prc(i,2) = prctile(gas_cold(i,1,:),84);
    gas_cold_prc(i,3) = prctile(gas_cold(i,1,:),16);
    gas_cold_BE(i,1)  = gas_cold(i,1,In);

    gas_true_prc(i,1) = prctile(gas_true(i,1,:),50);
    gas_true_prc(i,2) = prctile(gas_true(i,1,:),84);
    gas_true_prc(i,3) = prctile(gas_true(i,1,:),16);
    gas_true_BE(i,1)  = gas_true(i,1,In);

    inv_avg_prc(i,1) = prctile(inv_avg(i,1,:),50);
    inv_avg_prc(i,2) = prctile(inv_avg(i,1,:),84);
    inv_avg_prc(i,3) = prctile(inv_avg(i,1,:),16);
    inv_avg_BE(i,1)  = inv_avg(i,1,In);

    inv_cold_prc(i,1) = prctile(inv_cold(i,1,:),50);
    inv_cold_prc(i,2) = prctile(inv_cold(i,1,:),84);
    inv_cold_prc(i,3) = prctile(inv_cold(i,1,:),16);
    inv_cold_BE(i,1)  = inv_cold(i,1,In);

    inv_true_prc(i,1) = prctile(inv_true(i,1,:),50);
    inv_true_prc(i,2) = prctile(inv_true(i,1,:),84);
    inv_true_prc(i,3) = prctile(inv_true(i,1,:),16);
    inv_true_BE(i,1)  = inv_true(i,1,In);

    price_avg_prc(i,1) = prctile(price_avg(i,1,:),50);
    price_avg_prc(i,2) = prctile(price_avg(i,1,:),84);
    price_avg_prc(i,3) = prctile(price_avg(i,1,:),16);
    price_avg_BE(i,1)  = price_avg(i,1,In);

    price_cold_prc(i,1) = prctile(price_cold(i,1,:),50);
    price_cold_prc(i,2) = prctile(price_cold(i,1,:),84);
    price_cold_prc(i,3) = prctile(price_cold(i,1,:),16);
    price_cold_BE(i,1)  = price_cold(i,1,In);

    price_true_prc(i,1) = prctile(price_true(i,1,:),50);
    price_true_prc(i,2) = prctile(price_true(i,1,:),84);
    price_true_prc(i,3) = prctile(price_true(i,1,:),16);
    price_true_BE(i,1)  = price_true(i,1,In);
end
gas_avg_prc(h+1,1) = prctile(gas_avg(h+1,1,:),50);
gas_avg_prc(h+1,2) = prctile(gas_avg(h+1,1,:),84);
gas_avg_prc(h+1,3) = prctile(gas_avg(h+1,1,:),16);
gas_avg_BE(h+1,1)  = gas_avg(h+1,1,In);

gas_cold_prc(h+1,1) = prctile(gas_cold(h+1,1,:),50);
gas_cold_prc(h+1,2) = prctile(gas_cold(h+1,1,:),84);
gas_cold_prc(h+1,3) = prctile(gas_cold(h+1,1,:),16);
gas_cold_BE(h+1,1)  = gas_cold(h+1,1,In);

gas_true_prc(h+1,1) = prctile(gas_true(h+1,1,:),50);
gas_true_prc(h+1,2) = prctile(gas_true(h+1,1,:),84);
gas_true_prc(h+1,3) = prctile(gas_true(h+1,1,:),16);
gas_true_BE(h+1,1)  = gas_true(h+1,1,In);

%% plotting

% Effect of temperature on variable
dates_plot = dates_f(idxStart+1:idxEnd+h);
tt         = datetime(datevec(dates_plot));


ytemp = [y(idxEnd+1:ix_last,:); nan(idxEnd+h-ix_last,n)];
xband = [tt(1):calmonths(1):tt(end) tt(end):-calmonths(1):tt(1)];

f  = figure('Units','normalized','Position',[.1 .1 .5 .3]);
tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile;
plot(datetime(datevec(dates_plot)),Xtemp(end-h+1:end,2,4),'b-','LineWidth',2); hold on
plot(datetime(datevec(dates_plot)),Xtemp(end-h+1:end,2,3),'r:','LineWidth',2);
ax = gca; ax.Box = 'off'; ax.FontSize = 10; ax.YLabel.String = '°C';
ax.Title.String = 'Temperature'; ax.Title.FontSize = 12; ax.Title.FontName = 'Times';
grid on; 

nexttile;
percentiles1 = gas_true_prc(2:end,2:3)/1000;
percentiles2 = gas_cold_prc(2:end,2:3)/1000;
yband1       = [percentiles1(:,2)' percentiles1(end:-1:1,1)'];
yband2       = [percentiles2(:,2)' percentiles2(end:-1:1,1)'];
p1 = fill(xband,yband1,LightBlue,'FaceAlpha',.3,'EdgeColor',LightBlue,'EdgeAlpha',.3); hold on
p2 = fill(xband,yband2,LightRed,'FaceAlpha',.3,'EdgeColor',LightRed,'EdgeAlpha',.3); 
p3 = plot(datetime(datevec(dates_plot)),gas_true_prc(2:end,1)/1000,'b-','LineWidth',2); hold on
p4 = plot(datetime(datevec(dates_plot)),gas_cold_prc(2:end,1)/1000,'r:','LineWidth',2); 
% add data
gas_data(1,1) = (1+ytemp(1,1)/100)*netgas_level(idxEnd,1);
for i = 2:h
    gas_data(i,1) = (1 + ytemp(i,1)/100).*gas_data(i-1,1);
end
p5 = plot(datetime(datevec(dates_plot)),gas_data/1000,'k--','linewidth',2,'DisplayName','Data');
ax = gca; ax.Box = 'off'; ax.FontSize = 10; ax.YLabel.String = '1000 TJ';
ax.Title.String = varNames(1); ax.Title.FontSize = 12; ax.Title.FontName = 'Times';
grid on; 

nexttile;
percentiles1 = price_true_prc(:,2:3);
percentiles2 = price_cold_prc(:,2:3);
yband1       = [percentiles1(:,2)' percentiles1(end:-1:1,1)'];
yband2       = [percentiles2(:,2)' percentiles2(end:-1:1,1)'];
p1 = fill(xband,yband1,LightBlue,'FaceAlpha',.3,'EdgeColor',LightBlue,'EdgeAlpha',.3); hold on
p2 = fill(xband,yband2,LightRed,'FaceAlpha',.3,'EdgeColor',LightRed,'EdgeAlpha',.3); 
p3 = plot(datetime(datevec(dates_plot)),price_true_prc(:,1),'b-','LineWidth',2); hold on
p4 = plot(datetime(datevec(dates_plot)),price_cold_prc(:,1),'r:','LineWidth',2);
p5 = plot(datetime(datevec(dates_plot)),cumsum(ytemp(:,3)),'k--','linewidth',2,'DisplayName','Data');
ax = gca; ax.Box = 'off'; ax.FontSize = 10; ax.YLabel.String = 'cumulated %';
ax.Title.String = varNames(3); ax.Title.FontSize = 12; ax.Title.FontName = 'Times';
grid on; 

nexttile;
percentiles1 = inv_true_prc(:,2:3);
percentiles2 = inv_cold_prc(:,2:3);
yband1       = [percentiles1(:,2)' percentiles1(end:-1:1,1)'];
yband2       = [percentiles2(:,2)' percentiles2(end:-1:1,1)'];
p1 = fill(xband,yband1,LightBlue,'FaceAlpha',.3,'EdgeColor',LightBlue,'EdgeAlpha',.3); hold on
p2 = fill(xband,yband2,LightRed,'FaceAlpha',.3,'EdgeColor',LightRed,'EdgeAlpha',.3); 
p3 = plot(datetime(datevec(dates_plot)),inv_true_prc(:,1),'b-','LineWidth',2); hold on
p4 = plot(datetime(datevec(dates_plot)),inv_cold_prc(:,1),'r:','LineWidth',2); 
l = legend([p3,p4],'Actual','Cold','Location','southwest','Box','off','Orientation','horizontal','FontName','Times','FontSize',12);
p5 = plot(datetime(datevec(dates_plot)),stand(idxStand+1:idxStand+h),'k--','linewidth',2,'DisplayName','Data');
ax = gca; ax.Box = 'off'; ax.FontSize = 10; ax.YLabel.String = '% of capacity';
ax.Title.String = varNames(4); ax.Title.FontSize = 12; ax.Title.FontName = 'Times';
grid on; 
