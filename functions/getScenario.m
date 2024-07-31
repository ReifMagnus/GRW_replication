function Scenario = getScenario(type,h,n,exog,trend,month_end,month_scen)
%% define scenario for historical decomposition
% input:
% type : either shock or temp_warm,temp_mean,temp_cold
% h    : horizon
% n    : variables

%% housekeeping
release   = version('-release');


%% define dummies
years  = floor(h/12);
months = rem(h,12);

dummies     = [repmat([zeros(1,11); eye(11)],years,1); zeros(months,1) eye(months) zeros(months,10-months)];    % auxiliary variable for monthly dummies during forecast horizon
dummies     = [dummies(month_scen+1:end,:); dummies(1:month_scen,:)];                             % start monthly dummies in month when the sample ends

%% define a forecast scenario in terms of structural shock realizations in
% the future, i.e. x months after August 2022, h denotes forecast horizon

scenario        = zeros(h,n);   % 0 = NO structural gas supply and demand shocks
scenario(1:h,1) = 0;            % negative gas supply shock of ~1 SD 
scenario(1:h,2) = 0;            % positive aggregate demand shock of ~1 SD 
scenario(1:h,3) = 0;            % positive other gas demand shock of ~1 SD 

%% define a forecast scenario in terms of the temperature in the upcoming
% h months. The coefficient on monthly average temperature for a certain
% draw is given by Beta_save(2,:,draw).

all_temp = reshape([NaN; exog(:,2+trend); NaN(12-month_end,1)],12,24);      % months x years
if str2double(release(1:4)) < 2021
    mean_temp = nanmean(all_temp,2);
    std_temp  = nanstd(all_temp,2);
    max_temp  = nanmax(all_temp,2);
    min_temp  = -nanmax(-all_temp,2);
else
    mean_temp = mean(all_temp,2,'omitnan');
    std_temp  = std(all_temp,[],2,'omitnan');
    max_temp  = max(all_temp,[],2,'omitnan');
    min_temp  = -max(-all_temp,[],2,'omitnan');
end

if h > 0
    temp_scen(1:h,1) = mean_temp([month_scen+1:min(month_scen-1+h,12),1:min(max(h-(12-month_scen),0),12)],1);    % sample MEAN monthly temperatures during OOS period
    temp_scen(1:h,2) = max_temp([month_scen+1:min(month_scen-1+h,12),1:min(max(h-(12-month_scen),0),12)],1);     % sample MAX monthly temperatures during OOS period
    temp_scen(1:h,3) = min_temp([month_scen+1:min(month_scen-1+h,12),1:min(max(h-(12-month_scen),0),12)],1);     % sample MIN monthly temperatures during OOS period   
end

%% output
Scenario.temp        = temp_scen;
Scenario.shocks      = scenario;
Scenario.temperature = temp_scen(:,1);
Scenario.dummies     = dummies;
Scenario.type        = type;