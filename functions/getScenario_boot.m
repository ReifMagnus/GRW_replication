function Scenario = getScenario_boot(type,h,n,exog,trend,month_end,month_scen)

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

all_temp = reshape(exog(month_scen:end-12+month_scen,2+trend),12,size(exog(month_scen:end-12+month_scen,2+trend),1)/12);      % months x years starting with month_scen

if h > 0
    temp_scen(1:h,1) = all_temp(:,mean(all_temp,1,'omitnan')==median(mean(all_temp,1,'omitnan')));    % sample MEDIAN year of average monthly temperatures during OOS period
    temp_scen(1:h,2) = all_temp(:,mean(all_temp,1,'omitnan')==max(mean(all_temp,1,'omitnan')));     % sample WARMEST year of average monthly temperatures during OOS period
    temp_scen(1:h,3) = all_temp(:,mean(all_temp,1,'omitnan')==min(mean(all_temp,1,'omitnan')));     % sample COLDEST year of average monthly temperatures during OOS period   
end

%% output
Scenario.temp        = temp_scen;
Scenario.shocks      = scenario;
Scenario.temperature = temp_scen(:,1);
Scenario.dummies     = dummies;
Scenario.type        = type;