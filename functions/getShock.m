function [temp_diff,temp_max,temp_min,temp_mean] = getShock(h,exog,trend,month_end,month_scen)
%% define temperature shock for impulse response functions (IRFs)
% input:
% h    : horizon
% n    : variables

%% housekeeping
release   = version('-release');

%% define a temperature shock: monthly average temperature is 1 STD higher
%  in h future months. The coefficient on monthly average temperature for a
%  certain draw is given by Beta_save(2,:,draw).

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

years  = floor(h/12);
months = rem(h,12);

if h > 0
    if years > 0
        temp_diff = [std_temp(month_scen+1:min(month_scen-1+h,12),1);...
            repmat(std_temp,years,1); std_temp(1:min(max(h-(12-month_scen),0),12))];            % sample STD monthly temperatures over IRF horizon
        temp_max = [max_temp(month_scen+1:min(month_scen-1+h,12),1);...
            repmat(max_temp,years,1); max_temp(1:min(max(h-(12-month_scen),0),12))]...          % sample MAX monthly temperatures over IRF period
            - [mean_temp(month_scen+1:min(month_scen-1+h,12),1);...
            repmat(mean_temp,years,1); mean_temp(1:min(max(h-(12-month_scen),0),12))];
        temp_min =[min_temp(month_scen+1:min(month_scen-1+h,12),1);...
            repmat(min_temp,years,1); min_temp(1:min(max(h-(12-month_scen),0),12))]...          % sample MIN monthly temperatures over IRF period
            - [mean_temp(month_scen+1:min(month_scen-1+h,12),1);...
            repmat(mean_temp,years,1); mean_temp(1:min(max(h-(12-month_scen),0),12))];
        temp_mean = [mean_temp(month_scen+1:min(month_scen-1+h,12),1);...
            repmat(mean_temp,years,1); mean_temp(1:min(max(h-(12-month_scen),0),12))];
    else
        temp_diff = std_temp([month_scen+1:min(month_scen-1+h,12),1:min(max(h-(12-month_scen),0),12)],1);    % sample STD monthly temperatures over IRF horizon
        temp_max = max_temp([month_scen+1:min(month_scen-1+h,12),1:min(max(h-(12-month_scen),0),12)],1)...   % sample MAX monthly temperatures over IRF period
            - mean_temp([month_scen+1:min(month_scen-1+h,12),1:min(max(h-(12-month_scen),0),12)],1);
        temp_min = min_temp([month_scen+1:min(month_scen-1+h,12),1:min(max(h-(12-month_scen),0),12)],1)...   % sample MIN monthly temperatures over IRF period
            - mean_temp([month_scen+1:min(month_scen-1+h,12),1:min(max(h-(12-month_scen),0),12)],1);
        temp_mean = mean_temp([month_scen+1:min(month_scen-1+h,12),1:min(max(h-(12-month_scen),0),12)],1);
    end
end