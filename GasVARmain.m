%% Code for Structuctural VAR %%
% Supports Sign, Zero and Narrative Sign Restrictions, as well as Conjugate Priors on the Reduced Form

clear
close all
clc

addpath('functions')

%% manual settings

% Model Specification
constant           = 1;                           % add constant in VAR
trend              = 0;                           % add a linear time trend
endo               = 1:4;                         % specify Position of endogenous Variables
detrendIP          = 0;                           % detrend IP or not
GLP                = 1;                           % optimize hyperparameters following Giannone et al. and Lenza/Primiceri
lagged_temperature = 1;

p = 12;                                 % lag order of VAR
h = 12;                                 % Desired forecast horizon after end of sample

bands  = [16,50,84]; 

cumulateWhich = [1,2,3]; % Compute Cumulated IRFs and HDs for plot

varNames   = {'Net Gas Imports','Industrial Production','Real Gas Price','Gas Inventories'};        % enter those by-hand to match names in shock definitions

%% Read and transform data
load data_GRW.mat

Tyears     = floor(length(dates)/12);
Tmonths    = length(dates) - Tyears*12;
sa_dummies = [repmat([eye(11); zeros(1,11)],Tyears,1); eye(Tmonths) zeros(Tmonths,11-Tmonths)];     % dummies for seasonal adjustment
 
Year_start  = 1999;                       % Choose start year
Year_end    = 2022;
month_start = 2;
month_end   = 12;

dates_vec = datevec(dates);

ix_first = find(dates_vec(:,1) == Year_start & dates_vec(:,2) == month_start);
ix_last  = find(dates_vec(:,1) == Year_end & dates_vec(:,2) == month_end);

if lagged_temperature == 1
    p_temp      = 11;
    temperature = [temp1998; data(:,end)];
    for l = 1:11
        temp_lag(:,l) = temperature(12+1-l:end-l);
    end
    data = [data temp_lag sa_dummies];
else 
    p_temp = 0;
    data   = [data sa_dummies];
end

data  = data(ix_first:ix_last,:);
dates = dates(ix_first:ix_last,:);

exo = length(endo) + 1:size(data,2);

%% Reduced Form Priors

prior_settings.prior_family = 'conjugate';
prior_settings.prior        = 'Minnesota';                % Select 'flat' or 'Minnesota'

%% Structural Identification Settings

StructuralIdentification = 'Signs/Zeros';           % Chose 'None' or 'Signs/Zeros' or 'Choleski';
agnostic                 = 'irfs';                  % select: 'structural' or 'irfs';

% Set up Sign Restrictions
% SR{r} = {Shockname,{VariableNames}, Horizon, Sign (1 or -1),};
SR{1} = {'Flow Supply',{'Net Gas Imports','Industrial Production'}                      ,0,  -1};
SR{2} = {'Flow Supply',{'Real Gas Price'}                                               ,0,   1};
SR{3} = {'Flow Demand',{'Net Gas Imports','Industrial Production','Real Gas Price'}    ,0,   1};
SR{4} = {'Storage Demand',{'Net Gas Imports','Real Gas Price','Gas Inventories'}       ,0,   1};
SR{5} = {'Storage Demand',{'Industrial Production'}                                    ,0,  -1};
SR{6} = {'Gas Preference',{'Net Gas Imports','Real Gas Price'}                         ,0,   1};
SR{7} = {'Gas Preference',{'Industrial Production','Gas Inventories'}                  ,0,  -1};
 
% Set up NarrativeSign Restrictions
% NSR{r} = {'shockname',typeof restriction ('sign of shock' or 'contribution'), date, end date (for contributions only), variable
% (for contributionsonly), sign, 'strong' or 'weak' for contributions.
 
NSR{1} = {'Flow Supply','sign of shock',datenum(2009,01,01),1};
NSR{2} = {'Flow Supply','sign of shock',datenum(2022,06,01),1};          % Achtung: Dies muss ein positiver "1" Schock sein, da dieser so definiert ist, dass Supply und IP sinken und der Preis steigt
NSR{3} = {'Flow Supply','sign of shock',datenum(2022,07,01),1};
NSR{4} = {'Flow Supply','contribution',datenum(2022,06,01),datenum(2022,07,01),{'Net Gas Imports'},1,'strong'};
NSR{5} = {'Flow Demand','contribution',datenum(2020,04,01),datenum(2020,04,01),{'Industrial Production'},1,'strong'};

checkVarNames

shockNames = (unique(cellfun(@(v) v(1), SR(1,:)),'stable'));
nshocks    = length(shockNames);

if nshocks < length(endo) 
    shockNames_for_HD = [shockNames {'Residual'}];
else
    shockNames_for_HD = shockNames;
end

%% Gibss Sampler Settings
numDesiredDraws  = 5000;
BetaSigmaTries   = 100;
Qs_per_BetaSigma = 100;
nRepsWeights     = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Code that is Edited Frequently %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bayesian_settings;                      

rng('default');
rng(1234);                        % set seed to ensure replicability

tvp         = 0;                  % no time-varying coefficients
computeIRFs = 1;                  % =0 -> only estimation, =1 -> structural identification
estimate    = 1;                  % 0: only load results, 1: estimate
startYear   = Year_start;
endYear     = Year_end;

Run_SVAR

%% Some Results
save_figs = 0;

%% Impulse Responses
Plot_IRFs

%% Forecast Error Variance Decompositions
Plot_FEVDs

%% Historical Decompositions 
PlotHD

%% Scenario analysis

%%  temperature scenarios
dates_f(1) = []; % Kill January 1999
startDate  = datenum(2022,03,01);
endDate    = datenum(2022,03,01);
month_scen = month(endDate);
Scenario   = getScenario_boot('temp',h,n,exog,trend,month_end,month_scen);

%% Effect of temperature on CFs of gas-market variables
PlotHDTemp_bands

if save_figs == 1 
    exportgraphics(f,'Temp_dlIP_tlag.pdf')
end

%% Russian gas embargo starting in April 2022
startDate            = datenum(2022,03,01);
endDate              = datenum(2022,03,01);
month_scen           = month(endDate);
Scenario             = getScenario('shocks',h,n,exog,trend,month_end,month_scen);
add_data             = 1;
Scenario.shocks(1,1) = 6;                       % negative gas supply shock in first OOS period

PlotHDShocks_bands

if save_figs == 1
    exportgraphics(f,'Russia_dlIP_tlag.pdf')
end

%% Norway pipeline breaks down in December 2022
startDate            = datenum(2022,12,01);
endDate              = datenum(2022,12,01);
month_scen           = month(endDate);
Scenario             = getScenario('shocks',h,n,exog,trend,month_end,month_scen);
Scenario.shocks(1,1) = 4;                         % negative gas supply shock in first OOS period
add_data             = 0;

%% Effect of Norway supply disruption on CFs of gas-market variables
PlotHDShocks_bands

if save_figs == 1
    exportgraphics(f,'Norway_dlIP_tlag.pdf')
end

%% Additional figures
FigureGasPrices
FiguresAppendix
