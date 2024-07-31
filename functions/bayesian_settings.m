% Additional Settings of the Minnesota Prior 
prior_settings.stationary = [1 4];                  % select stationary variables that get a white-noise prior          
 
prior_settings.hyperparameters.lambda = .2;         % Tightness of the Minnesota Prior; 
prior_settings.hyperparameters.alpha  = 2;          % Decay of the Minnesota Prior; 
prior_settings.var_exogenous          = 10E3;       % Exogenous variables except lags of temperature
prior_settings.lags_temp              = 1;          % lags of temperature                       

prior_settings.NoCointegration  = 0;                % Also known as Sum of Coefficients
prior_settings.SingleUnitRoot   = 0;                % Also  known as Dummy Initial Observation
prior_settings.LongRunPrior     = 0;
prior_settings.LongRunPrior2    = 0;
prior_settings.GrowthRatePrior  = 0;

prior_settings.hyperparameters. mu   = 1;           % Tightness of the No Cointegration / Sum of Coefficients Prior
prior_settings.hyperparameters.theta = 1;           % Tightness of the Dummy Initial Observation Priors

