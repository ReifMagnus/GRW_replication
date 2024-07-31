function [S,ee] = stability(beta,n,p,n_exo)
%% check stability of VAR 
%beta   (1+n*l) x n matrix with VAR coefficients
%p      number of lags
%n      number of endog variables
%F      VAR companion matrix 
%S      dummy var: if equal zero->stability


F                    = zeros(n*p,n*p);
F(n+1:n*p,1:n*(p-1)) = eye(n*(p-1),n*(p-1));

temp         = beta(n_exo+1:end,1:n)';         % cut exogenous variables (be careful they have to be ordered first!)
F(1:n,1:n*p) = temp;

ee = max(abs(eig(F)));
S  = ee > 1;


