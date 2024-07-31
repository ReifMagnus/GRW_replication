function [b0,Vb0,nu0,S0,Psi,pos,eta0] = get_prior_nc(y,n,p,k,prior,Toutlier,n_exo,dummies,exo)
%% construct normal-inverse Wishart prior
% b0  = prior mean of VAR coefficients
% Vb0 = prior variance of VAR coefficents
% nu0 = prior scale for covariance
% S0  = prior shape for covariance
% c1  = hyperparameter for lag coefficients
% c2  = hyperparameter for intercept coefficients
% pos = position of variables that enter in first differences, and for which the prior mean on first own
%       lag and on sum-of-coefficients prior mean of are set to 0
Y = y;

pt = n_exo - size(dummies,2) - 1;

if ~isempty(Toutlier)
    y = y(1:Toutlier-1,:);
    T = Toutlier - 1;
end

% std. deviation of AR(p)'s 
sig2 = zeros(n,1);
for i = 1:n
    [~,sig2(i)] = ar_p(y(:,i),p,1,dummies(1:T,:));
end

% AR(1) on the variables to get persistence
delta = zeros(n,1);
for i = 1:n
    [beta(i,:),~] = ar_p(y(:,i),1,0,dummies(1:T,:));
    if beta(i,1) >= 0.9999
        delta(i) = 1;
    else
        delta(i) = 0;
    end  
end


% get position of white-noise priors (for sum-of-coefficients prior)
if sum(delta == 0) ~=0
    pos = find(delta == 0)';        
else
    pos = [];
end

b0                    = zeros(k,n);
b0(n_exo+1:n_exo+n,:) = diag(delta);

Vb0           = zeros(1,k);
Vb0(1:n_exo)  = 10e3;
Vb0(1+(1:pt)) = 10./((1:pt).^2*std(exo(:,2)));

for i = 1:k - n_exo
    if i <= n
        l = 1;
    else
        l = ceil((i)/n);
    end
    j = mod(i+n,n); % variable index
    if j == 0
        j = n;
    end
    Vb0(i + n_exo) = prior.lambda^2/(l^prior.alpha*sig2(j));
end

nu0 = n + 2;  % prior scale of IW distribution
S0  = diag(sig2);
Psi = sig2;


%% Outlier part

aux = mean(abs(Y(Toutlier:max([Toutlier+1 T]),:) - Y(Toutlier - 1:max([Toutlier+1 T])-1,:)),2)./mean(abs(Y(2:Toutlier-1,:) - Y(1:Toutlier-2,:)),'all');
if isempty(aux)
    eta0 = [];
elseif length(aux) == 2
    eta0 = [aux; aux(1); .7];               % volatility hyperparameters (taken from Lenza & Primiceri 2020)
elseif length(aux) >= 3
    eta0 = [aux(1:3); .7];                  % volatility hyperparameters 
end


