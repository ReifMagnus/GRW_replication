function [results,MIN,MAX] = find_post_mode(prior,y,X,p,T,n,k,b0,nu0,Psi,S0,eta0,pos,Toutlier,dummies)
%% find posterior mode of log ML of VAR


n_exo = k-n*p;
pt    = n_exo - size(dummies,2) - 1;
%% bounds for maximization
MIN.lambda1  = 0.0001;
MIN.lambda4  = 0.1;
MIN.lambda5  = 0.0001;
MIN.lambda6  = 0.0001;
MAX.lambda1  = 5;
MAX.lambda4  = 5;
MAX.lambda5  = 50;
MAX.lambda6  = 50;
MIN.eta(1:3) = [1;1;1];
MAX.eta(1:3) = [500;500;500];
MIN.eta(4)   = .005;
MAX.eta(4)   = .995;

%% starting values for minimization
lambda1_0 = .2;     % std of MN prior

a = mean(abs(y(Toutlier:max([Toutlier+1 T]),:) - y(Toutlier - 1:max([Toutlier+1 T])-1,:)),2);
b = mean(abs(y(2:Toutlier-1,:) - y(1:Toutlier-2,:)),'all');

aux = a./b;
aux = [2 2 2]';
if isempty(aux)
    eta0 = [];
elseif length(aux) == 2
    eta0 = [aux; aux(1); .7];               % volatility hyperparameters (taken from Lenza & Primiceri 2020)
elseif length(aux) >= 3
    eta0 = [aux(1:3); .7];                  % volatility hyperparameters 
end


in_lambda1  = -log((MAX.lambda1 - lambda1_0)/(lambda1_0 - MIN.lambda1));
inH_lambda1 = (1/(MAX.lambda1 - lambda1_0)+1/(lambda1_0 - MIN.lambda1))^2*(abs(lambda1_0)/1)^2;

ineta  = -log((MAX.eta' - eta0)./(eta0 - MIN.eta'));
inHeta = (1./(MAX.eta' - eta0)+1./(eta0 - MIN.eta')).^2.*(abs(eta0)/1).^2;


x0 = [in_lambda1;ineta];
H0 = diag([inH_lambda1;inHeta]); % initial guess for the inverse Hessian

% check for appropriatness of initial conditions
if ~isreal(x0) 
    error('Initial conditions includes complex number; make sure that data is on the same scale for Lenze-Primiceri approach')
end

% moments of prior distribution (taken from Giannone et al. 2015, ReStat)
mode.lambda1         =.2;     % hyperpriors modes
sd.lambda1           =.4;       % hyperpriors std
r.priorcoef.lambda1  = GammaCoef(mode.lambda1,sd.lambda1,0);  % coefficients of hyperpriors

if ~isempty(Toutlier)
    mode.eta(4) = .8;
    sd.eta(4)   = .2;
    mosd        = [mode.eta(4) sd.eta(4)];
    fun         = @(x) BetaCoef(x,mosd);
    albet       = fsolve(fun,[2,2])';
    r.priorcoef.eta4.alpha = albet(1);    
    r.priorcoef.eta4.beta  = albet(2);
end

r.hyperpriors = 1;


% get posterior mode
[~,post_mode,~,~,~,~,~] = csminwel('logMLVAR_formin',x0,H0,[],1e-16,10000,y,X,p,T,n,k,b0,nu0,MIN,MAX,Psi,S0,pos,prior,...
                                        r.hyperpriors,r.priorcoef,Toutlier,dummies);
% VAR coefficients and residuals at the posterior mode
[fh,r.postmax.betahat,r.postmax.sigmahat] = logMLVAR_formin(post_mode,y,X,p,T,n,k,b0,nu0,MIN,MAX,Psi,S0,pos,prior,...
                                        r.hyperpriors,r.priorcoef,Toutlier,dummies);
                                               
%% output of the maximization

r.lags            = p;                       %  lags
r.postmax.SSar1   = Psi;                     % residual variance of AR(1) for each variable
r.postmax.logPost = -fh;                     % value of the posterior of the hyperparameters at the peak
r.postmax.lambda1 = MIN.lambda1+(MAX.lambda1-MIN.lambda1)/(1 + exp(-post_mode(1)));    % std of MN prior at the peak
r.postmax.eta = MIN.eta'+(MAX.eta'-MIN.eta')./(1 + exp(-post_mode(2:end)));


omega    = zeros(k,1);
if n_exo ~= 0
    omega(1:n_exo) = prior.var_exogenous;
end

omega(1+(1:pt)) = prior.lags_temp./((1:pt).^2*std(X(:,2)));

for i = 1:p
    omega(n_exo + (i-1)*n + 1:n_exo + i*n) = (r.postmax.lambda1^2)./((i^2))./Psi;            % prior covariance
end


invweights           = ones(T,1);                                      % vector of s_t
if ~isempty(Toutlier)    
    invweights(Toutlier)   = r.postmax.eta(1);
    invweights(Toutlier+1) = r.postmax.eta(2);
    if T > Toutlier+1
        invweights(Toutlier+2:T) = 1+(r.postmax.eta(3)-1)*r.postmax.eta(4).^[0:T-Toutlier-2];
    end
end
YY = diag(1./invweights)*y;                                             % remove heteroscedasticity
XX = diag(1./invweights)*X;                                             % remove heteroscedasticity


KA      = sparse(1:k,1:k,1./omega) + XX'*XX;
betahat = KA\(sparse(1:k,1:k,omega)\b0 + XX'*YY);                                  % posterior mean

% check stability of posterior mean
if stability(betahat,n,p,n_exo)  > 0, warning('Posterior mean is not stable'); end


%% compute inverse of hessian at posterior mode
    
hyperpriors = r.hyperpriors;
priorcoef   = r.priorcoef;
r.postmode  = [r.postmax.lambda1; r.postmax.eta];

% new computation of the inverse Hessian
fun = @(par) logMLVAR_formcmc(par,prior,y,X,p,T,n,k,b0,nu0,MIN,MAX,Psi,hyperpriors,priorcoef,Toutlier,0,0,0,0,dummies);
Hess = hessian(fun,r.postmode);
HH   = -inv(Hess);
if ~isempty(Toutlier) && T <= Toutlier+1
    HessNew      = Hess;
    HessNew(4,:) = 0; 
    HessNew(:,4) = 0; 
    HessNew(4,4) = -1;    
    HH     = -inv(HessNew); 
    HH(4,4)= HH(2,2);
end
r.HH = HH;

%% starting value of the Metropolis algorithm
const    = 1;                                                              % scaling factor for acceptance ratio
logMLold = -10e15;
while logMLold == -10e15
    P         = mvnrnd(r.postmode,HH*const^2,1);
    logMLold  = logMLVAR_formcmc(P(1,:)',prior,y,X,p,T,n,k,b0,nu0,MIN,MAX,Psi,hyperpriors,priorcoef,Toutlier,0,0,0,0,dummies);
end
LOGML(1) = logMLold;

%% output
r.P     = P;
r.LogML = LOGML;
results = r;


