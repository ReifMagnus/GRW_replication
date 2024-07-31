function [logML,betahat,sigmahat]=logMLVAR_formin(par,y,x,lags,T,n,k,b,d,MIN,MAX,psi,PSI,pos,prior,hyperpriors,priorcoef,Toutlier,dummies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the log-posterior (or the logML if hyperpriors=0), the posterior mode of the coefficients and the covariance matrix of the 
% residuals of the BVAR of Giannone, Lenza and Primiceri (2015), augmented  with a change in volatility at the time of Covid (March 2020).
% Last modified: 06/02/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = MIN.lambda1+(MAX.lambda1 - MIN.lambda1)/(1 + exp(-par(1)));
eta    = MIN.eta'+(MAX.eta'-MIN.eta')./(1 + exp(-par(2:end)));

invweights             = ones(T,1);                                      % vector of s_t
invweights(Toutlier)   = eta(1);
invweights(Toutlier+1) = eta(2);
if T > Toutlier+1
    invweights(Toutlier+2:T) = 1+(eta(3)-1)*eta(4).^[0:T-Toutlier-2];
end
y = diag(1./invweights)*y;                                             % remove heteroscedasticity
x = diag(1./invweights)*x;                                             % remove heteroscedasticity

    
%% setting up the priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_exo    = k-n*lags;
pt       = n_exo - size(dummies,2) - 1;

omega    = zeros(k,1);
if n_exo ~= 0
    omega(1:n_exo) = prior.var_exogenous;
end

omega(1+(1:pt)) = prior.lags_temp./((1:pt).^2*std(x(:,2)));


for i = 1:lags
    omega(n_exo + (i-1)*n + 1:n_exo + i*n) = (lambda^2)./((i^2))./psi;            % prior covariance
end

%% output
XX = x'*x;
XY = x'*y;
iPsi = diag(1./sqrt(psi));

betahat  = (XX+diag(1./omega))\(XY+diag(1./omega)*b);                                    % posterior mode of the VAR coefficients
epshat   = y-x*betahat;                                                                   % VAR residuals
sigmahat = (epshat'*epshat + PSI + (betahat-b)'*diag(1./omega)*(betahat-b))/(T+d+n+1);    % Posterior mode of the covariance matrix

% -----------------------------------------------------
% logML
aaa = diag(sqrt(omega))*XX*diag(sqrt(omega));
bbb = iPsi*(epshat'*epshat + (betahat-b)'*diag(1./omega)*(betahat-b))*iPsi;

eigaaa = real(eig(aaa)); eigaaa(eigaaa<1e-12) = 0; eigaaa = eigaaa+1;
eigbbb = real(eig(bbb)); eigbbb(eigbbb<1e-12) = 0; eigbbb = eigbbb+1;

logML = - n*T*log(pi)/2 + sum(gammaln((T+d-[0:n-1])/2)-gammaln((d-[0:n-1])/2)) - T*sum(log(psi))/2 - n*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;

% log determinant to account for the re-weighting
if ~isempty(Toutlier)
    logML = logML - n*sum(log(invweights));
end

if hyperpriors==1
    logML = logML+logGammapdf(lambda,priorcoef.lambda1.k,priorcoef.lambda1.theta);
    if ~isempty(Toutlier)
         logML = logML - 2*log(eta(1)) - 2*log(eta(2)) - 2*log(eta(3)) + logBetapdf(eta(4),priorcoef.eta4.alpha,priorcoef.eta4.beta);
    end
end
logML = -logML;


function r=logGammapdf(x,k,theta)
r=(k-1)*log(x)-x/theta-k*log(theta)-gammaln(k);

function r=logBetapdf(x,al,bet)
r=(al-1)*log(x)+(bet-1)*log(1-x)-betaln(al,bet);

function r=logIG2pdf(x,alpha,beta)
r=alpha*log(beta)-(alpha+1)*log(x)-beta./x-gammaln(alpha);