function [logML,Beta_draw,Sigma_draw,CSig,count_unstable,betahat,maxRoot] = logMLVAR_formcmc(par,prior,y,x,p,T,n,k,b,d,MIN,MAX,psi,hyperpriors,priorcoef,Toutlier,check_stable,count_unstable,maxtrys,draw_coeffs,dummies)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the log-posterior (or the logML if hyperpriors=0),and draws from the posterior distribution of the coefficients and of the 
% covariance matrix of the residuals of the BVAR of Giannone, Lenza and % Primiceri (2015), augmented with a change in volatility at the time of 
% Covid (March 2020).
% Last modified: 06/02/2020

%% hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda1 = par(1);

if ~isempty(Toutlier)
    eta = par(2:end);
    
    invweights           = ones(T,1); 
    invweights(Toutlier)   = eta(1);
    invweights(Toutlier+1) = eta(2);
    if T>Toutlier+1
        invweights(Toutlier+2:T)=1+(eta(3)-1)*eta(4).^[0:T-Toutlier-2];
    end
    y = diag(1./invweights)*y;
    x = diag(1./invweights)*x;
else
    eta   = MIN.eta';
end

n_exo = k - n*p;
pt    = n_exo - size(dummies,2) - 1;


%% return a very low value of the posterior if the parameters are outside the bounds
if sum([lambda1;eta]<[MIN.lambda1;MIN.eta'])>0 || sum([lambda1;eta]>[MAX.lambda1;MAX.eta'])>0
    logML      = -10e15;
    Beta_draw  = [];
    Sigma_draw = [];
    CSig       = [];
    return
else
    
    %% priors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    omega    = zeros(k,1);    
    if n_exo ~= 0
        omega(1:n_exo) = prior.var_exogenous;
    end   

    omega(1+(1:pt)) = prior.lags_temp./((1:pt).^2*std(x(:,2)));

    for i = 1:p
        omega(n_exo + (i-1)*n + 1:n_exo + i*n) = (lambda1^2)./((i^2))./psi;            % prior covariance
    end

    
    % prior scale matrix for the covariance of the shocks
    PSI = diag(psi);
        
    %% output
    XX = x'*x;
    XY = x'*y;
    YY = y'*y;
    
    KA      = sparse(1:k,1:k,1./omega) + XX;
    CK      = chol(KA,'lower')';
    betahat = KA\(sparse(1:k,1:k,omega)\b + XY);                                   % posterior mean
    S       = PSI + b'*sparse(1:k,1:k,1./omega)*b + YY - betahat'*KA*betahat;    % posterior scale of covariance
    S       = (S+S')/2;                                                     % adjust for rounding errors

    epshat  = y - x*betahat;                                                  % VAR residuals

    % -----------------------------------------------------
    % logML    
    aaa = diag(sqrt(omega))*XX*diag(sqrt(omega));
    bbb = diag(1./sqrt(psi))*(epshat'*epshat + (betahat-b)'*diag(1./omega)*(betahat-b))*diag(1./sqrt(psi));
    
    eigaaa = real(eig(aaa)); eigaaa(eigaaa<1e-12) = 0; eigaaa = eigaaa+1;
    eigbbb = real(eig(bbb)); eigbbb(eigbbb<1e-12) = 0; eigbbb = eigbbb+1;
    
    logML = - n*T*log(pi)/2 + sum(gammaln((T+d-(0:n-1))/2)-gammaln((d-[0:n-1])/2)) - T*sum(log(psi))/2 - n*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;
    
     
    % log determinant to account for the re-weighting
    if ~isempty(Toutlier)
        logML = logML - n*sum(log(invweights));
    end
    
    if hyperpriors == 1
        logML = logML + logGammapdf(lambda1,priorcoef.lambda1.k,priorcoef.lambda1.theta);
        if ~isempty(Toutlier)            
            logML = logML - 2*log(eta(1)) - 2*log(eta(2)) - 2*log(eta(3)) + logBetapdf(eta(4),priorcoef.eta4.alpha,priorcoef.eta4.beta);
        end
    end
    
    % takes a draw from the posterior of SIGMA and beta, if draw is on   
    if draw_coeffs == 1
        Sigma_draw = iwishrnd(S,d+T);                                                    % draw of covariance matrix
        CSig       = chol(Sigma_draw,'lower');
        check      = 1;
        trys       = 1;
        if check_stable == 1
            while check > 0 %&& trys <= maxtrys
                Beta_draw       = betahat + (CK\randn(k,n))*CSig';                                     % draw of VAR coefficients
                [check,maxRoot] = stability(Beta_draw,n,p,n_exo);                                          % if equal zero->stability
                trys            = trys + 1;
            end
            if trys == maxtrys
                count_unstable = count_unstable + 1;
                warning('Posterior mean is not stable')          
            end        
        else
            Beta_draw = betahat + (CK\randn(k,n))*CSig';                                         % draw of VAR coefficients
        end     
    end
end

function r=logGammapdf(x,k,theta)
r=(k-1)*log(x)-x/theta-k*log(theta)-gammaln(k);

function r=logBetapdf(x,al,bet)
r=(al-1)*log(x)+(bet-1)*log(1-x)-betaln(al,bet);

function r=logIG2pdf(x,alpha,beta)
r=alpha*log(beta)-(alpha+1)*log(x)-beta./x-gammaln(alpha);