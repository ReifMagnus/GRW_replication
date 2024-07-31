function [Beta_draw,Sigma_draw,maxRoot,count_unstable] = draw_beta_sigma_GLP(par,y,x,p,T,n,k,b,d,psi,Toutlier,check_stable,count_unstable,maxtrys,dummies)
%% This function draws Beta and Sigma given a draw of hyperparameters


%% hyperparameters
lambda1 = par(1);
eta     = par(2:end);

invweights             = ones(T,1);
invweights(Toutlier)   = eta(1);
invweights(Toutlier+1) = eta(2);
if T > Toutlier+1
    invweights(Toutlier+2:T)=1+(eta(3)-1)*eta(4).^(0:T-Toutlier-2);
end
y = diag(1./invweights)*y;
x = diag(1./invweights)*x;


n_exo = k - n*p;
pt    = n_exo - size(dummies,2) - 1;

%% priors
omega    = zeros(k,1);
if n_exo ~= 0
    omega(1:n_exo) = 10E3;
end

omega(1+(1:pt)) = 1./((1:pt).^2*std(x(:,2)));

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


% takes a draw from the posterior of SIGMA and beta, if draw is on
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
