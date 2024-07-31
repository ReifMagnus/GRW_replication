
%% bounds for maximization
MIN.lambda   = 0.0001;
MAX.lambda   = 5;
MIN.eta(1:3) = [1;1;1];
MAX.eta(1:3) = [500;500;500];
MIN.eta(4)   = .005;
MAX.eta(4)   = .995;




%% starting values for the minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda0 =.2;       % std of minesota prior

aux=mean(abs(y(Tenergycrisis:max([Tenergycrisis+1 T]),:)-y(Tenergycrisis-1:max([Tenergycrisis+1 T])-1,:))',1)./...
    mean(mean(abs(y(2:Tenergycrisis-1,:)-y(1:Tenergycrisis-2,:))));

if isempty(aux)
    eta0=[];
elseif length(aux)==2
    eta0=[aux';aux(1);.8];      % volatility hyperparameters 
elseif length(aux)>=3
    eta0=[aux(1:3)';.8];             % volatility hyperparameters 
end

% residual variance of AR(1) for each variable
SS=zeros(n,1);
for i=1:n
    Tend = T; 
    if ~isempty(Tenergycrisis) 
        Tend=Tenergycrisis-1; 
    end
    ar1   = ols1(y(2:Tend,i),[ones(Tend-1,1),y(1:Tend-1,i)]);
    SS(i) = ar1.sig2hatols;
end

inlambda  = -log((MAX.lambda-lambda0)/(lambda0-MIN.lambda));
inHlambda = (1/(MAX.lambda-lambda0)+1/(lambda0-MIN.lambda))^2*(abs(lambda0)/1)^2;


if ~isempty(Tenergycrisis)
    ncp    = length(eta0);           % # "covid" hyperparameters
    ineta  = -log((MAX.eta'-eta0)./(eta0-MIN.eta'));
    inHeta = (1./(MAX.eta'-eta0)+1./(eta0-MIN.eta')).^2.*(abs(eta0)/1).^2;
else
    ineta  = [];
    inHeta = [];
end

x0=[inlambda;ineta];
H0=diag([inHlambda;inHeta]); % initial guess for the inverse Hessian

mode.lambda      = .2;     % hyperpriors modes
sd.lambda        = .4;       % hyperpriors std
priorcoef.lambda = GammaCoef(mode.lambda,sd.lambda,0);  % coefficients of hyperpriors

mode.eta(4)          =.8;
sd.eta(4)            =.2;
mosd                 = [mode.eta(4) sd.eta(4)];
fun                  = @(x) BetaCoef(x,mosd);
albet                = fsolve(fun,[2,2])';
priorcoef.eta4.alpha = albet(1);
priorcoef.eta4.beta  = albet(2);


pos = 2;

