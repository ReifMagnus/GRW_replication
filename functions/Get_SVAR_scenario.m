function [eps_big,innovations,HDinit,HDconst,HDexo,HDtrend,HDseasonal,HDshock,lrm] = Get_SVAR_scenario(Y,exog,A0,c,B,p,h,scenario,time_trend,temperature,dummies)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[T, n] = size(Y);

if isempty(c)
    c = zeros(n,1);
end
B_companion = [B; [eye(n*(p-1)) zeros(n*(p-1),n)]];  % 

if isempty(time_trend) == 0
    trend = 1;
else 
    trend = 0;
end

J = [eye(n); zeros(n*(p-1),n)];

% long-run means
aux =  (eye(n*p) - B_companion)\[c(:,1); zeros(n*(p-1),1)];    
lrm = aux(1:n);

%% Deterministic Components (now including temperature and seasonality)
% pre-allocation
[HDinit_big,HDconst_big,HDexo_big,HDtrend_big,HDseaonal_big] = deal(zeros(p*n,T+h));
[CC,DD,TT,ZZ]                                                = deal(zeros(p*n,T+h));       % holds constants, exogenous and trend

aux = Y(p:-1:1,:)';
y0  = aux(:);    

CC(1:n,:) = repmat(c,1,T+h);
DD(1:n,:) = temperature;
TT(1:n,:) = trend;
ZZ(1:n,:) = dummies;

HDinit_big(:,p) = y0;
for i = p+1:T+h
    % Contribution of initial values
    HDinit_big(:,i) = B_companion*HDinit_big(:,i-1);
    % Contribution of constant and exogenous
    HDconst_big(:,i) = CC(:,i)  + B_companion*HDconst_big(:,i-1);
%     % Contribution of exogenous
    HDexo_big(:,i)   = DD(:,i)  + B_companion*HDexo_big(:,i-1);
%     % Contribution of trend
    HDtrend_big(:,i) = TT(:,i)  + B_companion*HDtrend_big(:,i-1);
%     % Contribution of seasonality
    HDseaonal_big(:,i) = ZZ(:,i)  + B_companion*HDseaonal_big(:,i-1);    
end

HDinit     = HDinit_big(1:n,p:end)';
HDconst    = HDconst_big(1:n,p:end)';
HDexo      = HDexo_big(1:n,p:end)';
HDtrend    = HDtrend_big(1:n,p:end)';
HDseasonal = HDseaonal_big(1:n,p:end)';

% unconditional_means = (eye(n*p) - B_companion)\c_companion;

% yWithLags                    = [Y mlag2(Y,p-1)];
% y0                           = yWithLags(p,:);
% deterministic_component      = nan(T+h,n*p);
% deterministic_component(p,:) = y0;
% 
% for t=p+1:T+h
%     %     deterministic_component(t,:) = (c_companion + B_companion*(deterministic_component(t-1,:))')'; % original deterministic component (intercept only)
%     if trend == 1
%         deterministic_component(t,:) = (c_companion + tt_companion(:,t) + t_companion(:,t) + d_companion(:,t) + B_companion*(deterministic_component(t-1,:))')'; % extended deterministic component (intercept, temperature, dummies)
%     else
%         deterministic_component(t,:) = (c_companion + t_companion(:,t) + d_companion(:,t) + B_companion*(deterministic_component(t-1,:))')'; % extended deterministic component (intercept, temperature, dummies)
%     end
% end

% lr = repmat(unconditional_means(1:n)',T+h,1);
% dc = deterministic_component(:,1:n);

%% Compute the Shocks

Ylag = mlag2(Y,p); % Y is [T x M]. ylag is [T x (Mp)]

if isempty(exog)    
    X     = Ylag(p+1:T,:);
    ALPHA = B';
else  
    X     = [exog(p+1:T,:) Ylag(p+1:T,:)];
    ALPHA = [c';B'];    
end

T = T-p+1;
k = n*p;            % number of lagged endogenous per equation

if trend == 1
    innovations = (Y(p+1:end,:) - X*ALPHA - time_trend(:,p+1:end-h)' - temperature(:,p+1:end-h)' - dummies(:,p+1:end-h)');     % temperature and seasonal effects must also be subtracted from the data
else
    innovations = (Y(p+1:end,:) - X*ALPHA - temperature(:,p+1:end-h)' - dummies(:,p+1:end-h)');
end
shocks = (innovations*(A0));

%% Contribution of each shock

invA0_big        = zeros(k,n);
invA0_big(1:n,:) = inv2(A0)';
Icomp            = J';
HDshock_big      = zeros(p*n,T+h,n);
HDshock          = zeros(n,T+h,n);           % variable x periods x shocks

invA0_big   = sparse(invA0_big);
B_companion = sparse(B_companion);
Icomp       = sparse(Icomp);

for shock = 1:n % Loop over shocks    
    eps_big            = zeros(n,T+h);
    eps_big(shock,2:T) = shocks(:,shock)';
    if isempty(scenario) == 0
        eps_big(shock,T+1:T+h) = scenario(:,shock)';
    end
    for t = 2:T+h
        HDshock_big(:,t,shock) = invA0_big*eps_big(:,t) + B_companion*HDshock_big(:,t-1,shock);
        HDshock(:,t,shock)     = Icomp*HDshock_big(:,t,shock);
    end
end

% All decompositions must add up to the original data
HDendo = (HDinit + HDconst + HDexo + HDtrend + HDseasonal + sum(HDshock,3)');

if h == 0
    if max(HDendo(2:end,:) - Y(p+1:end,:),[],'all') > 1e-5
        warning('There seems to be a problem in the HDs')
    end
end
% 
% % end
% %
% HDshock = permute(HDshock,[2 1 3]);
% HDshock = HDshock(:,2:end,:);
