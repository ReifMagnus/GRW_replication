function [beta,sigma2] = ar_p(y,p,const,dummies)
%% estimate AR(p)-model

    
xtmp = [];
for j = 1:p
    xtmp = [xtmp lag0(y,j)];
end
xtmp   = [xtmp ones(size(xtmp,1),const) ]; %#ok<AGROW>
y      = y(p+1:end);
xtmp   = xtmp(p+1:end,:);
beta   = (xtmp'*xtmp)\(xtmp'*y);
e      = y - xtmp*beta;
sigma2 = (e'*e)/(size(xtmp,1)-size(xtmp,2));
