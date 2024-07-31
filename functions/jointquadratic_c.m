function [irfquad,credibleset,model,models]=jointquadratic_c(IRF,cl)
%% Additively separable quadratric loss function:
% Bayes estimator and joint credible set for any M x nirf matrix IRF, where
% M is the number of posterior draws and nirf is the number of impulse
% responses, when Bayes estimator is not unique

M              = size(IRF,1);
c              = log(M)/sqrt(M);
quadratic_loss = zeros(M,1);
for i = 1:M
    for j = 1:M
        quadratic_loss(i,1) = quadratic_loss(i,1) + (IRF(j,:)-IRF(i,:))*(IRF(j,:)-IRF(i,:))';
    end
    if mod(i,100) == 0        
        fprintf('Computing absolute loss for %d of %d\n',i,M)
    end    
end
quadratic_loss = quadratic_loss/M;
[a,I]          = sort(quadratic_loss);
quadratic_loss = 100*(quadratic_loss-a(1))/a(1);
irfquad        = IRF(quadratic_loss<c,:);
credibleset    = IRF(I(1:ceil(cl*M)),:);
models         = I(1:ceil(cl*M));
model          = find(quadratic_loss<c,1,'first');





