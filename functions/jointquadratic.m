% JOINTQUADRATIC.M
%
% Additively separable quadratic loss function:
% Bayes estimator and joint credible set for any M x nirf matrix IRF, where 
% M is the number of posterior draws and nirf is the number of impulse 
% responses, when Bayes estimator is unique

function [irfquad,credibleset]=jointquadratic(IRF)

M=size(IRF,1);
quadratic_loss = zeros(M,1);
for i=1:M
   for j=1:M
       quadratic_loss(i,1) = quadratic_loss(i,1) + (IRF(j,:)-IRF(i,:))*(IRF(j,:)-IRF(i,:))';
   end
end
[~,Imin]=min(quadratic_loss);
irfquad = IRF(Imin,:);
[a,I] = sort(quadratic_loss);
credibleset = IRF(I(1:ceil(0.68*M)),:);