function [irfabs,credibleset,model,models]=jointabsolute_c(IRF,cl)
%% Additively separable absolute loss function:
% Bayes estimator and joint credible set for any M x nirf matrix IRF, where 
% M is the number of posterior draws and nirf is the number of impulse responses, when Bayes estimator is not unique


M = size(IRF,1);
c = log(M)/sqrt(M);
absolute_loss = zeros(M,1);
for i = 1:M
    absolute_loss(i,1) =  sum(abs(IRF'-IRF(i,:)'),'all');
%     for j = 1:M
%         absolute_loss(i,1) = absolute_loss(i,1) + sum(abs(IRF(j,:)'-IRF(i,:)'));
%     end
    if mod(i,100) == 0        
        fprintf('Computing absolute loss for %d of %d\n',i,M)
    end
end
absolute_loss = absolute_loss/M;
[a,I]         = sort(absolute_loss);
absolute_loss = 100*(absolute_loss-a(1))/a(1);
model         = find(absolute_loss<c,1,'first');
irfabs        = IRF(model,:);
credibleset   = IRF(I(1:ceil(cl*M)),:);
models        = I(1:ceil(cl*M));



