function [irfabs,credibleset,model,models]=jointabsolute(IRF,cl)
%% Additively separable absolute loss function:
% Bayes estimator and joint credible set for any M x nirf matrix IRF, where
% M is the number of posterior draws and nirf is the number of impulse responses, when Bayes estimator is unique


M=size(IRF,1);
absolute_loss = zeros(M,1);
for i=1:M
    for j=1:M
        absolute_loss(i,1) = absolute_loss(i,1)+sum(abs(IRF(j,:)-IRF(i,:)));
    end
    if mod(i,100) == 0
        fprintf('Computing absolute loss for %d of %d\n',i,M)
    end
end
[~,model]    = min(absolute_loss);
irfabs      = IRF(model,:);
[a,I]       = sort(absolute_loss);
models      = I(1:ceil(cl*M));
credibleset = IRF(models,:);




