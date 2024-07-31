% JOINTANGULAR.M

function [irfangular,credibleset]=jointangular(IRF,hmax,n)

M=size(IRF,1);
angular_distance = zeros(M,1);
for i=1:M
    IRFi = reshape(IRF(i,:),hmax,n*n);
    for j=1:M
        IRFj = reshape(IRF(j,:),hmax,n*n);
        for k=1:n
            for ell=1:n
                irfi = IRFi(:,n*(k-1)+ell);
                irfi = irfi/norm(irfi);
                irfj = IRFj(:,n*(k-1)+ell);
                irfj = irfj/norm(irfj); 
                angular_distance(i,1) = angular_distance(i,1)+acos(dot(irfi,irfj))/(pi*n*n);
            end
        end
    end
end
[~,Iangular]=min(angular_distance);
irfangular = IRF(Iangular,:);
[a,I] = sort(angular_distance);
credibleset = IRF(I(1:floor(0.68*M)+1),:);