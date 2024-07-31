function [B_0,V_0] = SetMinnesotaPrior(n,p,alpha,Vc,lambda,psi,stationary,exo)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(exo)
    n_exo = size(exo,2);
else
    n_exo = 0;
end

K = n*p + n_exo;


lag_1_prior             = ones(1,n);         % Trending variables centered at one
lag_1_prior(stationary) = 0;                 % Stationary variables centered instead at zero

B_0 = [zeros(n_exo,n),; diag(lag_1_prior); zeros((n)*(p-1),n)];         % prior means
d   = n+2;                                                                                  % prior degrees of freedom    


Omega = zeros(K,1);

if ~isempty(exo)
    Omega(1:n_exo) = Vc;
end


for i = 1:p
    Omega(n_exo + (i-1)*n + 1:n_exo + i*n) = (d-n-1)*(lambda^2)*(1/(i^alpha))./psi;
end

V_0 = diag(1./Omega);

end

