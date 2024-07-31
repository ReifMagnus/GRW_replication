function Df = NumericalDerivative_g(x,zr,n,k,p,A0,Aplus,epsilon)
%%
%  x - a vector of length n
%  f - function from R^n into R^m
%
%  Df - m x n matrix of partial derivatives evaluated at x
%    Df(i,j) = df_i/dx_j
%

nx = numel(x);
if nargin < 8
    epsilon=1.0e-6;
end
z = x;
for j = 1:nx
    z(j) = x(j) + epsilon;
    y1   = IRFRestrictions(x,zr,n,k,p,0,A0,Aplus);
    z(j) = x(j) - epsilon;
    y0   = IRFRestrictions(x,zr,n,k,p,0,A0,Aplus);
    z(j) = x(j);
    if j == 1
        Df = zeros(numel(y0),nx);
    end
    Df(:,j) = (y1-y0)/(2*epsilon);
end