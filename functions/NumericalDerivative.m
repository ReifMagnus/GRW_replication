function Df = NumericalDerivative(x,f,epsilon)
%%
%  x - a vector of length n
%  f - function from R^n into R^m
%
%  Df - m x n matrix of partial derivatives evaluated at x
%    Df(i,j) = df_i/dx_j
%

n       =numel(x);
epsilon =1.0e-6;

z = x;
for j = 1:n
    z(j) = x(j) + epsilon;
    y1   = f(z);
    z(j) = x(j) - epsilon;
    y0   = f(z);
    z(j) = x(j);
    if j == 1
        Df = zeros(numel(y0),n);
    end
    Df(:,j) = (y1(:)-y0(:))/(2*epsilon);
end