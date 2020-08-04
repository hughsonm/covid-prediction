function [alpha,beta,mu] = ExponentialFit(x,y)
%ExponentialFit Fits exponential model to data by linearizing the model
%   Fits a model of y = alpha*exp(beta*x) to the given data
u = x;
v = log(y);
[coeffs,err_struct,mu] = polyfit(u,v,1);
alpha = exp(coeffs(2));
beta = coeffs(1);
end

