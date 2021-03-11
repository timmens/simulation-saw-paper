function [Y,X]=simu_pls1(N,T,alpha,regime,sigma)
beta0 = alpha2beta(alpha,regime);
e = randn(N*T,1);
X = randn(N*T,1);
u = mean(reshape(X,T,N))';
S = X.*repmat(beta0,N,1);
Y = S+kron(u,ones(T,1))+sigma*e;
