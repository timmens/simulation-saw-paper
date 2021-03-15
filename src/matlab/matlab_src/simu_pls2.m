function [Y,X]=simu_pls2(N,T,alpha,regime,sigma)
beta0 = alpha2beta(alpha,regime);
a = 0.5; 
e = reshape(ar1rnd(a,sqrt(1-a^2),T,N),N*T,1);
X = randn(N*T,1);
u = mean(reshape(X,T,N))';
S = X.*repmat(beta0,N,1);
Y = S+kron(u,ones(T,1))+sigma*e;
