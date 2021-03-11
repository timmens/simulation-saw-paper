function [Y,X,Z]=simu_pgmm1(N,T,alpha,regime,sigma)
beta0 = alpha2beta(alpha,regime);
e = randn(T,N);
xi = randn(T,N);
X = sqrt(2/3)*xi + sqrt(1/3)*e;
Z0 = sqrt(2/3)*xi + sqrt(1/3)*randn(T,N);
u = mean(X)';

X = reshape(X,N*T,1);
S = X.*repmat(beta0,N,1);
Y = S+kron(u,ones(T,1))+sigma*reshape(e,N*T,1);
Z = [reshape(Z0(2:T,:)-Z0(1:T-1,:),N*(T-1),1) reshape(Z0(2:T,:),N*(T-1),1)];

