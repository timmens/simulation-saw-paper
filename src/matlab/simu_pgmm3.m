function [Y,X,Z]=simu_pgmm3(N,T,alpha,regime,sigma)
beta0 = alpha2beta(alpha,regime);
e = garch11_simu(T,0.05,0.05,0.9,N);
xi = randn(T,N);
%xi = ar1rnd(0.5,sqrt(3/4),T,N);
X = sqrt(2/3)*xi + sqrt(1/3)*e;
Z0 = sqrt(2/3)*xi + sqrt(1/3)*randn(T,N);
Z = [Z0(1:T-1,:) Z0(2:T,:)];
u = mean(X)';

X = reshape(X,N*T,1);
Z = reshape(Z,N*(T-1),2);
S = X.*repmat(beta0,N,1);
Y = S+kron(u,ones(T,1))+sigma*reshape(e,N*T,1);

