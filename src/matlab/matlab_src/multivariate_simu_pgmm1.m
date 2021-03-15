function [Y,X,Z]=multivariate_simu_pgmm1(N,T,alpha,regime,sigma)
% two dimensional case, i.e. alpha must be mbar x 2  
beta0 = alpha2beta(alpha,regime);
e     = randn(T, N);
x1    = randn(T, N);
X1    = sqrt(2/3) * x1 + sqrt(1/3) * e;
x2    = randn(T, N);
X2    = 10 * x2 + e;
Z2    = 10 * x2 + randn(T, N);
Z1    = sqrt(2/3) * x1 + sqrt(1/3) * randn(T, N);
u     = mean(X1 + X2)';

X1    = reshape(X1, N * T, 1);
X2    = reshape(X2, N * T, 1);
S1    = X1.*repmat(beta0(:, 1), N, 1); 
S2    = X2.*repmat(beta0(:, 2), N, 1); 
Y     = S1 + S2 + kron(u,ones(T,1)) + sigma * reshape(e, N * T, 1);
X     = [X1 X2];
Z1    = [reshape(Z1(2:T,:)-Z1(1:(T-1),:),N*(T-1),1) reshape(Z1(2:T,:),N*(T-1),1)];
Z2    = [reshape(Z2(2:T,:)-Z2(1:(T-1),:),N*(T-1),1) reshape(Z2(2:T,:),N*(T-1),1)];
Z     = [Z1 Z2];
end 

