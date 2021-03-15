function [regime,alpha,se,V,IC_opt,lambda_opt,IC,K,lambda] = panelpgmm(Y,X,Z,N,option)
% [regime,alpha,se,V,IC,K,lambda] = panelpgmm(Y,X,N,Z,option)
% Estimate linear panel data model with multiple structural breaks
% Inputs:
%   Y: dependent variable, Y=(y1',...,yT')', NT by 1
%   X: independent variable, X=(x1',...,xT')', NT by p, p is num of
%           regressors.
%   Z: instruments, Z=(z2',...,zT')', N(T-1) by q, q is num of instruments
%   N: number of cross-sections
%   option: a construct of options
%       option.maxLambda: maximum lambda
%       option.minLambda: minimum lambda
%       option.nGrid: number of grids on the range of lambda
%
% Outputs:  
%   regime: set of break points, regime = {1,tau_1,tau_2,...,tau_m,T+1}, m
%           is num of break points. 
%   alpha:  estimated coefficients for each regime, (m+1) by p
%   se:     estimated standard error for each parameter, (m+1) by p
%   V:      value of the quadratic term
%   IC,K,lambda: for debug use only
%  
% Reference:
% Qian & Su, 2013, Shrinkage estimation of common breaks in panel data models via adaptive
%   group fused lasso.
%   Nov 1, 2013

[n,p]=size(X);
T = n/N;

lam_max = option.maxLambda;
lam_min = option.minLambda;
S = option.nGrid;

R = log(lam_max/lam_min)/(S-1);
lambda = lam_min*exp(R*(0:S-1)');

IC = zeros(S,1);
K = zeros(S,1);
THT = zeros(T,S*p);
XTol = 1e-4; 

[b,tmp2,W0] = postpgmm(Y,X,Z,N,(1:T+1)');
weight = norms(diff(b),2,2).^(-2);
for i=1:S
    %tht = pbcdgmm(Y,X,Z,N,lambda(i),weight,W0(2:T,:)); % slower
    tht = pgmmbcd(Y,X,Z,N,lambda(i),weight,W0(2:T,:)); % mex, faster
    THT(:,p*(i-1)+1:p*i) = tht; 
    regime = getregime(tht,XTol);
    [tmp1,tmp2,W] = postpgmm(Y,X,Z,N,regime,W0); % obtain weighting matrices 
    [tmp1,V] = postpgmm(Y,X,Z,N,regime,W);
    K(i) = length(regime)-2;
%    IC(i) = V/(T-1) + (K(i)+1)*p/sqrt(n-N);
    IC(i) = V/(T-1) + 0.05*(K(i)+1)*p*log(N*T)/(N*T)^(1/2);
end
%disp([IC K lambda])
[tmp,i]=min(IC);
IC_opt = IC(i);
lambda_opt = lambda(i);
theta = THT(:,p*(i-1)+1:p*i);

regime = getregime(theta,XTol);
%[alpha,V,tmp,se] = postpgmm(Y,X,Z,N,regime);
[tmp1,tmp2,W] = postpgmm(Y,X,Z,N,regime,W0);
[alpha,V,tmp,se] = postpgmm(Y,X,Z,N,regime,W);



