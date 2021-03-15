function [regime,alpha,se,V] = panelpgmm1(Y,X,Z,N,lambda)
% [regime,alpha,se,V,IC,K,lambda] = panelpgmm(Y,X,N,Z,option)
% Estimate linear panel data model with multiple structural breaks
% Inputs:
%   Y: dependent variable, Y=(y1',...,yT')', NT by 1
%   X: independent variable, X=(x1',...,xT')', NT by p, p is num of
%           regressors.
%   Z: instruments, Z=(z2',...,zT')', N(T-1) by q, q is num of instruments
%   N: number of cross-sections
%
% Outputs:  
%   regime: set of break points, regime = {1,tau_1,tau_2,...,tau_m,T+1}, m
%           is num of break points. 
%   alpha:  estimated coefficients for each regime, (m+1) by p
%   se:     estimated standard error for each parameter, (m+1) by p
%   V:      value of the quadratic term
%  
% Reference:
% Qian & Su, 2013, Shrinkage estimation of common breaks in panel data models via adaptive
%   group fused lasso.
%   Nov 1, 2013

[n,p]=size(X);
T = n/N;

XTol = 1e-4; 

[b,tmp2,W0] = postpgmm(Y,X,Z,N,(1:T+1)');
weight = norms(diff(b),2,2).^(-2);
%tht = pgmmbcd(Y,X,Z,N,lambda,weight,W0(2:T,:));
tht = pbcdgmm(Y,X,Z,N,lambda,weight,W0(2:T,:));
regime = getregime(tht,XTol);
[tmp1,tmp2,W] = postpgmm(Y,X,Z,N,regime,W0); % obtain weighting matrices 
[alpha,V,tmp1,se] = postpgmm(Y,X,Z,N,regime,W);


