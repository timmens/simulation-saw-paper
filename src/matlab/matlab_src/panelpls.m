function [regime,alpha,se,ssr,IC,K,lambda,DX,DY,Q] = panelpls(Y,X,N,option, mex)
% [regime,alpha,se,ssr,IC,K,lambda] = panelpls(Y,X,N,option)
% Estimate linear panel data model with multiple structural breaks
% Inputs:
%   Y: dependent variable, Y=(y1',...,yT')', NT by 1
%   X: independent variable, X=(x1',...,xT')', NT by p, p is num of
%           regressors.
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
%   ssr:    sum of squared residuals
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
XTol = 1e-3;%1e-4; 
maxIter = 200;%400;
%    b = prest_beta(Y,X,N);
b = postpls(Y,X,N,(1:T+1)');
weight = norms(diff(b),2,2).^(-2);
%B = [diag(-ones(T-1,1))+diag(ones(T-2,1),1) [zeros(T-2,1);1]]; 
%weight = norms(B*b,2,2).^(-2);
%weight = ones(T-1,1);
if mex == 1
    for i=1:S
    %    tht = pbcd(Y,X,N,lambda(i),weight,XTol,maxIter); % slower
        tht = plsbcd(Y,X,N,lambda(i),weight,XTol,maxIter); % mex, faster
        THT(:,p*(i-1)+1:p*i) = tht; 
        regime = getregime(tht,XTol);
        [~,ssr] = postpls(Y,X,N,regime);
        K(i) = length(regime)-2;
    %    IC(i) = ssr/(n-N) + (K(i)+1)*p*log(n-N)/(n-N);
    %    IC(i) = ssr/(n-N) + (K(i)+1)*p/sqrt(n-N);
        IC(i) = ssr/(n-N) + 0.05*(K(i)+1)*p*log(N*T)/(N*T)^(1/2);
    end
else 
    for i=1:S
        tht = pbcd(Y,X,N,lambda(i),weight,XTol,maxIter); % slower
    %    tht = plsbcd(Y,X,N,lambda(i),weight,XTol,maxIter); % mex, faster
        THT(:,p*(i-1)+1:p*i) = tht; 
        regime = getregime(tht,XTol);
        [~,ssr] = postpls(Y,X,N,regime);
        K(i) = length(regime)-2;
    %    IC(i) = ssr/(n-N) + (K(i)+1)*p*log(n-N)/(n-N);
    %    IC(i) = ssr/(n-N) + (K(i)+1)*p/sqrt(n-N);
        IC(i) = ssr/(n-N) + 0.05*(K(i)+1)*p*log(N*T)/(N*T)^(1/2);
    end
end
%disp([IC K lambda])
[~,i]=min(IC);
theta = THT(:,p*(i-1)+1:p*i);

regime = getregime(theta,XTol);
[alpha,ssr,se,DX,DY,Q] = postpls(Y,X,N,regime);


