clear all;

% **** PLS          *****
sigma=0.5;
N=100; T=6;
option.maxLambda = 100; option.minLambda = 0.01; option.nGrid = 40;
alpha0 = [0.5; 1.5; 0.5]; regime0 = [1 T/3+1 T*2/3+1 T+1]';
%alpha0 = [0.3; 0.7]; regime0 = [1 T/2+1 T+1]';
%alpha0 = 1; regime0 = [1 T+1]';
[Y,X]=simu_pls1(N,T,alpha0,regime0,sigma);
%[regime,alpha,se,ssr,IC,K,lambda] = panelpls(Y,X,N,option);
[regime,alpha,se,ssr,IC,K,lambda] = panelpls(Y,X,N,option, 1);
%[regime,alpha,se,ssr,IC,K,lambda] = panelpls_slow(Y,X,N,option);

% **** End of PLS   *****


beta0 = alpha2beta(alpha0,regime0);
beta = alpha2beta(alpha,regime);

%disp('Estimated Coefficients')
%disp([beta0(:,1) beta(:,1)])
%disp('IC')
%disp([IC K])
figure(1), plot([beta0(:,1) beta(:,1)]), legend('True','Est')
figure(2), subplot(211), plot(lambda,IC), set(gca,'xscale','log')
subplot(212), plot(lambda,K),set(gca,'xscale','log')
disp('True Set of Break Points')
disp(regime0')
disp('Estimated Set of Break Points')
disp(regime')

