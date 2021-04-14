clear all;

N=500; T=6; sigma = 0.3;
option.nGrid = 40; option.maxLambda = 100; option.minLambda = 0.001; 

alpha0 = [0.3; 0.7; 0.3]; regime0 = [1 T/3+1 T*2/3+1 T+1]';
%alpha0 = [0.3; 0.7]; regime0 = [1 T/2+1 T+1]';
%alpha0 = [0.5]; regime0 = [1 T+1]';

%[Y,X,Z]=simu_dp_a(N,T,alpha0,regime0,sigma);
[Y,X,Z]=simu_pgmm1(N,T,alpha0,regime0,sigma);

[regime,alpha,se,V,IC_opt,lambda_opt,IC,K,lambda] = panelpgmm(Y,X,Z,N,option);


beta0 = alpha2beta(alpha0,regime0);
beta = alpha2beta(alpha,regime);

disp('True Set of Break Points')
disp(regime0')
disp('Estimated Set of Break Points')
disp(regime')
disp('Estimated Coefficients')
disp([beta0(:,1) beta(:,1)])
%disp('IC')
%disp([IC K])
figure(1), plot([beta0(:,1) beta(:,1)]), legend('True','Est')
figure(2), subplot(211), plot(lambda,IC), set(gca,'xscale','log')
subplot(212), plot(lambda,K),set(gca,'xscale','log')
