profile on
N=100; T=20; sigma = 0.5;
% using option."whatever" creates a structure (struct in C lingo) option 
% with corresponding attributes (here nGrid, etc.)
option.nGrid = 40; option.maxLambda = 100; option.minLambda = 0.001; 

alpha0 = [0.3 1; 0.7 5]; regime0 = [1; 10; (T+1)];

%[Y,X,Z]=simu_dp_a(N,T,alpha0,regime0,sigma);
[Y,X,Z] = multivariate_simu_pgmm1(N,T,alpha0,regime0,sigma);

[regime,alpha,se,V,IC_opt,lambda_opt,IC,K,lambda] = panelpgmm(Y,X,Z,N,option);


beta0 = alpha2beta(alpha0,regime0); % true values 
beta = alpha2beta(alpha,regime);    % estimated values 

[beta0,beta]

figure(1), plot([beta0(:,1) beta(:,1)]), legend('True','Est')
figure(3), plot([beta0(:,2) beta(:,2)]), legend('True','Est')

profile viewer
%disp('True Set of Break Points')
%disp(regime0')
%disp('Estimated Set of Break Points')
%disp(regime')
%disp('Estimated Coefficients')
%disp([beta0(:,1) beta(:,1)])
%disp('IC')
%disp([IC K]) 
%figure(2), subplot(211), plot(lambda,IC), set(gca,'xscale','log')
%subplot(212), plot(lambda,K),set(gca,'xscale','log')

%disp([beta0(:,2) beta(:,2)])
%disp('IC')
%disp([IC K]) 

%figure(4), subplot(211), plot(lambda,IC), set(gca,'xscale','log')
%subplot(212), plot(lambda,K),set(gca,'xscale','log')