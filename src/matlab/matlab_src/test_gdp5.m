clear all;

rgdp_num = xlsread('rgdp.xls');
[N1,T1] = size(rgdp_num);

fdi_num = xlsread('fdi.xls');
fdi_num = fdi_num(:,1:T1);

gdp_num = xlsread('gdp.xls');
gdp_num = gdp_num(:,1:T1);

gdppc_num = xlsread('gdppc.xls');
gdppc_num = gdppc_num(:,1:T1);

ivs_num = xlsread('fixcap.xls');

f_t = (max(fdi_num,[],1)~=99999999999999900000);
f_i = (max(fdi_num,[],2)~=99999999999999900000);

%g_t = (max(rgdp_num,[],1)~=99999999999999900000);
%g_i = (max(rgdp_num,[],2)~=99999999999999900000);

Y = rgdp_num(f_i,:)';
X1 = (fdi_num(f_i,:)./gdp_num(f_i,:))'*100;
X2 = (ivs_num(f_i,:)./gdp_num(f_i,:))'*100;
X3 = log(gdppc_num(f_i,:))';

[T2,N] = size(Y);
outlier = []; %US 84; China 13;
index = 1:N;
Y = Y(:,setdiff(index,outlier));
X1 = X1(:,setdiff(index,outlier));
X2 = X2(:,setdiff(index,outlier));
X3 = X3(:,setdiff(index,outlier));
N = N - length(outlier);

tau = fix(T2/5);
tau0 = T2 - tau*5;
A = [ zeros(tau,tau0) kron(eye(tau),ones(1,5)/5)];
Y = A*Y;
X1 = A*X1;
X2 = A*X2;

A3 = [zeros(tau,tau0-1) kron(eye(tau),ones(1,5)/5) zeros(tau,1)];
%X3 = X3(3:5:T2,:);
X3 = A3*X3;

T0 = 1;
T = tau - T0;

%*** No exogenous variable
%option.maxLambda = 1e7;
%option.minLambda = 1e4;
%option.nGrid = 20;
%y = reshape(Y(T0+1:T+T0,:),N*T,1);
%x = reshape(Y(T0:T+T0-1,:),N*T,1);
%z = [reshape(X2(T0+1:T+T0-1,:),N*(T-1),1) reshape(Y(T0:T+T0-2,:),N*(T-1),1)];
% ***

% *** One exogenous variable
option.maxLambda = 5;
option.minLambda = 0.005;
option.nGrid = 50;
y = reshape(Y(T0+1:T+T0,:),N*T,1);
x = [reshape(Y(T0:T+T0-1,:),N*T,1) reshape(X1(T0+1:T+T0,:),N*T,1) reshape(X3(T0+1:T+T0,:),N*T,1)];
%z = [reshape(Y(T0:T+T0-2,:),N*(T-1),1) reshape(X1(T0+2:T+T0,:)-X1(T0+1:T+T0-1,:),N*(T-1),1) reshape(X1(T0+2:T+T0,:),N*(T-1),1) reshape(X3(T0+1:T+T0-1,:)-X3(T0:T+T0-2,:),N*(T-1),1) reshape(X3(T0:T+T0-2,:),N*(T-1),1)];
z = [reshape(Y(T0:T+T0-2,:),N*(T-1),1) reshape(X1(T0+2:T+T0,:)-X1(T0+1:T+T0-1,:),N*(T-1),1) reshape(X1(T0+2:T+T0,:),N*(T-1),1) reshape(X3(T0+2:T+T0,:),N*(T-1),1) reshape(X3(T0+1:T+T0-1,:),N*(T-1),1)];

% *** Estimate the model with the optimal tuning parameter
[regime,alpha,se,V,lambda_opt,IC_opt,IC,K,lambda] = panelpgmm(y,x,z,N,option);
pvalue = 2*(1-normcdf(abs(alpha./se)));
beta = alpha2beta(alpha,regime);
m = length(regime)-2;
year = 1977:5:2010;
disp('Estimated Set of Break Points')
disp(year(regime(2:m+1)))
disp('IC   K    Lambda')
disp([IC K lambda])
disp('Parameter estimates       S.E.      p-value');
disp([alpha se pvalue]);

figure(1), plot(year,beta), legend('AR coef');

figure(2), [ax, h1, h2] = plotyy(lambda,IC,lambda,K); set(ax(1),'xscale','log','FontSize',9), set(ax(2),'xscale','log','FontSize',9), h_legend=legend('IC_2', 'Number of Breaks');set(h_legend,'FontSize',9);
xlabel('Tuning Parameter'), ylabel(ax(1), 'IC_2'), ylabel(ax(2), 'Number of Breaks'), box off, legend boxoff; 
set(ax(2),'ytick',(0:6)), set(ax(1),'ytick',(22:2:32)/100);
% ***





