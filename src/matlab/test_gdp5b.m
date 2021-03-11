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
A = [zeros(tau,tau0) kron(eye(tau),ones(1,5)/5) ];
Y = A*Y;
X1 = A*X1;
X2 = A*X2;

A3 = [zeros(tau,tau0-1) kron(eye(tau),ones(1,5)/5) zeros(tau,1)];
%X3 = X3(3:5:T2,:);
X3 = A3*X3;

T0 = 1;
T = tau - T0;

% *** Estimate the model with pre-specified tuning parameter
y = reshape(Y(T0+1:T+T0,:),N*T,1);
x = [reshape(Y(T0:T+T0-1,:),N*T,1) reshape(X1(T0+1:T+T0,:),N*T,1) reshape(X3(T0+1:T+T0,:),N*T,1)];
%z = [reshape(Y(T0:T+T0-2,:),N*(T-1),1) reshape(X1(T0+2:T+T0,:)-X1(T0+1:T+T0-1,:),N*(T-1),1) reshape(X1(T0+2:T+T0,:),N*(T-1),1) reshape(X3(T0+1:T+T0-1,:)-X3(T0:T+T0-2,:),N*(T-1),1) reshape(X3(T0:T+T0-2,:),N*(T-1),1)];
z = [reshape(Y(T0:T+T0-2,:),N*(T-1),1) reshape(X1(T0+2:T+T0,:)-X1(T0+1:T+T0-1,:),N*(T-1),1) reshape(X1(T0+2:T+T0,:),N*(T-1),1) reshape(X3(T0+2:T+T0,:),N*(T-1),1) reshape(X3(T0+1:T+T0-1,:),N*(T-1),1)];
    
%lambda = [0.02; 0.04; 0.08; 1; 10];
lambda = [0.005; 0.05; 0.15; 0.3; 0.5; 1; 10];
S = length(lambda);
year = 1980:5:2010;

for i=1:S
    [regime,alpha,se,V] = panelpgmm1(y,x,z,N,lambda(i));
    m = length(regime)-2;
    disp('**********************************')
    disp(['lambda= ' num2str(lambda(i))]);
    disp(['Breakdates: ' num2str(year(regime(2:m+1)))]);
    disp('Parameter estimates       S.E.       p-value');
    disp([alpha se 2*(1-normcdf(abs(alpha./se)))]);
end
    

