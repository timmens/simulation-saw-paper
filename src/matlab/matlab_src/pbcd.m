function d = pbcd(y,x,N,lambda,weight,XTol,maxIter)
% The block coordinated descent algorithm for panel data model

[n,p] = size(x);
T = n/N;

if nargin<5 || isempty(weight)
    weight = ones(T,1);
end

if nargin<6 || isempty(XTol)
    XTol = 1e-4; 
end

if nargin<7 || isempty(maxIter)
    maxIter = 400;
end

A = kron(tril(ones(T)),eye(p)); % Tp by Tp
X = zeros(n-N,T*p); % (T-1)N by Tp
dy = zeros(n-N,1);
for i=1:T-1
    X((i-1)*N+1:i*N,(i-1)*p+1:i*p)=-x(i:T:n,:);
    X((i-1)*N+1:i*N,i*p+1:(i+1)*p)=x(i+1:T:n,:);
    dy((i-1)*N+1:i*N) = y(i+1:T:n)-y(i:T:n);
end
Z = X*A; % (T-1)N by Tp
R = Z'*Z/N; % Tp by Tp
r = Z'*dy/N; 

d0 = zeros(T*p,1);
d1 = zeros(T*p,1);

e = 1e10;
i = 0;
while e>XTol && i<maxIter
    i = i+1;
    for t = 1:T
        if t == 1
            g = r(1:p)-R(1:p,p+1:T*p)*d0(p+1:p*T);
            d1(1:p) = R(1:p,1:p)\g;
        else
            g = R((t-1)*p+1:t*p,1:(t-1)*p)*d1(1:(t-1)*p) + R((t-1)*p+1:t*p,t*p+1:T*p)*d0(t*p+1:T*p) - r((t-1)*p+1:t*p);
            if norm(g)<=lambda*weight(t-1)
                d1((t-1)*p+1:t*p) = 0;
            else
                [gam,fval,exitflag] = fminsearch(@myfun2,1/N,[],R((t-1)*p+1:t*p,(t-1)*p+1:t*p),lambda*weight(t-1),g);
%                [gam,fval,exitflag] = fminbnd(@myfun2,0,10,[],R((t-1)*p+1:t*p,(t-1)*p+1:t*p),lambda,g);
%                [gam,fval,exitflag] = fzero(@myfun,[0;1],[],R((t-1)*p+1:t*p,(t-1)*p+1:t*p),lambda,g);
                if exitflag==1
                    d1((t-1)*p+1:t*p) = -gam*((gam*R((t-1)*p+1:t*p,(t-1)*p+1:t*p) + (lambda*weight(t-1))^2/2*eye(p))\g);
                else
                    disp('Line search failed');
                    return;
                end
            end
        end
    end
    e = norm(d1-d0);
    d0 = d1;
%    disp([i e])
end
d = reshape(d1',p,T)';

function y = myfun(gam,R,lam,g)
p = length(g);
tmp=lam/2 * ((gam*R+lam^2/2*eye(p))\g);
y =  sum(tmp.^2) - 1;

function y = myfun2(gam,R,lam,g)
if gam<0
    y=1e10;
    return;
end
p = length(g);
y = gam*(1-1/2*g'*((gam*R+lam^2/2*eye(p))\g));


