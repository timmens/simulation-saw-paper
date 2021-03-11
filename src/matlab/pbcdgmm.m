function d = pbcdgmm(y,x,z,N,lambda,weight,W,XTol,maxIter)
% The block coordinated descent algorithm for GMM estimation of 
% panel data model with structural breaks. 
% y: NT by 1 (y1,...,yN)
% x: NT by p (x1',...,xN')
% z: N(T-1) by q (z1',...,zN')
% N: scalar, number of crosssections
% lambda: scalar, penalty term
% W: (T-1) by q*q, weighting matrices, identity matrices in default


[n,p] = size(x);
q = size(z,2);
T = n/N;

options = optimset('Display','off');

if nargin<6 || isempty(weight)
    weight = ones(T,1);
end

if nargin<7 || isempty(W)
    W = repmat(reshape(eye(q),1,q*q),T-1,1);
end

if nargin<8 || isempty(XTol)
    XTol = 1e-4; 
end

if nargin<9 || isempty(maxIter)
    maxIter = 400;
end

dy = zeros(n-N,1);
dx = zeros(n-N,p);
tmp = zeros(n-N,q);
tmp2 = zeros(n-N,p);
for i=1:T-1
    dx((i-1)*N+1:i*N,:) = x(i+1:T:n,:)-x(i:T:n,:);
    dy((i-1)*N+1:i*N) = y(i+1:T:n)-y(i:T:n);
    tmp((i-1)*N+1:i*N,:) = z(i:(T-1):n-N,:);
    tmp2((i-1)*N+1:i*N,:) = x(i+1:T:n,:);
end
z = tmp;
x = tmp2;

delta=zeros((T-1),q);
Q = zeros((T-1),q*p);
P = zeros((T-1),q*p);
for i=1:T-1
    delta(i,:) = dy((i-1)*N+1:i*N)'*z((i-1)*N+1:i*N,:)/N;
    Q(i,:) = reshape(x((i-1)*N+1:i*N,:)'*z((i-1)*N+1:i*N,:)/N,1,p*q);
    P(i,:) = reshape(dx((i-1)*N+1:i*N,:)'*z((i-1)*N+1:i*N,:)/N,1,p*q);
end


R = zeros(T-1,p*p);
U = zeros(T-1,p*p);
V = zeros(T-1,p);
M = zeros(T-1,p*p);
L = zeros(T-1,p);

for i=1:T-1
    P_ = reshape(P(i,:),q,p);
    W_ = reshape(W(i,:),q,q);
    Q_ = reshape(Q(i,:),q,p);
    R(i,:) = reshape(P_' * W_ * P_,1,p*p);
    U(i,:) = reshape(P_' * W_ * Q_,1,p*p);
    M(i,:) = reshape(Q_' * W_ * Q_,1,p*p);
    V(i,:) = (P_' * W_ * delta(i,:)')';
    L(i,:) = (Q_' * W_ * delta(i,:)')';
end

Rc = cumsum(R(T-1:-1:1,:));
Vc = cumsum(V(T-1:-1:1,:));

d0 = zeros(T,p);
d1 = zeros(T,p);

e = 1e10;
i = 0;
while e>XTol && i<maxIter
    i = i+1;
    for t = 1:T
        if t == 1
            d0c = cumsum(d0(2:T-1,:));
            g = reshape(U(1,:),p,p)*d0(2,:)';
            for s=1:T-2
                g = g + reshape(U(s+1,:),p,p)*d0(s+2,:)' + reshape(R(s+1,:),p,p)*d0c(s,:)';
            end
            g = g - Vc(T-1,:)';
            d1(1,:) = (reshape(Rc(T-1,:),p,p)\g)';
        elseif t < T
            d1s = sum(d1(1:t-1,:),1);
            g = -Vc(T-t,:)' + reshape(U(t-1,:),p,p)'*d1s' - L(t-1,:)';
            for s=t+1:T
                g = g + reshape(R(s-1,:),p,p)'*(d1s + sum(d0(t+1:s-1,:),1))' + reshape(U(s-1,:),p,p)'*d0(s,:)';
            end
            A = reshape(M(t-1,:) + Rc(T-t+1,:),p,p);
            if norm(g)<=lambda*weight(t-1)
                d1(t,:) = 0;
            else
                [tmp,p1] = chol(A);
                if p1==0 % check whether A is negative definite
                    [gam,fval,exitflag] = fminsearch(@myfun2,1/N,options,A,lambda*weight(t-1),g);
%                   [gam,fval,exitflag] = fminbnd(@myfun2,0,10,[],R((t-1)*p+1:t*p,(t-1)*p+1:t*p),lambda,g);
%                   [gam,fval,exitflag] = fzero(@myfun,[0;1],[],R((t-1)*p+1:t*p,(t-1)*p+1:t*p),lambda,g);
                    if exitflag==1
                        d1(t,:) = -gam*((gam*A + (lambda*weight(t-1))^2/2*eye(p))\g)';
                    else
                        d = d1;
                        return;
                    end
                end
                    
            end
        else
            d1s = sum(d1(1:T-1,:),1);
            g = reshape(U(T-1,:),p,p)'*d1s' - L(T-1,:)';
            A = reshape(M(T-1,:),p,p);
            if norm(g)<=lambda*weight(T-1)
                d1(T,:) = 0;
            else
                [tmp,p1] = chol(A);
                if p1==0
                    [gam,fval,exitflag] = fminsearch(@myfun2,1/N,options,A,lambda*weight(T-1),g);
%                    [gam,fval,exitflag] = fminbnd(@myfun2,0,10,[],R((t-1)*p+1:t*p,(t-1)*p+1:t*p),lambda,g);
%                    [gam,fval,exitflag] = fzero(@myfun,[0;1],[],R((t-1)*p+1:t*p,(t-1)*p+1:t*p),lambda,g);
                    if exitflag==1
                        d1(T,:) = -gam*((gam*A + (lambda*weight(T-1))^2/2*eye(p))\g)';
                    else
                        d = d1;
                        return;
                    end
                end
            end
            
        end
    end
    e = norm(d1-d0);
    d0 = d1;
end
d = d1;

function y = myfun(gam,R,lam,g)
p = length(g);
tmp=lam/2 * ((gam*R+lam^2/2*eye(p))\g);
y =  sum(tmp.^2) - 1;

function y = myfun2(gam,R,lam,g)
if gam<0
    y=1e20;
    return;
end
p = length(g);
y = gam*(1-1/2*g'*((gam*R+lam^2/2*eye(p))\g));


