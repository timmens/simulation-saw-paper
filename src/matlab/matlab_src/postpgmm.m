function [alpha,V,W1,se] = postpgmm(Y,X,Z,N,regime,W)

[n,p] = size(X);
T = n/N; 
n1 = (T-1)*N;
[n2,q] = size(Z);
if n2==n1
    Z = reshape([ones(1,N*q);reshape(Z,T-1,N*q)],T*N,q);
end
if nargin<6 || isempty(W)
    W = repmat(reshape(eye(q),1,q*q),T,1);
end
m = length(regime)-2;
DY = diff(reshape(Y,T,N)); % (T-1) by N
DX = reshape(diff(reshape(X,T,N*p)),n1,p); %(T-1)*N by p

Q = zeros((m+1)*p,(m+1)*p); 
R = zeros((m+1)*p,1);

if m==0
    for i=regime(1)+1:regime(2)-1
        phi_zdx_0 = Z(i:T:n,:)'*DX(i-1:(T-1):n1,:)/N;
        Q(1:p,1:p) = Q(1:p,1:p) + phi_zdx_0' * reshape(W(i,:),q,q) * phi_zdx_0;
        phi_zdy_0 = Z(i:T:n,:)'*DY(i-1,:)'/N;
        R(1:p) = R(1:p) + phi_zdx_0' * reshape(W(i,:),q,q) * phi_zdy_0;
    end
else

phi_zx_1 = Z(regime(2):T:n,:)'*X(regime(2)-1:T:n,:)/N; % q by p
Q(1:p,1:p)= phi_zx_1'*reshape(W(regime(2),:),q,q)*phi_zx_1;

phi_zdy_1 = Z(regime(2):T:n,:)'*DY(regime(2)-1,:)'/N; % q by 1
R(1:p) =-phi_zx_1'*reshape(W(regime(2),:),q,q)*phi_zdy_1;
for i=regime(1)+1:regime(2)-1
    phi_zdx_0 = Z(i:T:n,:)'*DX(i-1:(T-1):n1,:)/N;
    Q(1:p,1:p) = Q(1:p,1:p) + phi_zdx_0' * reshape(W(i,:),q,q) * phi_zdx_0;
    phi_zdy_0 = Z(i:T:n,:)'*DY(i-1,:)'/N;
    R(1:p) = R(1:p) + phi_zdx_0' * reshape(W(i,:),q,q) * phi_zdy_0;
end


for j=2:m
    phi_zx_0 = Z(regime(j):T:n,:)'*X(regime(j):T:n,:)/N;
    phi_zx_1a = Z(regime(j):T:n,:)'*X(regime(j)-1:T:n,:)/N;
    phi_zx_1b = Z(regime(j+1):T:n,:)'*X(regime(j+1)-1:T:n,:)/N;
    Q((j-1)*p+1:j*p,(j-2)*p+1:(j-1)*p) = -phi_zx_0' * reshape(W(regime(j),:),q,q) * phi_zx_1a; % lower diagnonal
    Q((j-2)*p+1:(j-1)*p,(j-1)*p+1:j*p) = Q((j-1)*p+1:j*p,(j-2)*p+1:(j-1)*p)'; % higher diagonal
    Q((j-1)*p+1:j*p,(j-1)*p+1:j*p) = phi_zx_0' * reshape(W(regime(j),:),q,q) * phi_zx_0 + phi_zx_1b' * reshape(W(regime(j+1),:),q,q) * phi_zx_1b;
    phi_zdy_0 = Z(regime(j):T:n,:)'*DY(regime(j)-1,:)'/N; 
    phi_zdy_1 = Z(regime(j+1):T:n,:)'*DY(regime(j+1)-1,:)'/N; 
    R((j-1)*p+1:j*p) = phi_zx_0' * reshape(W(regime(j),:),q,q) * phi_zdy_0 - phi_zx_1b' * reshape(W(regime(j+1),:),q,q) * phi_zdy_1;
    for i=regime(j)+1:regime(j+1)-1
        phi_zdx_0 = Z(i:T:n,:)'*DX(i-1:(T-1):n1,:)/N;
        Q((j-1)*p+1:j*p,(j-1)*p+1:j*p) = Q((j-1)*p+1:j*p,(j-1)*p+1:j*p) + phi_zdx_0' * reshape(W(i,:),q,q) * phi_zdx_0;
        phi_zdy_0 = Z(i:T:n,:)'*DY(i-1,:)'/N;
        R((j-1)*p+1:j*p) = R((j-1)*p+1:j*p) + phi_zdx_0' * reshape(W(i,:),q,q) * phi_zdy_0;
    end
end

phi_zx_0 = Z(regime(m+1):T:n,:)'*X(regime(m+1):T:n,:)/N;
phi_zx_1a = Z(regime(m+1):T:n,:)'*X(regime(m+1)-1:T:n,:)/N;
Q(m*p+1:(m+1)*p,m*p+1:(m+1)*p) = phi_zx_0' * reshape(W(regime(m+1),:),q,q) * phi_zx_0;
Q(m*p+1:(m+1)*p,(m-1)*p+1:m*p) = -phi_zx_0' * reshape(W(regime(m+1),:),q,q) * phi_zx_1a;
Q((m-1)*p+1:m*p, m*p+1:(m+1)*p) = Q(m*p+1:(m+1)*p,(m-1)*p+1:m*p)';

phi_zdy_0 = Z(regime(m+1):T:n,:)'*DY(regime(m+1)-1,:)'/N;
R(m*p+1:(m+1)*p) = phi_zx_0' * reshape(W(regime(m+1),:),q,q) * phi_zdy_0;

for i=regime(m+1)+1:T
    phi_zdx_0 = Z(i:T:n,:)'*DX(i-1:(T-1):n1,:)/N;
    Q(m*p+1:(m+1)*p,m*p+1:(m+1)*p) = Q(m*p+1:(m+1)*p,m*p+1:(m+1)*p) + phi_zdx_0' * reshape(W(i,:),q,q) * phi_zdx_0;
    phi_zdy_0 = Z(i:T:n,:)'*DY(i-1,:)'/N;
    R(m*p+1:(m+1)*p) = R(m*p+1:(m+1)*p) + phi_zdx_0' * reshape(W(i,:),q,q) * phi_zdy_0;
end

end

alpha = reshape(Q\R,p,m+1)';
beta = alpha2beta(alpha,regime);
y1 = reshape(DY,n1,1);
res = y1 - sum(reshape(diff(reshape(X .* repmat(beta,N,1),T,N*p)),n1,p),2);
V = 0;
for t=2:T
    tmp = mean(Z(t:T:n,:).*repmat(res(t-1:(T-1):n1),1,q))'; % N by q
    V = V + tmp'*reshape(W(t,:),q,q)*tmp;
end
%R2 = 1-(res'*res)/(y1'*y1);

if nargout > 2
    W1 = zeros(T,q*q);
    for t=2:T
        tmp = Z(t:T:n,:).*repmat(res(t-1:(T-1):n1),1,q); % N by q
        W1(t,:) = reshape(inv(tmp'*tmp/N),1,q*q);
    end
end
    

if nargin > 3
    DU = reshape(res,T-1,N);
    A = zeros(p*(m+1),N);

    if m==0
        tmp = zeros(N,p); 
        for i=regime(1)+1:regime(2)-1
            phi_zdx_0 = Z(i:T:n,:)'*DX(i-1:(T-1):n1,:)/N; % p by q
            tmp = tmp + (Z(i:T:n,:).*repmat(DU(i-1,:)',1,q))*reshape(W(i,:),q,q)*phi_zdx_0; % 
        end
        A = tmp';
    else
        phi_zx_1 = Z(regime(2):T:n,:)'*X(regime(2)-1:T:n,:)/N; % p by q
        tmp = (-Z(regime(2):T:n,:).*repmat(DU(regime(2)-1,:)',1,q))*reshape(W(regime(2),:),q,q)*phi_zx_1;
        for i=regime(1)+1:regime(2)-1
           phi_zdx_0 = Z(i:T:n,:)'*DX(i-1:(T-1):n1,:)/N;
           tmp = tmp + (Z(i:T:n,:).*repmat(DU(i-1,:)',1,q))*reshape(W(i,:),q,q)*phi_zdx_0; % 
        end
        A(1:p,:) = tmp';
        for j=2:m
            phi_zx_0 = Z(regime(j):T:n,:)'*X(regime(j):T:n,:)/N;
            tmp = (Z(regime(j):T:n,:).*repmat(DU(regime(j)-1,:)',1,q))*reshape(W(regime(j),:),q,q)*phi_zx_0;
            phi_zx_1 = Z(regime(j+1):T:n,:)'*X(regime(j+1)-1:T:n,:)/N;
            tmp = tmp + (-Z(regime(j+1):T:n,:).*repmat(DU(regime(j+1)-1,:)',1,q))*reshape(W(regime(j+1),:),q,q)*phi_zx_1;
            for i=regime(j)+1:regime(j+1)-1
                phi_zdx_0 = Z(i:T:n,:)'*DX(i-1:(T-1):n1,:)/N;
                tmp = tmp + (Z(i:T:n,:).*repmat(DU(i-1,:)',1,q))*reshape(W(i,:),q,q)*phi_zdx_0;
            end
            A((j-1)*p+1:j*p,:) = tmp';
        end
        
        phi_zx_0 = Z(regime(m+1):T:n,:)'*X(regime(m+1):T:n,:)/N;
        tmp = (Z(regime(m+1):T:n,:).*repmat(DU(regime(m+1)-1,:)',1,q))*reshape(W(regime(m+1),:),q,q)*phi_zx_0;
        for i=regime(m+1)+1:T
            phi_zdx_0 = Z(i:T:n,:)'*DX(i-1:(T-1):n1,:)/N;
            tmp = tmp + (Z(i:T:n,:).*repmat(DU(i-1,:)',1,q))*reshape(W(i,:),q,q)*phi_zdx_0; % 
        end
        A(m*p+1:(m+1)*p,:) = tmp';
    end
    M = kron(diag(diff(regime).^(1/2)),eye(p));
    Omega = M\(A*A')/M/N;
    Phi = M\Q/M;
    S = M\(Phi\Omega/Phi)/M / N;
    se = reshape(sqrt(diag(S)),p,m+1)';

end

