function [alpha,ssr,se, DX, DY, Q] = postpls(Y,X,N,regime)
[n,p] = size(X);
T = n/N;
m=size(regime,1)-2;
Q = zeros((m+1)*p,(m+1)*p);
n1 = (T-1)*N;
DX=reshape(diff(reshape(X,T,N*p)),n1,p);
R = zeros((m+1)*p,1);

DY = diff(reshape(Y,T,N));
if m==0
    
%      alpha = ((DX'*DX)\(DX'*reshape(DY,n1,1)))';
    Q = DX'*DX/N;
    R = DX'*reshape(DY,n1,1)/N;
    
else
    
    % first block 
   Q(1:p,1:p)=(X(regime(2)-1:T:(T*N),:)'*X(regime(2)-1:T:n,:))/N;
   % T_1 >= 3 ??? 
   for l=0:regime(2)-regime(1)-2
     Q(1:p,1:p)=Q(1:p,1:p)+(DX(regime(1)+l:T-1:n-N,:)'*DX(regime(1)+l:T-1:n-N,:))/N;
   end
   
   for t=2:m
    Q((t-1)*p+1:t*p,(t-1)*p+1:t*p)= X(regime(t+1)-1:T:n,:)'*X(regime(t+1)-1:T:n,:)/N + X(regime(t):T:n,:)'*X(regime(t):T:n,:)/N;
    for l=0:regime(t+1)-regime(t)-2
      Q((t-1)*p+1:t*p,(t-1)*p+1:t*p)= Q((t-1)*p+1:t*p,(t-1)*p+1:t*p)+DX(regime(t)+l:T-1:n-N,:)'*DX(regime(t)+l:T-1:n-N,:)/N;
    end
   end
   
   % last block 
   Q(m*p+1:(m+1)*p ,m*p+1:(m+1)*p)= X(regime(m+1):T:n,:)'*X(regime(m+1):T:n,:)/N;
   for l=0:T-regime(m+1)-1
     Q(m*p+1:(m+1)*p,m*p+1:(m+1)*p)= Q(m*p+1:(m+1)*p,m*p+1:(m+1)*p)+(DX(regime(m+1)+l:T-1:n-N,:)'*DX(regime(m+1)+l:T-1:n-N,:))/N;                                
   end
   
   for t=2:m+1
    Q((t-1)*p+1:t*p,(t-2)*p+1:(t-1)*p) = -X(regime(t):T:n,:)'*X(regime(t)-1:T:n,:)/N;
    Q((t-2)*p+1:(t-1)*p,(t-1)*p+1:t*p) = -X(regime(t)-1:T:n,:)'*X(regime(t):T:n,:)/N;
   end
   
   
   R(1:p) = -(X(regime(2)-1:T:n,:)'*DY(regime(2)-1,:)')/N ;
   for l=0:regime(2)-regime(1)-2
       R(1:p) = R(1:p)+(DX(regime(1)+l:T-1:n-N,:)'*DY(regime(1)+l,:)')/N;
   end
   for t=2:m
          R((t-1)*p+1:t*p) = -(X(regime(t+1)-1:T:n,:)'*DY(regime(t+1)-1,:)')/N +(X(regime(t):T:n,:)'*DY(regime(t)-1,:)')/N ;
      for l=0:regime(t+1)-regime(t)-2
         R((t-1)*p+1:t*p) =  R((t-1)*p+1:t*p)+DX(regime(t)+l:T-1:n-N,:)'*DY(regime(t)+l,:)'/N;
      end
   end
    R(m*p+1:(m+1)*p)=X(regime(m+1):T:n,:)'*DY(regime(m+1)-1,:)'/N;
   for l=0:T-regime(m+1)-1
    R(m*p+1:(m+1)*p)=R(m*p+1:(m+1)*p)+DX(regime(m+1)+l:T-1:n-N,:)'*DY(regime(m+1)+l,:)'/N;
   end

end

alpha =reshape(Q\R,p,m+1)';


beta = alpha2beta(alpha,regime);
res = reshape(DY,n1,1) - sum(reshape(diff(reshape(X .* repmat(beta,N,1),T,N*p)),n1,p),2);
ssr = res'*res;

if nargout > 2
    DU = reshape(res,T-1,N);
    A = zeros(p*(m+1),N);
if m==0
    tmp = zeros(N,p);
%    S = inv(DX'*DX)*ssr/(n1-(m+1)*p);
   for i=regime(1)+1:regime(2)-1
       tmp = tmp + DX(i-1:T-1:n-N,:).*repmat(DU(i-1,:)',1,p);
%       tmp(l*N+1:(l+1)*N,:) = DX(regime(1)+l:T-1:n-N,:).*repmat(DU(regime(1)+l,:)',1,p);
   end
   A = tmp';

else
   tmp = -X(regime(2)-1:T:n,:).*repmat(DU(regime(2)-1,:)',1,p);
   for i=regime(1)+1:regime(2)-1
       tmp = tmp + DX(i-1:T-1:n-N,:).*repmat(DU(i-1,:)',1,p);
%       tmp(l*N+1:(l+1)*N,:) = DX(regime(1)+l:T-1:n-N,:).*repmat(DU(regime(1)+l,:)',1,p);
   end
   A(1:p,:) = tmp';
   
   for t=2:m
      tmp = X(regime(t):T:n,:).*repmat(DU(regime(t)-1,:)',1,p);
      tmp = tmp - X(regime(t+1)-1:T:n,:).*repmat(DU(regime(t+1)-1,:)',1,p);
      for i=regime(t)+1:regime(t+1)-1
          tmp = tmp + DX(i-1:T-1:n-N,:).*repmat(DU(i-1,:)',1,p);
%          tmp(l*N+1:(l+1)*N,:) = DX(regime(t)+l:T-1:n-N,:).*repmat(DU(regime(t)+l,:)',1,p);
      end
      A((t-1)*p+1:t*p,:) = tmp';
   end
   
   tmp = X(regime(m+1):T:n,:).*repmat(DU(regime(m+1)-1,:)',1,p);
   for i=regime(m+1)+1:T
        tmp = tmp + DX(i-1:T-1:n-N,:).*repmat(DU(i-1,:)',1,p);
%       tmp(l*N+1:(l+1)*N,:) = DX(regime(m+1)+l:T-1:n-N,:).*repmat(DU(regime(m+1)+l,:)',1,p);
   end
   A(m*p+1:(m+1)*p,:)=tmp';
end
    M = kron(diag(diff(regime).^(1/2)),eye(p));
    Omega = M\(A*A')/M/N;
    Phi = M\Q/M;
    S = M\(Phi\Omega/Phi)/M / N;
   se = reshape(sqrt(diag(S)),p,m+1)';
end


