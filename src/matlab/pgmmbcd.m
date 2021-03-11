function tht = pgmmbcd(y,x,z,N,lambda,weight,W,XTol,maxIter)
% This function wraps up mpgmmbcd, the mex implementation of the block 
%   coordinated descent algorithm for the GMM estimation of panel model 
% 
%  tht = pgmmbcd(y,x,z,N,lambda,weight,W,XTol,maxIter)
%  Inputs:
%      y: NT-by-1 column vector of lhs variable, y = (y1,...,yN)
%      x: NT-by-p matrix of rhs variable, each column corresponds to a
%           regressor, x=(x1,...,xN)
%      z: N(T-1)-by-q matrix of IV's, z=(z1',...,zN')
%      N: number of cross-sections in the panel, an integer 
%      lambda: tuning parameter, scalar
%      weight: a (T-1)-by-1 weighting vector (default: ones(T-1,p))
%      XTol: a scalar error tolerance (default: 1e-4)
%      maxIter: an integer of the max num of iterations (default: 400)
%  Outputs:
%      tht: a T-by-p matrix of theta. cumsum(tht) gives an estimate of 
%               (beta_t) in the panel data model. See Qian and Su (2014)  
%  Note: the mex program mpgmmbcd is compiled on 64-bit Windows 7. To make it work 
%   on other platforms, you need to recompile the code using Matlab.To 
%   compile mex program, first make sure you have a compiler with which 
%   your Matlab is compatible. Then enter the following command:
%      >> mex mpgmmbcd.c nrutil.c mnbrak.c brent.c gaussj.c
%
%  Reference: Qian, J., Su, L., 2014, Shrinkage Estimation of Common Breaks
%       in Panel Data Models via Adaptive Group Fused Lasso, Working paper.
%
%  Junhui Qian
%  Aug 22, 2014
%  junhuiq@gmail.com


[n,p] = size(x);
[nz,q] = size(z);
T = n/N;
if nz ~= N*(T-1)
    tht = [];
    msgbox('Dimension of z must be N(T-1)-by-q, q>p.');
    return;
end
if q<=p 
    tht = [];
    msgbox('Dimension of z must be N(T-1)-by-q, q>p.');
    return;
end
    
if nargin<6 || isempty(weight)
    weight = ones(T-1,1);
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

if nargin>=7
    cw = zeros(T-1,1);
    for i=1:T-1;
        w = reshape(W(i,:),q,q);
        cw = cond(w);
    end
    if max(cw)>1e8
        W = repmat(reshape(eye(q),1,q*q),T-1,1);
    end
end
tht = mpgmmbcd(y,x,z,lambda,weight,W,XTol,maxIter);
tht = reshape(tht,T,p);