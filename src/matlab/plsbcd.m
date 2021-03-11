function tht = plsbcd(y,x,N,lambda,weight,XTol,maxIter)
% This function wraps up mplsbcd, the mex implementation of the block 
%   coordinated descent algorithm for the LS estimation of panel model
% 
%  tht = plsbcd(y,x,N,lambda,weight,XTol,maxIter)
%  Inputs:
%      y: lhs variable, an NT by 1 column vector
%      x: rhs variables, an NT by p matrix, each column corresponds to a
%           regressor
%      lambda: tuning parameter, scalar
%      weight: a (T-1)-by-1 weighting vector (default: ones(T-1,p))
%      XTol: a scalar error tolerance (default: 1e-4)
%      maxIter: an integer of the max num of iterations (default: 400)
%  Outputs:
%      tht: a T-by-p matrix of theta. cumsum(tht) gives an estimate of 
%               (beta_t) in the panel data model. See Qian and Su (2014)  
%  Note: the mex program is compiled on 64-bit Windows 7. To make it work 
%   on other platforms, you need to recompile the code using Matlab.To 
%   compile mex program, first make sure you have a compiler with which 
%   your Matlab is compatible. Then enter the following command:
%      >> mex mplsbcd.c nrutil.c mnbrak.c brent.c gaussj.c
%
%  Reference: Qian, J., Su, L., 2014, Shrinkage Estimation of Common Breaks
%       in Panel Data Models via Adaptive Group Fused Lasso, Working paper.
%
%  Junhui Qian
%  Aug 22, 2014
%  junhuiq@gmail.com

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

tht = mplsbcd(y,x,lambda,weight,XTol,maxIter);
tht = reshape(tht,T,p);
