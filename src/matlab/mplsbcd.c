/*********************************************************************
 * mplsbcd.c
 * This program implements the BCD (Block-Coordinate Descent) algorithm of 
 * PLS estimation of panel data model.
 * tht = mplsbcd(y,x,lambda,weight,XTol,maxIter);
 *  Output: 
 *      y: lhs variable, an n by 1 column vector
 *  Inputs:
 *      x: rhs variables, an n by m matrix
 *      lambda: tuning parameter, scalar
 *      weight: weighting vector, (n-1) by 1
 *      XTol: error tolerance, scalar
 *      maxIter: max num of iterations, scalar integer
 *  Subroutines needed:
 *      useful.c useful.h
 *  To compile mex program, enter the following statement in Maltab:
 *      >> mex mplsbcd.c useful.c
 *
 *  Junhui Qian
 *  Shanghai Jiao Tong University
 *  Aug 22, 2014
 *  junhuiq@gmail.com
 *
 ********************************************************************/
#include <matrix.h>
#include <mex.h>   
#include <math.h>
#include "useful.h"

double func(double gam, const double *R, mwSize p, const double lam, const double *g){
    mwSize i;
    double *A, *q, sum;
    
	if(gam<0)
		return 1e+6;

	A   = (double *)malloc(p*p*sizeof(double));
	q   = (double *)malloc(p*sizeof(double));
    for(i=0;i<p*p;i++)
        A[i] = gam*R[i];
    for(i=0;i<p;i++)
        A[i*p+i] = A[i*p+i]+lam*lam/2;
    for(i=0;i<p;i++)
        q[i] = g[i];
    gaussj(A, p, q, 1);
    sum=0;
    for(i=0;i<p;i++)
        sum = sum + g[i]*q[i];
    free(A);
    free(q);
    
    return gam*(1-0.5*sum);
}

void wmplsbcd(double *y,double *x,int N,int T,int p,double *lambda,double *weight,double *XTol,double *maxIter,double *b)
{
    double *dx,*dy,*S_xx,*S_dxdx,*S_xdx,*S_dxdy,*S_xdy,*CS_dxdx,*CS_dxdy;
    double *b0,*g,*vTmp1,*vTmp2,*mTmp,*A,*r,*theta,*err,*Rn,*rv,*rvn,*R; 
    double ax,bx,cx,fa,fb,fc,xmin,tmp1,tmp2,tmp3,errn,dTmp,lam;
    int n1,n,i,j,k,l,m,t,exitflag;
    
    ax = 1e-16; bx = 1e+6;

    n1 = (T-1)*N; n = N*T;

	err    = (double *)malloc(T*p*sizeof(double));
	dx     = (double *)malloc(n1*p*sizeof(double));
	dy     = (double *)malloc(n1*sizeof(double));
	S_xx   = (double *)malloc((T-1)*p*p*sizeof(double));
	S_dxdx = (double *)malloc((T-1)*p*p*sizeof(double));
	S_xdx  = (double *)malloc((T-1)*p*p*sizeof(double));
	S_dxdy = (double *)malloc((T-1)*p*sizeof(double));
	S_xdy  = (double *)malloc((T-1)*p*sizeof(double));
	CS_dxdx= (double *)malloc((T-1)*p*p*sizeof(double));
	CS_dxdy= (double *)malloc((T-1)*p*sizeof(double));

    
	A     = (double *)malloc(p*p*sizeof(double));
	r     = (double *)malloc(p*sizeof(double));
	theta     = (double *)malloc(p*sizeof(double));
    
    
	Rn  = (double *)malloc((T-1)*p*p*sizeof(double));
	R  = (double *)malloc((T-1)*p*p*sizeof(double));
	rv  = (double *)malloc((T-1)*p*sizeof(double));
	rvn  = (double *)malloc((T-1)*p*sizeof(double));

	b0     = (double *)malloc(T*p*sizeof(double));
	g     = (double *)malloc(p*sizeof(double));
	vTmp1     = (double *)malloc(p*sizeof(double));
	vTmp2     = (double *)malloc(p*sizeof(double));
	mTmp     = (double *)malloc(p*p*sizeof(double));
    
    for(k=0;k<p;k++){
        for(i=0;i<N;i++){
            for(j=0;j<T-1;j++){
                dx[k*n1+(T-1)*i+j] = x[k*n+T*i+j+1]-x[k*n+T*i+j];
            }
        }
    }
    for(i=0;i<N;i++){
        for(j=0;j<T-1;j++){
            dy[(T-1)*i+j] = y[T*i+j+1]-y[T*i+j];
        }
    }
    
    //R[n*(l*p+m)+row]
    for(i=0;i<T-1;i++){
        for(j=0;j<p;j++){
            for(k=0;k<p;k++){
                tmp1=0; tmp2=0; tmp3=0;
                for(l=0;l<N;l++){
                    tmp1 = tmp1 + x[j*n+l*T+i+1]*x[k*n+l*T+i+1];
                    tmp2 = tmp2 + dx[j*n1+l*(T-1)+i]*dx[k*n1+l*(T-1)+i];
                    tmp3 = tmp3 + x[j*n+l*T+i+1]*dx[k*n1+l*(T-1)+i];
                }
                S_xx[(T-1)*(k*p+j)+i] = tmp1/N;
                S_dxdx[(T-1)*(k*p+j)+i] = tmp2/N;
                S_xdx[(T-1)*(k*p+j)+i] = tmp3/N;
            }
        }
    }

    
    for(i=0;i<T-1;i++){
        for(j=0;j<p;j++){
            tmp1=0; tmp2=0;
            for(l=0;l<N;l++){
                tmp1 = tmp1 + dx[j*n1+l*(T-1)+i]*dy[l*(T-1)+i];
                tmp2 = tmp2 + x[j*n+l*T+i+1]*dy[l*(T-1)+i];
            }
            S_dxdy[(T-1)*j+i] = tmp1/N;
            S_xdy[(T-1)*j+i] = tmp2/N;
        }
    }

    invshift(S_dxdx,T-1,p*p,Rn);
    cumsum(Rn,T-1,p*p,R);
    invshift(R,T-1,p*p,CS_dxdx);

    invshift(S_dxdy,T-1,p,rvn);
    cumsum(rvn,T-1,p,rv);
    invshift(rv,T-1,p,CS_dxdy);

    
	fillzero(b0,T*p);
    errn=1e+10; 
    i = 0;
    while(errn>*XTol && i<*maxIter){
        i=i+1;
        for(t=1;t<=T;t++){
            if(t==1){
                /* 
                theta = reshape(d0,p,T); % T by p
                r = CS_dxdy(1,:)' - reshape(S_xdx(1,:),p,p)'*theta(:,2);
                for j = 3:T
                    r = r - reshape(S_dxdx(j-1,:),p,p)*sum(theta(:,2:j-1),2) - reshape(S_xdx(j-1,:),p,p)'*theta(:,j);
                end
                d1(1:p) = reshape(CS_dxdx(1,:),p,p)\r; 
                */
                getrow(CS_dxdy,T-1,p,1,0,r);
                getrow(S_xdx,T-1,p,p,0,mTmp);
                mtrans(mTmp,p,p,A);
                getrow(b0,T,p,1,1,vTmp1);
                mtimesv(A,p,p,vTmp1,vTmp2);
                vsub(r,p,vTmp2);
                for(j=3;j<=T;j++){
                    getrow(S_dxdx,T-1,p,p,j-2,A);
                    fillzero(theta,p);
                    for(k=2;k<=j-1;k++){
                        getrow(b0,T,p,1,k-1,vTmp1);
                        vadd(theta,p,vTmp1);
                    }
                    mtimesv(A,p,p,theta,vTmp2);
                    vsub(r,p,vTmp2);
                    
                    getrow(S_xdx,T-1,p,p,j-2,mTmp);
                    mtrans(mTmp,p,p,A);
                    getrow(b0,T,p,1,j-1,vTmp1);
                    mtimesv(A,p,p,vTmp1,vTmp2);
                    vsub(r,p,vTmp2);
                }
                                
                getrow(CS_dxdx,T-1,p,p,0,A);
                gaussj(A, p, r, 1);
                for(k=0;k<p;k++)
                    b[k*T] = r[k];
                
            }else{
                if(t==T){
                    /* 
                    g = S_xdy(T-1,:)' - reshape(S_xdx(T-1,:),p,p)*sum(theta(:,1:T-1),2);
                    R = reshape(S_xx(T-1,:),p,p);
                    */
                    getrow(S_xdx,T-1,p,p,T-2,A);
                    fillzero(theta,p);
                    for(k=1;k<=T-1;k++){
                        getrow(b,T,p,1,k-1,vTmp1);
                        vadd(theta,p,vTmp1);
                    }
                    mtimesv(A,p,p,theta,vTmp1);
                    getrow(S_xdy,T-1,p,1,T-2,g);
                    vsub(g,p,vTmp1); //g
                    getrow(S_xx,T-1,p,p,T-2,A);//R
                    
                }else{
                    /*
                    g = S_xdy(t-1,:)' - reshape(S_xdx(t-1,:),p,p)*sum(theta(:,1:t-1),2) + CS_dxdy(t,:)';
                    for j=t+1:T
                        g = g  - reshape(S_dxdx(j-1,:),p,p)*(sum(theta(:,1:t-1),2)+sum(theta(:,t+1:j-1),2)) ...
                            - reshape(S_xdx(j-1,:),p,p)'*theta(:,j); 
                    end */
                    getrow(S_xdy,T-1,p,1,t-2,g);
                    fillzero(theta,p);
                    for(k=1;k<=t-1;k++){       //sum(theta(:,1:t-1),2)
                        getrow(b,T,p,1,k-1,vTmp1);
                        vadd(theta,p,vTmp1);
                    }
                    getrow(S_xdx,T-1,p,p,t-2,A); //S_xdx(t-1,:)
                    mtimesv(A,p,p,theta,vTmp1); //
                    vsub(g,p,vTmp1);
                    getrow(CS_dxdy,T-1,p,1,t-1,vTmp1);
                    vadd(g,p,vTmp1);
                    for(m=t+1;m<=T;m++){
                        fillzero(theta,p);
                        for(k=1;k<=t-1;k++){
                            getrow(b,T,p,1,k-1,vTmp1);
                            vadd(theta,p,vTmp1);
                        }
                        if(m>t+1){
                            for(k=t+1;k<=m-1;k++){
                                getrow(b,T,p,1,k-1,vTmp1);
                                vadd(theta,p,vTmp1);
                            }
                        }
                        getrow(S_dxdx,T-1,p,p,m-2,mTmp);
                        mtimesv(mTmp,p,p,theta,vTmp1);
                        vsub(g,p,vTmp1);
                        
                        getrow(S_xdx,T-1,p,p,m-2,mTmp);
                        mtrans(mTmp,p,p,A);
                        getrow(b,T,p,1,m-1,theta);
                        mtimesv(A,p,p,theta,vTmp1);
                        vsub(g,p,vTmp1);
                    }
                    /* R = reshape(S_xx(t-1,:)+ CS_dxdx(t,:),p,p);*/    
                    getrow(S_xx,T-1,p,p,t-2,mTmp);
                    getrow(CS_dxdx,T-1,p,p,t-1,A);
                    madd(A,p,p,mTmp);
                }
                dTmp = norm2(g,p);
                if(dTmp <= (*lambda)*weight[t-2]){
                    for(k=0;k<p;k++)
                        b[k*T+t-1]=0;
                }else{
                    for(k=0;k<p;k++)
                        g[k]=-g[k];
                /* Equivalent Matlab code:
                 * [gam,fval,exitflag] = fminbnd(@myfun2,0,1e10,[],reshape(Rn(nn,:),p,p),lambda,g(nn,:)'); */
                    lam = (*lambda)*weight[t-2];
                    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,A,p,lam,g,func);
                    exitflag = brent(ax,bx,cx,A,p,lam,g,func,*XTol,&xmin);
                /* Equivalent Matlab code:
                 * d1(nn,:) = -gam*((gam*reshape(Rn(nn,:),p,p) + lambda^2/2*eye(p))\g(nn,:)')'; */    
                    if(exitflag==1){
                        for(k=0;k<p*p;k++)
                            A[k] = xmin*A[k];
                        for(k=0;k<p;k++)
                            A[k*p+k] = A[k*p+k]+lam*lam/2;
                        gaussj(A,p,g,1);
                        for(k=0;k<p;k++)
                            b[k*T+t-1]=-xmin*g[k];
                    }else{
                        for(k=0;k<p;k++)
                            b[k*T+t-1]=0;
                    }
                }
            }
        }
        for(k=0;k<T*p;k++){
            err[k]=b[k]-b0[k];
            b0[k]=b[k];
        }
        errn = norm2(err,T*p);
    }
    
    free(err);
    free(dx);
    free(dy);
    free(S_xx);
    free(S_dxdx);
    free(S_xdx);
    free(S_dxdy);
    free(S_xdy);
    free(CS_dxdx);
    free(CS_dxdy);
    free(A);
    free(r);
    free(theta);
    free(Rn);
    free(R);
    free(rv);
    free(rvn);

    free(b0);
    free(g);
    free(vTmp1);
    free(vTmp2);
    free(mTmp);

}
        
void mexFunction(mwSize nlhs, mxArray *plhs[], mwSize nrhs, const mxArray *prhs[])
{
    #define b_out       plhs[0]
    #define y_in        prhs[0]
    #define x_in        prhs[1]
    #define lambda_in   prhs[2]
    #define weight_in   prhs[3]
    #define XTol_in     prhs[4]
    #define maxIter_in  prhs[5]
    
//declare variables
    const mwSize *dims;
    double *y,*x,*lambda,*weight,*XTol,*maxIter,*b; // inputs and outputs
    int n,p,N,T;

    if (nrhs < 6) {
        mexErrMsgIdAndTxt( "MATLAB:mxsetdimensions:minrhs",
                "At least 6 input arguments required.");
    }
    if(nlhs > 1){
        mexErrMsgIdAndTxt( "MATLAB:mxsetdimensions:maxlhs",
                "Too many output arguments.");
    }

//figure out dimensions
    dims = mxGetDimensions(x_in);
    n = (mwSize)dims[0]; p = (mwSize)dims[1];
    dims = mxGetDimensions(weight_in);
    T = (mwSize)dims[0]; T = T+1;
    N = n/T; 
    
    //mexPrintf("Hello,N=%d, T=%d, p=%d\n",N,T,p);

    b_out = mxCreateDoubleMatrix(T,p,mxREAL);
    
//associate pointers
    y       = mxGetPr(y_in);
    x       = mxGetPr(x_in);
    lambda  = mxGetPr(lambda_in);
    weight  = mxGetPr(weight_in);
    XTol    = mxGetPr(XTol_in);
    maxIter = mxGetPr(maxIter_in);
//associate outputs
    b = mxGetPr(b_out);

    wmplsbcd(y,x,N,T,p,lambda,weight,XTol,maxIter,b);

    return;
}
