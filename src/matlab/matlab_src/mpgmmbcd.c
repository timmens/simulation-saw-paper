/*********************************************************************
 * mpgmmbcd.c
 * This program implements the BCD (Block-Coordinate Descent) algorithm of 
 * PGMM estimation of panel data model.
 * tht = mpgmmbcd(y,x,z,lambda,weight,XTol,maxIter);
 *  Output: 
 *      y: lhs variable, an NT by 1 column vector
 *  Inputs:
 *      x: rhs variables, an NT by p matrix
 *      z: IV, N(T-1) by q matrix, q>p
 *      lambda: tuning parameter, scalar
 *      weight: weighting vector, (n-1) by 1
 *      XTol: error tolerance, scalar
 *      maxIter: max num of iterations, scalar integer
 *  Subroutines needed:
 *      useful.c useful.h
 *  To compile mex program, enter the following statement in Maltab:
 *      >> mex mpgmmbcd.c useful.c
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

void wmpgmmbcd(double *y,double *x,double *z,int N,int T,int p,int q,double *lambda,double *weight,double *W,double *XTol,double *maxIter,double *b)
{   //wmpgmmbcd(y,x,z,N,T,p,q,lambda,weight,W,XTol,maxIter,b);
    double *dx,*dy,*S_zx,*S_zdx,*S_zdy,*M,*R,*U,*V,*L,*Q,*P,*D,*CR,*CV,*wt,*mt,*rt,*lt,*ut,*vt,*mWQ,*mWP;
    double *b0,*g,*vTmp1,*vTmp2,*mTmp1,*mTmp2,*A,*theta,*err,*Rn,*Rnc,*rn,*rnc; 
    double ax,bx,cx,fa,fb,fc,xmin,tmp1,tmp2,tmp3,errn,dTmp,lam;
    int n1,n,i,j,k,l,m,t,exitflag;
    
    ax = 1e-16; bx = 1e+6;

    n1 = (T-1)*N; n = N*T;

	err   = (double *)malloc(T*p*sizeof(double));
	dx   = (double *)malloc(n1*p*sizeof(double));
	dy   = (double *)malloc(n1*sizeof(double));
    
	S_zx   = (double *)malloc((T-1)*p*q*sizeof(double));
	S_zdx  = (double *)malloc((T-1)*p*q*sizeof(double));
	S_zdy  = (double *)malloc((T-1)*q*sizeof(double));
    
	M   = (double *)malloc((T-1)*p*p*sizeof(double));
	R   = (double *)malloc((T-1)*p*p*sizeof(double));
	U   = (double *)malloc((T-1)*p*p*sizeof(double));
	V   = (double *)malloc((T-1)*p*sizeof(double));
	L   = (double *)malloc((T-1)*p*sizeof(double));

	Q   = (double *)malloc(p*q*sizeof(double));
	P   = (double *)malloc(p*q*sizeof(double));
	D   = (double *)malloc(q*sizeof(double));
	CV   = (double *)malloc(p*(T-1)*sizeof(double));
	CR   = (double *)malloc((T-1)*p*p*sizeof(double));

	Rn   = (double *)malloc((T-1)*p*p*sizeof(double));
	Rnc  = (double *)malloc((T-1)*p*p*sizeof(double));
	rn   = (double *)malloc((T-1)*p*sizeof(double));
	rnc  = (double *)malloc((T-1)*p*sizeof(double));
    

	ut   = (double *)malloc(p*p*sizeof(double));
	wt   = (double *)malloc(q*q*sizeof(double));
	mt   = (double *)malloc(p*p*sizeof(double));
	rt   = (double *)malloc(p*p*sizeof(double));
	lt   = (double *)malloc(p*sizeof(double));
	vt   = (double *)malloc(p*sizeof(double));
   
	mWQ   = (double *)malloc(p*q*sizeof(double));
	mWP   = (double *)malloc(p*q*sizeof(double));
	A     = (double *)malloc(p*p*sizeof(double));
	theta = (double *)malloc(p*sizeof(double));
    
	b0   = (double *)malloc(p*T*sizeof(double));
	g   = (double *)malloc(p*sizeof(double));
	vTmp1   = (double *)malloc(p*sizeof(double));
	vTmp2   = (double *)malloc(p*sizeof(double));
	mTmp1   = (double *)malloc(p*q*sizeof(double));
	mTmp2   = (double *)malloc(p*p*sizeof(double));
    
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
    
    for(i=0;i<T-1;i++){
        for(j=0;j<q;j++){
            for(k=0;k<p;k++){
                tmp1=0; tmp2=0;
                for(l=0;l<N;l++){
                    tmp1 = tmp1 + z[j*n1+l*(T-1)+i]*x[k*n+l*T+i+1];
                    tmp2 = tmp2 + z[j*n1+l*(T-1)+i]*dx[k*n1+l*(T-1)+i];
                }
                S_zx[(T-1)*(k*q+j)+i] = tmp1/N;
                S_zdx[(T-1)*(k*q+j)+i] = tmp2/N;
            }
        }
    }

    
    for(i=0;i<T-1;i++){
        for(j=0;j<q;j++){
            tmp1=0;
            for(l=0;l<N;l++){
                tmp1 = tmp1 + z[j*n1+l*(T-1)+i]*dy[l*(T-1)+i];
            }
            S_zdy[(T-1)*j+i] = tmp1/N;
        }
    }
    
    for(i=0;i<T-1;i++){
        getrow(W,T-1,q,q,i,wt);
        getrow(S_zx,T-1,q,p,i,Q);
        getrow(S_zdx,T-1,q,p,i,P);
        getrow(S_zdy,T-1,q,1,i,D);
        // M = Q'WQ (p by p)
        mtrans(Q,q,p,mTmp1);
        mtimes(wt,q,q,Q,p,mWQ); //
        mtimes(mTmp1,p,q,mWQ,p,mt);
        insert_row(M,T-1,p,p,i,mt);
        
        // R = P'WP (p by p)
        mtrans(P,q,p,mTmp1);
        mtimes(wt,q,q,P,p,mWP);
        mtimes(mTmp1,p,q,mWP,p,rt);
        insert_row(R,T-1,p,p,i,rt);
        
        // U = P'WQ (p by p)
        mtimes(mTmp1,p,q,mWQ,p,ut);
        insert_row(U,T-1,p,p,i,ut);
        
        // V = P'WD (p by 1)
		mtrans(mWP,q,p,mTmp1);
        mtimes(mTmp1,p,q,D,1,vt);
        insert_row(V,T-1,p,1,i,vt);
        
        // L = Q'WD (p by 1)
		mtrans(mWQ,q,p,mTmp1);
        mtimes(mTmp1,p,q,D,1,lt);
        insert_row(L,T-1,p,1,i,lt);
    }

    invshift(R,T-1,p*p,Rn);
    cumsum(Rn,T-1,p*p,Rnc);
    invshift(Rnc,T-1,p*p,CR);

    invshift(V,T-1,p,rn);
    cumsum(rn,T-1,p,rnc);
    invshift(rnc,T-1,p,CV);

	fillzero(b0,T*p); fillzero(b,T*p);
    errn=1e+10; 
    i = 0;
    while(errn>*XTol && i<*maxIter){
        i=i+1;
        for(t=1;t<=T;t++){
            if(t==1){
                getrow(CV,T-1,p,1,t-1,g);
                getrow(U,T-1,p,p,t-1,ut);
                getrow(b0,T,p,1,t,theta);
                mtimesv(ut,p,p,theta,vTmp1);
                vsub(g,p,vTmp1);
                for(j=3;j<=T;j++){
                    getrow(R,T-1,p,p,j-2,rt);
                    fillzero(theta,p);
                    for(k=2;k<=j-1;k++){
                        getrow(b0,T,p,1,k-1,vTmp1);
                        vadd(theta,p,vTmp1);
                    }
                    mtimesv(rt,p,p,theta,vTmp1);
                    vsub(g,p,vTmp1);
                    
                    getrow(U,T-1,p,p,j-2,ut);
                    getrow(b0,T,p,1,j-1,vTmp1);
                    mtimesv(ut,p,p,vTmp1,vTmp2);
                    vsub(g,p,vTmp2);
                }
                                
                getrow(CR,T-1,p,p,0,A);
                gaussj(A, p, g, 1);
                for(k=0;k<p;k++)
                    b[k*T] = g[k];
                
            }else{
                if(t==T){
                    getrow(U,T-1,p,p,t-2,ut);
                    mtrans(ut,p,p,A);
                    fillzero(theta,p);
                    for(k=1;k<=T-1;k++){
                        getrow(b,T,p,1,k-1,vTmp1);
                        vadd(theta,p,vTmp1);
                    }
                    mtimesv(A,p,p,theta,vTmp1);
                    getrow(L,T-1,p,1,t-2,g);
                    vsub(g,p,vTmp1); //g
                    getrow(M,T-1,p,p,t-2,A);//R
                    
                }else{
                    getrow(L,T-1,p,1,t-2,g);
                    fillzero(theta,p);
                    for(k=1;k<=t-1;k++){       //sum(theta(:,1:t-1),2)
                        getrow(b,T,p,1,k-1,vTmp1);
                        vadd(theta,p,vTmp1);
                    }
                    getrow(U,T-1,p,p,t-2,ut); 
                    mtrans(ut,p,p,A);
                    mtimesv(A,p,p,theta,vTmp1); //
                    vsub(g,p,vTmp1);
                    getrow(CV,T-1,p,1,t-1,vTmp1);
                    vadd(g,p,vTmp1);
                    
                    for(m=t+1;m<=T;m++){
                        fillzero(theta,p);
                        for(k=1;k<=t-1;k++){
                            getrow(b,T,p,1,k-1,vTmp1);
                            vadd(theta,p,vTmp1);
                        }
                        if(m>t+1){
                            for(k=t+1;k<=m-1;k++){
                                getrow(b0,T,p,1,k-1,vTmp1);
                                vadd(theta,p,vTmp1);
                            }
                        }
                        getrow(R,T-1,p,p,m-2,rt);
                        mtimesv(rt,p,p,theta,vTmp1);
                        vsub(g,p,vTmp1);
                        
                        getrow(U,T-1,p,p,m-2,ut);
                        getrow(b0,T,p,1,m-1,theta);
                        mtimesv(ut,p,p,theta,vTmp1);
                        vsub(g,p,vTmp1);
                    }
                    /* R = reshape(S_xx(t-1,:)+ CS_dxdx(t,:),p,p);*/    
                    getrow(M,T-1,p,p,t-2,mTmp2);
                    getrow(CR,T-1,p,p,t-1,A);
                    madd(A,p,p,mTmp2);
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
                    exitflag=brent(ax,bx,cx,A,p,lam,g,func,*XTol,&xmin);
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
    
    free(S_zx);
    free(S_zdx);
    free(S_zdy);
    
    free(M);
    free(R);
    free(U);
    free(V);
    free(L);
    
    free(Q);
    free(P);
    free(D);
    free(CV);
    free(CR);

    free(ut);
    free(vt);
    free(wt);
    free(mt);
    free(rt);
    free(lt);

    free(mWQ);
    free(mWP);
    free(A);
    free(theta);

	free(Rn);
    free(Rnc);
	free(rn);
    free(rnc);

    free(b0);
    free(g);
    free(vTmp1);
    free(vTmp2);
    free(mTmp1);
    free(mTmp2);
}
        
void mexFunction(mwSize nlhs, mxArray *plhs[], mwSize nrhs, const mxArray *prhs[])
{
    //b = mpbcdgmm(Y,X,Z,lambda,weight,W,XTol,maxIter);
    #define b_out       plhs[0]
    #define y_in        prhs[0]
    #define x_in        prhs[1]
    #define z_in        prhs[2]
    #define lambda_in   prhs[3]
    #define weight_in   prhs[4]
    #define W_in        prhs[5]
    #define XTol_in     prhs[6]
    #define maxIter_in  prhs[7]

//declare variables
    const mwSize *dims;
    double *y,*x,*z,*lambda,*weight,*W,*XTol,*maxIter,*b; // inputs and outputs
    int n,p,N,T,q;
    
    if (nrhs < 8) {
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
    dims = mxGetDimensions(z_in);
    q = (mwSize)dims[1];
    dims = mxGetDimensions(weight_in);
    T = (mwSize)dims[0]; T = T+1;
    N = n/T; 
    
    b_out = mxCreateDoubleMatrix(T,p,mxREAL);
    
//associate pointers
    y       = mxGetPr(y_in);        // n by 1
    x       = mxGetPr(x_in);        // n by p
    z       = mxGetPr(z_in);        // n by q
    lambda  = mxGetPr(lambda_in);   // scalar
    weight  = mxGetPr(weight_in);   // (T-1) by 1
    W       = mxGetPr(W_in);        // (T-1) by (q*q)
    XTol    = mxGetPr(XTol_in);     // scalar
    maxIter = mxGetPr(maxIter_in);  // integer
//associate outputs
    b = mxGetPr(b_out);

    wmpgmmbcd(y,x,z,N,T,p,q,lambda,weight,W,XTol,maxIter,b);

    return;
}
