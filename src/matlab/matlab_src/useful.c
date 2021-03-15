#include <math.h>
#include <stdlib.h>

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(double *a, int n, double *b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	indxc = (int *)malloc(n*sizeof(int));
	indxr = (int *)malloc(n*sizeof(int));
	ipiv  = (int *)malloc(n*sizeof(int));

	for (j=0;j<n;j++) ipiv[j]=0; // change to 0-based indexing
	for (i=0;i<n;i++) {         //change to 0-based indexing
		big=0.0;
		for (j=0;j<n;j++)//change to 0-based indexing
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {//change to 0-based indexing
					if (ipiv[k] == 0) {
						if (fabs(a[j*n+k]) >= big) { //change a[j][k] to a[j*n+k]
							big=fabs(a[j*n+k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow*n+l],a[icol*n+l])//change to 0-based indexing
			for (l=0;l<m;l++) SWAP(b[irow*m+l],b[icol*m+l])//change to 0-based indexing
		}
		indxr[i]=irow;
		indxc[i]=icol;
//		if (a[icol*n+icol] == 0.0) nrerror("gaussj: Singular Matrix");
		pivinv=1.0/a[icol*n+icol];
		a[icol*n+icol]=1.0;
		for (l=0;l<n;l++) a[icol*n+l] *= pivinv;//change to 0-based indexing
		for (l=0;l<m;l++) b[icol*m+l] *= pivinv;//change to 0-based indexing
		for (ll=0;ll<n;ll++)//change to 0-based indexing
			if (ll != icol) {
				dum=a[ll*n+icol];
				a[ll*n+icol]=0.0;
				for (l=0;l<n;l++) a[ll*n+l] -= a[icol*n+l]*dum;//change to 0-based indexing
				for (l=0;l<m;l++) b[ll*m+l] -= b[icol*m+l]*dum;//change to 0-based indexing
			}
	}
	for (l=n-1;l>=0;l--) {//change to 0-based indexing
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)//change to 0-based indexing
				SWAP(a[k*n+indxr[l]],a[k*n+indxc[l]]);
	}
	free(ipiv);
	free(indxr);
	free(indxc);
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, const double *R, const int p, const double lam, const double *g,
	double (*func)(double,const double*, const int, const double, const double*))
{
	double ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax,R,p,lam,g);
	*fb=(*func)(*bx,R,p,lam,g);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx,R,p,lam,g);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u,R,p,lam,g);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u,R,p,lam,g);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u,R,p,lam,g);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u,R,p,lam,g))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u,R,p,lam,g);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u,R,p,lam,g);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT

#define ITMAX 400
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

int brent(double ax, double bx, double cx, const double *R, const int p1, const double lam, const double *g, double (*f)(double,const double*, const int, const double, const double*), double tol,
	double *xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x,R,p1,lam,g);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return 1;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u,R,p1,lam,g);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
//	nrerror("Too many iterations in brent");
	*xmin=x;
	return 0;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT

void vcopy(const double *y,int n, double *x)
{
    int i;
    for(i=0;i<n;i++)
        x[i]=y[i];
}

void sumr(const double *x, int n, int p, int istart, int length, double *s)
{
    int i,j;
    for(j=0;j<p;j++){
        s[j]=0;
        for(i=istart;i<istart+length;i++)
            s[j]=s[j]+x[j*n+i];
    }
}
    
void sumc(const double *x, int n, int p, int jstart, int length, double *s)
{
    int i,j;
    for(i=0;i<n;i++){
        s[i]=0;
        for(j=jstart;j<(jstart+length);j++)
            s[i]=s[i]+x[j*n+i];
    }
}

void summ(const double *x, int n, int p, int dim, double *s)
{
    int i,j;
    if(dim==1){
        for(j=0;j<p;j++){
            s[j]=0;
            for(i=0;i<n;i++)
                s[j]=s[j]+x[n*j+i];
        }
    }
    if(dim==2){
        for(i=0;i<n;i++){
            s[i]=0;
            for(j=0;j<p;j++)
                s[i]=s[i]+x[n*j+i];
        }
    }
}

void cumsum(double *x, int n, int p, double *ans)
{
    double sum;
    int i, j;

    for(j=0; j<p; j++) { 
        sum = 0.0;
        for( i=0 ; i<n ; i++) {
        	sum = sum+x[j*n+i];
            ans[j*n+i] = sum;
        }
    }
}

void invshift(double *x, int n, int p, double *y)
{
    int i, j;
    for(j=0;j<p;j++){
        for(i=0;i<n;i++){
            y[j*n+i] = x[j*n+(n-i-1)];
        }
    }
}

void fillzero(double *x,int n)
{
    int i;
    for(i=0;i<n;i++)
        x[i]=0;
}

void fillv(double *x,int n,double v)
{
    int i;
    for(i=0;i<n;i++)
        x[i]=v;
}

void dotmult(double *A, int n, int p, double *B)
{
    int i;
    for(i=0;i<n*p;i++)
        B[i]=A[i]*B[i];
}

double dotprod(double *x, double n, double *y)
{
    int i;
    double s;
    s = 0;
    for(i=0;i<n;i++)
        s=s+x[i]*y[i];
    return s;
}
/*
void getrow(double *b,int n, int p, int row, double *q)
{
    int i;
    for(i=0;i<p;i++)
        q[i] = b[n*i+row];
}*/
void getrow(double *R,int n,int p,int q,int k,double *A)
{// copy the k-th row of R, which is n-by-pq, to a q-by-p matrix A
    // set q=1 if copy to a vector.
    //A: q by p
    int l,m;
    for(l=0;l<p;l++){
        for(m=0;m<q;m++)
            A[l*q+m] = R[n*(l*q+m)+k];
    }
}
void insert_row(double *R,int n,int p,int q,int k,const double *A)
{// insert the matrix A, which is q-by-p, to the k-th row of R, 
 //   which is n-by-pq;
    int l,m;
    for(l=0;l<p;l++){
        for(m=0;m<q;m++)
            R[n*(l*q+m)+k] = A[l*q+m];
    }
}

void getrow2(double *Rn,int n,int p,int row,double *A)
{
    int l,m;
    for(l=0;l<p;l++){
        for(m=0;m<p;m++)
            A[l*p+m] = Rn[n*(l*p+m)+row];
    }
}

double norm2(double *x, int n)
{
    int i;
    double sum=0;
    for(i=0;i<n;i++)
        sum = sum + pow(x[i],2);
    return sqrt(sum);
}

void mtimesv(const double *a, int n, int p, const double *b, double *c)
{
    int i,j;
    double sum;
    for(i=0;i<n;i++){
        sum = 0;
        for(j=0;j<p;j++)
            sum = sum + a[j*n+i]*b[j];
        c[i]=sum;
    }
    
}

void mtimes(const double *a, int n, int p, const double *b, int q, double *c)
{
    int i,j,k;
    double sum;
    for(i=0;i<n;i++){
        for(j=0;j<q;j++){
            sum = 0;
            for(k=0;k<p;k++)
                sum = sum + a[k*n+i]*b[j*p+k];
            c[j*n+i]=sum;
        }
    }
}

void mtrans(const double *a, int n, int p, double *b)
{
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<p;j++)
            b[i*p+j] = a[j*n+i];
    }
    
}

void madd(double *a, int n, int p, const double *b)
{
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<p;j++)
        a[i*p+j]=a[i*p+j]+b[i*p+j];
    }
}


void vadd(double *a, int n, const double *b)
{
    int i;
    for(i=0;i<n;i++){
        a[i]=a[i]+b[i];
    }
}

void vsub(double *a, int n, const double *b)
{
    int i;
    for(i=0;i<n;i++){
        a[i]=a[i]-b[i];
    }
    
}

void sandwitch(double *a, int n, int p, double *b, double *c, int q, double *d)
{
	// calculate d=a'bc, where a is n-by-p matrix, b is n-by-n matrix, c is n-by-q, d is p-by-q. 
    double *ab,*at;

	ab   = (double *)malloc(n*p*sizeof(double));
	at   = (double *)malloc(n*p*sizeof(double));

	mtrans(a,n,p,at);
	mtimes(at,p,n,b,n,ab);
	mtimes(ab,p,n,c,q,d);

	free(ab);
	free(at);
    
}

void getblock(double *R,int M,int N,int m,int n, int nrow, int ncol,double *A)
{
	int i,j;
	for(i=0;i<nrow;i++){
		for(j=0;j<ncol;j++){
			A[j*nrow+i] = R[(n+j)*M+m+i];
		}
	}
}

void putblock(double *R,int M,int N,int m,int n, int nrow, int ncol,double *A)
{
	int i,j;
	for(i=0;i<nrow;i++){
		for(j=0;j<ncol;j++){
			R[(n+j)*M+m+i] = A[j*nrow+i];
		}
	}
}

