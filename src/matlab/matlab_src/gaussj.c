
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(double *a, int n, double *b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
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
		if (a[icol*n+icol] == 0.0) nrerror("gaussj: Singular Matrix");
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
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
#undef SWAP
#undef NRANSI
