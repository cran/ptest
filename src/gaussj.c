#include "stdafx.h"
#include "nrutil.h"

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

// Linear equation solution by Gauss-Jordan elimination (Numerical Recipes) a[1..n][1..n]
// is the input matrix. b[1..n][1..m] is input containing the m right-hand side vectors. On
// output, a is replaced by its matrix inverse, and b is replaced by the corresponding set 
// of solution vectors.

void gaussj(double a[3][3],double b[3][1]){
//	int *indxc,*indxr,*ipiv;
	int i,icol=1,irow=1,j,k,l,ll;
	double big,dum,pivinv,temp;
	int p,m;

	int indxc[3];
	int indxr[3];
	int ipiv[3];
	p=3;
	m=1;
	
//	indxc=ivector(1,p);
//	indxr=ivector(1,p);
//	ipiv=ivector(1,p);

	for (j=0;j<p;j++) ipiv[j]=0;
	for (i=0;i<p;i++) {
		big=0.0;
		for (j=0;j<p;j++)
			if (ipiv[j] != 1)
				for (k=0;k<p;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<p;l++) SWAP(a[irow][l],a[icol][l])
			for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		// 		if (a[icol][icol] == 0.0)
		// 			printf("gaussj: Singular Matrix");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<p;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<p;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<p;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=p-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<p;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
//	free_ivector(ipiv,1,p);
//	free_ivector(indxr,1,p);
//	free_ivector(indxc,1,p);
}
#undef SWAP
//#undef NRANSI


