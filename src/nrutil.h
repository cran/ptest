#ifndef _UTILS_H_
#define _UTILS_H_
   
typedef double*  VECTOR;
typedef float*  FVECTOR;
typedef double** MATRIX;
typedef int** IMATRIX;
typedef float** FMATRIX;

typedef int*  IVECTOR;


#define SQR(a) (a == 0.0 ? 0.0 : (a)*(a))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

VECTOR Vector( long n );
FVECTOR FVector(long n);
IVECTOR IVector( long n );

MATRIX Matrix( long n, long m );
IMATRIX IMatrix( long n, long m );
FMATRIX FMatrix( long n, long m );

void free_matrix( MATRIX hMatrix );
void free_imatrix(IMATRIX IMatrix );
void free_fmatrix( FMATRIX FMatrix );
void free_vector( VECTOR hVector );
void free_fvector(FVECTOR FVector );
void free_ivector(IVECTOR IVector );

#endif /* _UTILS_H_ */

