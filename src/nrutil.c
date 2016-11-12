#include "stdafx.h"

VECTOR Vector( long n ){
	// allocate a double vector and set default values to 0
	VECTOR vector = (VECTOR) malloc( n * sizeof(double) );
	if( !vector ) return NULL;
	else memset( vector, 0, n * sizeof(double) );

	return vector;
}

FVECTOR FVector( long n ){
	// allocate a double vector and set default values to 0
	FVECTOR fvector = (FVECTOR) malloc( n * sizeof(float) );
	if( !fvector ) return NULL;
	else memset( fvector, 0, n * sizeof(float) );

	return fvector;
}
IVECTOR IVector( long n ){
	// allocate a double vector and set default values to 0
	IVECTOR ivector = (IVECTOR) malloc( n * sizeof(int) );
	if( !ivector ) return NULL;
	else memset( ivector, 0, n * sizeof(int) );

	return ivector;
}

MATRIX Matrix( long n, long m ){
	int i;
	// allocate a double matrix with subscript range m[n x m]
	MATRIX matrix;
   
	/* allocate pointers to rows */
	matrix = (MATRIX) malloc( n * sizeof(double*) );
	if( !matrix ) return NULL;
   
	/* allocate rows and set pointers to them */
	matrix[0] = (double*) malloc( n * m * sizeof(double) );
	if( !matrix[0] ) return NULL; 

   	memset( matrix[0], 0, n * m * sizeof(double) );
   
	for( i = 0; i < n; i++ ) matrix[i] = matrix[0] + m*i;

	/* return pointer to array of pointers to rows */
	return matrix;
}

IMATRIX IMatrix( long n,  long m ){
	int i;
	// allocate a double matrix with subscript range m[n x m]
	IMATRIX imatrix;
   
	/* allocate pointers to rows */
	imatrix = (IMATRIX) malloc( n * sizeof(int*) );
	if( ! imatrix ) return NULL;
   
	/* allocate rows and set pointers to them */
	imatrix[0] = (int*) malloc( n * m * sizeof(int) );
	if( ! imatrix[0] ) return NULL; 

   	memset( imatrix[0], 0, n * m * sizeof(int) );
   
	for( i = 0; i < n; i++ ) imatrix[i] = imatrix[0] + m*i;

	/* return pointer to array of pointers to rows */
	return imatrix;
}

FMATRIX FMatrix( long n, long m ){
	int i;
	// allocate a double matrix with subscript range m[n x m]
	FMATRIX fmatrix;
   
	/* allocate pointers to rows */
	fmatrix = (FMATRIX) malloc( n * sizeof(float*) );
	if( !fmatrix ) return NULL;
   
	/* allocate rows and set pointers to them */
	fmatrix[0] = (float*) malloc( n * m * sizeof(float) );
	if( !fmatrix[0] ) return NULL; 

   	memset( fmatrix[0], 0, n * m * sizeof(float) );
   
	for( i = 0; i < n; i++ ) fmatrix[i] = fmatrix[0] + m*i;

	/* return pointer to array of pointers to rows */
	return fmatrix;
}

void free_matrix( MATRIX matrix ){
	free( matrix[0] );
	free( matrix );
} 

void free_imatrix( IMATRIX imatrix ){
	free( imatrix[0] );
	free( imatrix );
} 

void free_fmatrix( FMATRIX fmatrix ){
	free( fmatrix[0] );
	free( fmatrix );
} 

void free_vector( VECTOR vector ){
	free( vector );
}
 
void free_fvector( FVECTOR fvector ){
	free( fvector );
}  

void free_ivector( IVECTOR ivector ){
	free( ivector );
} 
 
