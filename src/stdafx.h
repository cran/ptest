// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#if !defined(AFX_STDAFX_H__2D7C85E9_6DCE_11D4_A6DC_00500447FC27__INCLUDED_)
#define AFX_STDAFX_H__2D7C85E9_6DCE_11D4_A6DC_00500447FC27__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

#define	TWOPI 6.283185307179586

#include <stdio.h>
#include "nrutil.h"
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h> 
#include <string.h>


//void GetX(double, double *,  double *, double [][3]);

double RegSS( double, double * ,  long n, double *, MATRIX);

void gaussj(double [3][3], double [3][1]);

//void LikelihoodRatio(double *,  double *, double *, double *);

// TODO: reference additional headers your program requires here    

//{{AFX_INSERT_LOCATION}}

// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__2D7C85E9_6DCE_11D4_A6DC_00500447FC27__INCLUDED_)

