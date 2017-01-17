#include <stdio.h>
#include <stdlib.h>

float*** make3DFloatArrays(float *array1D, int nx, int ny, int nz);
/************************************************************
Return 3D array that points to the proper locations in the
1D array compatible with FORTRAN.  z is the fastest index,
y is the next fastest, and x is the slowest.  If memory allocation
error occurs, NULL is returned.
************************************************************/

void free3DFloatArrays(float ***array3D, int nx);
/************************************************************
Free extra memory allocated by make3DFloat().
************************************************************/


int*** make3DIntArrays(int *array1D, int nx, int ny, int nz);
/************************************************************
Return 3D array that points to the proper locations in the
1D array compatible with FORTRAN.  z is the fastest index,
y is the next fastest, and x is the slowest.  If memory allocation
error occurs, NULL is returned.
************************************************************/

void free3DIntArrays(int ***array3D, int nx);
/************************************************************
Free extra memory allocated by make3DInt().
************************************************************/
