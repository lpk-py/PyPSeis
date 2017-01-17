#include "array3d.h"

float*** make3DFloatArrays(float *array1D, int nx, int ny, int nz)
/************************************************************
Return 3D array that points to the proper locations in the
1D array compatible with FORTRAN.  z is the fastest index,
y is the next fastest, and x is the slowest.
************************************************************/
{

    float ***array3D;
    int iy, ix;

    array3D=  (float***) malloc( nx*sizeof(float**) );
    if (array3D==NULL) {
        fprintf(stderr,"Cannot allocated %d bytes!\n", nx*sizeof(float**) );
        return NULL;
    }

    for (ix=0;ix<nx;ix++) {
        array3D[ix]=   (float**) malloc(ny*sizeof(float*) );
        if (array3D[ix]==NULL) {
            fprintf(stderr,"Cannot allocated %d bytes!\n", ny*sizeof(float*) );
            return NULL;
        }
        for (iy=0;iy<ny;iy++)
            array3D[ix][iy]=  &(array1D[ix*ny*nz + iy*nz]) ;
    }

    return array3D;
}

void free3DFloatArrays(float ***array3D, int nx)
/************************************************************
Free extra memory allocated by make3DFloat().
************************************************************/
{
    int ix;

    for (ix=0;ix<nx;ix++)
        free(array3D[ix]);

    free(array3D);
}


int*** make3DIntArrays(int *array1D, int nx, int ny, int nz)
/************************************************************
Return 3D array that points to the proper locations in the
1D array compatible with FORTRAN.  z is the fastest index,
y is the next fastest, and x is the slowest.
************************************************************/
{

    int ***array3D;
    int iy, ix;

    array3D= (int***) malloc( nx*sizeof(int**) );
    if (array3D==NULL) {
        fprintf(stderr,"Cannot allocate %d bytes!\n", nx*sizeof(int**) );
        return NULL;
    }

    for (ix=0;ix<nx;ix++) {
        array3D[ix]=  (int**) malloc(ny*sizeof(int*) );
        if (array3D[ix]==NULL) {
            fprintf(stderr,"Cannot allocate %d bytes!\n", ny*sizeof(int*) );
            return NULL;
        }
        for (iy=0;iy<ny;iy++)
            array3D[ix][iy]=  &(array1D[ix*ny*nz + iy*nz]) ;
    }

    return array3D;
}

void free3DIntArrays(int ***array3D, int nx)
/************************************************************
Free extra memory allocated by make3DFloat().
************************************************************/
{
    int ix;

    for (ix=0;ix<nx;ix++)
        free(array3D[ix]);

    free(array3D);
}
