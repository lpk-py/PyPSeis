#ifndef FDTIMES
#define FDTIMES

/* documentation : see heading comments in source codes */

int time_2d(float *hs, float *t, int nx, int ny,
            float xs, float ys, float eps_init, int messages);

int time_3d(float *hs, float *t, int nx, int ny, int nz,
            float xs, float ys, float zs, float eps_init, int messages);

/* Public routines time_2d_() and time_3d_() are not prototyped, as    */
/* they are intended to be called only by programs written in Fortran. */

#endif /*!FDTIMES*/
