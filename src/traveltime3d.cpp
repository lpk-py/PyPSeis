//FDRIMES package	

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>

#include <unistd.h>
#include <time.h>
#include "fdtimes.h"
//#include <omp.h>


#define STORE_STDMODEL 1
//static float gen_model(int nx, int ny, int nz, float h, float *buf);
//static int ABORT(char *m) { printf("%s\n",m); return 1; }

float ttbufpw1[35*45*55*10];//Traveltime Table of P-wave

float sp[35][45][55];//Slowness Model of P, take 121*121*55 grids as an example
//float ss[nx][ny][nz];//Slowness Model of S
int main()
/* args unused */
{
    int nx,ny,nz,nr;//nr is the receiver number
    int sx,sy,sz;
    float h,eps_init,pp;//sx,sy,sz,
    float *hsbufp,*tbufp,*ptr;//*tbufs,*hsbufs are temporal parameters for subroutine

    int msg,status,i,ii,iii;
    int bytes,n;//fd,
    FILE *fd;
    float vel0,vel1,vel2,vel3,vel4,vel5,vel6,vel7;
    int nz1,nz2,nz3,nz4,nz5,nz6,nz7,nz8,x,y,z;

    int indw1[]={20,23,26,29,32,35,38,41,44,47};//receiver positions in the grid model

    char filename[64];

    //printf("\nUsage: if model file entered below does not exist, this\n");
    //printf("program will create a simple two-layer model for you.\n\n");

    /* dimensioning the problem */
    nx=35;//grid number in x-axis
    ny=45;//grid number in y-axis
    nz=55;//grid number in z-axis
    nr=10;//receiver number

    bytes=nx*ny*nz*sizeof(float);
    //hsbuf=(float *)malloc(bytes);
    tbufp=(float *)malloc(bytes);
    //tbufs=(float *)malloc(bytes);
    //printf("Grid spacing h in meters : ");

    h=5;//spacing in x, y and z directions
    //if(h<=0.0) return ABORT("nonpositive spacing.");

    /* building slowness model */
    printf("Building a velocity model.\n");
    //gen_model(nx,ny,nz,h,hsbuf);
    vel0=2308;// velocity values
    vel1=2517;
    vel2=2798;
    vel3=3008;
    vel4=3230;
    vel5=3329;
    vel6=3403;
    vel7=3342;

    nz1=5;//layer depth in the grid model
    nz2=6;
    nz3=8;
    nz4=9;
    nz5=17;
    nz6=23;
    nz7=47;


    for(x=0;x<nx;x++)
    {	for(y=0;y<ny;y++)
        {
            for(z=0;z<nz1;z++) 	{sp[x][y][z]=1.0/vel0;}//ss[x][y][z]=1.0/(vel0/1.67);
            for(z=nz1;z<nz2;z++) 	{sp[x][y][z]=1.0/vel1;}//ss[x][y][z]=1.0/(vel1/1.67);
            for(z=nz2;z<nz3;z++) {sp[x][y][z]=1.0/vel2;}//ss[x][y][z]=1.0/(vel2/1.67);
            for(z=nz3;z<nz4;z++) {sp[x][y][z]=1.0/vel3;}
            for(z=nz4;z<nz5;z++) {sp[x][y][z]=1.0/vel4;}
            for(z=nz5;z<nz6;z++) {sp[x][y][z]=1.0/vel5;}
            for(z=nz6;z<nz7;z++) {sp[x][y][z]=1.0/vel6;}
            for(z=nz7;z<nz;z++) {sp[x][y][z]=1.0/vel7;}
        }
    }


    hsbufp=sp[0][0];
    //hsbufs=ss[0][0];

    sprintf(filename,"smodelp_%dx%dx%d_h%g.txt",nx,ny,nz,h);
    fd=fopen(filename,"w+");
    for(n=0;n<nx*ny*nz;n++) fprintf(fd,"%f\n",hsbufp[n]);

    printf("INFO: the p model was saved as file './%s'\n",filename);
    fclose(fd);

    //printf("nx=%d\n",nx);
    //for(n=0;n<nx*ny;n++) printf("velocity model=%f\n",hsbuf[n]);

    /* the product h*s must be passed to time_2d */
    for(n=0,ptr=hsbufp;n<nx*ny*nz;n++) *ptr++*=h;
    //for(n=0,ptr=hsbufs;n<nx*ny*nz;n++) *ptr++*=h;

    /* source coordinates
    //printf("Source coordinates (in meters) sx sy : ");
    sx=50;
    sy=100;
    sz=50;*/
    /* src coords. must be passed in h unit (indices)
    sx/=h;
    sy/=h;
    sz/=h;
*/
    /* timefield initialization */
    /* This is not required if source is located within model boundaries */
    /* Yet, if this is not true, tbuf will be treated as an initialized  */
    /* time map ("multiple source" case). In this case, time_2d() will   */
    /* fail and complain if you provide it with a uniform map (whatever  */
    /* the time value it contains), because it would be left unchanged.  */
    /* The line below thus ensures that an error will be diagnosed if a  */
    /* "wrong" source location is inadvertently entered.                 */
    /* You would not be informed of this if tbuf was left uninitialized  */
    /* (and thus potentially not uniform).                               */
    for(n=0,ptr=tbufp;n<nx*ny*nz;n++) *ptr++=0.0;
    //for(n=0,ptr=tbufs;n<nx*ny*nz;n++) *ptr++=0.0;

    /* other params */
    msg=2;          /* smalltalk mode */
    eps_init=0.001; /* my advised value: simply always use this... */

    /* calling time_3d() */
    ii=0;
    iii=0;
    pp=0;
    int t0=time(NULL);

    //#pragma omp parallel for
    for (sx=0;sx<5;sx++)//source positions in the grid model
    {
        //#pragma omp parallel for
        for (sy=0;sy<5;sy++)
        {

            //omp_set_num_threads(6);
            //#pragma omp parallel for
            for (sz=0;sz<5;sz++)
            {

                /* source coordinates */
                //printf("Source coordinates (in meters) sx sy : ");
                //sx=50;
                //sy=100;
                //sz=50;
                /* src coords. must be passed in h unit (indices)
    sx/=h;
    sy/=h;
    sz/=h;*/
                //#pragma omp critical
                status=time_3d(hsbufp,tbufp,nx,ny,nz,sx,sy,sz,eps_init,msg);//P
                //status=time_3d(hsbufs,tbufs,nx,ny,nz,sx,sy,sz,eps_init,msg);//S
                //if(status!=0) return ABORT("something wrong occurred in time_3d().");

                //omp_set_num_threads(6);
                //#pragma omp parallel for

                for (i=0;i<nr;i++)

                {
                    ttbufpw1[iii+i]=tbufp[indw1[i]];//

                }
                //
                //ttbuf is the traveltime to the receivers
                //tbuf is the traveltime of a certain source to all grids



                ii=ii+1;
                pp=ii/nx*ny*nz*100.0;//
                iii=nr*ii;

                printf("%f calculation is completed.\n ",pp);//calculation percentage
                //printf("ii is %d and iii is %d \n",ii,iii);


            }
        }
    }
    int t1=time(NULL);
    /* storing results */


    printf("Store timefield (values in seconds) of P-wave of w1 on file : ");
    //scanf("%s",filename);
    sprintf(filename,"traveltime_test.txt");
    fd=fopen(filename,"w+");
    for(n=0;n<nx*ny*nz*nr;n++) fprintf(fd,"%f\n",ttbufpw1[n]);
    fclose(fd);
    //nx*ny*nz
    /*
    printf("Store timefield (values in seconds) of S-wave on file : ");
    scanf("%s",filename);
    fd=fopen(filename,"w+");
    for(n=0;n<nx*ny*nz*nr;n++) fprintf(fd,"%f\n",ttbufs[n]);
    fclose(fd);
    */
    printf("total calculation time:%d s\n",t1-t0);//total calculation time
    /* ending */
    free(tbufp);
    //free(tbufs);
    //free(hsbuf);
    //
    free(ttbufpw1);
    return 0;

}




