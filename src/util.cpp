#include "util.h"
#pragma hdrstop
#include <memory>
#include <cstring>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define swapb(x,y) {byte tmp; tmp=(x); (x)=(y); (y)=tmp;}
#define PI 3.1415926535897

static int _swap=1;

int setswap(int swap)
{
    _swap=swap; return swap;
    
}


int swapi4(int x)
{
    
    byte* cbuf;
    byte  tem;
    
    cbuf=(byte*)&x;        /* assign address of input to char array */
    
    tem=cbuf[0]; cbuf[0]=cbuf[3]; cbuf[3]=tem;
    tem=cbuf[2]; cbuf[2]=cbuf[1]; cbuf[1]=tem;
    
    return x;
}

float swapf4(float x)
{
    char *sz = (char*)&x;
    char t=sz[0];sz[0]=sz[3],sz[3]=t;
    t=sz[1];sz[1]=sz[2],sz[2]=t;
    
    return x;
}

short swapi2(short x)
{
    
    byte* cbuf;
    byte  tem;
    
    cbuf=(byte*)&x;        /* assign address of input to char array */
    
    tem=cbuf[0]; cbuf[0]=cbuf[1]; cbuf[1]=tem;
    
    
    return x;
}



void ibm2ieee(float* input, int swap)
{
    typedef unsigned char byte;
    typedef unsigned long ulng;
    
    byte  *cbuf,expp,tem,sign;
    ulng  *umantis,expll,signl;
    long *mantis;
    int  shift;
    
    cbuf=(byte*)&input[0];        /* assign address of input to char array */
    umantis=(ulng*)&input[0];     /* two differnt points to the same spot  */
    mantis =(long*)&input[0];     /* signned & unsigned                    */
    
    if(swap)
    {
        /* now byte reverce for PC use if swap true */
        tem=cbuf[0]; cbuf[0]=cbuf[3]; cbuf[3]=tem;
        tem=cbuf[2]; cbuf[2]=cbuf[1]; cbuf[1]=tem;
    }
    
    /* start extraction information from number */
    
    expp=*mantis>>24;     /* get expo fro upper byte      */
    *mantis=(*mantis)<<8; /* shift off upper byte         */
    shift=1;              /* set a counter to 1           */
    while(*mantis>0 && shift<23) /* start of shifting data*/
    {
        *mantis=*mantis<<1;
        shift++;
    } /* shift until a 1 in msb */
    
    *mantis=*mantis<<1; /* need one more shift to get implied one bit */
    sign=expp & 0x80;   /* set sign to msb of exponent            */
    expp=expp & 0x7F;   /* kill sign bit                          */
    
    if(expp!=0)        /* don't do anymore if zero exponent       */
    {
        expp=expp-64;   /* compute what shift was (old exponent)*/
        *umantis=*umantis>>9; /* MOST IMPORTANT an UNSIGNED shift back down */
        expll=0x7F+(expp*4-shift); /* add in excess 172 */
        
        /* now mantissa is correctly aligned, now create the other two pairs */
        /* needed the extended sign word and the exponent word               */
        
        expll=expll<<23;        /* shift exponent up */
        
        /* combine them into a floating point IEEE format !     */
        
        if(sign) *umantis=expll | *mantis | 0x80000000;
        else     *umantis=expll | *mantis; /* set or don't set sign bit */
    }
}

void toibm(long *addr, int ll)
{
    unsigned long real;
    
    unsigned long lng;
    unsigned char *ch;
    
    unsigned int gain,sign,t;
    
    //    ch= (long *)&lng;
    ch= (unsigned char*)&lng;
    
    while(ll--)
    {
        if(real == *addr)
        {
            gain = (real>>23);
            sign = (gain & 0x100)>>1;
            gain = (gain & 0xff) - 0x80 +2;
            real = ((real& 0x7FFFFFL)| 0x800000L);
            
            real = (real& 0x7FFFFFL);
            real = (real | 0x800000L);
            
            t = gain & 3;
            if(t)
            {
                real >>= 4-t;
                gain +=  4-t;
            }
            gain = ((gain>>2) + 0x40) & 0x7f;
            gain |= sign;
            lng = real;
            ch[3] = ch[0];
            t = ch[1];
            ch[1] = ch[2];
            ch[2] = t;
            ch[0] = gain;
            real = lng;
        }
        *addr++ = real;
    }
}


void si4(char* buf, int nbyte, int i)
{
    union { byte h[4]; int x; } u;
    
    u.x=i;
    
    if(_swap) { swapb(u.h[0], u.h[3]); swapb(u.h[1], u.h[2]);}
    
    memcpy(buf+nbyte-1, &u.x, sizeof(int));
}

void si2(char* buf,int nbyte, short i)
{
    union { byte h[2]; short x; } u;
    
    u.x=i;
    
    if(_swap) swapb(u.h[0], u.h[1]);
    
    memcpy(buf+nbyte-1, &u.x, sizeof(short));
}

void si1(char* buf, int nbyte, char i)
{
    memcpy(buf+nbyte-1,&i,sizeof(char));
}


int i4(char* buf, int nbyte)
{
    int i;
    union { byte h[4]; long x; } u;
    
    memcpy(&i,buf+nbyte-1,sizeof(int));
    u.x=i;
    
    if(_swap)
    swapb(u.h[0], u.h[3]); swapb(u.h[1], u.h[2]);
    return u.x;
}

short i2(char* buf, int nbyte)
{
    short i;
    union { byte h[2]; int x; } u;
    
    memcpy(&i,buf+nbyte-1,sizeof(short));
    u.x=i;
    if(_swap)
    swapb(u.h[0], u.h[1]);
    
    return u.x;
}

char i1(char* buf, int nbyte)
{
    char i;
    
    memcpy(&i,buf+nbyte-1,sizeof(char));
    
    return i;
}

float s4(char *buf,int nbyte)
{
    float f;
    union { byte h[4]; float x; } u;
    
    memcpy(&f,buf+nbyte-1,sizeof(float));
    return f;
}

void swapCh4(char *ch)
{
    char szTemp;
    szTemp = ch[0]; ch[0] = ch[3]; ch[3] = szTemp;
    szTemp = ch[1]; ch[1] = ch[2]; ch[2] = szTemp;
}

void swapCh2(char *ch)
{
    char szTemp;
    szTemp = ch[0]; ch[0] = ch[1]; ch[1] = szTemp;
}

unsigned char EbcdicToAscii(unsigned char in_ascii_int)
{
    short ascii[96] =
    {
        32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
        48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
        64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
        80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
        96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
        112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127
    };
    
    short ebcdic[96] =
    {
        64, 79,127,123, 91,108, 80,125, 77, 93, 92, 78,107, 96, 75, 97,
        240,241,242,243,244,245,246,247,248,249,122, 94, 76,126,110,111,
        124,193,194,195,196,197,198,199,200,201,209,210,211,212,213,214,
        215,216,217,226,227,228,229,230,231,232,233, 74,224, 90, 95,108,
        64,129,130,131,132,133,134,135,136,137,145,146,147,148,149,150,
        151,152,153,162,163,164,165,166,167,168,169,192,106,208,161, 64
    };
    
    int j;
    for(j=0;j<96;j++){
        if (ebcdic[j]==in_ascii_int)
        {
            return ascii[j];
        }
    }
    return 32;
}

unsigned char AsciiToEbcdic(unsigned char in_ebcdic_int)
{
    int ascii[96] =
    {
        32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
        48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
        64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
        80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
        96, 97, 98, 99, 100,101,102,103,104,105,106,107,108,109,110,111,
        112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127
    };
    
    int ebcdic[96] =
    {
        64, 79,127,123, 91,108, 80,125, 77, 93, 92, 78,107, 96, 75, 97,
        240,241,242,243,244,245,246,247,248,249,122, 94, 76,126,110,111,
        124,193,194,195,196,197,198,199,200,201,209,210,211,212,213,214,
        215,216,217,226,227,228,229,230,231,232,233, 74,224, 90, 95,108,
        64,129,130,131,132,133,134,135,136,137,145,146,147,148,149,150,
        151,152,153,162,163,164,165,166,167,168,169,192,106,208,161, 64
    };
    
    int j;
    for(j=0;j<96;j++){
        if (ascii[j]==in_ebcdic_int){
            return ebcdic[j];
        }
    }
    return 32;
}


void float_to_ibm(int from[], int to[], int n, int endian)
/**********************************************************************
 float_to_ibm - convert between 32 bit IBM and IEEE floating numbers
 ***********************************************************************
 Input:
 from	   input vector
 n	   number of floats in vectors
 endian	   =0 for little endian machine, =1 for big endian machines
 
 Output:
 to	   output vector, can be same as input vector
 
 ***********************************************************************
 Notes:
 Up to 3 bits lost on IEEE -> IBM
 
 IBM -> IEEE may overflow or underflow, taken care of by
 substituting large number or zero
 
 Only integer shifting and masking are used.
 ***********************************************************************
 Credits:     CWP: Brian Sumner
 ***********************************************************************/
{
    register int fconv, fmant, i, t;
    
    for (i=0;i<n;++i) {
        fconv = from[i];
        if (fconv) {
            fmant = (0x007fffff & fconv) | 0x00800000;
            t = (int) ((0x7f800000 & fconv) >> 23) - 126;
            while (t & 0x3) { ++t; fmant >>= 1; }
            fconv = (0x80000000 & fconv) | (((t>>2) + 64) << 24) | fmant;
        }
        if(endian==0)
        fconv = (fconv<<24) | ((fconv>>24)&0xff) |
        ((fconv&0xff00)<<8) | ((fconv&0xff0000)>>8);
        
        to[i] = fconv;
    }
    return;
}

void ebasd( unsigned char* ascii, unsigned char* ebcd )
{
    return;
}




int swap(int x, int type)
{
    if (type == 2)
    return swapi2(x);
    else return swapi4(x);
}


void GetVectorMIN(vector<double> a, double &mindata)
{
    //求向量最小值
    mindata=a[0];
    int len=a.size();
    for(int i=1;i<len;i++)
    {
        if (a[i]<mindata)
            mindata=a[i];
    }
}


void findMaxMin(float a[], int n, float &max, float &min)
{
    int i;
    max = min = a[0];
    
    for (i = 1; i < n/2; ++i) {
        if (a[i*2] > a[i*2+1]) {
            if (a[i*2] > max)
            max = a[i*2];
            if (a[i*2+1] < min)
            min = a[i*2+1];
        } else {
            if (a[i*2+1] > max)
            max = a[i*2+1];
            if (a[i*2] < min)
            min = a[i*2];
        }
    }
    
    
    if (n%2 != 0) {
        max = (max >= a[n-1]) ? max : a[n-1];
        min = (min <= a[n-1]) ? min : a[n-1];
    }
    
    
}


int GetIndexOfMin(int *arry, int n)
{
    int i = 0, t = 0, temp;
    int *p;
    p = arry;
    
    for (i = 0; i <n - 1; i++)
    {
        if (*(p + i) >*(p + i + 1))
        t = i + 1;
        else
        {
            temp = *(p + i);
            *(p + i) = *(p + i + 1);
            *(p + i + 1) = temp;
        }
    }
    return t;
}



int GetIndexOfMin(double *arry, int n)
{
    int i=0,t=0,temp;//一定要把t初始化为0。要不如果第一个是最大的时候就会输出一个地址。
    double *p;
    p=arry;
    
    for(i=0;i <n-1;i++)/*两两比较*/              //这里数组越界了，改成了i<n-1
    {
        if(*(p+i) >*(p+i+1))/*如果前一个下标元素大于于后一个下标元素，
                             把后一个下标值赋给t*/
        t=i+1;
        else/*如果前一个下标元素小于后一个下标元素*/  //这里也改了一下
        {
            temp=*(p+i);
            *(p+i)=*(p+i+1);
            *(p+i+1)=temp;
        }
    }
    
    
    return t;/*返回t*/
}

void FindMaxMin(vector<double> pArr, int nLength, double &nMax, double &nMin)
{
    //    if (pArr == NULL || nLength <= 0)
    //        {
    //            cout << "输入有误！" << endl;
    //            return;
    //        }
    
    
    if (nLength == 2)//数组只有两个元素
    {
        if (pArr[0] > pArr[1])
        {
            nMax = pArr[0];
            nMin = pArr[1];
        }
        else
        {
            nMax = pArr[1];
            nMin = pArr[0];
        }
        return;
    }
    
    //初始赋值
    if (pArr[0] > pArr[1])
    {
        nMax = pArr[0];
        nMin = pArr[1];
    }
    else
    {
        nMax = pArr[1];
        nMin = pArr[0];
    }
    
    for (int i=2; i<nLength; i+=2)
    {
        if (i+1 < nLength && pArr[i] < pArr[i+1])
        {
            if (nMax < pArr[i+1])//将较大值和nMax比较
            {
                nMax = pArr[i+1];
            }
            
            if (nMin > pArr[i])//将较小值和nMin比较
            {
                nMin = pArr[i];
            }
        }
        else if (i+1 < nLength && pArr[i] > pArr[i+1])
        {
            if (nMax < pArr[i])//将较大值和nMax比较
            {
                nMax = pArr[i];
            }
            
            if (nMin > pArr[i+1])//将较小值和nMin比较
            {
                nMin = pArr[i+1];
            }
        }
        else//最后剩下一个元素
        {
            if (nMax < pArr[i])
            {
                nMax = pArr[i];
            }
            
            if (nMin > pArr[i])
            {
                nMin = pArr[i];
            }
        }
    }
}

void FindMaxMin(double *pArr, int nLength, double &nMax, double &nMin)
{
    //    if (pArr == NULL || nLength <= 0)
    //        {
    //            cout << "输入有误！" << endl;
    //            return;
    //        }
    
    
    if (nLength == 2)//数组只有两个元素
    {
        if (pArr[0] > pArr[1])
        {
            nMax = pArr[0];
            nMin = pArr[1];
        }
        else
        {
            nMax = pArr[1];
            nMin = pArr[0];
        }
        return;
    }
    
    //初始赋值
    if (pArr[0] > pArr[1])
    {
        nMax = pArr[0];
        nMin = pArr[1];
    }
    else
    {
        nMax = pArr[1];
        nMin = pArr[0];
    }
    
    for (int i=2; i<nLength; i+=2)
    {
        if (i+1 < nLength && pArr[i] < pArr[i+1])
        {
            if (nMax < pArr[i+1])//将较大值和nMax比较
            {
                nMax = pArr[i+1];
            }
            
            if (nMin > pArr[i])//将较小值和nMin比较
            {
                nMin = pArr[i];
            }
        }
        else if (i+1 < nLength && pArr[i] > pArr[i+1])
        {
            if (nMax < pArr[i])//将较大值和nMax比较
            {
                nMax = pArr[i];
            }
            
            if (nMin > pArr[i+1])//将较小值和nMin比较
            {
                nMin = pArr[i+1];
            }
        }
        else//最后剩下一个元素
        {
            if (nMax < pArr[i])
            {
                nMax = pArr[i];
            }
            
            if (nMin > pArr[i])
            {
                nMin = pArr[i];
            }
        }
    }
}

void FindMaxMin(float *pArr, int nLength, float &nMax, float &nMin)
{
    
    
    
    if (nLength == 2)//数组只有两个元素
    {
        if (fabs(pArr[0]) > fabs(pArr[1]))
        {
            nMax = fabs(pArr[0]);
            nMin = fabs(pArr[1]);
        }
        else
        {
            nMax = fabs(pArr[1]);
            nMin = fabs(pArr[0]);
        }
        return;
    }
    
    //初始赋值
    if (fabs(pArr[0]) > fabs(pArr[1]))
    {
        nMax = fabs(pArr[0]);
        nMin = fabs(pArr[1]);
    }
    else
    {
        nMax = fabs(pArr[1]);
        nMin = fabs(pArr[0]);
    }
    
    for (int i=2; i<nLength; i+=2)
    {
        if (i+1 < nLength && fabs(pArr[i]) < fabs(pArr[i+1]))
        {
            if (nMax < fabs(pArr[i+1]))//将较大值和nMax比较
            {
                nMax = fabs(pArr[i+1]);
            }
            
            if (nMin > fabs(pArr[i]))//将较小值和nMin比较
            {
                nMin = fabs(pArr[i]);
            }
        }
        else if (i+1 < nLength && fabs(pArr[i]) > fabs(pArr[i+1]))
        {
            if (nMax < fabs(pArr[i]))//将较大值和nMax比较
            {
                nMax = fabs(pArr[i]);
            }
            
            if (nMin > fabs(pArr[i+1]))//将较小值和nMin比较
            {
                nMin = fabs(pArr[i+1]);
            }
        }
        else//最后剩下一个元素
        {
            if (nMax < fabs(pArr[i]))
            {
                nMax = fabs(pArr[i]);
            }
            
            if (nMin > fabs(pArr[i]))
            {
                nMin = fabs(pArr[i]);
            }
        }
    }
}

void FindMaxMin(int *pArr, int nLength, int &nMax, int &nMin)
{
    
    if (nLength == 2)//数组只有两个元素
    {
        if (abs(pArr[0]) > abs(pArr[1]))
        {
            nMax = abs(pArr[0]);
            nMin = abs(pArr[1]);
        }
        else
        {
            nMax = abs(pArr[1]);
            nMin = abs(pArr[0]);
        }
        return;
    }
    
    //初始赋值
    if (abs(pArr[0]) > abs(pArr[1]))
    {
        nMax = abs(pArr[0]);
        nMin = abs(pArr[1]);
    }
    else
    {
        nMax = abs(pArr[1]);
        nMin = abs(pArr[0]);
    }
    
    for (int i=2; i<nLength; i+=2)
    {
        if (i+1 < nLength && abs(pArr[i]) < abs(pArr[i+1]))
        {
            if (nMax < abs(pArr[i+1]))//将较大值和nMax比较
            {
                nMax = abs(pArr[i+1]);
            }
            
            if (nMin > abs(pArr[i]))//将较小值和nMin比较
            {
                nMin = abs(pArr[i]);
            }
        }
        else if (i+1 < nLength && abs(pArr[i]) > abs(pArr[i+1]))
        {
            if (nMax < abs(pArr[i]))//将较大值和nMax比较
            {
                nMax = abs(pArr[i]);
            }
            
            if (nMin > abs(pArr[i+1]))//将较小值和nMin比较
            {
                nMin = abs(pArr[i+1]);
            }
        }
        else//最后剩下一个元素
        {
            if (nMax < abs(pArr[i]))
            {
                nMax = abs(pArr[i]);
            }
            
            if (nMin > abs(pArr[i]))
            {
                nMin = abs(pArr[i]);
            }
        }
    }
}


