#ifndef UTIL_H
#define UTIL_H

#include<string>
#include<vector>
#include<list>


using namespace std;

typedef unsigned char byte;


void ibm2ieee(float* input, int swap);
void toibm   (long *addr, int ll);

unsigned char ebasc(unsigned char ascii);
void          ebasd(unsigned char* ascii, unsigned char* ebcd);
void          asebd(char* ebcd, char* ascii);
char          asebc(char ascii);

int setswap(int swap);

int   i4(char* buf, int nbyte);
short i2(char* buf, int nbyte);
char  i1(char* buf, int nbyte);
void f4(char *buf,int nbute);

void si4(char* buf, int nbyte,  int i);
void si2(char* buf,int nbyte, short i);
void si1(char* buf, int nbyte, char i);


int   swapi4(int   x);
float swapf4(float x);
short swapi2(short x);
int swap (int x, int type);

float s4(char *buf,int nbyte);

void swapCh4(char *);
void swapCh2(char *);

unsigned char AsciiToEbcdic(unsigned char in_ebcdic_int);
unsigned char EbcdicToAscii(unsigned char in_ascii_int);
void float_to_ibm(int from[], int to[], int n, int endian);


void GetVectorMIN(vector<double> a, double &mindata);

void findMaxMin(float a[], int n, float &max, float &min);

int GetIndexOfMin(int *arry, int n);
int GetIndexOfMin(double *arry,int n);

void FindMaxMin(vector<double> pArr, int nLength, double &nMax, double &nMin);
void FindMaxMin(double *pArr, int nLength, double &nMax, double &nMin);
void FindMaxMin(float *pArr, int nLength, float &nMax, float &nMin);
void FindMaxMin(int *pArr, int nLength, int &nMax, int &nMin);



#endif // UTIL_H
