#ifndef RAYTRACE_H
#define RAYTRACE_H

#include <string>
#include <vector>
#include "survey.h"


using namespace std;

typedef unsigned char byte;



void MDtoTVD(double MD, double &TVD, int nStyle, int TreatWellX, int TreatWellY, int MonitorWellX, int MonitorWellY, int PerfWellX, int PerfWellY, vector<WellDevPara*> devTreatWell, vector<WellDevPara*> devMonitorWell, vector<WellDevPara*> devPerfWell);
bool TVDtoZXY(double TVD, double &MD, double &X, double &Y, int nStyle, int TreatWellX, int TreatWellY, int MonitorWellX, int MonitorWellY, int PerfWellX, int PerfWellY, vector<WellDevPara*> devTreatWell, vector<WellDevPara*> devMonitorWell, vector<WellDevPara*> devPerfWell);
bool MDtoZXY(double MD, double &TVD, double &X, double &Y, int nStyle, int TreatWellX, int TreatWellY, int MonitorWellX, int MonitorWellY, int PerfWellX, int PerfWellY, vector<WellDevPara*> devTreatWell, vector<WellDevPara*> devMonitorWell, vector<WellDevPara*> devPerfWell);


void Plane(vector<double> &plane_a, vector<double> &plane_b, vector<double> &plane_c, vector<double> &plane_d, double twX, double twY, double mwX, double mwY, double pwX, double pwY, vector<WellDevPara*> twDev, vector<WellDevPara*> mwDev, vector<WellDevPara*> pwDev, vector<double> depthTreatWell, vector<double> depthMonitorWell, vector<double> vpTreatWell, vector<double>vsTreatWell, vector<double> vpMonitorWell, vector<double> vsMonitorWell);

void Swhere(double x, double y, double z, int &i11, int &i22, vector<double> plane_a, vector<double> plane_b, vector<double> plane_c, vector<double> plane_d);

double GetVelocity(int Layer, double X, int Type, double twX, double twY, double mwX, double mwY, double pwX, double pwY, vector<WellDevPara*> twDev, vector<WellDevPara*> mwDev, vector<WellDevPara*> pwDev, vector<double> depthTreatWell, vector<double> depthMonitorWell, vector<double> vpTreatWell, vector<double>vsTreatWell, vector<double> vpMonitorWell, vector<double> vsMonitorWell);

void CalculateTime(vector<double> &pTime, vector<double> &sTime, double twX, double twY, double mwX, double mwY, double pwX, double pwY, vector<WellDevPara*> twDev, vector<WellDevPara*> mwDev, vector<WellDevPara*> pwDev, double sourceX, double sourceY, double sourceZ, vector<double> depthTreatWell, vector<double> depthMonitorWell, vector<double> vpTreatWell, vector<double>vsTreatWell, vector<double> vpMonitorWell, vector<double> vsMonitorWell, vector<double> plane_a, vector<double> plane_b, vector<double> plane_c, vector<double> plane_d, vector<double> rx, vector<double> ry, vector<double> rz);

void calAngle(double sx, double sy, double sz, vector<double> rx, vector<double> ry, vector<double> rz, vector<double> &azimuth, vector<double> &dip);
void genWavelet1(vector<double> &wave, int Fre, double LeftRatio, double RightRatio, int Len, int sRate);
void genWavelet2(vector<double> &wave, int A1, int A2, int Fre, int k1, int k2, int Len, int sRate);

#endif // RAYTRACE_H
