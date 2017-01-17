#ifndef Calculator_H
#define Calculator_H

#include <string>
#include <vector>
#include <list>
#include "csegy.h"
#include "util.h"
#include "survey.h"

#define PI 3.1415926


class Calculator
{
public:
    Calculator();

private:

    // Wellheads
    double treatWellX, treatWellY, monitorWellX, monitorWellY, perfWellX, perfWellY;
    // Well deviation data
    vector<WellDevPara*> devTreatWell, devMonitorWell, devPerfWell;
    // Geophone coordinates
    vector<double> rX, rY, rZ;
    // Velocity Model
    vector<double> vTreatWellDepth, vTreatWellVp, vTreatWellVs, vMonitorWellDepth, vMonitorWellVp, vMonitorWellVs;

    //Source coordinates
    double sX,sY,sZ;
    vector<double> planeA, planeB, planeC, planeD;

    vector<double> vAzimuth, vDip; // aizmuth and dip angles from source to receives
    vector<double> vRotation; // geophone rotation angles
    vector<double> vPTime, vSTime; //travel time from source to receivers, unit: ms

    //定义速度及扫描网格
    double m_VGridX1, m_VGridY1, m_VGridX2, m_VGridY2; //速度模型边界
    double m_GridX1, m_GridY1, m_GridX2, m_GridY2; //扫描网格边界
    double m_Z1, m_Z2;
    int m_StepX, m_StepY,m_StepZ;
    //叠加参数
    int m_nWinStart, m_nWinEnd; //叠加时窗大小
    int m_nQ1, m_nQ2; //相似叠加Q1值，取1~4; 线性叠加权重Q2取0时为只计算semblance相似系数，取1时计算相似加权叠加
    int m_nWaveType;



private:
    //SEGY数据信息
    const char* m_szSegyDir;
    const char* m_szEventDir;

    int m_nGeo;
    int m_nTraces;
    int m_nSamples;
    float m_nSampleRate;

    vector<vector<double>*> m_TracesZ,m_TracesX,m_TracesY;


private:
    CSegy m_segy;

public:

    //获取Events存储目录
    void GetEventDir(const char* s);

    // Velocity Model Max Boundary Points(BottomLeft, TopRight)
    bool setPseudoWell(double &x1, double &y1, double &x2, double &y2, double &y3, double &x3);



    //获取网格参数
    void GetGrid();

    //获取SEGY数据
    void GetSegyDir(const char* s);
    void GetSegyInfo(SegyPara sp);
    void SetParameters();
    void ReadData(const char *szFileName);
    void PseudoTrace();
    void LoadData(const char* szFileName);
    void ClearData(vector<vector<double> *> &v);


    //获取检波器坐标
    void ReadGeophones(const char* strPath);
    void GetGeoCoor(list<GeoPara *> &geophone);

    //获取速度值
    void LoadVModel(const char* &strPath);


    // Time to samples
    void Time2Sample(vector<double> time, vector<int> &sample, float rate);
    void CalTTT();

    //叠加
    void MoveoutCorrection(vector<vector<double>*> vOrigin, vector<vector<double>*> &vCorrected, vector<int> sample, int min, int max, int nGeo);
    void LinearStack(vector<vector<double>*> &vCorrected, double &lStack, int min, int max);
    void CalculateSemblance(vector<vector<double> *> vCorrected, double &semblance,int min, int max );
    void CalSET();

};

#endif // Calculator_H
