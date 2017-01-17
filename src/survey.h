#ifndef SURVEY_H
#define SURVEY_H

#include <string>
#include <vector>

using namespace std;

struct ForwardPara
{
    int ng;
    int ntr;
    int ns;
    int type; // perfshot event:0, microseismic event: 1
    int sRate; // Sample rate, unit: us
    double sX, sY, sZ; // Source point coordinate
    double snr; // signal to noise ratio
    double psr; // p-wave to s-wave energy ratio
    int pF, sF; // p-wave and s-wave frequencies
    int pLen, sLen; // p-wave and s-wave length
    double pA, sA; // p-wavelet and s-wavelet attenuation coefficients
    double attn; // p-wave and s-wave attenuation coefficient

};
struct WellHeadPara
{
    string name;
    double x;
    double y;
    int dev;
};

struct VModelPara
{
    vector<double> fd;   //压裂井的速度分层
    vector<double> md;   //监测井的速度分层
    vector<double> fvp;  //压裂井各层纵波速度
    vector<double> fvs;  //压裂井各层横波速度
    vector<double> mvp;  //监测井各层纵波速度
    vector<double> mvs;  //监测井各层横波速度
};

struct MsetGridPara
{
    //速度模型边界
    double VGridX1; //左下角
    double VGridY1;
    double VGridX2; //右上角
    double VGridY2;

    //扫描网格边界
    double GridX1; //左下角
    double GridY1;
    double GridX2; //右上角
    double GridY2;
    double GridZ1; //深度
    double GridZ2;

    //扫描网格步长
    int StepX;
    int StepY;
    int StepZ;

    //叠加波类型
    int WaveType;

    //叠加时窗
    int WinStart;
    int WinEnd;

    //叠加算法参数
    int Q1;
    int Q2;
};

struct GeoPara
{
    string name;
    double x;
    double y;
    double z;
};

struct PumpPara
{
    double time;
    double surfp;
    double tubep;
    double slurry;
    double proppant;
};

struct PrjPara{
    string pd; //项目路径
    string pn; //项目名称(不带后缀)
    string pp; //项目全称
    string event;
    string seis;
    string seisRaw;
    string perf;
    string treat;
    string vel;
    string geo;
    string well;

};

struct SegyPara{

    int nt; //道数
    int ns; //采样数
    int sr; //采样率
    int fc; //格式format code
};

struct DevPara
{
    double md;//斜深
    double incl; //倾角
    double azim;//方位角
    double tvd;//垂深
    double dx;//东西位移
    double dy;//南北位移
    double x; //东西绝对坐标
    double y;// 南北绝对坐标
    double z; // 考虑kb的深度，一般为负值
};

struct MsEventPara
{
    string date;
    string time;
    int stage;
    double x;
    double y;
    double tvd;
    double mag;
    double confi;
    double dip;
    double azim;
};

struct WellDevPara
{
    double MD;//斜深
    double Dip; //倾角
    double Azimuth;//方位角
    double TVD;//垂深
    double EW;//东西位移
    double SN;//南北位移
};


#endif // SURVEY_H

