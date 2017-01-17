#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<math.h>


#include"util.h"
#include"raytrace.h"


// Velocity Model Max Boundary Points(BottomLeft, TopRight)
bool setPseudoWell(double &x1, double &y1, double &x2, double &y2, double &x3, double &y3)
{
    x1 = 0.0;
    y1 = 0.0;
    x2 = 6000.0;
    y2 = 6000.0;
    x3 = 6000.0;
    y3 = 6000.0;
    
    return true;
}

// Read ASCII file containing geophone information, i.e., ID, X, Y, Z
bool readGeophones(const char* fileName, int num, vector<double> &rx, vector<double> &ry, vector<double> &rz)
{
    rx.clear();
    ry.clear();
    rz.clear();
    
    vector<GeoPara*> vGeo;
    GeoPara *geo;
    
    
    std::ifstream fin(fileName);
    std::stringstream strStream;
    std::string strTemp;
    int i = 0;
    while(fin.good() && !fin.eof())
    {
        
        geo = new GeoPara;
        
        getline(fin, strTemp,'\n');
        strStream.clear();
        strStream.str(strTemp);
        strStream>>geo->name>>geo->x>>geo->y>>geo->z;
        
        //printf("[Geophone Info] ID:%s, X:%f, Y:%f, Z:%f\n", geo->name.c_str(), geo->x,geo->y,geo->z);
        rx.push_back(geo->x);
        ry.push_back(geo->y);
        rz.push_back(geo->z);
        //printf("[Geophone Coordinates] ID:%d, X:%f, Y:%f, Z:%f\n", i, rx.at(i),ry.at(i),rz.at(i));
        
        vGeo.push_back(geo);
        //printf("[Geophone Coordinates] ID:%d, X:%f, Y:%f, Z:%f\n", i, vGeo.at(i)->x, vGeo.at(i)->y, vGeo.at(i)->z);
        if(strStream.fail())
            continue;
        
        i++;
    }
    delete geo;
    fin.clear();
    fin.close();
    
    return true;
    
}

// Read ASCII file containing velocity information
bool readVelModel(const char* fileName, vector<double> &depthTreatWell, vector<double> &vpTreatWell, vector<double> &vsTreatWell, vector<double> &depthMonitorWell, vector<double> &vpMonitorWell, vector<double> &vsMonitorWell)
{
    
    depthTreatWell.clear();
    depthMonitorWell.clear();
    vpTreatWell.clear();
    vsTreatWell.clear();
    vpMonitorWell.clear();
    vsMonitorWell.clear();
    
    
    std::ifstream fin(fileName);
    std::stringstream strStream;
    std::string strTemp;
    
    int i = 0;
    int j = 0;
    
    double a,b,c,d,e,f,g,h;
    while(fin.good() && !fin.eof())
    {
        
        getline(fin, strTemp,'\n');
        strStream.clear();
        strStream.str(strTemp);
        strStream>>a>>b>>c>>d>>e>>f>>g>>h;
        
        depthTreatWell.push_back(a);
        depthTreatWell.push_back(b);
        //printf("[Treatment Well layer detpth] Layer:%d, Top:%f, Bottom:%f\n", j, depthTreatWell[i], depthTreatWell[i+1]);
        vpTreatWell.push_back(c);
        vsTreatWell.push_back(d);
        
        depthMonitorWell.push_back(e);
        depthMonitorWell.push_back(f);
        //printf("[Monitor Well layer detpth] Layer:%d, Top:%f, Bottom:%f\n", j, depthMonitorWell[i], depthMonitorWell[i+1]);
        vpMonitorWell.push_back(g);
        vsMonitorWell.push_back(h);
        
        if(strStream.fail())
            continue;
        
        i+=2;
        j++;
    }
    fin.clear();
    fin.close();
    
    return true;
    
}

// Read ASCII file containing monitoring well trajectory information
bool readWellDev(const char* fileName, vector<WellDevPara*> &devData)
{
    devData.clear();
    
    WellDevPara *dev;
    std::ifstream fin(fileName);
    std::stringstream strStream;
    std::string strTemp;
    
    int i = 0;
    double a,b,c,d,e,f;
    while(fin.good() && !fin.eof())
    {
        dev = new WellDevPara;
        
        getline(fin, strTemp,'\n');
        strStream.clear();
        strStream.str(strTemp);
        strStream>>a>>b>>c>>d>>e>>f;
        dev->MD = a;
        dev->Dip = b;
        dev->Azimuth = c;
        dev->TVD = d;
        dev->EW = e;
        dev->SN = f;
        
        if(strStream.fail())
            continue;
        devData.push_back(dev);
        printf("[Well Deviation Info %d]: MD:%f, Dip:%f, Azimuth:%f, TVD:%f, EW:%f, SN:%f\n", i, devData[i]->MD, devData[i]->Dip, devData[i]->Azimuth, devData[i]->TVD, devData[i]->EW, devData[i]->SN);
        i++;
    }
    delete dev;
    
    
    fin.clear();
    fin.close();
    
    return true;
    
}



bool setForwardParameters(ForwardPara &fp)
{
    fp.ng = 49;
    fp.ntr = 49;
    fp.ns = 4000;
    fp.type = 1;
    fp.sRate = 500; // unit: us
    fp.sX = 3000.0;
    fp.sY = 3000.0;
    fp.sZ = 2000.0;
    fp.snr = 20.0;
    fp.psr = 1.0;
    fp.pF = 60;
    fp.sF = 20;
    fp.pLen = 800;
    fp.sLen = 800;
    fp.pA = 0.05;
    fp.sA = 0.05;
    fp.attn = 0.5;
    
    
    return true;
}

// Set attenuation coefficient according to travel time table
void setAttn(vector<double> &attn, vector<double> time, double minAttn)
{
    attn.clear();
    double max = 0.0;
    double min = 0.0;
    int count = time.size();
    
    FindMaxMin(time, count, max, min);
    
    double width = max - min;
    
    for (int i = 0; i <  count; i++)
    {
        double coef = minAttn/((time[i] - min)/width * (1.0 - minAttn)+ minAttn);
        attn.push_back(coef);
        //printf ("Attenuation Coefficient[%d]: %f \n", i, attn[i]);
    }
    
}



int main(int argc, char* argv[])
{
    // Velocity model grid boundary
    //double bLeftX, bLeftY, tRightX, tRightY;
    
    // Wellheads
    double treatWellX, treatWellY, monitorWellX, monitorWellY, perfWellX, perfWellY;
    // Well deviation data
    vector<WellDevPara*> devTreatWell, devMonitorWell, devPerfWell;
    // Geophone coordinates
    vector<double> rX, rY, rZ;
    // Velocity Model
    vector<double> vTreatWellDepth, vTreatWellVp, vTreatWellVs, vMonitorWellDepth, vMonitorWellVp, vMonitorWellVs;
    // Forward parameters
    ForwardPara fwdPara;
    // Synthetic seismic data
    vector<vector<double>*> synData;
    
    
    
    //
    vector<double> planeA, planeB, planeC, planeD;
    
    vector<double> vAzimuth, vDip; // aizmuth and dip angles from source to receives
    vector<double> vRotation; // geophone rotation angles
    vector<double> vPTime, vSTime; //travel time from source to receivers, unit: ms
    
    
    bool bRead = false;
    bool bSet = false;
    
    bSet = setForwardParameters(fwdPara);
    if (bSet)
    {
        printf("\n");
        printf("[Foward Parameters]: \n");
        printf("[Source Coordinates]: X:%f, Y:%f, Z:%f\n", fwdPara.sX, fwdPara.sY, fwdPara.sZ);
        printf("\n");
        bSet = false;
    }
    
    int geoNum = fwdPara.ng;
    int traceNum = fwdPara.ntr;
    int sampleNum = fwdPara.ns;
    int eventType = fwdPara.type;
    int sampleRate = fwdPara.sRate;
    double sX = fwdPara.sX;
    double sY = fwdPara.sY;
    double sZ = fwdPara.sZ;
    double SNR = fwdPara.snr;
    double PSR = fwdPara.psr;
    int pFreq = fwdPara.pF;
    int sFreq = fwdPara.sF;
    int pLen = fwdPara.pLen;
    int sLen = fwdPara.sLen;
    double pA = fwdPara.pA;
    double sA = fwdPara.sA;
    double attn = fwdPara.attn;
    
    //
    bSet = setPseudoWell(treatWellX, treatWellY, monitorWellX, monitorWellY, perfWellX, perfWellY);
    if (bSet)
    {
        printf("Info: Wellheads successfully set!\n");
        bSet = false;
    }
    
    char buff1[] = "../input_data/wellpath_monitor_vertical.txt";
    const char* mwDevFile = buff1;
    bRead = readWellDev(mwDevFile, devMonitorWell);
    if (bRead)
    {
        printf("Info: Monitor well deviation data successfully read and loaded!\n");
        bRead = false;
    }
    
    char buff2[] = "../input_data/wellpath_SS103H_vertical.txt";
    const char* twDevFile = buff2;
    bRead = readWellDev(twDevFile, devTreatWell);
    if (bRead)
    {
        printf("Info: Treatment well deviation data successfully read and loaded!\n");
        bRead = false;
    }
    
    char buff3[] = "../input_data/wellpath_SS103H_vertical.txt";
    const char* pwDevFile = buff3;
    bRead = readWellDev(pwDevFile, devPerfWell);
    if (bRead)
    {
        printf("Info: Perforation well deviation data successfully read and loaded!\n");
        bRead = false;
    }
    
    //char buf1[] = "../output_data/stararray.dat";//"../input_data/demo4_geo_coor.prn";
    char buf1[] = "../output_data/gridarray_fwd4.dat";
    const char* geoFile = buf1;
    bRead = readGeophones(geoFile, geoNum, rX, rY, rZ);
    if (bRead)
    {
        printf("Info: Geophones data successfully read and loaded!\n");
        bRead = false;
    }
    
    
    char buf2[] = "../input_data/velocity_SS103H_true_new01.txt";
    const char* vmFile = buf2;
    bRead = readVelModel(vmFile, vTreatWellDepth, vTreatWellVp, vTreatWellVs, vMonitorWellDepth, vMonitorWellVp, vMonitorWellVs);
    if (bRead)
    {
        printf("Info: Velocity model data successfully read and loaded!\n");
        bRead = false;
    }
    /* debug codes
     int size1 = vTreatWellVp.size();
     int size2 = vTreatWellDepth.size();
     printf ("Vp size: %d; Depth size: %d\n", size1, size2);
     int size = devTreatWell.size();
     for (int i = 0; i < size; i++)
     {
     printf("devTreatWell[%d] TVD: %f\n", i, devTreatWell.at(i)->TVD);
     }
     */
    
    Plane(planeA, planeB, planeC, planeD, treatWellX, treatWellY, monitorWellX, monitorWellY, perfWellX, perfWellY, devTreatWell, devMonitorWell, devPerfWell, vTreatWellDepth, vMonitorWellDepth, vTreatWellVp, vTreatWellVs, vMonitorWellVp, vMonitorWellVs);
    
    
    
    CalculateTime(vPTime, vSTime, treatWellX, treatWellY, monitorWellX, monitorWellY, perfWellX, perfWellY, devTreatWell, devMonitorWell, devPerfWell, sX, sY, sZ, vTreatWellDepth, vMonitorWellDepth, vTreatWellVp, vTreatWellVs, vMonitorWellVp, vMonitorWellVs, planeA, planeB, planeC, planeD, rX, rY, rZ);
    double minPTime = 0.0;
    GetVectorMIN(vPTime, minPTime);
    vector<double> vAttn;
    setAttn(vAttn, vPTime, attn);
    
    vector<double> waveData;

    // Generate wavelet for P and/or S waves
    //genWavelet1(waveData, pFreq, pA, pA, pLen, sampleRate);
    genWavelet2(waveData, 3, 1, pFreq, 80, 50, pLen, sampleRate);
    /*debug codes
     int size = waveData.size();
     for (int i = 0; i < size; i++)
     {
     printf("[Wavelet Sample Value:%d]: %f\n", i, waveData[i]);
     }
     */
    
    
    
    // Calculate aizmuth and dip angles from source to receives
    calAngle(sX, sY, sZ, rX, rY, rZ, vAzimuth, vDip);
    /* debug codes
     int num = vAzimuth.size();
     printf("num:%d\n", num);
     for (int i = 0; i < num; i++)
     {
     printf("Azimuth[%d]: %f, Dip[%d]: %f\n", i,vAzimuth[i], i, vDip[i]);
     }
     */
    
    
    vector<double> *vTemp;
    for (int i = 0; i < traceNum; i++)
    {
        vTemp = new vector<double>;
        vTemp->clear();
        
        for (int j = 0; j < sampleNum; j++)
        {
            vTemp->push_back(0.0);
        }
        synData.push_back(vTemp);
    }
    
    double pRate, nRate;
    nRate = pow(10, SNR/10);
    
    pRate = nRate * 1;
    
    delete vTemp;
    
    for (int i = 0; i < traceNum; i++)
    {
        int n = 100 + (int)((vPTime[i]-minPTime)*(1000/sampleRate));
        
        for(int j = n; j< n + waveData.size(); j++)
        {
            (*synData[i])[j] = ((*synData[i])[j] - cos(vDip[i]) * pRate * waveData[j-n]) * vAttn[i];
            //printf("Synthetic Trace[%d] Sample[%d]: %f\n", i, j, (*synData[i])[j]);
        }
        
    }
    
    
    // Export synthetic seismic data into .dat file
    FILE *pFile = NULL;
    pFile = fopen("../output_data/gridarray_fwd4_syndata.dat", "wb");
    if (NULL == pFile)
    {
        printf("Error: Failed to open file for writing!\n");
        return 0;
    }
    
    float *buffer;
    for (int i = 0; i < traceNum; i++)
    {
        buffer = new float[sampleNum];
        for(int j = 0; j < sampleNum; j++)
        {
            buffer[j] = (float)(*synData[i])[j];
        }
        fwrite(buffer, sizeof(float), sampleNum, pFile);
        
    }
    
    fclose(pFile);
    pFile = NULL;
    
  
    
    return 0;
}
