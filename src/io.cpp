#include "io.h"


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
