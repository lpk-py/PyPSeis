#include "raytrace.h"
#pragma hdrstop
#include <memory>
#include <cstring>
#include <string.h>
#include <math.h>
#include <stdlib.h>


#define PI 3.1415926535897


void MDtoTVD(double MD, double &TVD, int nStyle, int TreatWellX, int TreatWellY, int MonitorWellX, int MonitorWellY, int PerfWellX, int PerfWellY, vector<WellDevPara*> devTreatWell, vector<WellDevPara*> devMonitorWell, vector<WellDevPara*> devPerfWell)
{
    
    int nCount;
    int i,j;
    TVD=MD;
    if (nStyle==0)
    {
        nCount=devTreatWell.size();
        for (i=0;i<nCount;i++)
        {
            if (MD<=devTreatWell[i]->MD)
            break;
        }
        if (i>=(nCount-1))
        {
            TVD=devTreatWell[nCount-1]->TVD;
        }
        else if(i==0)
        {
            return;
        }
        else if (i<(nCount-1)  && i!=0)
        {
            if (MD==devTreatWell[i]->MD)
            TVD=devTreatWell[i]->TVD;
            else
            TVD=devTreatWell[i-1]->TVD+(MD-devTreatWell[i-1]->MD)*(devTreatWell[i]->TVD-devTreatWell[i-1]->TVD)/(devTreatWell[i]->MD-devTreatWell[i-1]->MD);
        }
    }
    else if (nStyle==1)
    {
        nCount=devMonitorWell.size();
        for (i=0;i<nCount;i++)
        {
            if (MD<=devMonitorWell[i]->MD)
            break;
        }
        if (i>=(nCount-1))
        {
            TVD=devMonitorWell[nCount-1]->TVD;
        }
        else if(i==0)
        {
            return;
        }
        else if (i<(nCount-1) && i!=0)
        {
            if (MD==devMonitorWell[i]->MD)
            TVD=devMonitorWell[i]->TVD;
            else
            TVD=devMonitorWell[i-1]->TVD+(MD-devMonitorWell[i-1]->MD)*(devMonitorWell[i]->TVD-devMonitorWell[i-1]->TVD)/(devMonitorWell[i]->MD-devMonitorWell[i-1]->MD);
        }
    }
    
}

bool TVDtoZXY(double TVD, double &MD, double &X, double &Y, int nStyle, int TreatWellX, int TreatWellY, int MonitorWellX, int MonitorWellY, int PerfWellX, int PerfWellY, vector<WellDevPara*> devTreatWell, vector<WellDevPara*> devMonitorWell, vector<WellDevPara*> devPerfWell)
{
    
    int nCount;
    int i,j;
    bool bFind,bDraw;
    bFind=false;
    bDraw=true;
    //vertical well
    if(nStyle==0)
    {
        X=TreatWellX;
        Y=TreatWellY;
    }
    else if(nStyle==1)
    {
        X=MonitorWellX;
        Y=MonitorWellY;
    }
    
    
    if (nStyle==0)
    {
        nCount=devTreatWell.size();
        
        for (i=0;i<nCount;i++)
        {
            if (TVD<=devTreatWell[i]->TVD)
            break;
        }
        
        if (i>=(nCount-1))
        {
            X=devTreatWell[nCount-1]->EW;
            Y=devTreatWell[nCount-1]->SN;
            MD=devTreatWell[nCount-1]->MD;
            bDraw=false;
        }
        
        else if (i<(nCount-1))
        {
            if (TVD==devTreatWell[i]->TVD )
            {
                
                X=devTreatWell[i]->EW;
                Y=devTreatWell[i]->SN;
                MD=devTreatWell[i]->MD;
                
            }
            
            else
            {
                X=devTreatWell[i-1]->EW+(TVD-devTreatWell[i-1]->TVD)*(devTreatWell[i]->EW-devTreatWell[i-1]->EW)/(devTreatWell[i]->TVD-devTreatWell[i-1]->TVD);
                Y=devTreatWell[i-1]->SN+(TVD-devTreatWell[i-1]->TVD)*(devTreatWell[i]->SN-devTreatWell[i-1]->SN)/(devTreatWell[i]->TVD-devTreatWell[i-1]->TVD);
                MD=devTreatWell[i-1]->MD+(TVD-devTreatWell[i-1]->TVD)*(devTreatWell[i]->MD-devTreatWell[i-1]->MD)/(devTreatWell[i]->TVD-devTreatWell[i-1]->TVD);
            }
            
        }
        X=X+TreatWellX;
        Y=Y+TreatWellY;
        return bDraw;
    }
    else if (nStyle==1)
    {
        nCount=devMonitorWell.size();
        for (i=0;i<nCount;i++)
        {
            if (TVD<=devMonitorWell[i]->TVD)
            break;
        }
        if (i>=(nCount-1))
        {
            X=devMonitorWell[nCount-1]->EW;
            Y=devMonitorWell[nCount-1]->SN;
            MD=devMonitorWell[nCount-1]->MD;
            bDraw=false;
        }
        else if (i<(nCount-1))
        {
            if (TVD==devMonitorWell[i]->TVD)
            {
                X=devMonitorWell[i]->EW;
                Y=devMonitorWell[i]->SN;
                MD=devMonitorWell[i]->MD;
            }
            else
            {
                X=devMonitorWell[i-1]->EW+(TVD-devMonitorWell[i-1]->TVD)*(devMonitorWell[i]->EW-devMonitorWell[i-1]->EW)/(devMonitorWell[i]->TVD-devMonitorWell[i-1]->TVD);
                Y=devMonitorWell[i-1]->SN+(TVD-devMonitorWell[i-1]->TVD)*(devMonitorWell[i]->SN-devMonitorWell[i-1]->SN)/(devMonitorWell[i]->TVD-devMonitorWell[i-1]->TVD);
                MD=devMonitorWell[i-1]->MD+(TVD-devMonitorWell[i-1]->TVD)*(devMonitorWell[i]->MD-devMonitorWell[i-1]->MD)/(devMonitorWell[i]->TVD-devMonitorWell[i-1]->TVD);
            }
        }
        X=X+MonitorWellX;
        Y=Y+MonitorWellY;
        return bDraw;
    }
    return false;
}

bool MDtoZXY(double MD, double &TVD, double &X, double &Y, int nStyle, int TreatWellX, int TreatWellY, int MonitorWellX, int MonitorWellY, int PerfWellX, int PerfWellY, vector<WellDevPara*> devTreatWell, vector<WellDevPara*> devMonitorWell, vector<WellDevPara*> devPerfWell)
{
    int nCount;
    int i;
    bool bDraw;
    bDraw=true;
    //vertical well
    if(nStyle==0)
    {
        X=TreatWellX;
        Y=TreatWellY;
    }
    else if(nStyle==1)
    {
        X=MonitorWellX;
        Y=MonitorWellY;
    }
    else if(nStyle==2)
    {
        X=PerfWellX;
        Y=PerfWellY;
    }
    
    
    //for treatment well
    if (nStyle==0)
    {
        nCount=devTreatWell.size();
        for (i=0;i<nCount;i++)
        {
            if (MD<=devTreatWell[i]->MD)
            break;
        }
        if (i>=(nCount-1))
        {
            TVD=devTreatWell[nCount-1]->TVD;
            X=devTreatWell[nCount-1]->EW;
            Y=devTreatWell[nCount-1]->SN;
            bDraw=false;
        }
        else if (i<(nCount-1))
        {
            if (MD==devTreatWell[i]->MD)
            {        TVD=devTreatWell[i]->TVD;
                X=devTreatWell[i]->EW;
                Y=devTreatWell[i]->SN;
            }
            else
            {        TVD=devTreatWell[i-1]->TVD+(MD-devTreatWell[i-1]->MD)*(devTreatWell[i]->TVD-devTreatWell[i-1]->TVD)/(devTreatWell[i]->MD-devTreatWell[i-1]->MD);
                X=devTreatWell[i-1]->EW+(MD-devTreatWell[i-1]->MD)*(devTreatWell[i]->EW-devTreatWell[i-1]->EW)/(devTreatWell[i]->MD-devTreatWell[i-1]->MD);
                Y=devTreatWell[i-1]->SN+(MD-devTreatWell[i-1]->MD)*(devTreatWell[i]->SN-devTreatWell[i-1]->SN)/(devTreatWell[i]->MD-devTreatWell[i-1]->MD);
            }
        }
        X=X+TreatWellX;
        Y=Y+TreatWellY;
        return bDraw;
    }
    //for monitoring well
    else if (nStyle==1)
    {
        nCount=devMonitorWell.size();
        for (i=0;i<nCount;i++)
        {
            if (MD<=devMonitorWell[i]->MD)
            break;
        }
        if (i>=(nCount-1))
        {
            TVD=devMonitorWell[nCount-1]->TVD;
            X=devMonitorWell[nCount-1]->EW;
            Y=devMonitorWell[nCount-1]->SN;
            bDraw=false;
        }
        else if (i<(nCount-1))
        {
            if (MD==devMonitorWell[i]->MD)
            {   TVD=devMonitorWell[i]->TVD;
                X=devMonitorWell[i]->EW;
                Y=devMonitorWell[i]->SN;
            }
            else
            {    TVD=devMonitorWell[i-1]->TVD+(MD-devMonitorWell[i-1]->MD)*(devMonitorWell[i]->TVD-devMonitorWell[i-1]->TVD)/(devMonitorWell[i]->MD-devMonitorWell[i-1]->MD);
                X=devMonitorWell[i-1]->EW+(MD-devMonitorWell[i-1]->MD)*(devMonitorWell[i]->EW-devMonitorWell[i-1]->EW)/(devMonitorWell[i]->MD-devMonitorWell[i-1]->MD);
                Y=devMonitorWell[i-1]->SN+(MD-devMonitorWell[i-1]->MD)*(devMonitorWell[i]->SN-devMonitorWell[i-1]->SN)/(devMonitorWell[i]->MD-devMonitorWell[i-1]->MD);
            }
        }
        X=X+MonitorWellX;
        Y=Y+MonitorWellY;
        return bDraw;
    }
    
    //for perforation well
    else if (nStyle==2)
    {
        nCount=devPerfWell.size();
        for (i=0;i<nCount;i++)
        {
            if (MD<=devPerfWell[i]->MD)
            break;
        }
        if (i>=(nCount-1))
        {
            TVD=devPerfWell[nCount-1]->TVD;
            X=devPerfWell[nCount-1]->EW;
            Y=devPerfWell[nCount-1]->SN;
            bDraw=false;
        }
        else if (i<(nCount-1))
        {
            if (MD==devPerfWell[i]->MD)
            {   TVD=devPerfWell[i]->TVD;
                X=devPerfWell[i]->EW;
                Y=devPerfWell[i]->SN;
            }
            else
            {    TVD=devPerfWell[i-1]->TVD+(MD-devPerfWell[i-1]->MD)*(devPerfWell[i]->TVD-devPerfWell[i-1]->TVD)/(devPerfWell[i]->MD-devPerfWell[i-1]->MD);
                X=devPerfWell[i-1]->EW+(MD-devPerfWell[i-1]->MD)*(devPerfWell[i]->EW-devPerfWell[i-1]->EW)/(devPerfWell[i]->MD-devPerfWell[i-1]->MD);
                Y=devPerfWell[i-1]->SN+(MD-devPerfWell[i-1]->MD)*(devPerfWell[i]->SN-devPerfWell[i-1]->SN)/(devPerfWell[i]->MD-devPerfWell[i-1]->MD);
            }
        }
        X=X+PerfWellX;
        Y=Y+PerfWellY;
        return bDraw;
    }
    
    return false;
}


void Plane(vector<double> &plane_a, vector<double> &plane_b, vector<double> &plane_c, vector<double> &plane_d, double twX, double twY, double mwX, double mwY, double pwX, double pwY, vector<WellDevPara*> twDev, vector<WellDevPara*> mwDev, vector<WellDevPara*> pwDev, vector<double> depthTreatWell, vector<double> depthMonitorWell, vector<double> vpTreatWell, vector<double>vsTreatWell, vector<double> vpMonitorWell, vector<double> vsMonitorWell)
{
    
    double d2,d1,c1;
    double x,y,TVD,Depth;
    int n,i,j;
    vector<double> x1,x2,y1,y2,z1,z2; //（x1,y1,z1）压裂井井点坐标  （x2,y2,z2）监测井井点坐标
    n=vpTreatWell.size()-1;//界面数=速度层数-1
    plane_a.clear();
    plane_b.clear();
    plane_c.clear();
    plane_d.clear();
    x1.clear();  y1.clear();   z1.clear();
    x2.clear();  y2.clear();   z2.clear();
    
    n=0;
    
    
    for(i=1;i<depthTreatWell.size();i+=2)
    {
        
        TVD=depthTreatWell[i];       //压裂井的 速度 ：是用垂深来输入的
        
        TVDtoZXY(TVD,Depth,x,y,0,twX, twY, mwX, mwY, pwX, pwY,
                 twDev, mwDev, pwDev); //压裂井：利用井斜获取空间坐标 //监测井的 速度 ：是用垂深来输入的
        
        x1.push_back(x);  y1.push_back(y);  z1.push_back(TVD);
        
        TVD=depthMonitorWell[i];        //监测井的 速度 ：是用垂深来输入的
        TVDtoZXY(TVD,Depth,x,y,1,twX, twY, mwX, mwY, pwX, pwY,
                 twDev, mwDev, pwDev); //监测井：利用井斜获取空间坐标 //监测井的 速度 ：是用垂深来输入的
        x2.push_back(x);  y2.push_back(y);  z2.push_back(TVD);
        
        n++;
        plane_a.push_back(0.0);  plane_b.push_back(0.0);  plane_c.push_back(0.0); plane_d.push_back(0.0);
    }
    
    for(int i2=0;i2<n;i2++)
    {
        d2=(x2[i2]-x1[i2])*(x2[i2]-x1[i2])+(y2[i2]-y1[i2])*(y2[i2]-y1[i2]);
        d1=sqrt(d2);
        d2=d2+(z2[i2]-z1[i2])*(z2[i2]-z1[i2]);
        d2=sqrt(d2);
        plane_c[i2]=d1/d2;
        c1=d1*d2;
        d2=(x2[i2]-x1[i2])*(x2[i2]-x1[i2])+(y2[i2]-y1[i2])*(y2[i2]-y1[i2]);
        plane_a[i2]=(x2[i2]-x1[i2])*(z2[i2]-z1[i2])/c1;
        plane_a[i2]=-plane_a[i2];
        plane_b[i2]=(y2[i2]-y1[i2])*(z2[i2]-z1[i2])/c1;
        plane_b[i2]=-plane_b[i2];
        plane_d[i2]=plane_a[i2]*x1[i2]+plane_b[i2]*y1[i2]+plane_c[i2]*z1[i2];
        plane_d[i2]=-plane_d[i2];
    }
}


void Swhere(double x, double y, double z, int &i11, int &i22, vector<double> plane_a, vector<double> plane_b, vector<double> plane_c, vector<double> plane_d)
{
    double f,f1,f2;
    int j6;
    int n;
    n=plane_a.size(); //界面层数
    f=plane_a[0]*x+plane_b[0]*y+plane_c[0]*z+plane_d[0];
    if(f<0.0)
    {
        i11=-1;
        i22=0;
        return;
    }
    if(f==0.0)
    {
        i11=0;
        i22=0;
        return;
    }
    f=plane_a[n-1]*x+plane_b[n-1]*y+plane_c[n-1]*z+plane_d[n-1];
    if(f>0.0)
    {
        i11=n-1;
        i22=n;
        return;
    }
    if(f==0.0)
    {
        i11=n-1;
        i22=n-1;
        return;
    }
    for(j6=0;j6<n-1;j6++)
    {
        f1=plane_a[j6]*x+plane_b[j6]*y+plane_c[j6]*z+plane_d[j6];
        f2=plane_a[j6+1]*x+plane_b[j6+1]*y+plane_c[j6+1]*z+plane_d[j6+1];
        if(f1>0.0&&f2<0.0)
        {
            i11=j6;
            i22=j6+1;
            return;
        }
        if(f1==0.0)
        {
            i11=j6;
            i22=j6;
            return;
        }
        if(f2==0.0)
        {
            i11=j6+1;
            i22=j6+1;
            return;
        }
    }
}

double GetVelocity(int Layer, double X, int Type, double twX, double twY, double mwX, double mwY, double pwX, double pwY, vector<WellDevPara*> twDev, vector<WellDevPara*> mwDev, vector<WellDevPara*> pwDev, vector<double> depthTreatWell, vector<double> depthMonitorWell, vector<double> vpTreatWell, vector<double>vsTreatWell, vector<double> vpMonitorWell, vector<double> vsMonitorWell)
{
    //计算某个点的速度,Layer：点在速度模型中的层位；X：点的x轴的坐标;Type：速度类型（0：Vp；1：Vs）
    int i;
    double TVD,Depth,x,y,midX;
    vector<float> treatWell,monitorWell;
    treatWell.clear();  monitorWell.clear();
    //按分层建立中点X坐标数组
    for(i=0;i<depthTreatWell.size();i+=2)
    {
        midX=0;
        TVD=depthTreatWell[i];
        TVDtoZXY(TVD,Depth,x,y,0,twX, twY, mwX, mwY, pwX, pwY,
                 twDev, mwDev, pwDev); //压裂井：利用井斜获取空间坐标  //压裂井的 速度 ：是用垂深来输入的
        midX=x;
        TVD=depthTreatWell[i+1];
        TVDtoZXY(TVD,Depth,x,y,0,twX, twY, mwX, mwY, pwX, pwY,
                 twDev, mwDev, pwDev); //压裂井：利用井斜获取空间坐标  //压裂井的 速度 ：是用垂深来输入的
        midX+=x;
        treatWell.push_back(midX/2);
        midX=0;
        TVD=depthMonitorWell[i];
        TVDtoZXY(TVD,Depth,x,y,1,twX, twY, mwX, mwY, pwX, pwY,
                 twDev, mwDev, pwDev); //监测井：利用井斜获取空间坐标  //监测井的 速度 ：是用垂深来输入的
        midX=x;
        TVD=depthMonitorWell[i+1];
        TVDtoZXY(TVD,Depth,x,y,1,twX, twY, mwX, mwY, pwX, pwY,
                 twDev, mwDev, pwDev); //监测井：利用井斜获取空间坐标  //监测井的 速度 ：是用垂深来输入的
        midX+=x;
        monitorWell.push_back(midX/2);
    }
    //进行对速度线性插值
    double Vel1,Vel2;
    if(Type==0)  //Vp
    {
        Vel1=vpTreatWell[Layer];  Vel2=vpMonitorWell[Layer];
    }
    else   //Vs
    {
        Vel1=vsTreatWell[Layer]; Vel2=vsMonitorWell[Layer];
    }
    double X1,X2;
    X1=treatWell[Layer];  X2=monitorWell[Layer];
    if(X1==X2)
        return (Vel1+Vel2)/2;
    else
        return Vel1+(Vel2-Vel1)*(X-X1)/(X2-X1);
}



void CalculateTime(vector<double> &pTime, vector<double> &sTime, double twX, double twY, double mwX, double mwY, double pwX, double pwY, vector<WellDevPara*> twDev, vector<WellDevPara*> mwDev, vector<WellDevPara*> pwDev, double sourceX, double sourceY, double sourceZ, vector<double> depthTreatWell, vector<double> depthMonitorWell, vector<double> vpTreatWell, vector<double>vsTreatWell, vector<double> vpMonitorWell, vector<double> vsMonitorWell, vector<double> plane_a, vector<double> plane_b, vector<double> plane_c, vector<double> plane_d, vector<double> rx, vector<double> ry, vector<double> rz)
{
    int i;
    int GLayer1,GLayer2; //检波器位置的上下分界面编号
    int SLayer1,SLayer2; //震源的位置的上下分界面编号
    int Layer1,Layer2,LayerCur;   //循环控制量
    double x,y,z,len;
    double vx,vy,vz;   //检波器和震源连线的方向数
    double vp,vs;
    double TimeP,TimeS;
    pTime.clear();
    sTime.clear();
    
    Swhere(sourceX,sourceY,sourceZ,SLayer1,SLayer2, plane_a, plane_b, plane_c, plane_d);
    double L0,L1,x0,y0,z0;
    
    for(i=0;i<rx.size();i++)
    {
        TimeP=TimeS=0;
        x=sourceX-rx[i]; y=sourceY-ry[i];  z=sourceZ-rz[i];
        len=sqrt(x*x+y*y+z*z);
        
        vx=x/len;  vy=y/len;  vz=z/len;
        Swhere(rx[i],ry[i],rz[i],GLayer1,GLayer2,plane_a, plane_b, plane_c, plane_d);
        if(SLayer1==GLayer1)  //震源和检波器在同一层中
        {
            
            TimeP=1000*len/((vpMonitorWell[SLayer2]+vpTreatWell[SLayer2])/2);
            TimeS=1000*len/((vsMonitorWell[SLayer2]+vsTreatWell[SLayer2])/2);
            pTime.push_back(TimeP);
            sTime.push_back(TimeS);
            
        }
        else if(GLayer1>SLayer1) //震源在检波器上方，地震波朝下传播
        {
            //
            Layer1=SLayer2;
            Layer2=GLayer1;  //下行波将穿越Layer1-Layer2
            x=sourceX; y=sourceY; z=sourceZ; //从震源开始计算
            for(LayerCur=Layer1;LayerCur<=Layer2;LayerCur+=1)
            {
                //由空间直线与空间平面之间的关系，求得下列关系式(求空间直线与空间平面之间的交点)
                L0=plane_a[LayerCur]*sourceX+plane_b[LayerCur]*sourceY+plane_c[LayerCur]*sourceZ+plane_d[LayerCur];
                L1=plane_a[LayerCur]*vx+plane_b[LayerCur]*vy+plane_c[LayerCur]*vz;
                L0=-L0/L1;
                x0=sourceX+vx*L0; y0=sourceY+vy*L0; z0=sourceZ+vz*L0;
                len=sqrt((x0-x)*(x0-x)+(y0-y)*(y0-y)+(z0-z)*(z0-z));
                vp=(GetVelocity(LayerCur,x0,0,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell)+GetVelocity(LayerCur,x,0, twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell))/2.0;
                vs=(GetVelocity(LayerCur,x0,1, twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell)+GetVelocity(LayerCur,x,1,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell))/2.0;
                TimeP+=1000*len/vp;
                TimeS+=1000*len/vs;
                x=x0;y=y0;z=z0;
            }
            //最后一段路程的旅行时
            len=sqrt((x0-rx[i])*(x0-rx[i])+(y0-ry[i])*(y0-ry[i])+(z0-rz[i])*(z0-rz[i]));
            vp=(GetVelocity(Layer2+1,rx[i],0,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell)+GetVelocity(Layer2+1,x0,0,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell))/2.0;
            vs=(GetVelocity(Layer2+1,rx[i],1,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell)+GetVelocity(Layer2+1,x0,1,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell))/2.0;
            TimeP+=1000*len/vp;
            TimeS+=1000*len/vs;
            pTime.push_back(TimeP);
            sTime.push_back(TimeS);
        }
        else if(GLayer1<SLayer1)               //震源在检波器下方，地震波朝上传播
        {
            
            Layer1=SLayer1;
            Layer2=GLayer2; //上行波将穿越Layer1-Layer2
            x=sourceX; y=sourceY; z=sourceZ; //从震源开始计算
            
            for(LayerCur=Layer1;LayerCur>=Layer2;LayerCur-=1) //
            {
                
                //由空间直线与空间平面之间的关系，求得下列关系式(求空间直线与空间平面之间的交点)
                L0=plane_a[LayerCur]*sourceX+plane_b[LayerCur]*sourceY+plane_c[LayerCur]*sourceZ+plane_d[LayerCur];
                L1=plane_a[LayerCur]*vx+plane_b[LayerCur]*vy+plane_c[LayerCur]*vz;
                L0=-L0/L1;
                x0=sourceX+vx*L0; y0=sourceY+vy*L0; z0=sourceZ+vz*L0;
                len=sqrt((x0-x)*(x0-x)+(y0-y)*(y0-y)+(z0-z)*(z0-z));
                
                vp=GetVelocity(LayerCur,x0,0,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell); //此处有问题,需要先确定两端速度，再求均值
                
                vp+=GetVelocity(LayerCur,x,0,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell);
                vp=vp/2.0;
                vs=GetVelocity(LayerCur,x0,1,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell); //此处有问题,需要先确定两端速度，再求均值
                vs+=GetVelocity(LayerCur,x,1,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev,depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell);
                vs=vs/2.0;
                TimeP+=1000*len/vp;
                TimeS+=1000*len/vs;
                x=x0;y=y0;z=z0;
            }
            
            //最后一段路程的旅行时
            len=sqrt((x0-rx[i])*(x0-rx[i])+(y0-ry[i])*(y0-ry[i])+(z0-rz[i])*(z0-rz[i]));
            vp=(GetVelocity(Layer2-1,rx[i],0,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev, depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell)+GetVelocity(Layer2-1,x0,0,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev, depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell))/2.0;
            vs=(GetVelocity(Layer2-1,rx[i],1,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev, depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell)+GetVelocity(Layer2-1,x0,1,twX, twY, mwX, mwY, pwX, pwY, twDev, mwDev, pwDev, depthTreatWell, depthMonitorWell, vpTreatWell,vsTreatWell, vpMonitorWell, vsMonitorWell))/2.0;
            
            TimeP+=1000*len/vp;
            TimeS+=1000*len/vs;
            pTime.push_back(TimeP);
            sTime.push_back(TimeS);
        }
    }
}

void calAngle(double sx, double sy, double sz, vector<double> rx, vector<double> ry, vector<double> rz, vector<double> &azimuth, vector<double> &dip)
{
    
    int i;
    double xx,yy,zz,tempAngle;
    
    azimuth.clear();
    dip.clear();
    
    for(i=0;i<rx.size();i++)
    {
        xx=sx-rx[i];
        yy=sy-ry[i];
        zz=sz-rz[i];
        zz=zz>0?zz:-1*zz;
        if(xx==0)
        {
            tempAngle=(yy>0)?90.0:-90.0;
            azimuth.push_back(tempAngle);
        }
        else if(yy==0)
        {
            tempAngle=(xx>0)?0:180.0;
            azimuth.push_back(tempAngle);
        }
        else
        {
            //应该考虑-180°~180°
            if(yy*xx>0 && yy>0)
                azimuth.push_back(180*atan(yy/xx)/PI);
            else if(yy*xx>0 && yy<0)
                azimuth.push_back(180*atan(yy/xx)/PI-180.0);
            else if(yy*xx<0 && yy>0)
                azimuth.push_back(180+180*atan(yy/xx)/PI);
            else if(yy*xx<0 && yy<0)
                azimuth.push_back(180*atan(yy/xx)/PI);
            
        }
        xx=sqrt(xx*xx+yy*yy);
        if(xx==0)
            tempAngle=90;
        else
            tempAngle=180*atan(zz/xx)/PI;
        if(sz>rz[i])
            dip.push_back(-1*tempAngle);
        else
            dip.push_back(tempAngle);
        
    }
}



void genWavelet1(vector<double> &wave, int Fre, double LeftRatio, double RightRatio, int Len, int sRate)
{
    int LeftLen,RightLen,i,j;
    vector<double> left,right;
    double Temp;
    double N2T;  //点数变时间
    LeftLen=Len/4;
    RightLen=Len-LeftLen;
    wave.clear();
    for(i=0;i<Len;i++)
        wave.push_back(0.0);
    N2T=1000000.0/sRate;  //1s对应采样点数
    left.clear(); right.clear();
    //(double)i/N2T功能是：将采样点数变成时间(s); exp(-1*LeftRatio*i): 衰减量
    for(i=1;i<=LeftLen;i++)
    {
        Temp=-1*sin(2*PI*(double)Fre*((double)i/N2T))*exp(-1*LeftRatio*i);
        left.push_back(Temp);
    }
    for(i=0;i<RightLen;i++)
    {
        Temp=sin(2*PI*(double)Fre*((double)i/N2T))*exp(-1*RightRatio*i);
        right.push_back(Temp);
    }
    for(i=0;i<Len;i++)
    {
        if(i<LeftLen)
            wave[i] = left[LeftLen-i-1];
        else
            wave[i] = right[i-LeftLen];
    }
}



// Both the P and S wavelets have the form
//  s(t) = A0 * sin(2πft) * exp(–kt)
// This is the response of a damped harmonic oscillator. Specific
// values for the parameters f and k are arbitrary; we have chosen
// f = 300 Hz, k = 80/s for the P arrivals, and f = 200 Hz, k = 50/s
// for the S arrivals.
void genWavelet2(vector<double> &wave, int A1, int A2, int Fre, int k1, int k2, int Len, int sRate)
{
    int LeftLen,RightLen,i,j;
    vector<double> left,right;
    double Temp;
    double N2T;  //点数变时间
    LeftLen=Len/4;
    RightLen=Len-LeftLen;
    wave.clear();
    for(i=0;i<Len;i++)
        wave.push_back(0.0);
    N2T=1000000.0/sRate;  //1s对应采样点数
    left.clear(); right.clear();
    //(double)i/N2T功能是：将采样点数变成时间(s); exp(-1*k*((double)i/N2T)): 衰减量
    for(i=1;i<=LeftLen;i++)
    {
        Temp=-1*A1*sin(2*PI*(double)Fre*((double)i/N2T))*exp(-1*k1*((double)i/N2T));
        left.push_back(Temp);
    }
    for(i=0;i<RightLen;i++)
    {
        Temp=A2*sin(2*PI*(double)Fre*((double)i/N2T))*exp(-1*k2*((double)i/N2T));
        right.push_back(Temp);
    }
    for(i=0;i<Len;i++)
    {
        if(i<LeftLen)
            wave[i] = left[LeftLen-i-1];
        else
            wave[i] = right[i-LeftLen];
    }
}

