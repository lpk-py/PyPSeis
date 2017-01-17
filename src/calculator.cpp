#include "calculator.h"
#include <math.h>
#include <stdio.h>

#include "raytrace.h"



Calculator::Calculator()
{

}

bool Calculator::setPseudoWell(double &x1, double &y1, double &x2, double &y2, double &y3, double &x3)
{
    x1 = 0.0;
    y1 = 0.0;
    x2 = 6000.0;
    y2 = 6000.0;
    x3 = 6000.0;
    y3 = 6000.0;

    return true;

}



void Calculator::GetGrid()
{
    m_VGridX1 = 0.0;
    m_VGridY1 = 0.0;
    m_VGridX2 = 6000.0;
    m_VGridY2 = 6000.0;

    m_GridX1 = 2000.0;
    m_GridY1 = 2000.0;
    m_GridX2 = 4048.0;
    m_GridY2 = 4048.0;

    m_Z1 = 2900.0;
    m_Z2 = 3100.0;

    m_StepX = 2;
    m_StepY = 2;
    m_StepZ = 10;

    m_nWaveType = 1;
    m_nWinStart = 2500;
    m_nWinEnd = 6596;

    m_nQ1 = 1;
    m_nQ2 = 0;

    printf("Velocity Model Boundary: \nBottom Left[X:%6.2f   Y:%6.2f]\nTop Right[X:%6.2f   Y:%6.2f]\n", m_VGridX1, m_VGridY1,m_VGridX2, m_VGridY2);
    printf("Scanning Grid Boundary: \nBottom Left[X:%6.2f   Y:%6.2f]\nTop Right[X:%6.2f   Y:%6.2f]\nDepth Interval[Z1:%6.2f   Z2:%6.2f]\n", m_GridX1, m_GridY1,m_GridX2, m_GridY2, m_Z1, m_Z2);
    printf("Scanning Grid Offset: \nOffset X:%d   Y:%d   Z:%d\n", m_StepX, m_StepY, m_StepZ);
    printf("Wave Type:\n%d (Note: P-wave = 0, S-wave = 1)\n", m_nWaveType);
    //printf("Time Window:\nStart:%d   End:%d\n", m_nWinStart, m_nWinEnd);


}

void Calculator::GetSegyDir(const char* s)
{
    m_szSegyDir = s;
}

void Calculator::GetSegyInfo(SegyPara sp)
{
    m_nTraces = sp.nt;
    m_nSamples = sp.ns;
    m_nSampleRate = sp.sr / 1000;

}

void Calculator::SetParameters()
{
    m_nGeo = 49;
    printf("Total Geophone Number:\n%d\n", m_nGeo );

}

void Calculator::ReadData(const char *szFileName)
{
    m_segy.OpenFile(szFileName);

    m_nTraces = m_segy.getTotalTraceNumber();
    m_nSamples =m_segy.getSamplesNumber();
    m_nSampleRate = m_segy.getSamplesInterval() / 1000;
    if (m_nSampleRate == 0)
        m_nSampleRate = 0.25;
    m_nGeo = m_nTraces/3;

    printf("Total Geophone Number:\n%d\n", m_nGeo );
    printf("Sample Rate:\n%3.2f\n", m_nSampleRate );
    printf("Samples per Trace:\n%d\n", m_nSamples );


    for (int i=0; i<m_nTraces; i++)
    {

        vector<double>* m_ChangeData=new vector<double>;
        m_ChangeData->clear();

        float *pTrace = m_segy.GetTraceData(i+1);

        for (int k=0; k<m_nSamples; k++)
        {
            m_ChangeData->push_back((double)pTrace[k]);

        }

        if (i%3 == 0)  //地面数据，先Z后XY;井下及正演数据，先xy后z

            m_TracesZ.push_back(m_ChangeData);
        else if (i%3 == 1)
            m_TracesX.push_back(m_ChangeData);
        else if (i%3 == 2)
            m_TracesY.push_back(m_ChangeData);


        delete []pTrace;
        pTrace = NULL;
    }

    m_segy.closeFile();
}

void Calculator::PseudoTrace()
{
    double sampleValue= 0.0;

    for (int i=0; i<m_nTraces; i++)
    {

        vector<double>* m_ChangeData = new vector<double>;
        m_ChangeData->clear();

        for (int k=0; k<m_nSamples; k++)
        {
            m_ChangeData->push_back(sampleValue);
        }

        if (i%3 == 0)  //地面数据，先Z后XY;井下及正演数据，先xy后z

            m_TracesZ.push_back(m_ChangeData);
        else if (i%3 == 1)
            m_TracesX.push_back(m_ChangeData);
        else if (i%3 == 2)
            m_TracesY.push_back(m_ChangeData);
    }
}

void Calculator::LoadData(const char* szFileName)
{

}

void Calculator::ClearData(vector<vector<double> *> &v)
{
    for (int i=0;i<v.size();i++)
    {

        delete v[i];
    }
    v.clear();

}

void Calculator::Time2Sample(vector<double> time, vector<int> &sample, float rate)
{
    sample.clear();
    for (int i=0; i< time.size(); i++)
    {
        sample.push_back(time[i]/rate);
    }

}



void Calculator::CalTTT()
{
    int imax,jmax,kmax, i, j, k;
    imax = (m_GridX2-m_GridX1)/m_StepX;
    jmax = (m_GridY2-m_GridY1)/m_StepY;
    kmax = (m_Z2-m_Z1)/m_StepZ;

    //int size = imax*jmax;

    FILE * pFile = NULL;

    pFile = fopen ("ttt.dat" , "wb" );
    if (NULL == pFile)
    {
        printf("File open fail!\n");
        return ;
    }

    printf("Travel time table calculation in progress[Waiting...]\n");

    Plane(planeA, planeB, planeC, planeD, treatWellX, treatWellY, monitorWellX, monitorWellY, perfWellX, perfWellY, devTreatWell, devMonitorWell, devPerfWell, vTreatWellDepth, vMonitorWellDepth, vTreatWellVp, vTreatWellVs, vMonitorWellVp, vMonitorWellVs);


    for (k = 0; k<kmax;k++)
    {
        for (j = 0; j<jmax;j++)
        {
            for (i = 0; i<imax;i++)
            {
                sX = m_GridX1 + i*m_StepX;
                sY = m_GridY1 + j*m_StepY;
                sZ = m_Z1 + k*m_StepZ;

                CalculateTime(vPTime, vSTime, treatWellX, treatWellY, monitorWellX, monitorWellY, perfWellX, perfWellY, devTreatWell, devMonitorWell, devPerfWell, sX, sY, sZ, vTreatWellDepth, vMonitorWellDepth, vTreatWellVp, vTreatWellVs, vMonitorWellVp, vMonitorWellVs, planeA, planeB, planeC, planeD, rX, rY, rZ);

                float *buffer = new float[m_nGeo];

                switch (m_nWaveType)
                {

                case 0:

                    for (int table=0;table<m_nGeo;table++)
                        buffer[table] = (float)vPTime[table] ;
                    break;

                case 1:

                    for (int table=0;table<m_nGeo;table++)
                        buffer[table] =  (float)vSTime[table] ;
                    break;
                default:
                    break;
                }


                printf("TravelTime[k:%d, j:%d, i:%d]...\n", k, j, i);
                fwrite (buffer , sizeof(float), m_nGeo , pFile );
                delete []buffer;

            }
        }

    }
    printf("Travel time table calculation completed[OK]\n");
    fclose (pFile);
    pFile = NULL;

}

void Calculator::MoveoutCorrection(vector<vector<double> *> vOrigin, vector<vector<double> *> &vCorrected, vector<int> sample, int min, int max, int nGeo)
{
    int i, j,k,smin;
    double tempArray[nGeo];

    for ( i=0; i <nGeo; i++)
    {
        vector<double> * vTemp = new vector<double>;

        for ( j= 0; j<max-min; j++)
        {
            vTemp->push_back(0.0);
        }
        vCorrected.push_back(vTemp);
    }

    for ( i=0; i <nGeo; i++)
    {

        tempArray[i] = sample[i];

    }
    smin = GetIndexOfMin(tempArray,nGeo);


    for ( i=0; i <nGeo; i++)
    {

        k= sample[i]-sample[smin];
        //        for ( j= 0; j<m_nSamples-k; j++)
        //        {
        //            vCorrected[i]->replace(j,(*vOrigin[i])[j+k]);
        //        }

        for ( j=0; j<max-min; j++)
        {
            if (k<m_nSamples-(j+min))
                (*vCorrected[i])[j] = (*vOrigin[i])[j+min+k];
            else
                (*vCorrected[i])[j] = 0.0;

        }


    }
}

void Calculator::LinearStack(vector<vector<double> *> &vCorrected, double &lStack, int min, int max)
{
    for (int i=0; i<max-min; i++)
    {

        double temp=0.0;

        for (int j=0; j<m_nGeo; j++)
        {

            temp+=(*vCorrected[j])[i];

        }
        lStack+=temp;

    }
}

void Calculator::CalculateSemblance(vector<vector<double> *> vCorrected, double &semblance, int min, int max)
{
    double temp1, temp2, sum1, sum2;
    for (int i=0; i<max-min; i++)
    {

        for (int j=0; j<m_nGeo; j++)
        {
            temp1+=(*vCorrected[j])[i];
            temp2+=(*vCorrected[j])[i]*(*vCorrected[j])[i];

        }
        sum1+=temp1*temp1;
        sum2+=temp2;

    }
    semblance = sum1/((max-min)*m_nGeo*sum2);

}



void Calculator::CalSET()
{
    int imax,jmax,kmax, i, j, k;
    imax = (m_GridX2-m_GridX1)/m_StepX;
    jmax = (m_GridY2-m_GridY1)/m_StepY;
    kmax = (m_Z2-m_Z1)/m_StepZ;

    int size = imax*jmax;

    for (int ii=0; ii<1; ii++)
    {
        FILE * pFile = NULL;



        pFile = fopen ( "set.dat" , "wb" );
        if (NULL == pFile)
        {
            printf("File open fail!\n");
            return ;
        }

        for ( k = 0; k<kmax; k++)
        {

            double *buffer = new double[size];

            for ( j = 0; j<jmax; j++)
            {
                for ( i = 0; i<imax; i++)
                {
                    sX = m_GridX1 + i*m_StepX;
                    sY = m_GridY1 + j*m_StepY;
                    sZ = m_Z1 + k*m_StepZ;

                    vector<vector<double>*>  vCorrectedTrace;

                    CalculateTime(vPTime, vSTime, treatWellX, treatWellY, monitorWellX, monitorWellY, perfWellX, perfWellY, devTreatWell, devMonitorWell, devPerfWell, sX, sY, sZ, vTreatWellDepth, vMonitorWellDepth, vTreatWellVp, vTreatWellVs, vMonitorWellVp, vMonitorWellVs, planeA, planeB, planeC, planeD, rX, rY, rZ);

                    vector<int> sample;
                    if (m_nWaveType == 0)
                        Time2Sample(vPTime, sample, m_nSampleRate);
                    else if (m_nWaveType ==1)
                        Time2Sample(vSTime, sample, m_nSampleRate);


                    switch(ii)
                    {
                    case 0:
                        MoveoutCorrection(m_TracesZ,vCorrectedTrace, sample, m_nWinStart,m_nWinEnd, m_nGeo);
                        break;
                    case 1:
                        MoveoutCorrection(m_TracesX,vCorrectedTrace, sample, m_nWinStart,m_nWinEnd, m_nGeo);
                        break;
                    case 2:
                        MoveoutCorrection(m_TracesY,vCorrectedTrace, sample, m_nWinStart,m_nWinEnd, m_nGeo);
                        break;
                    default:
                        break;

                    }

                    double linearStackedValue = 0.0;
                    double semblanceValue=0.0;
                    double semStackedValue=0.0;

                    if (m_nQ2 == 1)
                    {
                        LinearStack(vCorrectedTrace, linearStackedValue,m_nWinStart,m_nWinEnd);
                        CalculateSemblance(vCorrectedTrace, semblanceValue, m_nWinStart,m_nWinEnd);
                        ClearData(vCorrectedTrace);
                        semStackedValue = pow(semblanceValue, double(m_nQ1))*linearStackedValue;
                        buffer[j*imax+i] = semStackedValue;
                    }
                    else if (m_nQ2 == 0)
                    {
                        CalculateSemblance(vCorrectedTrace, semblanceValue, m_nWinStart,m_nWinEnd);
                        ClearData(vCorrectedTrace);
                        buffer[j*imax+i] = semblanceValue;
                    }

                    printf("[k:%d, j:%d, i:%d]...\n", k, j, i);

                }
            }

            fwrite (buffer , sizeof(double), size , pFile );

            delete []buffer;
        }

        if (m_nQ1 == 0 && m_nQ2 ==1)
        {
            printf("Linear stacking calculation completed[OK]\n");
        }
        else if (m_nQ1 != 0 && m_nQ2 ==1)
        {
            printf("Semblance weighted stacking calculation completed[OK]\n");
        }
        else if (m_nQ1 ==1 && m_nQ2 ==0)
        {
            printf("Semblance calculation completed[OK]\n");
        }

        fclose (pFile);
        pFile = NULL;

    }
}


