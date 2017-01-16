
#include "csegy.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <cstring>
#include <string.h>
#include "util.h"


CSegy::CSegy()
{
    m_bDataIBMFmt = true;
    m_bDataWS = true;
    m_bEBText = true; //是否有道头Y = true;
    m_nDataLenth = 4;
    m_nXCoordLenth = 4;
    m_nXCoordSet = 73;
    m_nYCoordSet = 77;
    m_nInlineSet = 9; //inline位置
    m_nXlineSet = 21; //xl
    m_nInlineLenth = 4; //线号长度
    m_nXLineLenth = 4; //道号长度
    m_nXCoordLenth = 4; //X坐标长度
    m_nYCoordLenth = 4; //Y坐标长度
    m_nSampleInterval = 0;
    m_nSampleNum = 0;
    m_fmt = IBMFLOAT;
    bDataOrdination = false;
    m_nStrtTime = 0;
    m_nStrtXLine = 0;
    m_nStrtInLine = 0;
    m_nStopXline = 0;
    m_nStopTime = 0;
    m_nStopInLine = 0;
    m_nScaleoSet = 71;

    m_fXInterval = 0.f;
    m_fLineInterval = 0.f;

    m_bIsEDBCode = true;
    m_bIsBigEndian = true;
    m_bHasTextualHeader = true;
    m_bOpen = false;

    m_bLineSet = true;
    m_bCoordSet = true;
    segyRegular = false;

    m_pData = 0;
    m_pszData = 0;

    m_fpr = NULL;

}

CSegy::~CSegy()
{
    fclose(m_fpr);
    if (m_pData)
        delete []m_pData;
    if (m_pszData)
        delete []m_pszData;
}


const unsigned char *CSegy::ReadTextualHeader()
{
    fseek(m_fpr,0,0);
    if (m_bHasTextualHeader)
    {

        fread(m_szTextInfo,sizeof(unsigned char),TEXTUAL_HEADER_SIZE,m_fpr);

    }
    else
    {

        char* temTextInfo  = new char[3200];
        temTextInfo = "This SEG-Y file has no 3200-byte EBCDIC textual header!";
        memcpy(m_szTextInfo, temTextInfo,3200);
    }

    return m_szTextInfo;
}

const VOLUMEHEADER &CSegy::ReadVolumeHeader()
{
    fseek(m_fpr,0,0);

    unsigned long long lenth;

    if (m_bHasTextualHeader)
    {
        fseek(m_fpr,TEXTUAL_HEADER_SIZE,0);
        lenth = m_lFileLenth;

    }
    else lenth = m_lFileLenth-TEXTUAL_HEADER_SIZE;


    fread(&m_volHead,sizeof(VOLUMEHEADER),1,m_fpr);


    if (m_bDataWS)
        SwapVolumeHead(m_volHead);

    if (m_nSampleNum != 0 && m_nXLineNum != 0)
        m_nInLineNum = lenth/(m_nSampleNum*m_nXLineNum*m_nDataLenth);
    if (m_nInLineNum < 0)
        m_nInLineNum = 0;

    return m_volHead;
}

const TRACEHEADER &CSegy::ReadTraceHeader(int nIndex/*=0*/)
{
    long long lOffset = TEXTUAL_HEADER_SIZE+VOLUME_HEADER_SIZE;
    if (!m_bHasTextualHeader)
        lOffset = VOLUME_HEADER_SIZE;
    lOffset += (TRACE_HEADER_SIZE+m_nSampleNum*m_nDataLenth)*(nIndex-1);
    fseek(m_fpr,lOffset,SEEK_SET);
    fread(m_szTraceHead,sizeof(char),TRACE_HEADER_SIZE,m_fpr);
    memcpy((char*)&m_traHead,m_szTraceHead,sizeof(char)*TRACE_HEADER_SIZE);

    if (m_bDataWS)
        SwapTraceHead(m_traHead);

    return m_traHead;
}

float* CSegy::ReadTraceData(int nIndex)
{
    //if (nIndex < 0)
    //    nIndex = m_nCurrentTraceNum;
    long long lOffset = TEXTUAL_HEADER_SIZE+VOLUME_HEADER_SIZE;
    if (!m_bHasTextualHeader)
        lOffset = VOLUME_HEADER_SIZE;
    lOffset += (nIndex-1)*(TRACE_HEADER_SIZE+m_nSampleNum*m_nDataLenth)+TRACE_HEADER_SIZE;
    fseek(m_fpr,lOffset,0);

    if (m_pData != 0)
        delete []m_pData;
    if (m_pszData != 0)
        delete []m_pszData;

    m_pszData = new char[m_nSampleNum*m_nDataLenth];
    m_pData = new float[m_nSampleNum];
    fread(m_pszData,sizeof(char),m_nSampleNum*m_nDataLenth,m_fpr);


    for (int i=0;i<m_nSampleNum;i++)
    {
        memcpy((char*)m_pData+i*sizeof(float),m_pszData+i*m_nDataLenth,m_nDataLenth);
        if (m_bDataIBMFmt)
            ibm2ieee(&m_pData[i],m_bDataWS);

    }



    return m_pData;


}

float *formatTransform(const char *szBuff,int len,int format)
{
    float *pData = new float[len];
    switch (format)
    {
    case 1:
    {
        for (int i=0;i<len;i++)
        {
            float x;
            memcpy((char*)&x,szBuff+i*4,4);
            ibm2ieee(&x,1);

            pData[i] = x;
        }
    }
        break;
    case 2:
    {
        for (int i=0;i<len;i++)
        {
            int x;
            memcpy((char*)&x,szBuff+i*4,4);
            pData[i] = swapi4(x);
        }
    }
        break;
    case 3:
    {
        for (int i=0;i<len;i++)
        {
            int x;
            memcpy((char*)&x,szBuff+i*2,2);
            pData[i] = swapi2(x);
        }
    }


        break;
    case 5:
    {
        for (int i=0;i<len;i++)
        {
            int x;
            memcpy((char*)&x,szBuff+i*4,4);
            pData[i] = swapi4(x);
        }
    }
        break;
    }

    return pData;
}

float* CSegy::GetTraceData(int nIndex )
{

    int bufLen = m_nSampleNum*m_nDataLenth;

    char *szBuff = new char[bufLen];

    unsigned long long lOffset = TEXTUAL_HEADER_SIZE+VOLUME_HEADER_SIZE;

    if (!m_bHasTextualHeader)
        lOffset =VOLUME_HEADER_SIZE;
    lOffset +=(nIndex-1)*(TRACE_HEADER_SIZE+bufLen)+TRACE_HEADER_SIZE;
    fseek(m_fpr,lOffset,0);

    fread(szBuff,sizeof(char),bufLen,m_fpr);

    int len = m_nSampleNum;
    float *m_pData = new float[len];
    memset((char*)m_pData,0,len*sizeof(float));

    switch (m_nFormatCode)
    {
    case 1:
    {
        for (int i=0;i<len;i++)
        {
            float x;
            memcpy((char*)&x,szBuff+i*4,4);
            ibm2ieee(&x,1);
            m_pData[i] = x;
        }
    }
        break;
    case 2:
    {
        for (int i=0;i<len;i++)
        {
            int x;
            memcpy((char*)&x,szBuff+i*4,4);
            m_pData[i] = swapi4(x);
        }
    }
        break;
    case 3:
    {
        for (int i=0;i<len;i++)
        {
            int x;
            memcpy((char*)&x,szBuff+i*2,2);
            m_pData[i] = swapi2(x);
        }
    }


        break;
    case 5:
    {
        for (int i=0;i<len;i++)
        {
            float x;
            memcpy((char*)&x,szBuff+i*4,4);
            m_pData[i] = swapf4(x);

        }
    }
        break;

    }



    delete []szBuff;
    return (float*)m_pData;
}

void CSegy::GetAllTraceData(float *m_pData, int nIndex)
{
	int bufLen = m_nSampleNum*m_nDataLenth;

	char *szBuff = new char[bufLen];

	unsigned long long lOffset = TEXTUAL_HEADER_SIZE + VOLUME_HEADER_SIZE;

	if (!m_bHasTextualHeader)
		lOffset = VOLUME_HEADER_SIZE;
	lOffset += (nIndex - 1)*(TRACE_HEADER_SIZE + bufLen) + TRACE_HEADER_SIZE;
	fseek(m_fpr, lOffset, 0);

	fread(szBuff, sizeof(char), bufLen, m_fpr);

	int len = m_nSampleNum;
	
	switch (m_nFormatCode)
	{
	case 1:
	{
		for (int i = 0; i<len; i++)
		{
			float x;
			memcpy((char*)&x, szBuff + i * 4, 4);
			ibm2ieee(&x, 1);
			m_pData[i] = x;
		}
	}
	break;
	case 2:
	{
		for (int i = 0; i<len; i++)
		{
			int x;
			memcpy((char*)&x, szBuff + i * 4, 4);
			m_pData[i] = swapi4(x);
		}
	}
	break;
	case 3:
	{
		for (int i = 0; i<len; i++)
		{
			int x;
			memcpy((char*)&x, szBuff + i * 2, 2);
			m_pData[i] = swapi2(x);
		}
	}


	break;
	case 5:
	{
		for (int i = 0; i<len; i++)
		{
			float x;
			memcpy((char*)&x, szBuff + i * 4, 4);
			m_pData[i] = swapf4(x);

		}
	}
	break;

	}

	
	delete[]szBuff;
	
}

float CSegy::GetSampleValue(int idx)
{

    m_pszData = new char[m_nSampleNum*m_nDataLenth];
    m_pData = new float[m_nSampleNum];
    int*    xmp = (int*)&smp;
    float*  dat = (float*)(m_pszData+240);

    short*       int2ptr = (short*)       dat;
    int*         int4ptr = (int*)         dat;
    signed char* int1ptr = (signed char*) dat;

    if(m_nFormatCode==1)
    {
        smp =  dat[idx];
        ibm2ieee(&smp,1);
    }
    else if(m_nFormatCode==2)
    {
        smp = swapi4(int4ptr[idx]);
    }
    else if(m_nFormatCode==3)
    {
        smp = swapi2(int2ptr[idx]);
    }
    else if(m_nFormatCode==4)
    {
        smp = dat[idx];
    }
    else if(m_nFormatCode==5)
    {
        *xmp = swapi4(int4ptr[idx]);
    }
    else if(m_nFormatCode==6)
    {
        smp = int1ptr[idx];
    }
    else
        smp = 0;

    return smp;
}

/*
float CSegy::GetTrace(long long n)
{
    long long l;
    int trl = 240+m_nTotalTraceNum*m_nDataLenth;

    //if(_F == NULL) return false;

    if(n == 0) return false;

    if(n < 1 || n > m_nTotalTraceNum)
    {
        printf("Bad Trace Number");
        return false;
    }

    fseek(m_fpr, 3600+(n-1)*trl, SEEK_SET);
    l = fread(_INPTRC, 1, _TRL, _F );

    _CurrentTrace = n;

    for(n=0; n < _TotalSamples; n++)
        m_pData[n] =  GetSampleValue(n);

    return m_pData;
}
*/


bool CSegy::OpenFile(const char *csReadFile)
{

    closeFile();


    m_fpr = fopen(csReadFile,"rb");

    if (m_fpr == NULL)
    {
        prinfErr("文件打开失败");
        return false ;
    }

    //strcpy(m_csFilepath,csReadFile);
    fseek(m_fpr,0,2);

    m_lFileLenth = ftell(m_fpr);

    fseek(m_fpr,0,0);


    ReadTextualHeader();
    ReadVolumeHeader();
    ReadTraceHeader(1);
    SetInfo();

    //ReadTraceData(20);
    //ScanData();

    return true;
}

void CSegy::WriteSegy(const char *csReadFile, const char *csWriteFile)
{

    m_fpr = fopen(csReadFile,"rb");
    m_fprWrite = fopen(csWriteFile,"wb");

    if (m_fpr == NULL)
    {
        prinfErr("文件打开失败");
        return ;
    }

    //写入3200字节文本信息
    fseek(m_fprWrite,0,0);
    for (int i=0; i<3200; i++)
    {
        fputc(m_szNewHeader[i], m_fprWrite);
    }


    char ch;
    fseek(m_fpr,3200,0);
    fseek(m_fprWrite,3200,0);
    while(!feof(m_fpr))
    {
        ch = fgetc(m_fpr);
        fputc(ch, m_fprWrite);

    }

    fclose(m_fprWrite);


}

void CSegy::closeFile()
{
    if (m_fpr != NULL)
    {
        fclose(m_fpr);
        m_fpr = 0;
    }
}

void CSegy::SetNewTextualHeader()
{

    m_szNewHeader = "This is a SEG-Y file containing microseismic event signals"; //3200字节文本EBCDIC信息,用于更改文件头
}



void CSegy::SetInfo()
{
    m_nFormatCode = 5;//m_volHead.format; //设置数据类型
    m_nSegyRev = m_volHead.rev; //设置SEGY版本号

    m_nSampleNum = m_volHead.hns;
    m_nSampleInterval = m_volHead.hdt;

    //m_fmt = m_volHead.format == 1?IBMFLOAT:IEEEFLOAT;
    m_nStrtTime = m_traHead.delrt;
    m_nStopTime = m_traHead.delrt + m_nSampleInterval*m_nSampleNum / 1000;
    m_nXLineNum = m_volHead.ntrpr;

    //重置数据长度
    switch(m_nFormatCode)
    {
    case 1:
    {
        m_nDataLenth = 4;
        break;
    }

    case 2:
    {
        m_nDataLenth = 4;
        break;
    }
    case 3:
    {
        m_nDataLenth = 2;
        break;
    }
    case 5:
    {
        m_nDataLenth = 4;
        break;
    }
    case 8:
    {
        m_nDataLenth = 1;
        break;
    }
    }



    unsigned long long len = m_lFileLenth;
    if (m_bHasTextualHeader)
    {
        len = len-TEXTUAL_HEADER_SIZE-VOLUME_HEADER_SIZE;
    }
    else
    {
        len = len -VOLUME_HEADER_SIZE;
    }

    m_nTotalTraceNum = len/(TRACE_HEADER_SIZE+m_nSampleNum*m_nDataLenth);
}

void CSegy::ScanData()
{
    Scan3DData();
}

void CSegy::Scan3DData()
{
    m_nInLineNum = 0;
    m_nXLineNum = 0;

    unsigned long long len = m_lFileLenth;
    if (m_bHasTextualHeader)
        len = len-TEXTUAL_HEADER_SIZE-VOLUME_HEADER_SIZE;
    m_nTotalTraceNum = len/(TRACE_HEADER_SIZE+m_nSampleNum*m_nDataLenth);

    m_nInLineNum  = 1;

    P1 = ScanPerTrace(0);
    m_nStrtInLine = m_traHead.iTSTNum;
    m_nStrtXLine = m_traHead.iCdpNum;
    m_nStrtTime = m_traHead.delrt;
    m_nStopTime = m_nStrtTime + m_nSampleNum*m_nSampleInterval/1000;

    int flag = m_traHead.iOriNum;
    for (int i=1;i<m_nTotalTraceNum;i++)
    {
        POINTF pt = ScanPerTrace(i);
        if (flag != m_traHead.iOriNum) //扫描道数的规则 再仔细斟酌
        {
            m_nXLineNum = i;
            break;
        }
    }
    m_nInLineNum = m_nTotalTraceNum / m_nXLineNum;
    m_nStopInLine = m_nStrtInLine + m_nInLineNum;
    P2 = ScanPerTrace(m_nXLineNum-1);
    m_nStopXline = m_traHead.iCdpNum;
    P3 = ScanPerTrace(m_nTotalTraceNum-1);
    P4 = ScanPerTrace(m_nTotalTraceNum-m_nXLineNum);
}

POINTF CSegy::ScanPerTrace(int num)
{
    unsigned long long lOffset = TEXTUAL_HEADER_SIZE+VOLUME_HEADER_SIZE;
    if (!m_bHasTextualHeader)
        lOffset = 0;
    lOffset += (TRACE_HEADER_SIZE+m_nSampleNum*m_nDataLenth)*(unsigned long long)num;

    char sz[TRACE_HEADER_SIZE];
    int nLineFlag = -1;
    int scale;

    fseek(m_fpr,lOffset,0);
    fread(sz,TRACE_HEADER_SIZE,1,m_fpr);
    memcpy((char*)&m_traHead,sz,TRACE_HEADER_SIZE);

    memcpy((char*)&m_traHead.iTSTNum,sz+m_nInlineSet-1,m_nInlineLenth);
    memcpy((char*)&m_traHead.iCdpNum,sz+m_nXlineSet-1,m_nXLineLenth);
    memcpy((char*)&scale,sz+m_nScaleoSet-1,sizeof(short));
    memcpy((char*)&m_traHead.iGroupX,sz+m_nXCoordSet-1,m_nXCoordLenth);
    memcpy((char*)&m_traHead.iGroupY,sz+m_nYCoordSet-1,m_nYCoordLenth);

    if (m_bDataWS)
    {
        SwapTraceHead(m_traHead);
        scale = swapi2(scale);
    }

    scale = abs(scale);
    if (scale==0)
        scale = 1;

    POINTF pt;
    pt.x = (double)m_traHead.iGroupX/(double)scale;
    pt.y = (double)m_traHead.iGroupY/(double)scale;

    return pt;
}

int CSegy::getSamplesNumber() const
{
    return m_nSampleNum;
}

int CSegy::getSamplesInterval() const
{
    return m_nSampleInterval;
}

short CSegy::getFormat() const
{
    return m_nFormatCode;
}



int CSegy::getTotalTraceNumber() const
{
    return m_nTotalTraceNum;
}

int CSegy::getSegyRev() const
{
    return m_nSegyRev;
}



/*
const unsigned char *CSegy::GetEbcdicTextualHeader() const
{
    fseek(m_fpr,0,0);
    if (m_bHasTextualHeader)
    {
        fread(m_szTextInfo,sizeof(unsigned char),TEXTUAL_HEADER_SIZE,m_fpr);
        if (m_bIsEDBCode)
        {
            for (int i=0;i<TEXTUAL_HEADER_SIZE;i++)
                m_szTextInfo[i] = EbcdicToAscii(m_szTextInfo[i]);
        }
    }

    return m_szTextInfo;
}
*/


void CSegy::SwapVolumeHead(VOLUMEHEADER &m_volHead)
{
    m_volHead.jobid   = swapi4(m_volHead.jobid);
    m_volHead.lino    = swapi4(m_volHead.lino);
    m_volHead.reno    = swapi4(m_volHead.reno);
    m_volHead.ntrpr   = swapi2(m_volHead.ntrpr);
    m_volHead.nart    = swapi2(m_volHead.nart);
    m_volHead.hdt     = swapi2(m_volHead.hdt);
    m_volHead.hns     = swapi2(m_volHead.hns);
    m_volHead.nso     = swapi2(m_volHead.nso);
    m_volHead.format  = swapi2(m_volHead.format);
    m_volHead.fold    = swapi2(m_volHead.fold);
    m_volHead.tsort   = swapi2(m_volHead.tsort);
    m_volHead.vscode  = swapi2(m_volHead.vscode);
    m_volHead.hsfs    = swapi2(m_volHead.hsfs);
    m_volHead.hsfe    = swapi2(m_volHead.hsfe);
    m_volHead.hslen   = swapi2(m_volHead.hslen);
    m_volHead.hstyp   = swapi2(m_volHead.hstyp);
    m_volHead.schn    = swapi2(m_volHead.schn);
    m_volHead.hstas   = swapi2(m_volHead.hstas);
    m_volHead.hstae   = swapi2(m_volHead.hstae);
    m_volHead.htatyp  = swapi2(m_volHead.htatyp);
    m_volHead.hcorr   = swapi2(m_volHead.hcorr);
    m_volHead.bgrcv   = swapi2(m_volHead.bgrcv);
    m_volHead.rcvm    = swapi2(m_volHead.rcvm);
    m_volHead.mfeet   = swapi2(m_volHead.mfeet);
    m_volHead.polyt   = swapi2(m_volHead.polyt);
    m_volHead.vpol    = swapi2(m_volHead.vpol);
}

void CSegy::SwapTraceHead(TRACEHEADER &m_traHead)
{
    m_traHead.iSLNum    = swapi4(m_traHead.iSLNum);
    m_traHead.iTSTNum   = swapi4(m_traHead.iTSTNum);
    m_traHead.iOriNum   = swapi4(m_traHead.iOriNum);
    m_traHead.iTraNum   = swapi4(m_traHead.iTraNum);
    m_traHead.iSorceNum = swapi4(m_traHead.iSorceNum);
    m_traHead.iCdpNum   = swapi4(m_traHead.iCdpNum);
    m_traHead.iCdpTra   = swapi4(m_traHead.iCdpTra);
    m_traHead.nTraId    = swapi2(m_traHead.nTraId);
    m_traHead.nVsTra    = swapi2(m_traHead.nVsTra);
    m_traHead.nHsTra    = swapi2(m_traHead.nHsTra);
    m_traHead.nDuse     = swapi2(m_traHead.nDuse);
    m_traHead.iOffset   = swapi4(m_traHead.iOffset);
    m_traHead.gelev     = swapi4(m_traHead.gelev);
    m_traHead.selev     = swapi4(m_traHead.selev);
    m_traHead.sdepth    = swapi4(m_traHead.sdepth);
    m_traHead.gdel      = swapi4(m_traHead.gdel);
    m_traHead.sdel      = swapi4(m_traHead.sdel);
    m_traHead.swdep     = swapi4(m_traHead.swdep);
    m_traHead.gwdep     = swapi4(m_traHead.gwdep);
    m_traHead.scalel    = swapi2(m_traHead.scalel);
    m_traHead.scalco    = swapi2(m_traHead.scalco);
    m_traHead.iSourseX  = swapi4(m_traHead.iSourseX);
    m_traHead.iSourseY  = swapi4(m_traHead.iSourseY);
    m_traHead.iGroupX   = swapi4(m_traHead.iGroupX);
    m_traHead.iGroupY   = swapi4(m_traHead.iGroupY);
    m_traHead.counit    = swapi2(m_traHead.counit);
    m_traHead.wevel     = swapi2(m_traHead.wevel);
    m_traHead.swevel    = swapi2(m_traHead.swevel);
    m_traHead.sut       = swapi2(m_traHead.sut);
    m_traHead.gut       = swapi2(m_traHead.gut);
    m_traHead.sstat     = swapi2(m_traHead.sstat);
    m_traHead.gstat     = swapi2(m_traHead.gstat);
    m_traHead.tstat     = swapi2(m_traHead.tstat);
    m_traHead.laga      = swapi2(m_traHead.laga);
    m_traHead.lagb      = swapi2(m_traHead.lagb);
    m_traHead.delrt     = swapi2(m_traHead.delrt);
    m_traHead.muts      = swapi2(m_traHead.muts);
    m_traHead.mute      = swapi2(m_traHead.mute);
    m_traHead.uNs       = swapi2(m_traHead.uNs);
    m_traHead.uRate     = swapi2(m_traHead.uRate);
    m_traHead.gain      = swapi2(m_traHead.gain);
    m_traHead.igc       = swapi2(m_traHead.igc);
    m_traHead.igi       = swapi2(m_traHead.igi);
    m_traHead.corr      = swapi2(m_traHead.corr);
    m_traHead.sfs       = swapi2(m_traHead.sfs);
    m_traHead.sfe       = swapi2(m_traHead.sfe);
    m_traHead.slen      = swapi2(m_traHead.slen);
    m_traHead.styp      = swapi2(m_traHead.styp);
    m_traHead.stas      = swapi2(m_traHead.stas);
    m_traHead.stae      = swapi2(m_traHead.stae);
    m_traHead.tatyp     = swapi2(m_traHead.tatyp);
    m_traHead.afilf     = swapi2(m_traHead.afilf);
    m_traHead.nofilf    = swapi2(m_traHead.nofilf);
    m_traHead.nofils    = swapi2(m_traHead.nofils);
    m_traHead.lcf       = swapi2(m_traHead.lcf);
    m_traHead.hcf       = swapi2(m_traHead.hcf);
    m_traHead.lcs       = swapi2(m_traHead.lcs);
    m_traHead.hcs       = swapi2(m_traHead.hcs);
    m_traHead.year      = swapi2(m_traHead.year);
    m_traHead.day       = swapi2(m_traHead.day);
    m_traHead.hour      = swapi2(m_traHead.hour);
    m_traHead.minute    = swapi2(m_traHead.minute);
    m_traHead.sec       = swapi2(m_traHead.sec);
    m_traHead.timbas    = swapi2(m_traHead.timbas);
    m_traHead.trwf      = swapi2(m_traHead.trwf);
    m_traHead.grnors    = swapi2(m_traHead.grnors);
    m_traHead.grnofr    = swapi2(m_traHead.grnofr);
    m_traHead.grnlof    = swapi2(m_traHead.grnlof);
    m_traHead.gaps      = swapi2(m_traHead.gaps);
    m_traHead.otrav     = swapi2(m_traHead.otrav);
}
