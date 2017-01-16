#ifndef CSEGY_H
#define CSEGY_H


#include <stdio.h>



//SEGY卷头:3600字节
const int TEXTUAL_HEADER_SIZE = 3200; // ASCII 区域: 3200字节(40条记录 x 80 字节/每条记录)
const int VOLUME_HEADER_SIZE = 400; //二进制数区域: 400字节(3201~3600)

//SEGY道头字区: 含: 60个字/4字节整或120个字/2字节整，共240个字节，按二进制格式存放
const int TRACE_HEADER_SIZE = 240;



struct POINTF
{
    double x;
    double y;
};


//SEGY卷头，标准rev1,IEEE FLOATING POINT，400字节
struct VOLUMEHEADER
{
    int jobid;	/* 3201-3204                作业标识号  */
    int lino;	 	//*线号(每卷只有一个线号) 3205-3208
    int reno;		//*卷号 3209-3212
    short ntrpr;	//*每个道集的数据道数。叠前数据强制要求 3213-3214
    short nart;		//*每个道集的辅助道数。叠前数据强制要求 3215-3216
    unsigned short hdt; 	//*微秒（us）形式的采样间隔。叠前数据强制要求 3217-3218
    unsigned short dto; /* 3219-3220    微秒（us）形式的原始野外记录采样间隔*/
    unsigned short hns; //*数据道采样点数。叠前数据强制要求 3221-3222
    unsigned short nso; /* 原始野外记录每道采样点数 */
    short format;//*1=4字节IBM浮点； 2=4字节，两互补整数； 3=2字节，两互补整数； 4=4字节带增益定点数（过时，不再使用）
    //	5=4字节IEEE浮点数
    //	6=现在没有使用
    //	7=现在没有使用
    //	8=1字节，两互补整数
    short fold;		//*CDP覆盖次数 3227-3228
    short tsort;	//道分选码（即集合类型）：-1=其他（应在用户扩展文件头文本段中解释）
    //	0=未知
    //	1=同记录（未分选）
    //	2=CDP道集
    //	3=单次覆盖连续剖面
    //	4=水平叠加
    //	5=共炮点
    //	6=共接收点
    //	7=共偏移距
    //	8=共中心点
    //	9=共转换点
    //	强烈推荐所有类型的数据使用
    short vscode;	//3231-3232	垂直组合码：1-无组合，2-两道组合，······，N-N道组合
    short hsfs;	/* 3233-3234	起始扫描频率 */
    short hsfe;	/* 3235-3236	终止扫描频率 */
    short hslen;	/* 3237-3238	扫描长度（ms） */
    short hstyp; //3239-3240	扫描类型码：
    short schn;	/* 3241-3242	扫描通道道号 */
    short hstas; //3243-3244	如果使用斜坡，扫描道起始点用的斜坡长（ms）
    short hstae; //3245-3246	扫描结束时用的斜坡长度（ms）
    short htatyp; //3247-3248	斜坡函数类型码：
    short hcorr; //3249-3250	相关的数据道：1-不是；2-是
    short bgrcv; //3251-3252	二进制增益恢复：1-是；2-不是
    short rcvm;	 //3253-3254	振幅恢复方法：1-无 2-AGC 3-球面扩散 4-其他
    short mfeet; //  //计量单位=1 3255-3256
    short polyt; //3257-3258	脉冲信号极性：
    short vpol;	//3259-3260	振动极性码：地震信号延迟控制信号
    short hunass1[120];	/* 3261-3500 unassigned */
    short rev; //3501-3502 SEGY版本号， SEGY Rev1 为1975年标准,2002年修订
    short hunass2[49];	/* 3261-3500 unassigned */



};


//SEGY道头，标准rev1,IEEE FLOATING POINT，240字节
struct TRACEHEADER
{
    int iSLNum;	//1-4	测线内的道顺序号
    int iTSTNum; //5-8	卷内顺序道号，每盘磁带从1开始
    int iOriNum;		// 9-12 **原始的野外记录号（炮号）
    int iTraNum;  //13-16	原始记录内的道号
    int iSorceNum;	 //17-20	爆炸点号，当同一地表位置有多个记录时用此编号
    int iCdpNum;	/* 21-24	CDP号，当记录为CDP道集用此编号 CMP号 */
    int iCdpTra;    //25-28	CDP内道号，每个CDP从1开始
    short nTraId;	   //29-30	道标识码：-1=其他。0=未知。
    //1=地震数据。2=死道。3=哑道。4=时断。5=井口。6=扫描。7=定时。8=水断。9=近场枪信号。
    //10=远场枪信号。11=地震压力传感器。
    //12=多分量地震传感器――垂直分量。13=多分量地震传感器――横向分量。14=多分量地震传感器――纵向分量。
    //15=旋转多分量地震传感器――垂直分量。16=旋转多分量地震传感器――切向分量。17=旋转多分量地震传感器――径向分量。
    //18=可控源反应质量。19=可控源底盘。20=可控源估计地面力。21=可控源参考22=时间速度对。
    //23…N=选用，（最大N=32767）。 强烈推荐所有类型的数据使用。
    short nVsTra;	   // 31-32	产生此道的垂直组合道数（1-1道，2-2道）垂直叠加道数
    short nHsTra;	  //33-34	产生此道的水平叠加道数（1-1道，2-2道）水平叠加道数
    short nDuse;	  //35-36	数据用途：1-生产；2-试验
    int iOffset;  //37-40	从炮点到接收点的距离（如果测线为反方向放炮则为负值）偏移距
    int gelev; //41-44	接收点高程；海平面以上为正，海平面以上为负
    int selev;	//45-48	炮点高程
    int sdepth;	/* 49-52	炮点在地表下的深度（正数） */
    int gdel;	/* 53-56	接收点基准线高程 */
    int sdel;	/* 57-60                    炮点基准线高程 */
    int swdep;	/* 61-64	炮点水深*/
    int gwdep;	/* 65-68	接收点水深 */
    short scalel;//69-70	对字节41-68的所有高程和深度应用的比例
    short scalco;//1-72	对字节73-88内所有确定的坐标应用的比例，
    int  iSourseX;	/* 73-76	炮点坐标-x       如果坐标单位为弧的秒数，则x表示经度- X */
    int  iSourseY;	/* 77-80	炮点坐标-y       如果坐标单位为弧的秒数，则y表示经度，- Y */
    int  iGroupX;	/* 81-81	接收点坐标-x	东赤道以北，负值指示南或西的秒数 - X */
    int  iGroupY;	/* 85-88	接收点坐标- Y */
    short counit; //89-90	坐标单位：1-长度（米或尺）；2-秒
    short wevel;	/* 91-92	风化层速度*/
    short swevel;	/* 93-94	风化层下速度*/
    short sut;	/* 95-96	炮点井口时间 */
    short gut;	/*97-98	接收点井口时间 */
    short sstat;	/* 99-100	炮点静校正量 */
    short gstat;	/* 101-102	接收点静校正量*/
    short tstat; //103-104	应用的总静校正量
    short laga; //105-106	延迟时间A（ms），如果时断信号在头块结束之后出现则为正
    short lagb;	//107-108	延迟时间B（ms），可正可负，是时断信号和爆炸源起爆时间
    short delrt; //109-110	延迟记录时间（ms），是数据采样开始时间和爆炸源起爆时间
    short muts;	/* 111-112	切除开始时间 */
    short mute;	/* 113-114	切除结束时间 */
    unsigned short uNs; //**采样点数 115-116
    unsigned short uRate;//**采样率us 117-118
    short gain;	/* 119-120	野外仪器增益类型*/
    short igc;	/* 121-122	仪器增益常数 */
    short igi;	/* 123-124	仪器初始增益（d B） */
    short corr;	/* 125-126	相关的地震道：1-不是；2-是	*/
    short sfs;	/* 127-128	起始扫描频率*/
    short sfe;	/* 129-130	终止扫描频率 */
    short slen;	/* 131-132	扫描长度（ms） */
    short styp;	/* 133-134	扫描类型1 = linear 2 = cos-squared 3 = other */
    short stas;	/* 135-136	扫描道开始点用的斜坡长度（ms） */
    short stae;	/* 137-138	扫描道结束时用的斜坡长度（ms） */
    short tatyp;	/* 139-140	斜坡类型：1-线型；2-cos2；3-其他 */
    short afilf;	/* 141-142	假频滤波器频率 */
    short afils;	/* 143-144	假频滤波器斜坡 */
    short nofilf;	/* 145-146	陷波滤波器频率 */
    short nofils;	/* 147-148	陷波滤波器斜坡 */
    short lcf;	/*149-150	低截止频率 */
    short hcf;	/* 151-152	高截止频率 */
    short lcs;	/* 153-154	低截止斜坡*/
    short hcs;	/* 155-156	高截止斜坡*/
    short year;	//数据记录的年 157-158
    short day;	//年中的天数 159-160
    short hour;	//小时（24小时制） 161-162
    short minute;	//分 163-164
    short sec;	//秒 165-166
    short timbas;	/*167-168	时间标准码：1-本地时间，2-格林威治时间，3-其他 */
    short trwf;	/* 169-170	地震道加权因子-按2的-N次方所置，对于最低有效位（N=0，1， */
    short grnors;	/* 171-172	滚动开关位置1的检波器组号 */
    short grnofr;	/* 173-174	在原始记录内，道1的检波器组号 */
    short grnlof;	/* 175-176	在原始记录内，最后一道的检波器组号 */
    short gaps;	/* 177-178	间隔大小（去掉的检波器组的总数） */
    short otrav; //179-180	在测线开始或结束超越的斜坡

    int cdpX; //181-184 CDP道X坐标
    int cdpY;//185-188 CDP道Y坐标
    int Inline3D;//189-192  三维地震测线号
    int CrossLine3D;//193-196 三维地震道号
    int ShotPoint;//197-200 炮点号，可能只适用于二维地震
    short ShotPointScalar;//201-202 炮点标量，与197-200相关
    short TraceValueMeasurementUnit;
    /*
     203-204 测量单位：
     -1 = 其它（应在DataSampleMeasurementUnits Stanza中描述）
     0 = 未知
     1 = Pascal(Pa)
     2 = Volts(V)
     3 = Millivolts(mV)
     4 = Amperes(A)
     5 = Meters(m)
     6 = Meters per second(m/s)
     7 = Meters per second squared(m/s2)
     8 = Newton(N)
     9 = Watt(W)
     */

    int TransductionConstantMantissa; //205-208 转换常量
    short TransductionConstantPower;//209-210 转换常量权重
    short TransductionUnit; //211-212 转换单位
    /*
     -1 = 其它（应在DataSampleMeasurementUnits Stanza中描述）
     0 = 未知
     1 = Pascal(Pa)
     2 = Volts(V)
     3 = Millivolts(mV)
     4 = Amperes(A)
     5 = Meters(m)
     6 = Meters per second(m/s)
     7 = Meters per second squared(m/s2)
     8 = Newton(N)
     9 = Watt(W)
     */
    short TraceIdentifier;//213-214 设备/道识别号
    short ScalarTraceHeader; //215-216 道头标量
    short SourceType;//217-218 震源类型
    /*
     -1~-n = 其它（应在SourceType/Orientation stanza中描述）
     0 = 未知
     1 = 振动-垂直方向
     2 = 振动-联络测线方向
     3 = 振动-测线方向
     4 = 脉冲-垂直方向
     5 = 脉冲-联络测线方向
     6 = 脉冲-测线方向
     7 = 分布式脉冲-垂直方向
     8 = 分布式脉冲-联络测线方向
     9 = 分布式脉冲-测线方向
     */
    int SourceEnergyDirectionMantissa;//219-222 相对于炮点方位的炮点能量方向
    short SourceEnergyDirectionExponent;//223-225 炮点方向指数
    int SourceMeasurementMantissa;//225-228 炮点测量
    short SourceMeasurementExponent;//229-230  炮点测量指数
    short SourceMeasurementUnit; //2031-232 炮测量单位
    /*
     -1 = 其它（应在DataSampleMeasurementUnits Stanza中描述）
     0 = 未知
     1 = Joule(J)
     2 = Kilowatt(kW)
     3 = Pascal(Pa)
     4 = Bar(Bar)
     4 = Bar-meter(Bar-m)
     5 = Newton(N)
     6 = Kilograms(kg)
     */

    int UnassignedInt1; //未指定
    int UnassignedInt2;//未指定
};




class CSegy
{

public:
    CSegy();
    ~CSegy();

    bool m_bHasTextualHeader;//判断含文本卷头信息
    bool m_bIsEDBCode;//判断是否EDBC码
    bool m_bIsBigEndian;//判断是否大头

    enum DataFormat
    {
        IBMFLOAT,
        IEEEFLOAT,
        IEEEINT32,
        IEEEINT16,
        IEEEINT8
    };

    void SwapVolumeHead(VOLUMEHEADER &);
    void SwapTraceHead(TRACEHEADER &);

    bool OpenFile(const char *);//打开文件


    const unsigned char *ReadTextualHeader();    //读取3200字节文本信息
    const VOLUMEHEADER &ReadVolumeHeader(); //读取400字节二进制卷头信息
    const TRACEHEADER &ReadTraceHeader(int nIndex);  //读取240字节道头信息
    float * ReadTraceData(int nIndex);  //读取地震道数据

    //float * GetTraceData(){return m_pData;}  //获取地震道数据


    const unsigned char* GetEbcdic() const { return m_szTextInfo;}
    void SetHasTextualHeader() { m_bHasTextualHeader =false;}
    VOLUMEHEADER GetVolumeHeader(){ return m_volHead;}
    TRACEHEADER GetTraceHeader() { return m_traHead;}
    float* GetTraceData(int nIndex );
	void GetAllTraceData(float *m_pData, int nIndex);

    const char* GetTraceHeadAsBytes() { return m_szTraceHead;}

    int getSamplesNumber() const; //获取采样数
    int getSamplesInterval() const; // 获取采样率
    short getFormat() const; //获取数据格式类型
    int getTotalTraceNumber() const; //获取总道数
    int getSegyRev() const; //获取SEGY 版本信息，rev1 or rev0

    float GetSampleValue(int idx);//获取采样点值
    //float GetTrace(long long n);//获取单道数据


    void closeFile();
    void SetNewTextualHeader ();
    void WriteSegy(const char *csReadFile, const char *csWriteFile);



    //char m_csFilepath[300];//文件路径
public:
    char *m_szNewHeader;



private:
    char * csReadFile;  //文件路径
    DataFormat m_fmt;   //1-浮点数（4字节） 2-定点数（4字节） 3-带增益定点数（4字节）
    TRACEHEADER m_traHead;  //240字节道头信息
    VOLUMEHEADER m_volHead;     //400字节二进制卷头信息
    unsigned long long m_lFileLenth;         //文件大小

    unsigned  m_nTraNumPerInLine;//每个记录的道数
    int m_nTraceTolNum; //总道数
    int m_nInLineNum; //线数
    int m_nXLineNum; //道数
    //int m_nCurrentTraceNum;
    int m_nSmpNum;//采样点数

    long lStartPosition; //数据开始长度

    unsigned char m_szTextInfo[TEXTUAL_HEADER_SIZE]; //3200字节文本EBCDIC信息
    char m_szVolumeHead[VOLUME_HEADER_SIZE];//400字节Binary卷头信息
    char m_szTraceHead[TRACE_HEADER_SIZE]; //240字节道头信息

    short m_nFormatCode;//格式码
    short m_nSegyRev;//SEGY 版本号

    int m_nSampleInterval; //采样间隔
    int m_nSampleNum; //采样点数
    int m_nDataLenth; //数据长度
    int m_nStrtInLine; //起始线号
    int m_nStopInLine;//终止线号
    int m_nStrtXLine;//起始道号
    int m_nStopXline;//终止道号
    int m_nStrtTime;
    int m_nStopTime;
    int m_nScaleoSet; //比例因子

    int m_nXCoordLenth; //x坐标长度
    int m_nYCoordLenth;//Y坐标长度
    int m_nXCoordSet; //x坐标位置
    int m_nYCoordSet; //Y坐标位置
    int m_nInlineSet; //inline位置
    int m_nXlineSet; //xline位置
    int m_nInlineLenth; //线号长度
    int m_nXLineLenth; //道号长度
    int m_nTotalTraceNum; //总道数


    float m_fXInterval; //x方向道间距
    float m_fYInterval; //y方向道间距
    float m_fTraceInterval; //道间距
    float m_fLineInterval; //线间距

    bool bDataOrdination;
    bool m_bDataWS; //是否工作站格式
    bool m_bDataIBMFmt; //是否IBM浮点数
    bool m_bCoordUserSet;//是否用户自定义坐标
    bool m_bLineSet; //是否有线道号
    bool m_bCoordSet; //是否有x,y坐标
    bool segyRegular;//是否是规则的segy
    bool m_bEBText; //是否有道头Y
    bool m_bOpen;

    float smp;//采样点值
    float* m_pData; //
    char* m_pszData;


    POINTF P1,P2,P3,P4;

public:

    void SetInfo();
    void ScanData();
    void Scan3DData();


    POINTF ScanPerTrace(int num);

    FILE *m_fpr;
    FILE *m_fprWrite;







};

#endif // CSegy_H
#define prinfErr printf
