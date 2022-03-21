#ifndef _SMS_H
#define _SMS_H

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(WIN64) || defined(__WIN64__) || defined(_WIN64)
#define SMSWINDOWS
#define _CRT_SECURE_NO_DEPRECATE
#endif

#if defined (__linux) || defined (linux) || defined(__linux__)
#define SMSLINUX
#endif

#if defined (__MACH) || defined (MACH) || defined(__MACH__)
#define SMSMAC
#endif

#include <stdio.h>
#include <stdlib.h>
// added my self stdint.h
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <string.h>

#ifdef SMSWINDOWS
#include "direct.h"
#endif 

#ifdef SMSLINUX
#include <unistd.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <dlfcn.h>
#endif
 
#ifdef SMSMAC
#include <dlfcn.h>
#endif 

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b)) 
#endif

#ifdef SMSWINDOWS
#define TIMEB _timeb
#define ISNAN _isnan
#define FINITE _finite
#define DLLEXPORT __declspec(dllexport)
#define FTIME _ftime
#define GETCWD _getcwd
#define putw _putw
#define getw _getw
#endif

#ifdef SMSLINUX
#define GETCWD getcwd
#define _copysign copysign
#define FTIME ftime
#define TIMEB timeb
#define _MAX_PATH 4096
#define _isnan isnan
#define _finite finite
#define ISNAN isnan
#define FINITE finite
#define CALLBACK 
#define DLLEXPORT
#define _inline
#endif 

#ifdef SMSMAC
#define GETCWD getcwd
#define _copysign copysign
#define FTIME ftime
#define TIMEB timeb
#define _MAX_PATH 4096
#define _isnan isnan
#define _finite finite
#define ISNAN isnan
#define FINITE finite
#define FINITE finite
#define CALLBACK 
#define DLLEXPORT
#define _inline
#endif 

#define ID_IDataLength          0
#define ID_RDataLength          1
#define ID_IDataLast            2
#define ID_RDataLast            3
#define ID_feap1                4
#define ID_Iteration            5
#define ID_TotalIteration       6
#define ID_LinearEstimate       7
#define ID_ErrorStatus          8
#define ID_MaterialState        9
#define ID_NoSensParameters     10
#define ID_ElementShape         11
#define ID_SensIndex            12
#define ID_OutputFile           13
#define ID_MissingSubroutine    14
#define ID_SubDivergence        15
#define ID_ElementState         16
#define ID_feap2                17
#define ID_feap3                18
#define ID_NoNodes              19
#define ID_NoElements           20
#define ID_NoESpec              21
#define ID_MaxEquations         22 /*maximum number of equations - no constraints*/
#define ID_Debug                23
#define ID_NoDimensions         24
#define ID_SymmetricTangent     25 
#define ID_MinNoTmpData         26
#define ID_NoEquations          27
#define ID_NegativePivots       28
#define ID_Task                 29
#define ID_NoSubIterations      30
#define ID_CurrentElement       31  /* 1,2,..*/
#define ID_MaxPhysicalState     32
#define ID_ExtrapolationType    33
#define ID_TmpContents          34
#define ID_AssemblyNodeResidual 35
#define ID_NodeReordering       36
#define ID_SkipSolver           37
#define ID_NoNSpec              38
#define ID_SetSolver            39
#define ID_NoDesignVelocityFields    40 /* number of design velocity fields: user defined + 1 + 0 */
#define ID_GeometricTangentMatrix  41
#define ID_DataMemory           42
#define ID_SolverMemory         43
#define ID_Solver               44
#define ID_ErrorElement         45
#define ID_SkipTangent          46
#define ID_SkipResidual         47
#define ID_SubIterationMode     48
#define ID_Solver1              49
#define ID_Solver2              50
#define ID_Solver3              51
#define ID_Solver4              52
#define ID_Solver5              53
#define ID_ContactProblem       54
#define ID_Contact1             55
#define ID_Contact2             56
#define ID_Contact3             57
#define ID_Contact4             58
#define ID_Contact5             59
#define ID_PostIteration        60
#define ID_DummyNodes           61
#define ID_PostIterationCall    62
#define ID_Step                 63
#define ID_DebugElement         64
#define ID_ZeroPivots           65
#define ID_GlobalIterationMode  66
#define ID_NoDiscreteEvents     67
#define ID_LineSearchUpdate     68
#define ID_NoBackStep           69  /* number of failed iterative solution steps */
#define ID_SkipSubIteration     70
#define ID_HideIndeterminate    71
#define ID_NoMessages           72
#define ID_MaxMessages          73
#define ID_NoThreads            74
#define ID_NoMultipliers        75   /* number of boundary conditions multipliers*/
#define ID_NoBCFields       76   /* number of bounday conditions design sensitivity fields*/
#define ID_NoLinearConstraints  77   /* number of linear contrain equations */
#define ID_DOFScaling           78   /* weather DOF are scaled in the case of multifield problems*/
#define ID_SensIndexStart       79   /* starting sensitivity  index in a group calculated together*/
#define ID_SensIndexEnd         80   /* end sensitivity  index in a group calculated together*/
#define ID_SensMaxGroupLength   81   /* maximal number of sensitivity parameters calculated together - for dimensioning of working fields*/
#define ID_NoFirstOrderParameters 82
#define ID_NoSecondOrderParameters 83
#define ID_NoThirdOrderParameters 84
#define ID_TangentMatrixMemory  85
#define ID_ArcLengthStatus      86   /* 0 -arclength is not active */
#define ID_ArcLengthPredictorSignCriterion  87 /* 1-> Sign=Sign(INCRp*up) - by Feng 2-> Sign=Sign(det[Kt]) */
#define ID_ArcLengthMethod      88   /* 1 - Ordinary Arc Length analysis 2 - Arc Length analysis by Pohl */
#define ID_NoTrueNodes          89   /* number of true nodes NoNodes-number of dummy nodes
                                        node selectors and operations only operates on true nodes */
#define ID_Solver6              90
#define ID_SolverIDataLength    91  /*actual length of SolverIData is NoSolverGlobalConst+length of additional i parameters*/
#define ID_SolverRDataLength    92  /*actual length of SolverRData is NoSolverGlobalConst+length of additional a parameters*/
#define ID_NoSensProblems       93  /*number of solved sensitvity problems*/
#define ID_NoShapeFields        94
#define ID_NoParameterFields    95
#define ID_VVectorSize       96
#define IData_Last              97   /* used space in IData field*/
#define IData_Length            400  /* length of idata vector */

#define ERROR_NoError            0
#define ERROR_Warning            1
#define ERROR_StepCut            2
#define ERROR_Fatal              3

#define TASK_Quit                   -1
#define TASK_User                   0
#define TASK_TangentAndResidual     1
#define TASK_Residual               2
#define TASK_SensitivityPseudoLoad  3
#define TASK_DependentSensitivity   4
#define TASK_Postprocessing         5
#define TASK_Tangent                6
#define TASK_GlobalInitialisation   7
#define TASK_ElementInitialisation  8
#define TASK_NodePositionsNormals   9

#define TMPContents_ZERO            0
#define TMPContents_Residual        1
#define TMPContents_Post            2

#define SOLVER_FALSE         0
#define SOLVER_LU            1
#define SOLVER_LDL           2
#define SOLVER_SLU           3
#define SOLVER_UMF           4
#define SOLVER_PARDISO       5
#define SOLVER_ITERATIVE     6

#define RD_Multiplier            0
#define RD_ResidualError         1
#define RD_IncrementError        2
#define RD_MFlops                3
#define RD_SubMFlops             4
#define RD_Time                  5
#define RD_TimeIncrement         6
#define RD_MultiplierIncrement   7
#define RD_SubIterationTolerance 8
#define RD_LineSearchStepLength  9
#define RD_PostMaxValue          10
#define RD_PostMinValue          11
#define RD_TotalEnergy           12
#define RD_PotentialEnergy       13
#define RD_KineticEnergy         14
#define RD_Solver1               15
#define RD_Solver2               16
#define RD_Solver3               17
#define RD_Solver4               18
#define RD_Solver5               19
#define RD_Contact1              20
#define RD_Contact2              21
#define RD_Contact3              22
#define RD_Contact4              23
#define RD_Contact5              24
#define RD_ContactSearchTolerance 25
#define RD_MinX                  26
#define RD_MinY                  27
#define RD_MinZ                  28
#define RD_MaxX                  29
#define RD_MaxY                  30
#define RD_MaxZ                  31
#define RD_XYZRange              32
#define RD_KAndRTime             33
#define RD_SolverTime            34
#define RD_PreviousMultiplierIncrement 35
#define RD_SpatialTolerance      36
#define RD_Parameter             37
#define RD_ParameterIncrement    38
#define RD_ArcLengthPsi1          39 /*scaling for natuaral BC  norm*/
#define RD_ArcLengthPsi2          40 /*scaling for essential BC  norm*/
#define RD_ArcLengthDirectionCheck  41 /*chect for false direction if RData[RD_IncrementError]<tol*/
#define RD_ArcLengthDirectionTolerance  42 /*direction is false when cos fi < -1+ DirectionTolerance*/
#define RD_MinElementSize        43
#define RD_Solver6               44
#define RD_SensRTime             45  /*time to evaluate sensitivity pseudo-load vectors and resolve locally coupled sensistivity problems*/
#define RD_SensSolTime           46  /*time to solve right-hand sides of sensitvity problem*/
#define RD_ProfileTime           47  /* slot for profiling - non specific*/
#define RData_Last               48  /* used space in RData field*/
#define RData_Length             400 /* length of rdata vector */

/*comptibility with old C source codes before version 6.2*/
#define NoGroupData NoDomainData
#define GroupDataNames DomainDataNames

typedef struct{
    struct {
        int NodeIndex; /*0,1,2,...*/
        int NoDOF;
        int SpecIndex; /*0,1,...*/
        int NoElements; /* set by SetSolver */
    } id;
    int *DOF;  /* 0,1,... active DOF
                  -1 constrained by input data
                  -2 constrained by user after SMTAnalysis (for e.g. staggered schemes)
                  -3 constrained automatically due to the deactivation of NodeSpec - set by SMTSetSolver
               */
    int *Elements; /* pointer set by SetSolver 1,2,3... */
    double *X;
    double *Bt;   /* holds data fields Bt */
    double *Bp;   /* holds data fields Bp */
    double *dB;   /* holds data fields dB */
    double *at;
    double *ap;
    double *da;
    double *st;  /*sens in time t  dof/par*/
    double *sp;
    double *ht;
    double *hp;
    double *Data; /* holds data fields Data,ADVF length NoNodeData+NoDesignVelocityFields*/
    double *tmp;
}NodeData;

typedef struct{
    struct{
        int SpecIndex;      /*1 0,1,...*/
        int NoDOF;          /*2*/
        int NoNodeStorage;  /*3 length of node data ht and hp field*/
        int NoNodeData;     /*4 lenght of node data Data field */
        int NoData;         /*5 length of node spec Data field */
        int NoTmpData;      /*6 max(SMTIData["MinNoTmpData"],NoDOF)*/
        int Constrained;    /*7 */
        int Fictive;        /*8* not a topological node*/
        int Dummy;          /*9  True - specification is a dummy node specification*/
        int DummyNode;      /*10 index of dummy node if exists 1,2,...*/
        int Active;         /*11 0=>freeze all nodes with given spec */
        int dummy1;         /* compatibility */
        int dummy2;
        int dummy3;
        int dummy4;
    }id;
    double *Data;
    char *NodeID;
   /* zero Bt, Bp, dB in the case that there are no boundary conditions specified
    dimension (NoDOF+(1+NoSensParameters))NoMultipliers */ 
    void *dummy1;
    void *dummy2; 
    void *dummy3;
    double *DOFScaling;/* DOF scaling factor for multifield problems - NoDOF numbers*/
    int *SensBVFIndex;/* boundary conditions design velocity field index NoDOF*NoSensParameters 
                      index to ADVF node Data field */
} NodeSpec;

typedef struct{
    struct{
        int ElemIndex; /*0,1...*/
        int SpecIndex; /*0,1,..*/
        int Active;
    }id;
    double *Data;
    int *Nodes;        /*1,2,...  cn=Nodes[ce->Nodes[l]-1] !! */
    double *ht;
    double *hp;
} ElementData;

#define DEFAULT_USER_ARGUMENTS double*,void*,ElementData*,NodeSpec**,NodeData**,double*,int*
#ifndef USER1_ARGUMENTS 
#define USER1_ARGUMENTS DEFAULT_USER_ARGUMENTS
#endif
#ifndef USER2_ARGUMENTS 
#define USER2_ARGUMENTS DEFAULT_USER_ARGUMENTS
#endif
#ifndef USER3_ARGUMENTS 
#define USER3_ARGUMENTS DEFAULT_USER_ARGUMENTS
#endif
#ifndef USER4_ARGUMENTS 
#define USER4_ARGUMENTS DEFAULT_USER_ARGUMENTS
#endif
#ifndef USER5_ARGUMENTS 
#define USER5_ARGUMENTS DEFAULT_USER_ARGUMENTS
#endif
#ifndef USER6_ARGUMENTS 
#define USER6_ARGUMENTS DEFAULT_USER_ARGUMENTS
#endif
#ifndef USER7_ARGUMENTS 
#define USER7_ARGUMENTS DEFAULT_USER_ARGUMENTS
#endif

typedef struct ElementSpec{
    char *Code;
    struct{
        int (*SetElSpec)(struct ElementSpec *,int *,int ,int);
        void (*SKR)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double *,double **);
        void (*SRE)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double *);
        void (*SSE)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double **);
        void (*SHI)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *);
        void (*SPP)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double **,double **);
        void (*ELI)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *);
        void (*PAN)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double **);
        void (*User1)(USER1_ARGUMENTS);
        void (*User2)(USER2_ARGUMENTS);
        void (*User3)(USER3_ARGUMENTS);
        void (*User4)(USER4_ARGUMENTS);
        void (*User5)(USER5_ARGUMENTS);
        void(*SS2)(double *, struct ElementSpec *, ElementData *, NodeSpec **, NodeData **, double *, int *, double **);
        void (*User7)(USER7_ARGUMENTS);
        void (*Tasks)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double *,double **,
              int *,int *,int *,double *,int *, double *);
        void (*CondensationData)(struct ElementSpec *,ElementData *, NodeSpec **,NodeData **, double *rdata,
           int *idata,double **h, double **ha, double **hb,double **dh);
    }user;
    struct{
        int SpecIndex;         /*0,1,2...*/
        int NoDimensions;
        int NoDOFGlobal;
        int NoDOFCondense;
        int NoNodes;
        int NoDomainData;
        int NoSegmentPoints;
        int IntCode;
        int NoTimeStorage;
        int NoElementData;
        int NoIntPoints;
        int NoGPostData;
        int NoNPostData;
        int SymmetricTangent;
        int NoIntPointsA;
        int NoIntPointsB;
        int NoIntPointsC;
        int NoSensNames;
        int ShapeSensitivity;
        int NoIData;
        int NoRData;
        int DefaultIntegrationCode;
        int DemoLimitation;
        int NoAdditionalData;
        int NoCharSwitch;
        int NoIntSwitch;
        int NoDoubleSwitch;
        int CreateDummyNodes;
        int PostIterationCall;
        int Active;      /* 11 0=>freeze all elements with given spec */
        int DOFScaling;  /* weather DOF scaling is supported */
        int EBCSensitivity;
        int SensitivityOrder;
        int WorkingVectorSize;
    }id;
    char *Topology;
    char **DomainDataNames;
    char **GPostNames;
    char **NPostNames;
    int *Segments;
    int *DOFGlobal;
    int *dummy1;
    int *SensPVFIndex; /* */
    int *NoNodeStorage;
    int *NoNodeData;
    double *Data;
    double *IntPoints;
    double *ReferenceNodes;
    int    *NodeSpecs; /* index 0,1,2,.. */
    double *AdditionalData;
    char *OutFileName;
    char *AdditionalNodes;
    char **NodeID;
    char *AdditionalGraphics;
    char **CharSwitch;
    int  *IntSwitch;
    double *DoubleSwitch;
    char **SensitivityNames;
    char *MainTitle;
    char *SubTitle;
    char *SubSubTitle;
    int *SensLowerOrderIndex;/* fi_i, fi_j, fi_k, fi_ij, fi_ik, fi_jk ... for all sens par , all orders */
    char *MMANextStep;
    char *MMAStepBack;
    char *MMAPreIteration;
    char **IDataNames;
    int  *IDataIndex;
    char **RDataNames;
    int  *RDataIndex;
    void *dummy2;
    char *Bibliography;
    void *dummy3;
    void *dummy4;
    double *PostNodeWeights;
    int  *MMAInitialisation;
    double *Version;
    void *dummy5;
    void *dummy6;
    void *dummy7;
    void *dummy8;
    void *dummy9;
    void *dummy10;
    void *dummy11;
    void (*Execute)(char *f,char *m);
    void (*BreakPoint)(char *s, char *c,int m, int p);
    int (*BoundCheck)(int i,int j,int k,char *m,char *s,char *d,char *e);
} ElementSpec;

double *SMTIntPoints(int *icode,int *ngpo);

double *SMTMultiIntPoints(int *icode,int *idata, int *ngpo, int *na,int *nb,int *nc,int all);

double *SMTSetReferenceNodes(char *code);

double SMSDot(double a[],double b[],int n);

double SMSSum(double a[],int n);

_inline double SMSDeltaPart(double *a,int i, int j,int k);

void SMSMove(double a[],double b[],int n);

void SMSZero(double a[],int n);

void SMSRecover(double *hb, double *hg, double *du, int nb, int ne);

void SMSCondense(double **s,double *p,double *hb, double *hg,int nb, int ne);

double SMSIsDOFConstrained(int i,double a,double b);

double Time();

/*sensitivity related definitions*/
#define PD_MAXSensFirstOrder 200  /*max first order parameters in the case of second order analysis
                                  to speed up some algorithms*/
void SensDotProduct(ElementSpec *es,ElementData *ed,NodeSpec **ns,NodeData **nd,int *IData,double *KK,double **s,int symm);

#define SMSKDelta(i,j) (((i)==(j))?(1e0):(0e0))
#define Power SMSPower
#define SMSPower(x,y) (pow((double)(x),(double)(y)))
#define Sqrt(x)        (sqrt((double)(x)))
#define Cbrt(x) ( x<0 ? -pow(-(double)(x),1./3.) : pow((double)(x),1./3.) )
#define Abs(x)        (fabs((double)(x)))
#define Less(a,b,c)        (a<b && b<c)
#define LessEqual(a,b,c)        (a<=b && b<=c)
#define Greater(a,b,c)        (a>b && b>c)
#define GreaterEqual(a,b,c)        (a>=b && b>=c)

/* structures and function definitions for Steering */

#define COPT_LEN 4
#define MAX_DATA_LEN 100

#define OPT_DIVERGENCE "Divergence"
#define OPT_AUTOMATIC "Automatic"
#define OPT_TRUE "True"
#define OPT_FALSE "False"
#define OPT_LIST "List"
#define OPT_IGNORE "Ignore"

#define STAT_N_A "N/A"
#define STAT_NULL "Null"
#define STAT_MAX_ITERATIONS "MaxIterations"
#define STAT_CONVERGENCE "Convergence"
#define STAT_ERR_STATUS_2 "ErrorStatus-2"
#define STAT_ALTERNATE "Alternate"
#define STAT_GLOB_DIVERGENCE "Global divergence"
#define STAT_ERR_AND_DIVERGENCE "ErrorStatus & Divergence"
#define STAT_IGNORE_MINB_ALTERNATE "IgnoredMinBound-Alternate"
#define STAT_IGNORE_MINB_CONV "IgnoredMinBound-Convergence"
#define STAT_MIN_BOUND "MinBound"

#define T_Abort "Abort"
#define T_Interrupt "Interrupt"
#define T_Analyze "Analyze"
#define T_AdaptiveBC "Adaptive BC"
#define T_AdaptivePower "Adaptive Power"
#define T_AdaptiveTime "Adaptive Time"
#define T_MaxBound "MaxBound"
#define T_IgnoredMinBound "IgnoredMinBound"
#define T_MinBound "MinBound"

#define TYPE_INT "int"
#define TYPE_DOUBLE "double"
#define TYPE_LONG "long"

#define MAX_ITER 1000

typedef struct
{
  char Alternate[MAX_DATA_LEN];
  double AlternateList[6]; /* {{_,_},{_,_},{_,_}} */
  int PostIteration; /* T/F (-1:automatic) */
  int IgnoreMinBound; /* T/F */
  int AlternativeTarget; /* T/F */
  char Type[MAX_DATA_LEN];
  int TargetIterations;
  double MinIncrement; /* lambda or time */
  double MaxIncrement;
  double Target;
  int ConvFunctOpt; /* 1,2,3 */
  double LargeDrop;
} ConvergenceOptions;

typedef struct
{
  char Indeterminate[MAX_DATA_LEN];
  char MaxIterations[MAX_DATA_LEN];
  char Convergence[MAX_DATA_LEN];
  char ErrorStatus[MAX_DATA_LEN];
  char Alternate[MAX_DATA_LEN];
  char Divergence[MAX_DATA_LEN];
  char Other[MAX_DATA_LEN];
  int IsOtherDefined; /* T/F */
} ConvergenceStatus;
  
typedef struct
{
  int StepBack;
  int StepForward;
  double Increment;
  char Report[MAX_DATA_LEN];
  ConvergenceStatus StateDiagnostic;
} ConvergenceReturn;

typedef struct
{
  double CurrentDa;
  int OneMoreIter; /* T/F */
  int Converged; /* T/F */
  double* Da;
  int DaLen;
  double DLambda;
  int IsDLambdaNull;
  double Dt;
  int isDtNull;
  ConvergenceStatus StateDiagnostic;
  int FirstAlternateIter;
  int Level;
  int Error;
  int PostIteration; /* T/F */
  int ShortReturn; /* T/F */
  int SubIterMode;
  int GlobalIterMode;
  double* R;
  int RLen;
  int* NoDiscreteEvents;
  int NoDiscreteEventsLen;
  int ConvergedIteration;
  double ResidualError;
  ConvergenceReturn LastReturnValueBeforePostIter;
} IterationHistoryStruct;
/* 19 - last 10 steps for tail window */
/* default {} */

typedef struct
{
  double TimeFrequency;
  double MultiplierFrequency;
  double StepFrequency;
} PutOptions;

typedef struct
{
  double Tolerance;
  int Summation;
  int *IntegerInput;
  double *RealInput;
  int *Elements;
  int NoElements;
  int LastIntegerOutputSize;
  int LastRealOutputSize;
} TaskOptions;

typedef struct
{
  int Failed;
  int Type;
  int IntegerOutputLength;
  int *IntegerOutput;
  int RealOutputLength;
  double *RealOutput;
  int *IntegerOutputForElements;
  double *RealOutputForElements;
  double *PostVector;           /* postprocessing vector length NoNodes */
  double *MatrixGlobal;         /* global  matrix length NoDOF*NoDOF*/
  double *VectorGlobal;         /* global  vector length NoDOF*/
} TaskReturn;

typedef struct
{
  char *SimulationDllFile;
  char *H5InputFile;
  void (*InitializeConvergence)();
  int (*NextStep)(double, double, double);
  int (*StepBack)();
  int (*Convergence)(double, double,char *);
  double (*NewtonIteration)();
  void (*FreeConvergence)();
  void (*Put)(int);
  void (*DumpState)(char *dumpkey);
  void (*RestartState)(char *dumpkey);
  double (*TimeStamp)(char *msg);

  void (*ExportVectorToH5File)(char *, char *, void *, int, char *);
//  ConvergenceOptions ConvergenceOptions;
//  ConvergenceReturn ConvergenceReturn;
//  PutOptions PutOptions;

  int (*Task)(char *code);
//  TaskReturn TaskReturn;
//  TaskOptions TaskOptions;

  int (*Sensitivity)();

  /* nodes and elements selection functions */
  int (*FindNodes)(char *, int **);
  int (*FindNodesPoint)(double *, char *,double tol, int **);
  int (*FindNodesLine)(double **, int, char *,double, int **);
  int (*FindNodesPolygon)(double **, int, char *,double tol, int **);
  int (*FindNodesTetrahedron)(double **, char *,double tol, int **);
  int (*FindNodesHexahedron)(double **, char *,double tol, int **);
  int (*FindNodesFunction)(int (*)(double *, double), char *, int **);
  int (*FindElements)(int *, int ,int , int **);
  void (*Free)(void *);
} SMTStructure;

typedef int (*SimulationDll)(SMTStructure *SMT, ElementSpec** ElemSpecs, ElementData** Elements, NodeSpec** NodeSpecs, NodeData** Nodes, int *IData, double *RData);
int RunDllSimulation(char *dllFileName);


/*Matrix functions definitions*/
void MatrixExponentialOMMF(double m[9], int *r, double mf[9]);
void MatrixExponentialOMMFS(double m[6], int *r, double mf[6]);
void MatrixExponentialOMM2D(double m[5], int *r, double mf[5]);
void MatrixExponentialOMM2DS(double m[4], int *r, double mf[4]);
void MatrixExponentialOMDMF(double m[9], int *r, double mf[9], double dmf[45]);
void MatrixExponentialOMDMFS(double m[6], int *r, double mf[6], double dmf[21]);
void MatrixExponentialOMDM2D(double m[5], int *r, double mf[5], double dmf[15]);
void MatrixExponentialOMDM2DS(double m[4], int *r, double mf[4], double dmf[10]);
void MatrixExponentialOMDDMF(double m[9], int *r, double mf[9], double dmf[45], double ddmf[165]);
void MatrixExponentialOMDDMFS(double m[6], int *r, double mf[6], double dmf[21], double ddmf[74]);
void MatrixExponentialOMDDM2D(double m[5], int *r, double mf[5], double dmf[15], double ddmf[35]);
void MatrixExponentialOMDDM2DS(double m[4], int *r, double mf[4], double dmf[10], double ddmf[21]);
void MatrixLogarithmOMMF(double m[9], int *r, double mf[9]);
void MatrixLogarithmOMMFS(double m[6], int *r, double mf[6]);
void MatrixLogarithmOMM2D(double m[5], int *r, double mf[5]);
void MatrixLogarithmOMM2DS(double m[4], int *r, double mf[4]);
void MatrixLogarithmOMDMF(double m[9], int *r, double mf[9], double dmf[45]);
void MatrixLogarithmOMDMFS(double m[6], int *r, double mf[6], double dmf[21]);
void MatrixLogarithmOMDM2D(double m[5], int *r, double mf[5], double dmf[15]);
void MatrixLogarithmOMDM2DS(double m[4], int *r, double mf[4], double dmf[10]);
void MatrixLogarithmOMDDMF(double m[9], int *r, double mf[9], double dmf[45], double ddmf[165]);
void MatrixLogarithmOMDDMFS(double m[6], int *r, double mf[6], double dmf[21], double ddmf[74]);
void MatrixLogarithmOMDDM2D(double m[5], int *r, double mf[5], double dmf[15], double ddmf[35]);
void MatrixLogarithmOMDDM2DS(double m[4], int *r, double mf[4], double dmf[10], double ddmf[21]);
void MatrixSquareRootOMMF(double m[9], int *r, double mf[9]);
void MatrixSquareRootOMMFS(double m[6], int *r, double mf[6]);
void MatrixSquareRootOMM2D(double m[5], int *r, double mf[5]);
void MatrixSquareRootOMM2DS(double m[4], int *r, double mf[4]);
void MatrixSquareRootOMDMF(double m[9], int *r, double mf[9], double dmf[45]);
void MatrixSquareRootOMDMFS(double m[6], int *r, double mf[6], double dmf[21]);
void MatrixSquareRootOMDM2D(double m[5], int *r, double mf[5], double dmf[15]);
void MatrixSquareRootOMDM2DS(double m[4], int *r, double mf[4], double dmf[10]);
void MatrixSquareRootOMDDMF(double m[9], int *r, double mf[9], double dmf[45], double ddmf[165]);
void MatrixSquareRootOMDDMFS(double m[6], int *r, double mf[6], double dmf[21], double ddmf[74]);
void MatrixSquareRootOMDDM2D(double m[5], int *r, double mf[5], double dmf[15], double ddmf[35]);
void MatrixSquareRootOMDDM2DS(double m[4], int *r, double mf[4], double dmf[10], double ddmf[21]);
void MatrixPowerOMMF(double m[9], double *p, int *r, double mf[9]);
void MatrixPowerOMMFS(double m[6], double *p, int *r, double mf[6]);
void MatrixPowerOMM2D(double m[5], double *p, int *r, double mf[5]);
void MatrixPowerOMM2DS(double m[4], double *p, int *r, double mf[4]);
void MatrixPowerOMDMF(double m[9], double *p, int *r, double mf[9], double dmf[45]);
void MatrixPowerOMDMFS(double m[6], double *p, int *r, double mf[6], double dmf[21]);
void MatrixPowerOMDM2D(double m[5], double *p, int *r, double mf[5], double dmf[15]);
void MatrixPowerOMDM2DS(double m[4], double *p, int *r, double mf[4], double dmf[10]);
void MatrixPowerOMDDMF(double m[9], double *p, int *r, double mf[9], double dmf[45], double ddmf[165]);
void MatrixPowerOMDDMFS(double m[6], double *p, int *r, double mf[6], double dmf[21], double ddmf[74]);
void MatrixPowerOMDDM2D(double m[5], double *p, int *r, double mf[5], double dmf[15], double ddmf[35]);
void MatrixPowerOMDDM2DS(double m[4], double *p, int *r, double mf[4], double dmf[10], double ddmf[21]);
void MatrixPowerSeriesOMMF(double m[9], double *cp, int *n, int *r, double mf[9]);
void MatrixPowerSeriesOMMFS(double m[6], double *cp, int *n, int *r, double mf[6]);
void MatrixPowerSeriesOMM2D(double m[5], double *cp, int *n, int *r, double mf[5]);
void MatrixPowerSeriesOMM2DS(double m[4], double *cp, int *n, int *r, double mf[4]);
void MatrixPowerSeriesOMDMF(double m[9], double *cp, int *n, int *r, double mf[9], double dmf[45]);
void MatrixPowerSeriesOMDMFS(double m[6], double *cp, int *n, int *r, double mf[6], double dmf[21]);
void MatrixPowerSeriesOMDM2D(double m[5], double *cp, int *n, int *r, double mf[5], double dmf[15]);
void MatrixPowerSeriesOMDM2DS(double m[4], double *cp, int *n, int *r, double mf[4], double dmf[10]);
void MatrixPowerSeriesOMDDMF(double m[9], double *cp, int *n, int *r, double mf[9], double dmf[45], double ddmf[165]);
void MatrixPowerSeriesOMDDMFS(double m[6], double *cp, int *n, int *r, double mf[6], double dmf[21], double ddmf[74]);
void MatrixPowerSeriesOMDDM2D(double m[5], double *cp, int *n, int *r, double mf[5], double dmf[15], double ddmf[35]);
void MatrixPowerSeriesOMDDM2DS(double m[4], double *cp, int *n, int *r, double mf[4], double dmf[10], double ddmf[21]);

void MatrixFunctionCompressMMF(double m[3][3],double mc[9]);
void MatrixFunctionUncompressMFMF(double mf[9],double mffull[3][3]);
void MatrixFunctionUncompressDMFMF(double dmf[45],double dmffull[3][3][3][3]);
void MatrixFunctionUncompressDDMFMF(double ddmf[165],double ddmffull[3][3][3][3][3][3]);
void MatrixFunctionCompressMMFS(double m[3][3],double mc[6]);
void MatrixFunctionUncompressMFMFS(double mf[6],double mffull[3][3]);
void MatrixFunctionUncompressDMFMFS(double dmf[21],double dmffull[3][3][3][3]);
void MatrixFunctionUncompressDDMFMFS(double ddmf[74],double ddmffull[3][3][3][3][3][3]);
void MatrixFunctionCompressMM2D(double m[3][3],double mc[5]);
void MatrixFunctionUncompressMFM2D(double mf[5],double mffull[3][3]);
void MatrixFunctionUncompressDMFM2D(double dmf[15],double dmffull[3][3][3][3]);
void MatrixFunctionUncompressDDMFM2D(double ddmf[35],double ddmffull[3][3][3][3][3][3]);
void MatrixFunctionCompressMM2DS(double m[3][3],double mc[4]);
void MatrixFunctionUncompressMFM2DS(double mf[4],double mffull[3][3]);
void MatrixFunctionUncompressDMFM2DS(double dmf[10],double dmffull[3][3][3][3]);
void MatrixFunctionUncompressDDMFM2DS(double ddmf[21],double ddmffull[3][3][3][3][3][3]);



#endif
